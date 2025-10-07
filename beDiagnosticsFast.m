function Diag = beDiagnosticsFast(ydata, X, data, pointEst, Xfull, coeff, paramNames, opts)
    % Fast diagnostics: build J by forward differences around pointEst.
    % Keeps the same sim2fit signature you already use.
    
    if nargin < 8 || isempty(opts), opts = struct; end
    ds   = getfielddef(opts,'downsample',1);         % take every ds-th timepoint for J
    par  = getfielddef(opts,'useParallel',false);    % parfor over params
    hrel = getfielddef(opts,'relStep',1e-3);         % relative FD step
    habs = getfielddef(opts,'absStep',1e-4);         % min absolute step
    lb   = getfielddef(opts,'lb',[]);                % optional bounds
    ub   = getfielddef(opts,'ub',[]);
    plotRun = 0;
    
    % ---- Time grid (same as your full fit) ----
    hlag = -7; projection = 0;
    if projection==1
        tvec  = [1,2,61,93,134,141,148,155,162,169,176,186,200,211,218,223,227,236,250,258,266,271,279,294,310,322,330,338,349,370,384,397,407,418,433,445,463,468,474,491,504,517,540,561,567,575];
        xdata = 85:tvec(end);
    else
        tvec  = [1,2,61,91,[127,134,141,148,153,155,162,167,169,175,176,186,200,211,216,218,223,227,230,236,250,258,265,266,271,279,288,294,305,310,322,330,337,338,349,354,356,361,370,372,384,397,407,418,433,445,454,463,468,474,491,503,504,517,540,561,566,567,575]+hlag];
        xdata = 85:tvec(end-7);
    end
    lt    = length(tvec);
    X     = X(:,1:lt-1);
    Xfull = Xfull(:,1:lt-1);
    [~,lx2] = size(X);
    
    % ---- Data (England scaling) ----
    ydata  = ydata(1:length(xdata));
    ydata  = ydata * (sum(data.Npop)/56286961);
    ymean=mean(ydata)/5e3;
    idxJ   = 1:ds:numel(xdata);     % downsample for J if requested
    y_ds   = ydata(idxJ);
    
    % ---- Model handle ----
    fun = @(p) sim2fit(p, data, xdata, X, 1, Xfull, coeff, tvec, size(X,1), lx2, plotRun, ymean); % column out
    
    % ---- Base evaluation at pointEst ----
    p0   = pointEst(:).';
    f0   = fun(p0);
    f0   = f0(:);
    f0_ds = f0(idxJ);
    N    = numel(f0_ds);
    K    = numel(p0);
    
    % ---- Bound-aware forward differences for J ----
    J = zeros(N, K);
    if isempty(lb), lb = -inf(1,K); end
    if isempty(ub), ub =  inf(1,K); end
    
    steps = max(hrel*max(abs(p0),1), habs);     % per-parameter step
    % push inward if sitting at/near a bound
    for k=1:K
        if p0(k) + steps(k) > ub(k)
            steps(k) = -steps(k);               % forward step inward
        elseif p0(k) + steps(k) < lb(k)
            steps(k) = -steps(k);               % forward step inward
        end
    end
    
% ---- Bound-aware forward differences for J ----
J = zeros(N, K);
if isempty(lb), lb = -inf(1,K); end
if isempty(ub), ub =  inf(1,K); end

steps = max(hrel*max(abs(p0),1), habs);  % per-parameter step
for k=1:K
    % push inward if sitting near a bound
    if p0(k) + steps(k) > ub(k)
        steps(k) = -steps(k);
    elseif p0(k) + steps(k) < lb(k)
        steps(k) = -steps(k);
    end
end

if par
    % Parallel loop over parameters
    parfor k = 1:K
        pk = p0;
        pk(k) = min(max(pk(k) + steps(k), lb(k)), ub(k));
        fk = fun(pk); fk = fk(:);
        J(:,k) = (fk(idxJ) - f0_ds) / steps(k);
    end
else
    % Serial fallback
    for k = 1:K
        pk = p0;
        pk(k) = min(max(pk(k) + steps(k), lb(k)), ub(k));
        fk = fun(pk); fk = fk(:);
        J(:,k) = (fk(idxJ) - f0_ds) / steps(k);
    end
end

    
    % ---- Residuals on the DS grid (for s^2 and cov) ----
    res_ds = f0_ds - y_ds(:);
    rss    = sum(res_ds.^2);
    dof    = max(N - K, 1);
    s2     = rss / dof;
    
    % ---- SVD / conditioning / ridge ----
    [U,S,V]  = svd(J,'econ'); %#ok<ASGLU>
    s        = diag(S);
    condJ    = s(1)/max(s(end),eps);
    v_ridge  = V(:,end);
    frac_info = s(end)^2 / sum(s.^2);
    
% ---- Correlations & VIF ----
Z  = zscore(J,0,1);
R  = corrcoef(Z);
if rcond(R) < 1e-10
    VIF = diag(pinv(R));
else
    VIF = diag(inv(R));
end
maxAbsCorr = max(max(abs(R - diag(diag(R)))));

    
    % ---- Gauss–Newton covariance (DS grid) ----
    H     = J.'*J;
    pcov  = s2 * pinv(H);
    se    = sqrt(max(diag(pcov),0));
    tcrit = tinv(0.975, dof);
    CI    = [p0(:) - tcrit.*se(:),  p0(:) + tcrit.*se(:)];
    
    % ---- Pack results ----
    if nargin < 7 || isempty(paramNames)
        paramNames = arrayfun(@(i)sprintf('p%d',i), 1:K, 'uni',0);
    end
    Diag = struct();
    Diag.p_hat       = p0(:).';
    Diag.residuals   = (f0 - ydata(:));
    Diag.J           = J;
    Diag.rss         = rss;
    Diag.sigma2      = s2;
    Diag.S           = s;
    Diag.V           = V;
    Diag.H           = H;
    Diag.pcov        = pcov;
    Diag.se          = se(:);
    Diag.CI          = CI;
    Diag.paramNames  = paramNames(:);
    Diag.condJ       = condJ;
    Diag.v_ridge     = v_ridge(:);
    Diag.frac_info   = frac_info;
    Diag.R           = R;
    Diag.maxAbsCorr  = maxAbsCorr;
    Diag.VIF         = VIF(:);
    Diag.dof         = dof;
    Diag.tcrit       = tcrit;
    Diag.idxJ        = idxJ(:);
    Diag.f0          = f0;
    Diag.ydata       = ydata(:);
    Diag.lb          = lb(:).';
    Diag.ub          = ub(:).';
end

function val = getfielddef(s,fn,def)
if isfield(s,fn), val = s.(fn); else, val = def; end
end

function parfor_k(K,body) % serial fallback for “parfor”
for k=1:K, body(k); end
end

function [f,rhohat]=sim2fit(params,data,xdata,Xfit,intrinsic,Xfull,coeff,tvec,lx1,lx2,plotRun,ymean)
R0=2.8;%2.2;
tvec(1)=-80;%-59;
alpha=params([1,1,1]);
propIn=1;
%reducedParams=[1,params(2:end)];
reducedParams=[1,params(2),0,params(3),4.3792];
%BH
%Fitting link function:
[pr,be,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=bePrepCovid19(data,R0,ones(1,lx2-2),reducedParams,coeff,zeros(5,lx2),alpha,propIn);
pr.leak=0; pr.xfull=Xfull; be.BiFirstFit=1; pr.phi2=0;%.186;
%[pr,be,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=bePrepCovid19(data,R0,ones(1,lx2-2),[params(2:end),0.8036*params(3)-0.3232],coeff,zeros(5,lx2),alpha);
%Interaction term:
%{
delta = params(4);                 % behaviour lead/lag in days
idx_start = 7;                             % first window to shift
Xfull_shift = shift_driver_by_days(Xfull, tvec, delta, idx_start, 'previous');
pr.xfull = Xfull_shift;
%}
pr.xfull=Xfull;
pr.ymean=ymean;

%Fitting individual p's:
%[pr,be,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=bePrepCovid19(data,R0,ones(1,size(Xfull,2)-2),ones(1,3),1,zeros(5,lx2),alpha);%repmat([1,1,params(2:end)]
%pr.xfull=[1,1,1,params(2:end)];%Use xfull as the value of p


Wfit=Xfit.^(1/pr.a);
if intrinsic==1
    %Fit to ocupancy:
    %[simu,~,~]=heRunCovid19(pr,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,[ones(1,length(tvec)-1)],tvec(1:numInt+1),0,data);
    %Fit to admissions:
    %%BH
    %Fitting link function:
    [simu,simu2,~,rhohat]=beRunCovid19(pr,be,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,[ones(1,length(tvec)-1)],tvec(1:lx2+1),plotRun,data);
    %Fitting individual p's:
    %[simu,simu2,~,rhohat]=beRunCovid19(pr,be,vx,n,NN,NNbar,beta,[ones(1,length(tvec)-1)],tvec(1:lx2+1),0,data);

else
    %Fit to ocupancy:
    %[simu,~,~,rhohat]=beRunCovid19(pr,be,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,Wfit,tvec(1:lx2+1),0,data);
    %Fit to admissions:
    [simu,simu2,~,rhohat]=beRunCovid19(pr,be,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,Wfit,tvec(1:lx2+1),0,data);
end
t=simu(:,1)';

%Fit to ocupancy:
%h=simu(:,4)';
%Fit to admissions:
h=simu2';

f=interp1(t,h,xdata); 

%plot(simu(:,1),simu2);
%plot(xdata,f)

%f(isinf(f))=-1e6;
%f(isnan(f))=-1e6;

end