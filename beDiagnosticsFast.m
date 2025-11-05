function Diag = beDiagnosticsFast(ydata, data, pointEst, Xfull, paramNames, scheme, tvec, b0, opts)
    % Fast diagnostics: build J by forward differences around pointEst.
    % Keeps the same sim2fit signature you already use.

    % ---- Time grid (same as your full fit) ----
    hlag = 0; projection = 0;
    xdata=85:tvec(end);

    if nargin < 9 || isempty(opts), opts = struct; end
    ds   = getfielddef(opts,'downsample',1);         % take every ds-th timepoint for J
    par  = getfielddef(opts,'useParallel',false);    % parfor over params
    hrel = getfielddef(opts,'relStep',1e-3);         % relative FD step
    habs = getfielddef(opts,'absStep',1e-4);         % min absolute step
    lb   = getfielddef(opts,'lb',[]);                % optional bounds
    ub   = getfielddef(opts,'ub',[]);
    plotRun = 0;

    [~,lx2] = size(Xfull);
    
    % ---- Data (England scaling) ----
    ydata  = ydata(1:length(xdata));
    ydata  = ydata * (sum(data.Npop)/56286961);
    ymean=mean(ydata)/5e3;
    idxJ   = 1:ds:numel(xdata);     % downsample for J if requested
    y_ds   = ydata(idxJ);
    
    % ---- Model handle ----
    fun = @(p) sim2fit_global(p, data, xdata, 1, Xfull, tvec, lx2, plotRun, ymean, b0); % column out
    
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




Diag.scheme   = scheme;                  % 'ols' | 'wls1' | 'wls05' | 'poiss'
Diag.k        = numel(pointEst);
Diag.n        = numel(ydata);

% Base errors
yhat = f0;
e   = yhat(:) - ydata(:);
RSS = sum(e.^2);
WRSS_1  = sum( (1./(1+ydata(:))).*e.^2 );         % for 'wls1'
WRSS_05 = sum( (1./sqrt(1+ydata(:))).*e.^2 );     % for 'wls05'
DEV     = 2*sum( max(yhat(:),eps) - ydata(:) + ydata(:).*log( max(ydata(:),eps)./max(yhat(:),eps) ) );

% Store raw fit metrics (useful for diagnostics)
Diag.RSS      = RSS;
Diag.WRSS_1   = WRSS_1;
Diag.WRSS_05  = WRSS_05;
Diag.DEV      = DEV;

% Information criteria (use scheme-appropriate surrogate for -2 logL)
switch lower(scheme)
    case {'ols','wls1','wls05'}
        n   = Diag.n;
        ll2 = n*log(RSS/n);     % constants cancel for model comparison on same data
        Diag.AIC = ll2 + 2*Diag.k;
        Diag.BIC = ll2 + Diag.k*log(n);
    case 'poiss'
        Diag.AIC = DEV + 2*Diag.k;
        Diag.BIC = DEV + Diag.k*log(Diag.n);
end
end

function val = getfielddef(s,fn,def)
if isfield(s,fn), val = s.(fn); else, val = def; end
end

function parfor_k(K,body) % serial fallback for “parfor”
for k=1:K, body(k); end
end

function [f,rhohat]=sim2fit(params,data,xdata,Xfit,intrinsic,Xfull,coeff,tvec,lx1,lx2,plotRun,ymean,arg)
R0=2.8;
tvec(1)=-80;
alpha=params([1,1,1]);
propIn=1;
if length(params)==5
    %arg=b0
    k1=params(2);
    k2=params(3);
    k3=params(4);
    delta=params(5);%delta>eps_safe
    v0pca=-dot([k1,k2],arg)-delta;

    reducedParams=[1,k1,k2,k3,v0pca];
    %softplus = @(z) log1p(exp(-abs(z))) + max(z,0);
    %reducedParams = [1, params(2:4), -softplus(params(5))];
elseif length(params)==4
    pc=2;%PC1 OR PC2
    kstar=params(2);
    k3=params(3);
    delta=params(4);%delta>eps_safe
    v0pca=-kstar*arg(pc)-delta;
    if pc==1
        reducedParams=[1,kstar,0,k3,v0pca];   
    else
        reducedParams=[1,0,kstar,k3,v0pca];                      
    end
elseif length(params)==3
    pc=2;%PC1 OR PC2
    kstar=params(2);
    k3=params(3);
    delta=16.8943;
    v0pca=-kstar*arg(pc)-delta;
    if pc==1
        reducedParams=[1,kstar,0,k3,v0pca];   
    else
        reducedParams=[1,0,kstar,k3,v0pca];                      
    end
elseif length(params)==2
    alpha=0.5832*ones(1,3);
    pc=2;%PC1 OR PC2
    kstar=params(1);
    k3=params(2);
    delta=16.8943;
    v0pca=-kstar*arg(pc)-delta;
    if pc==1
        reducedParams=[1,kstar,0,k3,v0pca];   
    else
        reducedParams=[1,0,kstar,k3,v0pca];                      
    end
end
%BH
%Fitting link function:
[pr,be,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=bePrepCovid19(data,R0,ones(1,lx2-2),reducedParams,coeff,zeros(5,lx2),alpha,propIn);
pr.leak=0; pr.xfull=Xfull; be.BiFirstFit=1; pr.phi2=0;%.186;
pr.xfull=Xfull;
pr.ymean=ymean;

Wfit=Xfit.^(1/pr.a);
if intrinsic==1
    %Fit to ocupancy:
    %[simu,~,~]=heRunCovid19(pr,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,[ones(1,length(tvec)-1)],tvec(1:numInt+1),0,data);
    %Fit to admissions:
    %%BH
    %Fitting link function:
    [simu,simu2,~,rhohat]=beRunCovid19(pr,be,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,[ones(1,length(tvec)-1)],tvec(1:lx2+1),plotRun,data);
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
end