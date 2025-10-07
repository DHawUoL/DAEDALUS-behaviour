function plotmat=bePlotFrequentist(ydata,X,data,xsto,Xfull,coeff,pointEst,Diag)
%[xsto, outsto, history, accept_rate,covmat]=fitEpiBayesian(dataOcc,ones(1,19),dataUK1,[0.9058   -0.7800   -7.1241   15.4944    0.0000    0.0243],X2(4:22,:)',1);
hlag=0;
addmodifier=1;
intrinsic=1;
projection=0;

%%

if projection==1
    tvec=[1,2,61,93,134,141,148,155,162,169,176,186,200,211,218,223,227,236,250,258,266,271,279,294,310,322,330,338,349,370,384,397,407,418,433,445,463,468,474,491,504,517,540,561,567,575];
    xdata=85:tvec(end);%-2
else
    %tvec=[1,2,61,93,134,141,148,155,162,169,176,186,200,211,218,223,227,236,250,258,266,271,279,294,310,322,330,338,349,370,384,397,407,418,433,445,463,468,474,491,504,517,540,561,567,575];
    %Footfall:
    %tvec=[1,2,61,94,[134,141,148,155,162,169,176,186,200,211,218,223,227,236,250,258,266,271,279,294,310,322,330,338,349,370,384,397,407,418,433,445,463,468,474,491,504,517,540,561,567,575]+hlag];
    %Stringency:
    tvec=[1,2,61,91,[127,134,141,148,153,155,162,167,169,175,176,186,200,211,216,218,223,227,230,236,250,258,265,266,271,279,288,294,305,310,322,330,337,338,349,354,356,361,370,372,384,397,407,418,433,445,454,463,468,474,491,503,504,517,540,561,566,567,575]+hlag];

    xdata=85:tvec(end-7);%-2
end

lt=length(tvec);
X=X(:,1:lt-1);
Xfull=Xfull(:,1:lt-1);

[lx1,lx2]=size(X);

ydata=ydata(0+(1:length(xdata)));
%If data is just England:
%ydata=ydata*(sum(data.Npop)/56286961);%England, mid-2019 (ONS)
%If data is just England:
ydata=ydata*(sum(data.Npop)/56286961);%England, mid-2019 (ONS)

fun=@(params)sim2fit(params,data,xdata,X,intrinsic,Xfull,coeff,tvec,lx1,lx2,0);

%% Generate sample
l1        = 200;
%sample=drawParamSamplesx(pointEst, se, l1);
sample=drawParamSamples(Diag, l1, struct('lb',[0 0 -25 -25], 'ub',[1 1 25 25], 'truncate','none'));

%% Generate epidemics from sample etc.
y1=fun(sample(1,:));
l2=length(y1);
plotmat=zeros(l1,l2);
plotmat(1,:)=y1;
for i=2:l1
    plotmat(i,:)=fun(sample(i,:));
end
[ypoint,rhohat]=fun(pointEst);%fun(mean(xsto,1));%

ymax=max(max(ydata),max(max(plotmat)));
factor=5e3;%1e2*ceil(ymax/1e2);

if addmodifier==1
    %Periods, real H values:
    tvecPlus=[1,tvec(2:end)];
    value=zeros(1,tvec(end));%numPeriods);
    value2=pointEst(1)*ones(1,tvec(end));
    for i=1:lt-1%******** 3:
        ti=round(tvecPlus(i)):round(tvecPlus(i+1));
        value(ti)=rhohat(i);
    end
    value(1:round(tvec(3))-1)=0;%Kicks in at tvec(3)
    %}
    %factor=3e3;
    ydata=ydata/factor;
    plotmat=plotmat/factor;
    ypoint=ypoint/factor;
    ymax=1;
end

%% Plot:
fs=10; lw=2;
cmap=lines(7);
col1=cmap(1,:);
col2=.5*[1,1,1];
%xvec=[1,32,60,tvec(4:end)];
%xvec=xvec(1:2:end);
monthDur=[1,31,29,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31,31];
monthStart=cumsum(monthDur);
if projection==1
    xnames={'Jan 2020','Feb 2020','Mar 2020','Apr 2020','May 2020','Jun 2020','Jul 2020','Aug 2020','Sep 2020','Oct 2020','Nov 2020','Dec 2020',...
    'Jan 2021','Feb 2021','March 2021','Apr 2021','May 2021','Jun 2021','July2021','Aug 2021','Sep 2021','Oct 2021','Nov 2021','Dec 2021','Jan 2022'};
    %xnamesAdd={'Aug 2021','Sep 2021','Oct 2021','Nov 2021','Dec 2021','Jan 2022'};
else
    xnames={'Jan 2020','Feb 2020','Mar 2020','Apr 2020','May 2020','Jun 2020','Jul 2020','Aug 2020','Sep 2020','Oct 2020','Nov 2020','Dec 2020',...
    'Jan 2021','Feb 2021','March 2021','Apr 2021'};
end
xvec=monthStart(1:length(xnames));
%xnames=[xnames];%,xnamesAdd];
%xvec=monthStart(1:length(xnames));
figure
hold on
%plot(plotmat');
bar(xdata,ydata,'FaceColor',col2,'EdgeColor',col2,'LineWidth',.01);
if projection==1
    %plot(tvec(timeThresh)*[1,1],[0,factor],'k:','linewidth',2)
    plot((tvec(end)+1)*[1,1],[0,factor],'k:','linewidth',2)
end
plot_distribution_prctile(xdata,plotmat,'color',col1,'prctile',(0:25:100));%5
%plot(xdata,plotmat(1,:),'color',col1,'linewidth',2);%5

h5=plot(xdata,ypoint,'k-','linewidth',2);

%for i=[1,32,61,92,122,153,183,214,245,275,306,336,367,398,426,457,487,518,548,579,610,640,671,701,731]
%    plot(i*[1,1],[0,1.25*max(max(ydata),max(max(plotmat)))],'k-','linewidth',0.01);
%end
if addmodifier==1
    h1=plot([-1,-1],[-1-1],'linewidth',2,'color',col2);
    h2=plot([-1,-1],[-1-1],'linewidth',2,'color',col1);
    h3=plot(xdata,value(xdata),'k--','linewidth',2);
    h4=plot(xdata,value2(xdata),'--','linewidth',2,'color',[.5,0,0]);
    hleglines=[h1(1),h2(1),h5,h3(1),h4(1)];
    legend(hleglines,'Data','Model fit (Bayesian)','Model fit (point est.)','p (point est.)','\alpha','location','northeastoutside')
end
xlim([xdata(1),xdata(end)]);
ylim([0,1])%.25*max(ydata)]);
xticks(xvec)
xticklabels(xnames)
xtickangle(45)
box on;
grid on;
xlabel('Time');
ylabel('Hospital Admissions/5k');
%title('Model Fit');
end

function [f,rhohat]=sim2fit(params,data,xdata,Xfit,intrinsic,Xfull,coeff,tvec,lx1,lx2,plotRun)
R0=2.8;%2.2;
tvec(1)=-60;%-70;%-195;%-206;%-195;%Seasonal;-206;%-70;%-85;%-70;%params(2);
alpha=params([1,1,1]);
propIn=1;
%{
a1=-.814;
b1=8.0161;
a2=-.8067;
b2=-8.0887;
%}
%ks=params(2);
%reducedParams=[1,a1*ks+b1,a2*ks+b2,params(3),0];
%reducedParams=[1,ks,0,params(2),0];
reducedParams=[1,params(2:end)];%,0,params(3),0];
%BH
%Fitting link function:
[pr,be,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=bePrepCovid19(data,R0,ones(1,lx2-2),reducedParams,coeff,zeros(5,lx2),alpha,propIn);
pr.leak=0; pr.xfull=Xfull; be.BiFirstFit=1; pr.phi2=0.1668;
%[pr,be,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=bePrepCovid19(data,R0,ones(1,lx2-2),[params(2:end),0.8036*params(3)-0.3232],coeff,zeros(5,lx2),alpha);
%Interaction term:
pr.xfull=Xfull;
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

function sample = drawParamSamples(Diag, n, opts)
% drawParamSamples  Draw parameter samples respecting covariance.
%   sample = drawParamSamples(Diag, n) uses Diag.pcov (Gauss–Newton).
%   sample = drawParamSamples(Diag, n, opts) with fields:
%       .useRobust   (false)  -> use Diag.pcov_robust if available
%       .jitter      (1e-10)  -> small diagonal added if needed
%       .truncate    ('clip') -> 'clip' or 'none' to enforce bounds
%       .lb, .ub     ([] )    -> override bounds (defaults to Diag.lb/ub if present)
%       .seed        ([] )    -> rng seed (scalar). If empty, no seeding.
%
% Returns:
%   sample : n x K matrix, each row a parameter draw.

    if nargin < 3, opts = struct; end
    if ~isfield(opts,'useRobust'), opts.useRobust = false; end
    if ~isfield(opts,'jitter'),    opts.jitter    = 1e-10; end
    if ~isfield(opts,'truncate'),  opts.truncate  = 'clip'; end

    % Pick covariance
    C = [];
    if opts.useRobust && isfield(Diag,'pcov_robust') && ~isempty(Diag.pcov_robust)
        C = Diag.pcov_robust;
    elseif isfield(Diag,'pcov') && ~isempty(Diag.pcov)
        C = Diag.pcov;
    elseif isfield(Diag,'se') && ~isempty(Diag.se)
        C = diag(Diag.se(:).^2);  % fallback to diagonal only
    else
        error('No covariance or standard errors found in Diag.');
    end

    % Symmetrize and make PD (with graceful fallback)
    C  = (C + C.')/2;
    K  = numel(Diag.p_hat);
    L  = [];
    [L,p] = chol(C,'lower');
    if p ~= 0
        % Eigen clean + jitter
        [V,D] = eig(C);
        D     = diag(max(diag(D),0));
        Cfix  = V*D*V' + opts.jitter*eye(K);
        [L,p] = chol((Cfix+Cfix.')/2,'lower');
        if p ~= 0
            % Last resort: keep only diagonal variances (still correlated mean)
            warning('Covariance not PD after repair; falling back to diagonal.');
            Cdiag = diag(max(diag(C), 0)) + opts.jitter*eye(K);
            L = chol(Cdiag,'lower');
        end
    end

    % Seed if requested
    if isfield(opts,'seed') && ~isempty(opts.seed)
        rng(opts.seed);
    end

    % Draw
    Z      = randn(K, n);             % K x n
    mu     = Diag.p_hat(:);           % K x 1
    draws  = mu + L*Z;                % K x n
    sample = draws.';                 % n x K

    % ----- Optional truncation to bounds -----
    % Fill in opts.lb / opts.ub from Diag if not provided
    if ~isfield(opts,'lb') || isempty(opts.lb)
        if isfield(Diag,'lb') && ~isempty(Diag.lb)
            opts.lb = Diag.lb(:).';   % row vector 1×K
        else
            opts.lb = [];
        end
    end
    if ~isfield(opts,'ub') || isempty(opts.ub)
        if isfield(Diag,'ub') && ~isempty(Diag.ub)
            opts.ub = Diag.ub(:).';   % row vector 1×K
        else
            opts.ub = [];
        end
    end
    
    % Ensure row-shape for broadcasting
    if ~isempty(opts.lb), opts.lb = reshape(opts.lb,1,[]); end
    if ~isempty(opts.ub), opts.ub = reshape(opts.ub,1,[]); end
    
    % Clip if requested
    if strcmpi(opts.truncate,'clip')
        % If your MATLAB supports implicit expansion (R2016b+), these two lines are enough:
        if ~isempty(opts.lb), sample = max(sample, opts.lb); end
        if ~isempty(opts.ub), sample = min(sample, opts.ub); end
    
        % If you need compatibility with older MATLAB, uncomment these instead:
        % if ~isempty(opts.lb), sample = max(sample, repmat(opts.lb, n, 1)); end
        % if ~isempty(opts.ub), sample = min(sample, repmat(opts.ub, n, 1)); end
    end

end


function sample = drawParamSamplesx(pointEst, se, n)
    % pointEst : 1x4 vector of parameter estimates
    % se       : 1x4 vector of standard errors
    % n        : number of samples
    % sample   : n x 4 matrix of draws
    
    % replicate means and stds to size (n,4)
    muMat = repmat(pointEst(:).', n, 1);   % n x 4
    seMat = repmat(se(:).', n, 1);         % n x 4
    
    % independent normal draws
    sample = muMat + seMat .* randn(n, numel(pointEst));
end
