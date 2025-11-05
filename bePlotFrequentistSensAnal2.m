function f = bePlotFrequentistSensAnal2(ydata,X,data,Xfull,coeff,Diag,coeff_fit,mu_fit)
% Xfull is 2×T raw drivers [trust; stringency]; row 2 is ignored in the model
% coeff_fit (2×2) and mu_fit (1×2) are from the FIT (computed on cols 4:end)

hlag       = 0;
intrinsic  = 1;
projection = 0;
b0=[-0.5836   -0.7019];

% ----- timeline (unchanged) -----
if projection==1
    tvec  = [1,2,61,91,134,141,148,155,162,169,176,186,200,211,218,223,227,236,250,258,266,271,279,294,310,322,330,338,349,370,384,397,407,418,433,445,463,468,474,491,504,517,540,561,567,575];
    xdata = 85:tvec(end);
else
    tvec  = [1,2,61,91,127,134,141,148,153,155,162,167,169,175,176,186,200,211,216,218,223,227,230,236,[250,258,265,266,271,279,288,294,305,310,322,330,337,338,349,354,356,361,370,372,384,397,407,418,433,445,454,463,468,474,491,503,504,517,540,561,566,567,575]+hlag];
    xdata = 85:tvec(end-7);
end
lt    = numel(tvec);
X     = X(:,1:lt-1);
Xfull = Xfull(:,1:lt-1);                      % 2×(lt-1)
[lx1,lx2] = size(X);

% admissions + England scaling
ydata = ydata(1:numel(xdata));
ydata = ydata*(sum(data.Npop)/56286961);
ymean = 0;

% model wrapper (expects Xfull already “PCA’d” to 1×(lt-1) in row 1, and an ignored row 2)
fun = @(params,Xfull_drv) sim2fit(params, data, xdata, X, intrinsic, Xfull_drv, ...
                                  coeff, tvec, lx1, lx2, 0, ymean, b0);

% ===== helpers =====
project_pc1 = @(raw12) (raw12(1:2,:)' - mu_fit) * coeff_fit(:,1);   % returns [T×1] PC1
mk_drv      = @(pc1)   [pc1.' ; zeros(1,size(Xfull,2))];            % 2×T driver, row2 ignored

% Draw parameter sample around frequentist estimates
pointEst = Diag.p_hat(:).';
se       = Diag.se(:).';
nS       = 00;
ParamS   = repmat(pointEst, nS, 1) + randn(nS,numel(pointEst)).*repmat(se,nS,1);

% ===== scenarios =====
from_idx   = 22;                             % when to start trust changes
reduce_vec = 1 - (0.00:0.02:0.10);           % 0%,2%,…10% trust reduction
nSc        = numel(reduce_vec);
plotmatset = cell(1,nSc);
maxvals    = zeros(1,nSc+1);

% Baseline PC1 from FIT basis (no re-PCA)
pc1_base   = project_pc1(Xfull);
Xdrv_base  = mk_drv(pc1_base);

% Baseline uncertainty band (for legend/reference)
l2    = numel(xdata);
PM0   = zeros(nS,l2);
for j=1:nS
    PM0(j,:) = fun(ParamS(j,:), Xdrv_base);
end
maxvals(1) = max(PM0(:));

% Trust-only counterfactuals
for i=1:nSc
    Xraw_cf = Xfull;                               % 2×T raw drivers
    Xraw_cf(1,from_idx:end) = Xraw_cf(1,from_idx:end) * reduce_vec(i);  % scale trust only
    pc1_cf  = project_pc1(Xraw_cf);                % same coeff/mean as FIT
    Xdrv_cf = mk_drv(pc1_cf);
    PMi     = zeros(nS,l2);
    for j=1:nS
        PMi(j,:) = fun(ParamS(j,:), Xdrv_cf);
    end
    plotmatset{i} = PMi;
    maxvals(i+1)  = max(PMi(:));
end

% Optional: trust & stringency combined example (-2% both)
Xboth = Xfull;
Xboth(1,from_idx:end) = Xboth(1,from_idx:end)*0.98;
Xboth(2,from_idx:end) = Xboth(2,from_idx:end)*0.98;
pc1_b  = project_pc1(Xboth);
Xdrv_b = mk_drv(pc1_b);
PMboth = zeros(nS,l2);
for j=1:nS
    PMboth(j,:) = fun(ParamS(j,:), Xdrv_b);
end

% ===== Plot =====
factor = 5e3;
fs=10; lw=2;
cmap=lines(7); col2=.5*[1,1,1];
monthDur=[1,31,29,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31,31];
monthStart=cumsum(monthDur);
if projection==1
    xnames={'Jan 2020','Feb 2020','Mar 2020','Apr 2020','May 2020','Jun 2020','Jul 2020','Aug 2020','Sep 2020','Oct 2020','Nov 2020','Dec 2020', ...
            'Jan 2021','Feb 2021','March 2021','Apr 2021','May 2021','Jun 2021','July2021','Aug 2021','Sep 2021','Oct 2021','Nov 2021','Dec 2021','Jan 2022'};
else
    xnames={'Jan 2020','Feb 2020','Mar 2020','Apr 2020','May 2020','Jun 2020','Jul 2020','Aug 2020','Sep 2020','Oct 2020','Nov 2020','Dec 2020', ...
            'Jan 2021','Feb 2021','March 2021','Apr 2021'};
end
xvec=monthStart(1:numel(xnames));

figure; hold on
bar(xdata, ydata/factor, 'FaceColor',col2,'EdgeColor',col2,'LineWidth',.01);

% baseline band (black)
plot_distribution_prctile(xdata, PM0/factor, 'color', 0*[1,1,1], 'prctile', (0:25:100));

% trust & stringency (-2% both) band (light gray)
plot_distribution_prctile(xdata, PMboth/factor, 'color', .8*[1,1,1], 'prctile', (0:25:100));

% trust-only counterfactuals
parcmap = parula(nSc);
for i=1:nSc
    plot_distribution_prctile(xdata, plotmatset{i}/factor, 'color', parcmap(i,:), 'prctile', (0:25:100));
end

xlim([xdata(1), xdata(end)]);
ylim([0, 1.2]);
xticks(xvec); xticklabels(xnames); xtickangle(45);
box on; grid on;
xlabel('Time'); ylabel('Hospital Admissions/5k');

legend({'Data','Model fit','-2% trust & stringency', '-0%/-2%/-4%/... trust only'}, 'Location','northeastoutside');

f = 1;
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

function sample = drawParamSamples(pointEst, se, n)
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