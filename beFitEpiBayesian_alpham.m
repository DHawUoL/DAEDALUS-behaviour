function [xsto, outsto, history, accept_rate,covmat]=beFitEpiBayesian_alpham(ydata,X,data,thetaIn,Xfull,coeff)
%function [x_all, x_final, acc_rates]=beFitEpiBayesian_alpham(ydata,X,data,thetaIn,Xfull,coeff)
hlag=0;
%{
a=.6121;%Also hard-coded in sim2fit
b=.5987;
if a > 0
    alpha_min = -b / a;
    alpha_max = (1 - b) / a;
    alpha_max=min(alpha_max,1);
    alpha_min=max(alpha_min,.01);
else
    alpha_min = (1 - b) / a;
    alpha_max = -b / a;
    alpha_min=max(alpha_min,.01);
    alpha_max=min(alpha_max,1);
end
%}
alpha_min=.01;
alpha_max=.7;

%% Parameters to fit:
%R0, t0 - while no mitigation
%alpha, explicit p's (as previous deltas)
%alpha, logistic parameters

intrinsic=1;%=1 for behaviour as feed in fully open economy for "original" DAEDALUS
numPCA=size(coeff,1);%Number of x's in logistic regression, including H

%%
%Footfall:
%tvec=[1,2,61,94,[134,141,148,155,162,169,176,186,200,211,218,223,227,236,250,258,266,271,279,294,310,322,330,338,349,370,384,397,407,418,433,445,463,468,474,491,504,517,540,561,567,575]+hlag];
%Stringency:
tvec=[1,2,61,94,[127,134,141,148,153,155,162,167,169,175,176,186,200,211,216,218,223,227,230,236,250,258,265,266,271,279,288,294,305,310,322,330,337,338,349,354,356,361,370,372,384,397,407,418,433,445,454,463,468,474,491,503,504,517,540,561,566,567,575]+hlag];

xdata=85:tvec(end-15);%-2 -7
%}
lt=length(tvec);
X=X(:,1:lt-1);
Xfull=Xfull(:,1:lt-1);

[lx1,lx2]=size(X);

ydata=ydata(0+(1:length(xdata)));
%If data is just England:
ydata=ydata*(sum(data.Npop)/56286961);%England, mid-2019 (ONS)

x0=thetaIn;%[thetaIn([1,3:length(thetaIn)])];%Remove m

%Fitting link function:
%lb=[alpha_min,.1,-100,0];%a,a,m,k,k,k,h0 poptim5
%ub=[alpha_max,10,-1,.4];
%Fitting link function:
%lb=[alpha_min,.1,-5,0,0];%a,a,m,k,k,k,h0 poptim5
%ub=[alpha_max,5,-.01,2,20];
%Fitting link function:
lb=[.01,.3,1,-40,-50,.001];%a,a,m,k,k,k,h0 poptim5
ub=[.7,1,100,-1,50,.1];

%expandFactor=.25;
%lb(2)=x0(2)*(1-expandFactor); ub(2)=x0(2)*(1+expandFactor);%Was 3
%lb(5)=x0(5)*(1-expandFactor); ub(5)=x0(5)*(1+expandFactor);%Was 6

plim=[ub;lb];

n=1e5;
sigma=1;
fixinds=[];
blockind=[];
displ=false;

F=@(params)lhood1(params,data,xdata,ydata,X,intrinsic,Xfull,coeff,tvec,plim,lx1,lx2);
tStart = tic;
fixinds = [];%6;%[2,3,6; x0([2,3,6])];
[xsto, outsto, history, accept_rate,covmat] = MCMC_adaptive(F, x0, n, sigma, fixinds, blockind, displ);
%{
blockList = {
    [1];  % Most promising movers
    [1];    % Still stuck
    [3,4];    % Not moving, low risk
    [5]
};
sigmas = [1,1,.1,1];  % must match length(blockList)

% MCMC settings
sigma = 1;
nPerBlock = 500;
nCycles = 3;
displ = false;

% Your likelihood function F from before
%[x_all, x_final, acc_rates] = blockMCMC_wrapper(F, x0, sigma, nPerBlock, nCycles, blockList, displ);
[x_all, x_final, acc_rates] = blockMCMC_wrapper(F, x0, sigmas, nPerBlock, nCycles, blockList, displ);

%}
end

%% LIKELIHOOD %%

function f=lhood1(params,data,xdata,ydata,Xfit,intrinsic,Xfull,coeff,tvec,plim,lx1,lx2)
k=.003;
x0=50;
w=@(x)(1/(1+exp(-k*(x-x0))));
overdisp=params(end);

ymodel=sim2fit(params(1:end-1),data,xdata,Xfit,intrinsic,Xfull,coeff,tvec,lx1,lx2);%(1:end-1)
ymodel=max(ymodel,50);
%lhood=(1-arrayfun(w,ymodel)).*log(normpdf(ydata',ymodel,sqrt(ymodel+ymodel.^2*overdisp)))+arrayfun(w,ymodel).*log(normpdf(ydata',ymodel,sqrt(ymodel)));
wvec=arrayfun(w,ymodel);
logL_P=log(normpdf(ydata',ymodel,sqrt(ymodel)));
logL_NB=log(normpdf(ydata',ymodel,sqrt(ymodel+ymodel.^2*overdisp)));
lhood = log((1 - wvec) .* exp(logL_NB) + wvec .* exp(logL_P)); %NB low
%lhood = log((1 - wvec) .* exp(logL_P) + wvec .* exp(logL_NB)); %NB high
%{
df = 4;
sd=sqrt(ymodel);
resid = (ydata' - ymodel) ./ sd;
lhood = log(tpdf(resid, df)) - log(sd);
%}
f = .01*sum(lhood) + sum(log(unif(params, plim)));
end

%% SIMULATION %%

function [f,rhohat]=sim2fit(params,data,xdata,Xfit,intrinsic,Xfull,coeff,tvec,lx1,lx2)

a=-0.3771;%k3=ak2+b
b=-0.8084;
%a=0.3948;%v0=ak2+b
%b=2.1942;

R0=2.2;%2.75;%1.9;%2.75;%params(1);
tvec(1)=-145;%-84;%-70;%-195;%-206;%-195;%Seasonal;-206;%-70;%-85;%-70;%params(2);
alpha=params([1,1,1]);
%tvec(5:end)=tvec(5:end)+params(end);
%BH
%Fitting link function:
[pr,be,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=bePrepCovid19(data,R0,ones(1,lx2-2),[params(2:4),a*params(4)+b,params(5)],coeff,zeros(5,lx2),alpha);
pr.xfull=Xfull;

Wfit=Xfit.^(1/pr.a);
if intrinsic==1
    %Fit to admissions:
    %%BH
    %Fitting link function:
    [simu,simu2,~,rhohat]=beRunCovid19(pr,be,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,[ones(1,length(tvec)-1)],tvec(1:lx2+1),0,data);
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

%% PRIOR %%

function f=unif(x,plim)
val=1./(plim(1,:)-plim(2,:));
in=(x-plim(1,:)).*(x-plim(2,:));
in(in>0)=0;
in(in<0)=1;
f=val.*in;
end

function [all_xsto, final_x, accept_rates] = blockMCMC_wrapper(F, x0, sigma, nPerBlock, nCycles, blockList, displ)

% Initialize
x = x0;
d = length(x0);
all_xsto = [];
accept_rates = zeros(nCycles, length(blockList));

for cycle = 1:nCycles
    fprintf('\nCycle %d / %d\n', cycle, nCycles);
    
    for b = 1:length(blockList)
        block = blockList{b};
        fixed = setdiff(1:d, block);
        
        fixinds = [fixed; x(fixed)];
        
        % Run MCMC on this block
        fprintf('  Block %d: sampling params [%s]\n', b, num2str(block));

        %[xsto, ~, ~, accept_rate, ~] = MCMC_adaptive(F, x, nPerBlock, sigma, fixinds, [], displ);
        %
        if numel(sigma) == 1
            sigma_b = sigma;
        else
            sigma_b = sigma(b);
        end
        [xsto, ~, ~, accept_rate, ~] = MCMC_adaptive(F, x, nPerBlock, sigma_b, fixinds, [], displ);
        %}

        % Update current state
        x = xsto(end, :);
        all_xsto = [all_xsto; xsto];  % optionally store all samples
        accept_rates(cycle, b) = accept_rate;
    end
end

final_x = x;

end


%% Likelihood junk:
%
%ymodel=sim2fit(params(1:end-1),data,xdata,Xfit,intrinsic,Xfull,coeff,tvec);
%p=params(end);
%p=1/(1+2*psi);
%lhood=binopdf(round(ydata),round(ymodel/p),p);%Bin
%lhood=normpdf(ydata,ymodel,sqrt(ymodel/p/(1-p)));%Bin approx
%lhood=nbinpdf(ydata,ymodel*p/(1-p),p);%Neg Bin: r,p
%lhood=normpdf(ydata,ymodel,sqrt(ymodel/p));%Neg Bin approx
%{
ymodel=sim2fit(params(1:end-2),data,xdata,Xfit,intrinsic,Xfull,coeff,tvec);
p=params(end-1);
psi=params(end);
sd=ymodel*p*(1-p)+p^2*ymodel.^2*psi^2;
sd(sd<1)=1;
lhood=normpdf(ydata,ymodel*p,sd);%Neg Bin approx
%}