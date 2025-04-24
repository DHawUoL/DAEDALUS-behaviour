function [xsto, outsto, history, accept_rate,covmat]=beFitEpiBayesian(ydata,X,data,thetaIn,Xfull,coeff)
%function [x_all, x_final, acc_rates]=beFitEpiBayesian_alpham(ydata,X,data,thetaIn,Xfull,coeff)

a=.2293;%Also hard-coded in sim2fit
b=.6942;
if a > 0
    alpha_min = -b / a;
    alpha_max = (1 - b) / a;
    alpha_max=max(alpha_max,1);
else
    alpha_min = (1 - b) / a;
    alpha_max = -b / a;
    alpha_min=min(alpha_min,1);
end


%% Parameters to fit:
%R0, t0 - while no mitigation
%alpha, explicit p's (as previous deltas)
%alpha, logistic parameters

intrinsic=1;%=1 for behaviour as feed in fully open economy for "original" DAEDALUS
numPCA=size(coeff,1);%Number of x's in logistic regression, including H

%%
tvec=[1,2,61,92,134,141,148,155,162,169,176,186,200,211,218,223,227,236,250,258,266,271,279,294,310,322,330,338,349,370,384,397,407,418,433,445,463,468,474,491,504,517,540,561,567,575];
xdata=85:tvec(end-7);%-2
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
lb=[alpha_min,-10,-10,-10,-50,0];%a,a,m,k,k,k,h0 poptim5
ub=[alpha_max,20,20,20,50,.4];

%expandFactor=.25;
%lb(2)=x0(2)*(1-expandFactor); ub(2)=x0(2)*(1+expandFactor);%Was 3
%lb(5)=x0(5)*(1-expandFactor); ub(5)=x0(5)*(1+expandFactor);%Was 6

plim=[ub;lb];

n=5e3;
sigma=.2;
fixinds=[];
blockind=[];
displ=false;

F=@(params)lhood1(params,data,xdata,ydata,X,intrinsic,Xfull,coeff,tvec,plim,lx1,lx2);
tStart = tic;
fixinds = [];%6;%[2,3,6; x0([2,3,6])];
[xsto, outsto, history, accept_rate,covmat] = MCMC_adaptive(F, x0, n, sigma, fixinds, blockind, displ);
%{
blockList = {
    [2 3 6];  % Most promising movers
    [1];    % Still stuck
    [4 5];    % Not moving, low risk
};
sigmas = [.5, 0.5, 1.0];  % must match length(blockList)

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
k=.005;
x0=5e2;
w=@(x)1;%(1/(1+exp(-k*(x-x0))));
overdisp=params(end);

ymodel=sim2fit(params(1:end-1),data,xdata,Xfit,intrinsic,Xfull,coeff,tvec,lx1,lx2);%(1:end-1)
ymodel=max(ymodel,1);
lhood=(1-arrayfun(w,ymodel)).*log(normpdf(ydata',ymodel,sqrt(ymodel+ymodel.^2*overdisp)))+arrayfun(w,ymodel).*log(normpdf(ydata',ymodel,sqrt(ymodel)));

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

a=.2293;%Also hard-coded in sim2fit
b=.6942;

R0=2.75;%1.9;%2.75;%params(1);
tvec(1)=-84;%-70;%-195;%-206;%-195;%Seasonal;-206;%-70;%-85;%-70;%params(2);
alpha=params([1,1,1]);
%tvec(5:end)=tvec(5:end)+params(end);
%BH
%Fitting link function:
[pr,be,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=bePrepCovid19(data,R0,ones(1,lx2-2),[a*params(1)+b,params(2:end)],coeff,zeros(5,lx2),alpha);
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