function [xsto, outsto, history, accept_rate,covmat]=beFitEpiBayesian(ydata,X,data,x0,Xfull,coeff)
%function [x_all, x_final, acc_rates]=beFitEpiBayesian(ydata,X,data,thetaIn,Xfull,coeff)
hlag=0;
%% Parameters to fit:
%R0, t0 - while no mitigation
%alpha, explicit p's (as previous deltas)
%alpha, logistic parameters

%%
%Stringency:
tvec=[1,2,61,94,[127,134,141,148,153,155,162,167,169,175,176,186,200,211,216,218,223,227,230,236,250,258,265,266,271,279,288,294,305,310,322,330,337,338,349,354,356,361,370,372,384,397,407,418,433,445,454,463,468,474,491,503,504,517,540,561,566,567,575]+hlag];
xdata=85:tvec(end-15);%-2 -7

lt=length(tvec);
X=X(:,1:lt-1);
Xfull=Xfull(:,1:lt-1);

[lx1,lx2]=size(X);

ydata=ydata(0+(1:length(xdata)));
%If data is just England:
ydata=ydata*(sum(data.Npop)/56286961);%England, mid-2019 (ONS)

%Fitting link function:
%lb=[.01,.3,1,-40,1,-50,.001];%a,a,m,k,k,k,h0 poptim5
%ub=[.7,1,100,-1,40,50,.1];
lb=[.1,-15,-15,-15,-30];%a,a,m,k,k,k,h0 poptim5
ub=[.7,15,15,15,30];

plim=[ub;lb];

n=1e4;
sigma=1;
fixinds=[];
blockind=[];
displ=false;
flat=0.01;

F=@(params)lhood1(params,data,xdata,ydata,X,Xfull,coeff,tvec,plim,lx1,lx2,flat);
tStart = tic;
fixinds = [];%[2,3,6; x0([2,3,6])];
[xsto, outsto, history, accept_rate,covmat] = MCMC_adaptive(F, x0, n, sigma, fixinds, blockind, displ);
%{
blockList = {
    [1,2];  % Most promising movers
    [3,6];    % Still stuck
    [4,5];    % Not moving, low risk
    [7];
};
sigmas = [.1, 0.3, .5, 1];  % must match length(blockList)

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

function f=lhood1(params,data,xdata,ydata,Xfit,Xfull,coeff,tvec,plim,lx1,lx2,flat)
ymodel=sim2fit(params(1:end),data,xdata,Xfit,Xfull,coeff,tvec,lx1,lx2);%(1:end-1)
ymodel=max(ymodel,50);
lhood = -ymodel' + ydata .* log(ymodel') - gammaln(ydata+1);
% Beta(2,2) for a symmetric prior around 0.6
log_prior = log(betapdf((params(1) - 0.2) / 0.7, 2, 2));%(params(1) - 0.3) / 0.6
f = flat*sum(lhood) + sum(log(unif(params(2:end), plim(:,2:end)))) + log_prior;
end

%% SIMULATION %%

function [f,rhohat]=sim2fit(params,data,xdata,Xfit,Xfull,coeff,tvec,lx1,lx2)
R0=2.2;%1.9;%2.75;%params(1);
tvec(1)=-145;%-70;%-195;%-206;%-195;%Seasonal;-206;%-70;%-85;%-70;%params(2);
alpha=params([1,1,1]);
%tvec(5:end)=tvec(5:end)+params(end);
%BH
%Fitting link function:
[pr,be,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=bePrepCovid19(data,R0,ones(1,lx2-2),[1,params(2:end)],coeff,zeros(5,lx2),alpha);
%[pr,be,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=bePrepCovid19(data,R0,ones(1,lx2-2),[params(2:end),0.8036*params(3)-0.3232],coeff,zeros(5,lx2),alpha);
%Interaction term:
pr.xfull=Xfull;
%Fitting individual p's:
%[pr,be,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=bePrepCovid19(data,R0,ones(1,size(Xfull,2)-2),ones(1,3),1,zeros(5,lx2),alpha);%repmat([1,1,params(2:end)]
%pr.xfull=[1,1,1,params(2:end)];%Use xfull as the value of p


Wfit=Xfit.^(1/pr.a);

%Fit to ocupancy:
%[simu,~,~]=heRunCovid19(pr,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,[ones(1,length(tvec)-1)],tvec(1:numInt+1),0,data);
%Fit to admissions:
%%BH
%Fitting link function:
[simu,simu2,~,rhohat]=beRunCovid19(pr,be,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,[ones(1,length(tvec)-1)],tvec(1:lx2+1),1,data);
%Fitting individual p's:
%[simu,simu2,~,rhohat]=beRunCovid19(pr,be,vx,n,NN,NNbar,beta,[ones(1,length(tvec)-1)],tvec(1:lx2+1),0,data);

t=simu(:,1)';

%Fit to ocupancy:
%h=simu(:,4)';
%Fit to admissions:
h=simu2';
f=interp1(t,h,xdata); 

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
ac
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