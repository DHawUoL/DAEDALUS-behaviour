function chains = run_blocked_MCMC(start, nIter, Sig1_0, Sig2_0, flat)

%% Likelihoods:
F=@(params)lhood1(params,data,xdata,ydata,X,Xfull,coeff,tvec,plim,lx1,lx2,flat);

% F: log posterior in *natural* params [alpha,k1,k*] with simulator try/catch
% start: 1x3
% Sig1_0: 2x2 for [alpha,k1], Sig2_0: 1x1 for [k*]
% flat only to pass through / label runs

rng('shuffle');
d1=2; d2=1; d=3;
x = start(:)';  % [alpha,k1,k*]
X = zeros(nIter, d); X(1,:) = x;
L = F(x);       % must be finite (handle failures in F)

% book-keeping
acc1=0; acc2=0; accG=0;
C1 = Sig1_0; C2 = Sig2_0; CG = blkdiag(Sig1_0, Sig2_0);
m1 = x(1:2); m2 = x(3);   mG = x;

% adaptation schedule
burn = max(2000, round(0.2*nIter));
epsJ = 1e-9;               % jitter
scale1 = (2.38^2)/d1; 
scale2 = (2.38^2)/d2; 
scaleG = (2.38^2)/d;

for t = 2:nIter
    u = rand;
    xprop = x; Lprop = -Inf;

    if u < 0.45
        % --- block 1: [alpha,k1]
        z = mvnrnd([0 0], scale1*C1 + epsJ*eye(d1));
        xprop(1:2) = x(1:2) + z;
        Lprop = F(xprop);
        if log(rand) < Lprop - L
            x = xprop; L = Lprop; acc1 = acc1+1;
        end

        % online mean/cov (after burn) for block 1
        if t > burn
            m1 = m1 + (x(1:2)-m1)/(t-burn);
            C1 = C1 + ((x(1:2)-m1)'*(x(1:2)-m1) - C1)/(t-burn+1);
        end

    elseif u < 0.90
        % --- block 2: [k*]
        z = mvnrnd(0, scale2*C2 + epsJ);
        xprop(3) = x(3) + z;
        Lprop = F(xprop);
        if log(rand) < Lprop - L
            x = xprop; L = Lprop; acc2 = acc2+1;
        end

        if t > burn
            m2 = m2 + (x(3)-m2)/(t-burn);
            C2 = C2 + ((x(3)-m2)'*(x(3)-m2) - C2)/(t-burn+1);
        end

    else
        % --- global occasional move
        z = mvnrnd([0 0 0], scaleG*CG + epsJ*eye(d));
        xprop = x + z;
        Lprop = F(xprop);
        if log(rand) < Lprop - L
            x = xprop; L = Lprop; accG = accG+1;
        end

        if t > burn
            mG = mG + (x-mG)/(t-burn);
            CG = CG + ((x-mG)'*(x-mG) - CG)/(t-burn+1);
        end
    end

    X(t,:) = x;
end

chains.samples = X;
chains.acc = [acc1, acc2, accG]/nIter;
chains.flat = flat;
chains.Sig1 = C1; chains.Sig2 = C2; chains.SigG = CG;
end

function f=lhood1(params,data,xdata,ydata,Xfit,Xfull,coeff,tvec,plim,lx1,lx2,flat)
k=.003;
x0=50;
w=@(x)(1/(1+exp(-k*(x-x0))));
ymodel=sim2fit(params(1:end),data,xdata,Xfit,Xfull,coeff,tvec,lx1,lx2);%(1:end-1)
ymodel=max(ymodel,50);
lhood = -ymodel' + ydata .* log(ymodel') - gammaln(ydata+1);
log_prior = log(betapdf((params(1) - 0.2) / 0.7, 2, 2));%(params(1) - 0.3) / 0.6
f = flat*sum(lhood) + sum(log(unif(params(2:end), plim(:,2:end)))) + log_prior;
end

%% SIMULATION %%

function [f,rhohat]=sim2fit(params,data,xdata,Xfit,Xfull,coeff,tvec,lx1,lx2)
R0=2.2;
tvec(1)=-145;
alpha=params([1,1,1]);

a2=0.8072;
b2=-11.1656;
a3=-0.3567;
b3=4.3175;
av0=0.4695;
bv0=-5.5286;
ks=params(3);
reducedParams=[1,params(2),a2*ks+b2,a3*ks+b3,av0*ks+bv0];
%BH
%Fitting link function:
%[pr,be,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=bePrepCovid19(data,R0,ones(1,lx2-2),[1,params(2:end)],coeff,zeros(5,lx2),alpha);
[pr,be,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=bePrepCovid19(data,R0,ones(1,lx2-2),reducedParams,coeff,zeros(5,lx2),alpha);
pr.xfull=Xfull;

Wfit=Xfit.^(1/pr.a);
%Fit to ocupancy:
%[simu,~,~,rhohat]=beRunCovid19(pr,be,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,Wfit,tvec(1:lx2+1),0,data);
%Fit to admissions:
[simu,simu2,~,rhohat]=beRunCovid19(pr,be,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,Wfit,tvec(1:lx2+1),0,data);
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
