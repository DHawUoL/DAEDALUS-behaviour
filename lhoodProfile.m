function f=lhoodProfile(ydata,X,data,thetaIn,Xfull,coeff)
intrinsic=1;%=1 for behaviour as feed in fully open economy for "original" DAEDALUS
numPCA=size(coeff,1);%Number of x's in logistic regression, including H
tvec=[1,2,61,92,134,141,148,155,162,169,176,186,200,211,218,223,227,236,250,258,266,271,279,294,310,322,330,338,349,370,384,397,407,418,433,445,463,468,474,491,504,517,540,561];%,567,575];
xdata=85:tvec(end-7);
lt=length(tvec);
X=X(:,1:lt-1);
Xfull=Xfull(:,1:lt-1);
[lx1,lx2]=size(X);
ydata=ydata(0+(1:length(xdata)));

x0=[thetaIn];%,.3];
%Fitting link function:
lb=[.1,.5,-10,-10,-10,-50,0];%a,a,m,k,k,k,h0 poptim5
ub=[.7,1,20,20,20,50,.4];
plim=[ub;lb];

F=@(params)lhood1(params,data,xdata,ydata,X,intrinsic,Xfull,coeff,tvec,plim,lx1,lx2);
%{
i = 2;  % the parameter you're testing
x0_prof = x0;
vals = linspace(plim(2,i), plim(1,i), 100);
loglik = zeros(size(vals));

for j = 1:length(vals)
    x0_prof(i) = vals(j);
    try
        loglik(j) = F(x0_prof);
    catch ME
        fprintf('Failure at %f: %s\n', vals(j), ME.message);
        loglik(j) = NaN;
    end
end

figure;
plot(vals, loglik, 'LineWidth', 2);
xlabel(sprintf('Parameter %d', i));
ylabel('Log-posterior');
title(sprintf('Likelihood profile of parameter %d', i));
%}
%
%i = 1;  % index of the parameter you want to profile
for i=1:length(thetaIn)
x0_prof = x0;  % best starting point
vals = linspace(plim(2,i), plim(1,i), 100);  % from lower to upper bound
loglik = zeros(size(vals));
for j = 1:length(vals)
    x0_prof(i) = vals(j);
    loglik(j) = F(x0_prof);
end
figure;
plot(vals, loglik, 'LineWidth', 2);
xlabel(sprintf('Parameter %d', i));
ylabel('Log-posterior');
title(sprintf('Likelihood profile of parameter %d', i));
end
%}
end

function f=lhood1(params,data,xdata,ydata,Xfit,intrinsic,Xfull,coeff,tvec,plim,lx1,lx2)
k=.005;
x0=5e2;
w=@(x)(1/(1+exp(-k*(x-x0))));
overdisp=params(end);

ymodel=sim2fit(params(1:end-1),data,xdata,Xfit,intrinsic,Xfull,coeff,tvec,lx1,lx2);%(1:end-1)
lhood=(1-arrayfun(w,ymodel)).*log(normpdf(ydata',ymodel,sqrt(ymodel+ymodel.^2*overdisp)))+arrayfun(w,ymodel).*log(normpdf(ydata',ymodel,sqrt(ymodel)));

f = sum(lhood) + sum(log(unif(params, plim)));
end

%% SIMULATION %%

function [f,rhohat]=sim2fit(params,data,xdata,Xfit,intrinsic,Xfull,coeff,tvec,lx1,lx2)
R0=2.75;%1.9;%2.75;%params(1);
tvec(1)=-84;%-70;%-195;%-206;%-195;%Seasonal;-206;%-70;%-85;%-70;%params(2);
alpha=params([1,1,1]);
%tvec(5:end)=tvec(5:end)+params(end);
%BH
%Fitting link function:
[pr,be,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=bePrepCovid19(data,R0,ones(1,lx2-2),params(2:end),coeff,zeros(5,lx2),alpha);
%Interaction term:
pr.xfull=Xfull;
Wfit=Xfit.^(1/pr.a);
if intrinsic==1
    %Fitting link function:
    [simu,simu2,~,rhohat]=beRunCovid19(pr,be,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,[ones(1,length(tvec)-1)],tvec(1:lx2+1),0,data);
else
    [simu,simu2,~,rhohat]=beRunCovid19(pr,be,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,Wfit,tvec(1:lx2+1),0,data);
end
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