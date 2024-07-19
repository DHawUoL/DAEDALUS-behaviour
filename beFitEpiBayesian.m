function [xsto, outsto, history, accept_rate,covmat]=beFitEpiBayesian(ydata,X,data,thetaIn,Xfull,coeff)
%% Parameters to fit:
%R0, t0 - while no mitigation
%alpha, explicit p's (as previous deltas)
%alpha, logistic parameters

intrinsic=1;%=1 for behaviour as feed in fully open economy for "original" DAEDALUS
numPCA=size(coeff,1);%Number of x's in logistic regression, including H

%%
monthDur=[1,31,29,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31];
monthStart=cumsum(monthDur);
tvec=monthStart;
%
tvec=tvec(1:20);
lt=length(tvec);

X=X(:,1:lt-1);
Xfull=Xfull(:,1:lt-1);

[lx1,lx2]=size(X);

%xdata=85:578;%End of July 21
xdata=85:(tvec(17)-1);%End of April 21

ydata=ydata(0+(1:length(xdata)));
%If data is just England:
%ydata=ydata*(sum(data.Npop)/56286961);%England, mid-2019 (ONS)

x0=[thetaIn];

lb=[0,.15,0,-20,-200,0];
ub=[.6,1,1,0,0,.3];

plim=[ub;lb];

n=2e4;
sigma=1;
fixinds=[];
blockind=[];
displ=false;

F=@(params)lhood1(params,data,xdata,ydata,X,intrinsic,Xfull,coeff,tvec,plim,lx1,lx2);
tStart = tic;
[xsto, outsto, history, accept_rate,covmat] = MCMC_adaptive(F, x0, n, sigma, fixinds, blockind, displ);
tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));
end

%% LIKELIHOOD %%

function f=lhood1(params,data,xdata,ydata,Xfit,intrinsic,Xfull,coeff,tvec,plim,lx1,lx2)
ymodel=sim2fit(params(1:end-1),data,xdata,Xfit,intrinsic,Xfull,coeff,tvec,lx1,lx2);
%Deselect low occupancy:
thresh=0;
ymodel(ydata<thresh)=[];
ydata(ydata<thresh)=[];
%}
%
overdisp=params(end);
sd=sqrt(ymodel+ymodel.^2*overdisp);
%sd(sd<1)=1;
lhood=normpdf(ydata',ymodel,sd);%Neg Bin approx
%}
f=sum(log(lhood))+sum(log(unif(params,plim)));
end

%% SIMULATION %%

function [f,rhohat]=sim2fit(params,data,xdata,Xfit,intrinsic,Xfull,coeff,tvec,lx1,lx2)
R0=2.75;%params(1);
tvec(1)=-70;%params(2);
alpha=params([1,2,2]);
%BH
[pr,be,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=bePrepCovid19(data,R0,ones(1,lx2-2),params(3:end),coeff,zeros(5,lx2),alpha);

pr.xfull=Xfull;
Wfit=Xfit.^(1/pr.a);
if intrinsic==1
    %Fit to ocupancy:
    %[simu,~,~]=heRunCovid19(pr,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,[ones(1,length(tvec)-1)],tvec(1:numInt+1),0,data);
    %Fit to admissions:
    %%BH [simu,simu2,~,rhohat]=beRunCovid19(pr,be,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,[ones(1,length(tvec)-1)],tvec(1:lx2+1),0,data);
    [simu,simu2,~,rhohat]=beRunCovid19(pr,be,vx,n,NN,NNbar,beta,[ones(1,length(tvec)-1)],tvec(1:lx2+1),0,data);

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