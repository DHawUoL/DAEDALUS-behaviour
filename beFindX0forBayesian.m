function [f,g]=beFindX0forBayesian(ydata,X,data,thetaIn,Xfull,coeff)
%% Parameters to fit:
%R0, t0 - while no mitigation
%alpha, explicit p's (as previous deltas)
%alpha, logistic parameters

intrinsic=0;
numPCA=size(coeff,1);%Number of x's in logistic regression, including H

%%
monthDur=[1,31,29,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31];
monthStart=cumsum(monthDur);
tvec=[-68.7792,32,88.2788,monthStart(5:end)];
%
tvec=tvec(1:end);%[-40,365:368];
lt=length(tvec);
X=X(:,1:lt-1);
Xfull=Xfull(:,1:lt-1);

[lx1,lx2]=size(X);

xdata=70:tvec(end);%579;%245;%70:...
ydata=ydata(0+(1:length(xdata)));
ydata=ydata*(sum(data.Npop)/56286961);%England, mid-2019 (ONS)

x0=thetaIn;
%{
%R0, t0, t1, alpha, p1:
lb=[2.5,-150,61,0,0];%zeros(1,lx-2)];
ub=[3.5,-40,92,1,1];%zeros(1,lx-2)];
%}
%{
%alpha, p's:
lb=[.2,zeros(1,16),0];%zeros(1,lx-2)];
ub=[.5,ones(1,16),.2];%zeros(1,lx-2)];
%}
%{
%alpha, feedback params L1, k's, H0:
lb=[.1,0,-5*ones(1,numPCA),-20,0];%alpha, L, kvec, H0
ub=[.5,1,5*ones(1,numPCA),20,1];
%}
%
%2 alphas, feedback params L1, k's, H0:
lb=[.15,.15,0,-20*ones(1,numPCA),-60,0];%alpha, L, kvec, H0
ub=[.5,1,1,20*ones(1,numPCA),60,.2];
%}
%x0=[x0,.1];
plim=[ub;lb];

lb=thetaIn*.9; ub=thetaIn*1.1;

F=@(params)lhood1(params,data,xdata,ydata,X,intrinsic,Xfull,coeff,tvec,plim,lx1,lx2);

nparams=length(lb);%-1;%Discount nbinm variance - keep at .2
%phi=.1;
inc=1e4;%Increment
ssize=5e3;%Sample size
sampleMat=zeros(nparams,ssize);
sample=randsample(inc,ssize)/inc;%Between 0 1nd 1

for i=1:nparams
    sampleMat(i,:)=lb(i)+(ub(i)-lb(i))*sample(randperm(ssize));%randsample(inc,ssize)/inc;
end
lhoods=zeros(ssize,1);
for i=1:ssize
    lhoods(i)=F([sampleMat(:,i)']);%,phi]);
end
[~,indmax]=max(lhoods);
f=sampleMat(:,indmax)';
g=F([f]);%,phi]);
end

%% LIKELIHOOD %%

function f=lhood1(params,data,xdata,ydata,Xfit,intrinsic,Xfull,coeff,tvec,plim,lx1,lx2)
%{
ymodel=sim2fit(params,data,xdata,Xfit,intrinsic,Xfull,coeff,tvec);
%{
%Deselect low occupancy:
thresh=1e3;
ymodel(ydata<thresh)=[];
ydata(ydata<thresh)=[];
%}
%
%lhood=poisspdf(ydata,ymodel);%Poisson
lhood=normpdf(ydata,ymodel,sqrt(ymodel));%Bin approx
%}
%
ymodel=sim2fit(params(1:end-1),data,xdata,Xfit,intrinsic,Xfull,coeff,tvec,lx1,lx2);
%Deselect low occupancy:
thresh=0;
ymodel(ydata<thresh)=[];
ydata(ydata<thresh)=[];
%}
%
alpha=params(end);
sd=sqrt(ymodel+ymodel.^2*alpha);
%sd(sd<1)=1;
lhood=normpdf(ydata,ymodel,sd);%Neg Bin approx
%}
f=sum(log(lhood))+sum(log(unif(params,plim)));
end

%% SIMULATION %%

function [f,rhohat]=sim2fit(params,data,xdata,Xfit,intrinsic,Xfull,coeff,tvec,lx1,lx2)
R0=2.8025;%params(1);%2.0733;%
%t0=params(2);
%t1=params(3);
%tvec(1)=t0;
%tvec(3)=t1;
alpha=params(1:2);
beProps=[zeros(5,2),repmat(params(2:end),5,1)];
%[pr,be,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=bePrepCovid19(data,lx2,R0,[ones(1,lx2-4),1.25,1.5],0,coeff,beProps,alpha);
[pr,be,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=bePrepCovid19(data,lx2,R0,[ones(1,lx2-2)],params(3:end),coeff,ones(5,lx2),alpha);
pr.xfull=Xfull;
Wfit=Xfit.^(1/pr.a);
if intrinsic==1
    %Fit to ocupancy:
    %[simu,~,~]=heRunCovid19(pr,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,[ones(1,length(tvec)-1)],tvec(1:numInt+1),0,data);
    %Fit to admissions:
    [simu,simu2,~,rhohat]=beRunCovid19(pr,be,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,Xfit,tvec(1:lx2+1),0,data);
else
    [simu,~,~,rhohat]=beRunCovid19(pr,be,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,Wfit,tvec,0,data);
end

t=simu(:,1)';
%Fit to ocupancy:
h=simu(:,4)';
%Fit to admissions:
%h=simu2';

f=interp1(t,h,xdata); 

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