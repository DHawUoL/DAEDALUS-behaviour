function [poptim,Ypred,delta,resnorm]=beFitEpiStart(ydata,X,data,thetaIn,Xfull,coeff)%,lb,ub)
%% Parameters to fit:
%R0, t0 - while no mitigation
%alpha, explicit p's (as previous deltas)
%alpha, logistic parameters

intrinsic=1;%=1 for behaviour as feed in fully open economy for "original" DAEDALUS
numPCA=size(coeff,1);%Number of x's in logistic regression, including H

%%
monthDur=[1,31,29,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31];
monthStart=cumsum(monthDur);
%tvec=[-68.7792,32,88.2788,monthStart(5:end)];
tvec=[-68.7792,monthStart(2),93,monthStart(5:end)];
%
tvec=tvec(1:4);%[-40,365:368]; End=19
lt=length(tvec);
X=X(:,1:lt-1);
Xfull=Xfull(:,1:lt-1);

[lx1,lx2]=size(X);

xdata=85:tvec(end);
ydata=ydata(0+(1:length(xdata)));
%If data is just England:
%ydata=ydata*(sum(data.Npop)/56286961);%England, mid-2019 (ONS)

x0=thetaIn;
%
%R0, t0, t1, alpha, p1:
lb=[1.5,-200,85,0,0];%zeros(1,lx-2)];
ub=[3.5,0,120,1,1];%zeros(1,lx-2)];
%}
%
%alpha, p's:
%lb=[.2,zeros(1,16)];%zeros(1,lx-2)];
%ub=[.5,ones(1,16)];%zeros(1,lx-2)];
%}
%{
%alpha, feedback params L1, k's, H0:
lb=[.15,.15,0,-10*ones(1,numPCA),-10];%alpha, L, kvec, H0
ub=[.5,1,1,10*ones(1,numPCA),10];
%lb=[.15,0,0,-20,0,-50];%alpha, L, kvec, H0
%ub=[.5,1,20,0,20,50];
%}
%{
%alpha out:
lb=[0,-20*ones(1,numPCA),-100];%alpha, L, kvec, H0
ub=[1,20*ones(1,numPCA),100];
%lb=[.19,.8,-.4356,-1.1,4,-1];%alpha, L, kvec, H0
%ub=[.21,1,-.2356,-.9,6,3];
%}

%%
fun=@(params,xdata)sim2fit(params,data,xdata,X,intrinsic,Xfull,coeff,tvec,lx1,lx2);
%
tic
rng default;%for reproducibility
options=optimoptions(@lsqcurvefit,'MaxFunctionEvaluations',1e2,'MaxIterations',1e2);
problem=createOptimProblem('lsqcurvefit','x0',x0,'objective',fun,'xdata',xdata,'ydata',ydata','lb',lb,'ub',ub,'options',options);
ms=MultiStart;
[poptim,resnorm]=run(ms,problem,500);
toc
Ypred=1;%sim2fit(poptim,data,xdata,X,thetaIn,intrinsic,Xfull);
delta=1;
%
%%
%Plotting
ymod=fun(poptim,xdata);%sim2fit(poptim,data,xdata,X,thetaIn,intrinsic,Xfull,coeff,tvec,lx1,lx2);
f=figure('Units','centimeters','Position',[0 0 20 20]);
set(f,'defaulttextInterpreter','latex');
set(f,'defaultAxesTickLabelInterpreter','latex');
set(f,'defaultLegendInterpreter','latex');
set(f,'DefaultAxesFontSize',15);
hold on;
bar(xdata,ydata);
plot(xdata,ymod,'linewidth',2.5,'color','red');
for i=[1,32,61,92,122,153,183,214,245,275,306,336,367,398,426,457,487,518,548,579,610,640,671,701,731]
    plot(i*[1,1],[0,1.25*max(ydata)],'k-','linewidth',0.01);    
end
xlim([xdata(1),xdata(end)]);
ylim([0,1.25*max(ydata)]);
axis square;
box on;
grid on;
xlabel('Time');
ylabel('Hospital Occupancy');
title('Model Fit');
%}
end

function [f,rhohat]=sim2fit(params,data,xdata,Xfit,intrinsic,Xfull,coeff,tvec,lx1,lx2)
R0=params(1);
t0=params(2);
t1=params(3);
tvec(1)=t0;
tvec(3)=t1;
alpha=params(4)*ones(1,3);%0.3536;%params(1);

%%BH [pr,be,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=bePrepCovid19(data,lx2,R0,[ones(1,lx2-2)],zeros(1,5),coeff,ones(5,lx2),alpha);%5=numPCA+2
[pr,be,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=bePrepCovid19(data,R0,ones(1,lx2-2),zeros(1,5),ones(1,3)',zeros(5,lx2),alpha);

pr.xfull=Xfull;
Wfit=Xfit.^(1/pr.a);
if intrinsic==1
    %Fit to ocupancy:
    %[simu,~,~]=heRunCovid19(pr,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,[ones(1,length(tvec)-1)],tvec(1:numInt+1),0,data);
    %Fit to admissions:
    %%BH [simu,simu2,~,rhohat]=beRunCovid19(pr,be,vx,n,ntot,NN,NNbar,beta,[ones(1,length(tvec)-1)],tvec(1:numInt+1),0,data);
    be.BiFirstFit=params(5);
    [simu,simu2,~,rhohat]=beRunCovid19(pr,be,vx,n,NN,NNbar,beta,[ones(1,length(tvec)-1)],tvec,10,data);
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

%f(isinf(f))=-1e6;
%f(isnan(f))=-1e6;

end