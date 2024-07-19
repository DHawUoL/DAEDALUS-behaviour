function [poptim,Ypred,delta,resnorm]=beFitEpiFull(ydata,X,data,thetaIn,Xfull,coeff)%,lb,ub)
%% Parameters to fit:
%R0, t0 - while no mitigation
%alpha, explicit p's (as previous deltas)
%alpha, logistic parameters

intrinsic=1;%=1 for behaviour as feed in fully open economy for "original" DAEDALUS
%numPCA=size(coeff,1);%Number of x's in logistic regression, including H
projection=1;

%%
if projection==1
    monthDur=[1,31,29,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31,31];
    monthStart=cumsum(monthDur);
    tvec=monthStart;
    tvec=tvec(1:26);
    xdata=85:763;
else
    monthDur=[1,31,29,31,30,31,30,31,31,30,31,30,31,31,28,31,30];%,31,30,31];
    monthStart=cumsum(monthDur);
    tvec=monthStart;
    tvec=tvec(1:17);%20);
    xdata=85:tvec(17);%578;
end
lt=length(tvec);
X=X(:,1:lt-1);
Xfull=Xfull(:,1:lt-1);

[lx1,lx2]=size(X);

ydata=ydata(0+(1:length(xdata)));
%If data is just England:
%ydata=ydata*(sum(data.Npop)/56286961);%England, mid-2019 (ONS)

x0=thetaIn;

%Fitting link function:
%lb=[0,.15,0,-20,-200];
%ub=[.6,1,1,0,0];
%Fitting individual p's:
lb=zeros(1,lt-2);
ub=ones(1,lt-2);
ub([7,8,10])=.35;

%%
fun=@(params,xdata)sim2fit(params,data,xdata,X,intrinsic,Xfull,coeff,tvec,lx1,lx2);
plot(xdata,[fun(thetaIn,xdata);ydata'])
%
tic
rng default;%for reproducibility
options=optimoptions(@lsqcurvefit,'MaxFunctionEvaluations',1e3,'MaxIterations',1e3);
problem=createOptimProblem('lsqcurvefit','x0',x0,'objective',fun,'xdata',xdata,'ydata',ydata','lb',lb,'ub',ub,'options',options);
ms=MultiStart;
[poptim,resnorm]=run(ms,problem,10);
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
R0=2.75;%params(1);
tvec(1)=-70;%params(2);
alpha=params([1,1,1]);
%BH
%[pr,be,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=bePrepCovid19(data,R0,ones(1,lx2-2),params(3:end),coeff,zeros(5,lx2),alpha);
pr.xfull=Xfull;

[pr,be,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=bePrepCovid19(data,R0,ones(1,size(Xfull,2)-2),ones(1,3),1,zeros(5,lx2),alpha);%repmat([1,1,params(2:end)]
pr.xfull=[1,1,1,params(2:end)];%Use xfull as the value of p


Wfit=Xfit.^(1/pr.a);
if intrinsic==1
    %Fit to ocupancy:
    %[simu,~,~]=heRunCovid19(pr,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,[ones(1,length(tvec)-1)],tvec(1:numInt+1),0,data);
    %Fit to admissions:
    %%BH
    %Fitting link function:
    %[simu,simu2,~,rhohat]=beRunCovid19(pr,be,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,[ones(1,length(tvec)-1)],tvec(1:lx2+1),0,data);
    %Fitting individual p's:
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