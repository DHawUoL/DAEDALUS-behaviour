function f=beOnePlot(ydata,X,data,thetaIn,Xfull,coeff)%,lb,ub)

%beOnePlot(dataAdmissions',ones(1,19),dataUK1,xstoOneAlpha2(end,1:end-1),Xin20ld(:,[1,3:end]),coeff(:,1:3)');

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
tvec=[-177.4386,monthStart(2),92.5467,monthStart(5:end)];
%
tvec=tvec(1:19);%[-40,365:368]; End=19
lt=length(tvec);
X=X(:,1:lt-1);
Xfull=Xfull(:,1:lt-1);

[lx1,lx2]=size(X);

xdata=85:578;%70:tvec(19);%579;%245;%70:...
ydata=ydata(0+(1:length(xdata)));
%If data is just England:
%ydata=ydata*(sum(data.Npop)/56286961);%England, mid-2019 (ONS)

x0=thetaIn;

%%
fun=@(params,xdata)sim2fit(params,data,xdata,X,intrinsic,Xfull,coeff,tvec,lx1,lx2);
ymod=fun(x0,xdata);
%
%%
%Plotting
%
f=figure('Units','centimeters','Position',[0 0 20 20]);
set(f,'defaulttextInterpreter','latex');
set(f,'defaultAxesTickLabelInterpreter','latex');
set(f,'defaultLegendInterpreter','latex');
set(f,'DefaultAxesFontSize',15);
%}
cmap=lines(7);
hold on;
bar(xdata,ydata,'facecolor',cmap(1,:),'edgecolor',cmap(1,:));
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
R0=1.9928;%2.8025;%params(1);%2.0733;%
%t0=params(2);
%t1=params(3);
%tvec(1)=t0;
%tvec(3)=t1;
alpha=params(1)*ones(1,3);%[params(1:2),params(2)];%0.3536;%params(1); %params(1)*ones(1,3);%

%beProps=[zeros(5,2),repmat(params(2:end),5,1)];
%[pr,be,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=bePrepCovid19(data,lx2,R0,[ones(1,lx2-4),1.25,1.5],0,coeff,beProps,alpha);

[pr,be,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=bePrepCovid19(data,lx2,R0,[ones(1,lx2-2)],params(2:end),coeff,ones(5,lx2),alpha);
pr.xfull=Xfull;
Wfit=Xfit.^(1/pr.a);
if intrinsic==1
    %Fit to ocupancy:
    %[simu,~,~]=heRunCovid19(pr,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,[ones(1,length(tvec)-1)],tvec(1:numInt+1),0,data);
    %Fit to admissions:
    [simu,simu2,~,rhohat]=beRunCovid19(pr,be,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,[ones(1,length(tvec)-1)],tvec(1:lx2+1),0,data);
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

%plot(xdata,f)

%f(isinf(f))=-1e6;
%f(isnan(f))=-1e6;

end