function plotmat=bePlotSingleRun(ydata,X,data,params,Xfull,coeff,beProps,alphaIn)
%X used in contact rate file (including behaviour)
%Xfull used for PCA
%
%[xsto, outsto, history, accept_rate,covmat]=fitEpiBayesian(dataOcc,ones(1,19),dataUK1,[0.9058   -0.7800   -7.1241   15.4944    0.0000    0.0243],X2(4:22,:)',1);
intrinsic=1;
%nx=4;%Number of x's in logistic regression, including H
addmodifier=0;
%%
%
monthDur=[1,31,29,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31];
monthStart=cumsum(monthDur);
tvec=[-68.7792,32,88.2788,monthStart(5:end)];%-68.7792
%
tvec=tvec(1:19);%[-40,365:368];
lt=length(tvec);
X=X(1:lt-1);
Xfull=Xfull(:,1:lt-1);
beProps=beProps(:,1:lt-1);
%}
%{
weekDur=[1,repmat(7,1,51),9,repmat(7,1,51),8,repmat(7,1,51),8];
weekStart=cumsum(weekDur);
tvec=[-40.5541,weekStart(12:81)];%8 ~ last week of feb; 81 ~ Unum(195)
%}
X=X(:,1:lt-1);
[lx1,lx2]=size(X);

xdata=70:tvec(19);%365;%70:426;%245;%70:...
ydata=ydata(0+(1:length(xdata)));
ydata=ydata*(sum(data.Npop)/56286961);%England, mid-2019 (ONS)plotS

fun=@(params)sim2fit(params,data,xdata,X,intrinsic,Xfull,coeff,beProps,tvec,alphaIn,lx1,lx2);

[plotmat,rhohat]=fun(params);
%{
params.L1=pointEst(4);
params.k1=pointEst(5:7)';
params.H01=pointEst(8);
params.m1=pointEst(9);
%}
%Periods, real H values:
%
if addmodifier==1
    tvecPlus=[1,tvec(2:end)];
    value=ones(1,tvec(end));%numPeriods);
    for i=3:lt-1
        ti=floor(tvecPlus(i)):floor(tvecPlus(i+1));
        value(ti)=params(i-2);
    end
    value(1:floor(tvec(3))-1)=1;%Kicks in at tvec(3)
    %}
    factor=4e4;
    ydata=ydata/factor;
    plotmat=plotmat/factor;
end
%%

%

fs=10; lw=2;
cmap=lines(7);
col1=cmap(1,:);
col2=.5*[1,1,1];
xvec=monthStart;%[1,tvec(2:end)];
xnames={'Jan 2020','Feb 2020','Mar 2020','Apr 2020','May 2020','Jun 2020','Jul 2020','Aug 2020','Sep 2020','Oct 2020','Nov 2020','Dec 2020','Jan 2021','Feb 2021','March 2021','Apr 2021','May 2021','Jun 2021','July2021'};
figure
hold on
bar(xdata,ydata,'FaceColor',col2,'EdgeColor',col2,'LineWidth',.01);
plot(xdata,plotmat,'-','color',col1,'linewidth',lw);
%plot_distribution_prctile(xdata,plotmat,'color',col1,'prctile',(0:5:100));
xlim([xdata(1),xdata(end)]);
ylim([0,1.25*max(ydata)]);
%axis tight
if addmodifier==1
    plot(xdata,value(xdata),'k--','linewidth',2)
end
xticks(xvec)
xticklabels(xnames)
xtickangle(45)
box on;
grid on;
xlabel('Time');
ylabel('Hospital Occupancy/40k');
title('Model Fit');

%}

end

function [f,rhohat]=sim2fit(params,data,xdata,Xfit,intrinsic,Xfull,coeff,beProps,tvec,alpha,lx1,lx2)
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