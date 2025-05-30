function plotmat=bePlotBayesianFitEpiSensAnal(ydata,X,data,xsto,Xfull,coeff,pointEst)
%[xsto, outsto, history, accept_rate,covmat]=fitEpiBayesian(dataOcc,ones(1,19),dataUK1,[0.9058   -0.7800   -7.1241   15.4944    0.0000    0.0243],X2(4:22,:)',1);
addmodifier=0;
intrinsic=1;
nx=size(coeff,1);%Number of x's in logistic regression, including H
Xmu=mean(Xfull,2);
projection=0;
timeThresh=17;%17 end of April 21; 20 end of July 21
%%
if projection==1
    tvec=[1,2,61,91,134,141,148,155,162,169,176,186,200,211,218,223,227,236,250,258,266,271,279,294,310,322,330,338,349,370,384,397,407,418,433,445,463,468,474,491,504,517,540,561,567,575];
    xdata=85:tvec(end);%-2
else
    tvec=[1,2,61,91,134,141,148,155,162,169,176,186,200,211,218,223,227,236,250,258,266,271,279,294,310,322,330,338,349,370,384,397,407,418,433,445,463,468,474,491,504,517,540,561,567,575];
    xdata=85:tvec(end-7);%-2
end
lt=length(tvec);
X=X(:,1:lt-1);
Xfull=Xfull(:,1:lt-1);

[lx1,lx2]=size(X);

ydata=ydata(0+(1:length(xdata)));
%If data is just England:
%ydata=ydata*(sum(data.Npop)/56286961);%England, mid-2019 (ONS)
%If data is just England:
ydata=ydata*(sum(data.Npop)/56286961);%England, mid-2019 (ONS)
%%
burn=2.5e3;
numit=size(xsto,1);
int=floor((numit-burn)/5);
sample=xsto(burn:int:end,:);
l1=size(sample,1);

fun=@(params,Xfullin)sim2fit(params,data,xdata,X,intrinsic,Xfullin,coeff,tvec,lx1,lx2);
y1=fun(sample(1,1:end-1),Xfull);
l2=length(y1);
plotmat=zeros(l1,l2);
plotmat(1,:)=y1;
for i=2:l1
    plotmat(i,:)=fun(sample(i,1:end-1),Xfull);
end
%%
reduction=(.05:.05:.5);
reduction=1-reduction;
from=22;
lr=length(reduction);
plotmatset=cell(1,lr);
maxvals=zeros(1,lr+2);
maxvals(1)=max(max(plotmat));
for i=1:lr
    Xfulli=Xfull;
    Xfulli(1,from:end)=Xfulli(1,from:end)*reduction(i);
    plotmati=zeros(l1,l2);
    for j=1:l1
        plotmati(j,:)=fun(sample(j,1:end-1),Xfulli);
    end
    plotmatset{i}=plotmati;
    maxvals(i+1)=max(max(plotmati));
end
%%
plotmatboth=zeros(l1,l2);
Xfullboth=Xfull;
Xfullboth(:,from:end)=Xfullboth(:,from:end)*.95;
for i=1:l1
    plotmatboth(i,:)=fun(sample(i,1:end-1),Xfullboth);
end
maxvals(end)=max(max(plotmati));

ymax=max(max(ydata),max(maxvals));
factor=5e3;%1e2*ceil(ymax/1e2);
ymax=ymax/factor;

%% Plot:
fs=10; lw=2;
cmap=lines(7);
col1=cmap(1,:);
col2=.5*[1,1,1];
monthDur=[1,31,29,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31,31];
monthStart=cumsum(monthDur);
if projection==1
    xnames={'Jan 2020','Feb 2020','Mar 2020','Apr 2020','May 2020','Jun 2020','Jul 2020','Aug 2020','Sep 2020','Oct 2020','Nov 2020','Dec 2020',...
    'Jan 2021','Feb 2021','March 2021','Apr 2021','May 2021','Jun 2021','July2021','Aug 2021','Sep 2021','Oct 2021','Nov 2021','Dec 2021','Jan 2022'};
else
    xnames={'Jan 2020','Feb 2020','Mar 2020','Apr 2020','May 2020','Jun 2020','Jul 2020','Aug 2020','Sep 2020','Oct 2020','Nov 2020','Dec 2020',...
    'Jan 2021','Feb 2021','March 2021','Apr 2021'};
end
xvec=monthStart(1:length(xnames));

figure
cmap=parula(lr);
hleg=zeros(1,lr+2);
hold on
bar(xdata,ydata/factor,'FaceColor',col2,'EdgeColor',col2,'LineWidth',.01);
if projection==1
    %plot(tvec(timeThresh)*[1,1],[0,factor],'k:','linewidth',2)
    plot((tvec(end)+1)*[1,1],[0,factor],'k:','linewidth',2)
end
plot_distribution_prctile(xdata,plotmat/factor,'color',0*[1,1,1],'prctile',(0:25:100));
%hleg(1)=h(1);
plot_distribution_prctile(xdata,plotmatboth/factor,'color',.8*[1,1,1],'prctile',(0:25:100));
%hleg(end)=h(1);
for i=1:lr
    plot_distribution_prctile(xdata,plotmatset{i}/factor,'color',cmap(i,:),'prctile',(0:25:100));
    %hleg(i+1)=h(1);
end
plot(275*[1,1],[0,1],'-','linewidth',lw,'color',.2*[1,1,1])
%%
h=plot([-1,-1],[-1,-1],'linewidth',lw,'color',0*[1,1,1]);
hleg(1)=h(1);
h=plot([-1,-1],[-1,-1],'linewidth',lw,'color',.8*[1,1,1]);
hleg(end)=h(1);
for i=1:lr
    h=plot([-1,-1],[-1,-1],'linewidth',lw,'color',cmap(i,:));
    hleg(i+1)=h(1);
end
%legend(hleg,'model fit','-2% trust','-4% trust','-6% trust','-8% trust','-10%','-2% trust/mobility','location','north')
legend(hleg,'model fit','-5% trust','-10% trust','-15% trust','-20% trust','-25%','-30% trust','-35% trust','-40%','-45% trust','-50% trust','-5% trust/mobility','location','northeastoutside')
%%
%axis([xdata(1),xdata(end),0,1])
xlim([xdata(1),xdata(end)]);
ylim([0,1])
xticks(xvec)
xticklabels(xnames)
xtickangle(45)
box on;
grid on;
xlabel('Time');
ylabel('Hospital Admissions/5k');
end

function [f,rhohat]=sim2fit(params,data,xdata,Xfit,intrinsic,Xfull,coeff,tvec,lx1,lx2)

a=.6121;%Also hard-coded in sim2fit
b=.5987;

R0=2.75;%1.9;%2.75;%params(1);
tvec(1)=-84;%-70;%-195;%-206;%-195;%Seasonal;-206;%-70;%-85;%-70;%params(2);
alpha=params([1,1,1]);
%tvec(5:end)=tvec(5:end)+params(end);
%BH
%Fitting link function:
[pr,be,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=bePrepCovid19(data,R0,ones(1,lx2-2),[a*params(1)+b,params(2:end),0.8036*params(2)-0.3232],coeff,zeros(5,lx2),alpha);
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

function f=rho2ofh(H,pr)
%%**
%f=pr.m1+(1-pr.m1)*pr.L1/(1+exp(dot(pr.k1,H)-pr.H01'));
f=pr.L1;

%f=pr.m1+(1-pr.m1)*pr.L1/(1+exp(dot(pr.k1,H)-pr.H01));
%)
%f=min(1,1-pr.L1*pr.k1*(H-pr.H01)/(1+pr.L1*(H-pr.H01)));
%f=pr.L1*exp(-pr.k1*H);
%f=pr.L1+pr.k1*H;
end

function f=beFeedback(H,pr)
f=pr.L1/(1+exp(dot(pr.k1',H)-pr.H01'));
end

function f=beFeedback2(H,pr)
f=pr.L1/(1+exp(dot(pr.k1,H)-pr.H01'));
end