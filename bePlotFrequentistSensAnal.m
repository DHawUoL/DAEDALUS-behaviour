function plotmat=bePlotFrequentistSensAnal(ydata,X,data,xsto,Xfull,coeff,pointEst)
%[xsto, outsto, history, accept_rate,covmat]=fitEpiBayesian(dataOcc,ones(1,19),dataUK1,[0.9058   -0.7800   -7.1241   15.4944    0.0000    0.0243],X2(4:22,:)',1);
hlag=0;
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
    %Footfall:
    %tvec=[1,2,61,94,[134,141,148,155,162,169,176,186,200,211,218,223,227,236,250,258,266,271,279,294,310,322,330,338,349,370,384,397,407,418,433,445,463,468,474,491,504,517,540,561,567,575]+hlag];
    %Stringency:
    tvec=[1,2,61,94,[127,134,141,148,153,155,162,167,169,175,176,186,200,211,216,218,223,227,230,236,250,258,265,266,271,279,288,294,305,310,322,330,337,338,349,354,356,361,370,372,384,397,407,418,433,445,454,463,468,474,491,503,504,517,540,561,566,567,575]+hlag];

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

fun=@(params,Xfullin)sim2fit(params,data,xdata,X,intrinsic,Xfullin,coeff,tvec,lx1,lx2,0);

%% Generate sample
tab=[0.39381     0.0047523     0.38447      0.40316;
0.9954    0.00061135      0.9942       0.9966;
-0.11555      0.041358    -0.19685    -0.034255;
2.3062         0.107      2.0958       2.5165];
pointEst=tab(:,1)';
se=tab(:,2)';
l1        = 200;
sample=drawParamSamples(pointEst, se, l1);

y1=fun(sample(1,:),Xfull);
l2=length(y1);
plotmat=zeros(l1,l2);
plotmat(1,:)=y1;
for i=2:l1
    plotmat(i,:)=fun(sample(i,:),Xfull);
end
%%
reduction=(.02:.02:.1);
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
        plotmati(j,:)=fun(sample(j,:),Xfulli);
    end
    plotmatset{i}=plotmati;
    maxvals(i+1)=max(max(plotmati));
end
%%
plotmatboth=zeros(l1,l2);
Xfullboth=Xfull;
Xfullboth(:,from:end)=Xfullboth(:,from:end)*.98;
for i=1:l1
    plotmatboth(i,:)=fun(sample(i,:),Xfullboth);
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
plot(275*[1,1],[0,10],'-','linewidth',lw,'color',.2*[1,1,1])
%%
h=plot([-1,-1],[-1,-1],'linewidth',lw,'color',0*[1,1,1]);
hleg(1)=h(1);
h=plot([-1,-1],[-1,-1],'linewidth',lw,'color',.8*[1,1,1]);
hleg(end)=h(1);
for i=1:lr
    h=plot([-1,-1],[-1,-1],'linewidth',lw,'color',cmap(i,:));
    hleg(i+1)=h(1);
end
legend(hleg,'model fit','-2% trust','-4% trust','-6% trust','-8% trust','-10%','-2% trust/stringency','location','northeastoutside')
%legend(hleg,'model fit','-5% trust','-10% trust','-15% trust','-20% trust','-25%','-30% trust','-35% trust','-40%','-45% trust','-50% trust','-5% trust/stringency','location','northeastoutside')
%%
%axis([xdata(1),xdata(end),0,1])
xlim([xdata(1),xdata(end)]);
ylim([0,1.2])
xticks(xvec)
xticklabels(xnames)
xtickangle(45)
box on;
grid on;
xlabel('Time');
ylabel('Hospital Admissions/5k');
end

function [f,rhohat]=sim2fit(params,data,xdata,Xfit,intrinsic,Xfull,coeff,tvec,lx1,lx2,plotRun)
R0=2.2;%1.9;%2.75;%params(1);
tvec(1)=-145;%-70;%-195;%-206;%-195;%Seasonal;-206;%-70;%-85;%-70;%params(2);
alpha=params([1,1,1]);
%tvec(5:end)=tvec(5:end)+params(end);
a1=-.814;
b1=8.0161;
a2=-.8067;
b2=-8.0887;
ks=params(3);%ksrat0=0.1613
reducedParams=[1,a1*ks+b1,a2*ks+b2,params(4),0];
%BH
%Fitting link function:
[pr,be,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=bePrepCovid19(data,R0,ones(1,lx2-2),reducedParams,coeff,zeros(5,lx2),alpha);
pr.leak=params(2);
%[pr,be,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=bePrepCovid19(data,R0,ones(1,lx2-2),[params(2:end),0.8036*params(3)-0.3232],coeff,zeros(5,lx2),alpha);
%Interaction term:
pr.xfull=Xfull;
%Fitting individual p's:
%[pr,be,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=bePrepCovid19(data,R0,ones(1,size(Xfull,2)-2),ones(1,3),1,zeros(5,lx2),alpha);%repmat([1,1,params(2:end)]
%pr.xfull=[1,1,1,params(2:end)];%Use xfull as the value of p


Wfit=Xfit.^(1/pr.a);
if intrinsic==1
    %Fit to ocupancy:
    %[simu,~,~]=heRunCovid19(pr,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,[ones(1,length(tvec)-1)],tvec(1:numInt+1),0,data);
    %Fit to admissions:
    %%BH
    %Fitting link function:
    [simu,simu2,~,rhohat]=beRunCovid19(pr,be,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,[ones(1,length(tvec)-1)],tvec(1:lx2+1),plotRun,data);
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

function sample = drawParamSamples(pointEst, se, n)
    % pointEst : 1x4 vector of parameter estimates
    % se       : 1x4 vector of standard errors
    % n        : number of samples
    % sample   : n x 4 matrix of draws
    
    % replicate means and stds to size (n,4)
    muMat = repmat(pointEst(:).', n, 1);   % n x 4
    seMat = repmat(se(:).', n, 1);         % n x 4
    
    % independent normal draws
    sample = muMat + seMat .* randn(n, numel(pointEst));
end