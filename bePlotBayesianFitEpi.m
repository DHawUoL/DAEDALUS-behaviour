function plotmat=bePlotBayesianFitEpi(ydata,X,data,xsto,Xfull,coeff,pointEst)
%[xsto, outsto, history, accept_rate,covmat]=fitEpiBayesian(dataOcc,ones(1,19),dataUK1,[0.9058   -0.7800   -7.1241   15.4944    0.0000    0.0243],X2(4:22,:)',1);
addmodifier=1;
intrinsic=1;
nx=size(coeff,1);%Number of x's in logistic regression, including H
Xmu=mean(Xfull,2);
projection=0;
timeThresh=17;%17 end of April 21; 20 end of July 21
%%
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
%If data is just England:
%ydata=ydata*(sum(data.Npop)/56286961);%England, mid-2019 (ONS)
%%
burn=1e4;
numit=size(xsto,1);
int=floor((numit-burn)/50);%20;
sample=xsto(burn:int:end,1:end-1);
l1=size(sample,1);

fun=@(params)sim2fit(params,data,xdata,X,intrinsic,Xfull,coeff,tvec,lx1,lx2);
y1=fun(sample(1,:));
l2=length(y1);
plotmat=zeros(l1,l2);
plotmat(1,:)=y1;
for i=2:l1
    plotmat(i,:)=fun(sample(i,:));
end
[ypoint,rhohat]=fun(pointEst);%fun(mean(xsto,1));%

ymax=max(max(ydata),max(max(plotmat)));
factor=1e2*ceil(ymax/1e2);
if addmodifier==1
    xmean=xsto(end,3:end-1);%pointEst(3:end);%xsto(end:end,4:end-1);%mean(xsto(burn+1:end,4:end-1),1);%********
    Xmu=mean(Xfull,2);
    %
    params.L1=xmean(1);
    params.k1=xmean(2:end-1)';
    params.H01=xmean(end);
    %params.m1=xmean(end);
    %}
    %Periods, real H values:
    tvecPlus=[1,tvec(2:end)];
    value=zeros(1,tvec(end));%numPeriods);
    value2=ones(1,tvec(end));
    for i=4:lt-1%********
        ti=round(tvecPlus(i)):round(tvecPlus(i+1));
        %value(ti)=rhohat(i);
        value(ti)=beFeedback2(Xfull(i),params);%(coeff*(Xfull(:,i)-Xmu),params); (Xfull(:,i),params)
    end
    value(1:round(tvec(3))-1)=0;%Kicks in at tvec(3)
    %}
    %factor=3e3;
    ydata=ydata/factor;
    plotmat=plotmat/factor;
    ypoint=ypoint/factor;
    ymax=1;

    %value2(round(tvec(3)):tvec(5))=mean(xsto(burn+1:end,1));
    %value2(tvec(5)+1:tvec(8))=mean(xsto(burn+1:end,1));
    %value2(tvec(8)+1:end)=mean(xsto(burn+1:end,1));
    value2(round(tvec(4)):tvec(7))=xsto(end,1);%mean(xsto(burn+1:end,1));
    value2(tvec(7)+1:tvec(8))=xsto(end,2);%mean(xsto(burn+1:end,1));
    value2(tvec(8)+1:end)=xsto(end,2);%mean(xsto(burn+1:end,1));
end

%% Plot:
fs=10; lw=2;
cmap=lines(7);
col1=cmap(1,:);
col2=.5*[1,1,1];
%xvec=[1,32,60,tvec(4:end)];
%xvec=xvec(1:2:end);
if projection==1
    xnames={'Jan 2020','Feb 2020','Mar 2020','Apr 2020','May 2020','Jun 2020','Jul 2020','Aug 2020','Sep 2020','Oct 2020','Nov 2020','Dec 2020',...
    'Jan 2021','Feb 2021','March 2021','Apr 2021','May 2021','Jun 2021','July2021','Aug 2021','Sep 2021','Oct 2021','Nov 2021','Dec 2021','Jan 2022'};
    %xnamesAdd={'Aug 2021','Sep 2021','Oct 2021','Nov 2021','Dec 2021','Jan 2022'};
else
    xnames={'Jan 2020','Feb 2020','Mar 2020','Apr 2020','May 2020','Jun 2020','Jul 2020','Aug 2020','Sep 2020','Oct 2020','Nov 2020','Dec 2020',...
    'Jan 2021','Feb 2021','March 2021','Apr 2021'};
end
xnames=[xnames];%,xnamesAdd];
xvec=monthStart(1:length(xnames));
figure
hold on
%plot(plotmat');
bar(xdata,ydata,'FaceColor',col2,'EdgeColor',col2,'LineWidth',.01);
if projection==1
    plot(tvec(timeThresh)*[1,1],[0,factor],'k:','linewidth',2)
end
plot_distribution_prctile(xdata,plotmat,'color',col1,'prctile',(0:5:100));

%plot(xdata,ypoint,'k-','linewidth',2)

%for i=[1,32,61,92,122,153,183,214,245,275,306,336,367,398,426,457,487,518,548,579,610,640,671,701,731]
%    plot(i*[1,1],[0,1.25*max(max(ydata),max(max(plotmat)))],'k-','linewidth',0.01);
%end
if addmodifier==1
    h1=plot([-1,-1],[-1-1],'linewidth',2,'color',col2);
    h2=plot([-1,-1],[-1-1],'linewidth',2,'color',col1);
    h3=plot(xdata,value(xdata),'k--','linewidth',2);
    h4=plot(xdata,value2(xdata),'--','linewidth',2,'color',[.5,0,0]);
    hleglines=[h1(1),h2(1),h3(1),h4(1)];
    legend(hleglines,'Data','Model fit','p','\alpha','location','NW')
end
xlim([xdata(1),xdata(end)]);
ylim([0,ymax])%.25*max(ydata)]);
xticks(xvec)
xticklabels(xnames)
xtickangle(45)
box on;
grid on;
xlabel('Time');
ylabel('Hospital Admissions/50k');
%title('Model Fit');
end

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