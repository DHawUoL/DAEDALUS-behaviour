function zout=sensAnalGrid(ydata,X,data,xsto,Xfull,coeff,pointEst)
%[xsto, outsto, history, accept_rate,covmat]=fitEpiBayesian(dataOcc,ones(1,19),dataUK1,[0.9058   -0.7800   -7.1241   15.4944    0.0000    0.0243],X2(4:22,:)',1);
addmodifier=0;
intrinsic=1;
nx=size(coeff,1);%Number of x's in logistic regression, including H
Xmu=mean(Xfull,2);
projection=0;
timeThresh=17;%17 end of April 21; 20 end of July 21
%%
if projection==1
    tvec=[1,2,61,92,134,141,148,155,162,169,176,186,200,211,218,223,227,236,250,258,266,271,279,294,310,322,330,338,349,370,384,397,407,418,433,445,463,468,474,491,504,517,540,561,567,575];
    xdata=85:tvec(end);%-2
else
    tvec=[1,2,61,92,134,141,148,155,162,169,176,186,200,211,218,223,227,236,250,258,266,271,279,294,310,322,330,338,349,370,384,397,407,418,433,445,463,468,474,491,504,517,540,561,567,575];
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
%ydata=ydata*(sum(data.Npop)/56286961);%England, mid-2019 (ONS)
%%
fun=@(params,Xfullin)sim2fit(params,data,xdata,X,intrinsic,Xfullin,coeff,tvec,lx1,lx2);
%%
reduction=(0:.01:.3);
reduction=1-reduction;
lr=length(reduction);
reduction2=(0:.01:.5);%Check which index is which axis - 1 is y axis (up-down)/i axis in matrix
reduction2=1-reduction2;
lr2=length(reduction2);
%reduction2=reduction1;

from=22;
zout=zeros(lr,lr2);
for i=1:lr
    Xfulli=Xfull;
    Xfulli(1,from:end)=Xfulli(1,from:end)*reduction(i);
    parfor j=1:lr2
        Xfullj=Xfulli;
        Xfullj(2,from:end)=Xfullj(2,from:end)*reduction2(j);
        [ypointij,~]=fun(pointEst,Xfullj);
        zout(i,j)=sum(ypointij);
    end
end
zout=zout/zout(1,1)*100;
%% Plot:
fs=10; lw=2; ms=5;
figure
colormap parula
set(gca,'fontsize',fs)
imagesc(100*(1-reduction),100*(1-reduction),zout)
xlabel('% Reduction in mobility index')
ylabel('% Reduction in trust index')
%axis([xdata(1),xdata(end),0,1])
%xlim([xdata(1),xdata(end)]);
%ylim([0,1])
%xticks(xvec)
%xticklabels(xnames)
%xtickangle(45)
set(gca,'YDir','normal')
clim([0,100])%max(max(zout))])
colorbar
box on;
grid on;
end

function [f,rhohat]=sim2fit(params,data,xdata,Xfit,intrinsic,Xfull,coeff,tvec,lx1,lx2)
R0=2.75;%1.9;%2.75;%params(1);
tvec(1)=-84;%-70;%-206;%-195;%Seasonal;-206;%-70;%-85;%-70;%params(2);
alpha=params([1,1,1]);
%tvec(5:end)=tvec(5:end)+params(end);
%BH
%Fitting link function:
[pr,be,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=bePrepCovid19(data,R0,ones(1,lx2-2),params(2:end),coeff,zeros(5,lx2),alpha);
%Interaction term:
pr.xfull=Xfull;
%Fitting individual p's:
%[pr,be,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=bePrepCovid19(data,R0,ones(1,size(Xfull,2)-2),ones(1,3),1,zeros(5,lx2),alpha);%repmat([1,1,params(2:end)]
%pr.xfull=[1,1,1,params(2:end)];%Use xfull as the value of p


Wfit=Xfit.^(1/pr.a);
if intrinsic==1
    %Fitting link function:
    [simu,simu2,~,rhohat]=beRunCovid19(pr,be,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,[ones(1,length(tvec)-1)],tvec(1:lx2+1),0,data);
else
    [simu,simu2,~,rhohat]=beRunCovid19(pr,be,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,Wfit,tvec(1:lx2+1),0,data);
end
t=simu(:,1)';
h=simu2';
f=interp1(t,h,xdata); 
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