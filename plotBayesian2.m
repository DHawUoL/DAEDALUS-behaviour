function f=plotBayesian2(ydata,X,data,xsto,Xfull,coeff,pointEst)
intrinsic=1;%=1 for behaviour as feed in fully open economy for "original" DAEDALUS

burn=2500;
int=100;

%%
tvec=[1,2,61,93,134,141,148,155,162,169,176,186,200,211,218,223,227,236,250,258,266,271,279,294,310,322,330,338,349,370,384,397,407,418,433,445,463,468,474,491,504,517,540,561,567,575];
xdata=85:tvec(end-7);%-2
lt=length(tvec);
X=X(:,1:lt-1);
Xfull=Xfull(:,1:lt-1);
[lx1,lx2]=size(X);
ydata=ydata(0+(1:length(xdata)));
%If data is just England:
ydata=ydata*(sum(data.Npop)/56286961);%England, mid-2019 (ONS)
fun=@(params,xdata)sim2fit(params,data,xdata,X,intrinsic,Xfull,coeff,tvec,lx1,lx2);

hold on
for i=burn+1:int:length(xsto)
    plot(xdata,fun(xsto(i,1:end-1),xdata))
end


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