function [f,g]=bePlotMultiRun(ydata,X,data,params,Xfull,coeff)%,matOut1)%,matOut2)
%[xsto, outsto, history, accept_rate,covmat]=fitEpiBayesian(dataOcc,ones(1,19),dataUK1,[0.9058   -0.7800   -7.1241   15.4944    0.0000    0.0243],X2(4:22,:)',1);
intrinsic=0;
%nx=4;%Number of x's in logistic regression, including H
addmodifier=0;

be1=1*ones(5,size(Xfull,2));
be1(:,1:2)=0;
%%
%
monthDur=[1,31,29,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31];
monthStart=cumsum(monthDur);
tvec=[-68.7792,32,88.2788,monthStart(5:end)];
%
%Test with fewer periods:
%
%be1(:,1:2)=1;%Mitigate 1st 2 periods
tvec=tvec([1:12]);%[-40,365:368];
X=X(1:13);
Xfull=Xfull(:,1:13);
%}
%{
weekDur=[1,repmat(7,1,51),9,repmat(7,1,51),8,repmat(7,1,51),8];
weekStart=cumsum(weekDur);
tvec=[-40.5541,weekStart(12:81)];%8 ~ last week of feb; 81 ~ Unum(195)
%}
X=X(:,1:length(tvec)-1);
xdata=70:365;%426;%245;%70:...

pvec=(0:.02:1); lp=length(pvec);
alphavec=(0:.02:1); la=length(alphavec);

%% If need to generate data for plots:#
%
fun=@(beProps,alpha)sim2fit(params,data,xdata,X,intrinsic,Xfull,coeff,beProps,tvec,alpha);
matOut1=zeros(lp,la);
matOut2=matOut1;
vecOut=zeros(length(xdata),lp,la);
for i=1:lp
    bei=be1;
    bei(:,1:end)=be1(:,1:end)*pvec(i);
    for j=1:la
        alphaj=alphavec(j);
        [f1,f2,f3]=fun(bei,alphaj);
        matOut1(i,j)=f2;
        matOut2(i,j)=f3;
        vecOut(:,i,j)=f1;
    end
end
%}
f=matOut1;
%g=matOut2;

%% Plotting:
fs=12;
figure
imagesc(pvec,alphavec,matOut1);
hold on
plot(0.6,.3,'ko','markersize',7,'linewidth',2)
set(gca,'YDir','normal') 
set(gca,'fontsize',fs)
xlabel('Effect \alpha of behavioural change')
ylabel('Proportion p changing behaviour')
title('Maximum occupancy')
colorbar
caxis([0,max(max(matOut1))])
box on;
%{
fs=12;
figure
imagesc(pvec,alphavec,matOut2);
set(gca,'YDir','normal') 
set(gca,'fontsize',fs)
xlabel('Proportion p changing behaviour')
ylabel('Effectiveness \alpha of behavioural change')
title('Total bed days')
colorbar
caxis([0,max(max(matOut2))])
box on;
%}
%figure
%plot(xdata,vecOut(:,:,:))

end

function [f,g,h]=sim2fit(params,data,xdata,Xfit,intrinsic,Xfull,coeff,beProps,tvec,alpha)
R0=params(1);%2.0733;%
t0=params(2);
%t1=params(3);
tvec(1)=t0;
%tvec(3)=t1;
params=params(4:end);%"end-1" already accounted for in Lhood function
    [numSectors,numInt]=size(Xfit);

    [pr,be,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=bePrepCovid19(data,numInt,R0,[ones(1,length(tvec)-3)],params,coeff,beProps,alpha);%-5,,1.25,1.5
    %be.alphaB=alpha;
    pr.xfull=Xfull;
    Wfit=Xfit.^(1/pr.a);
    if intrinsic==1
        %Fit to ocupancy:
        %[simu,~,~]=heRunCovid19(pr,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,[ones(1,length(tvec)-1)],tvec(1:numInt+1),0,data);
        %Fit to admissions:
        [simu,simu2,~,rhohat]=beRunCovid19(pr,be,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,[ones(1,length(tvec)-1)],tvec(1:numInt+1),0,data);
    else
        [simu,~,~,rhohat]=beRunCovid19(pr,be,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,Wfit,tvec(1:numInt+1),0,data);
    end

    t=simu(:,1)';
    %Fit to ocupancy:
    h=simu(:,4)';
    %Fit to admissions:
    %h=simu2';

    f=interp1(t,h,xdata); 

    f(isinf(f))=-1e6;
    f(isnan(f))=-1e6;
    
    g=max(f);
    h=sum(f);
end