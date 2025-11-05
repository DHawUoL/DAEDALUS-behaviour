function [f,rhohat]=sim2fit_global(params,data,xdata,intrinsic,Xfull,tvec,lx2,plotRun,ymean,arg)
Xfit=ones(size(Xfull,2));
coeff=ones(1,3)';

R0=2.8;
tvec(1)=-81; tvec(2)=1;
alpha=params([1,1,1]);
propIn=1;
if length(params)==5
    %arg=b0
    k1=params(2);
    k2=params(3);
    k3=params(4);
    delta=params(5);%delta>eps_safe
    v0pca=-delta;%-dot([k1,k2],arg)-delta;

    reducedParams=[1,k1,k2,k3,v0pca];
    %softplus = @(z) log1p(exp(-abs(z))) + max(z,0);
    %reducedParams = [1, params(2:4), -softplus(params(5))];
elseif length(params)==4
    %{
    p0=.01;
    k1=params(2);
    k2=params(3);
    k3=params(4);
    delta=log((1-p0)/p0);
    v0pca=-dot([k1,k2],arg)-delta;
    reducedParams=[1,k1,k2,k3,v0pca];
    %}
    %
    pc=1;%PC1 OR PC2
    kstar=params(2);
    k3=params(3);
    delta=params(4);%delta>eps_safe
    v0pca=-kstar*arg(pc)-delta;
    %}
    if pc==1
        reducedParams=[1,kstar,0,k3,v0pca];   
    else
        reducedParams=[1,0,kstar,k3,v0pca];                      
    end
    %}
elseif length(params)==3
    %reducedParams=[1,0,0,params(2:3)];
    %
    pc=2;%PC1 OR PC2
    kstar=params(2);
    k3=params(3);
    delta=16.8943;
    v0pca=-kstar*arg(pc)-delta;
    if pc==1
        reducedParams=[1,kstar,0,k3,v0pca];   
    else
        reducedParams=[1,0,kstar,k3,v0pca];                      
    end
    %}
elseif length(params)==2
    alpha=0.5832*ones(1,3);
    pc=2;%PC1 OR PC2
    kstar=params(1);
    k3=params(2);
    delta=16.8943;
    v0pca=-kstar*arg(pc)-delta;
    if pc==1
        reducedParams=[1,kstar,0,k3,v0pca];   
    else
        reducedParams=[1,0,kstar,k3,v0pca];                      
    end
end
%BH
%Fitting link function:
[pr,be,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=bePrepCovid19(data,R0,ones(1,lx2-2),reducedParams,coeff,zeros(5,lx2),alpha,propIn);
pr.leak=0; pr.xfull=Xfull; be.BiFirstFit=1; pr.phi2=0;%.186;
pr.xfull=Xfull;
pr.ymean=ymean;
pr.t_alpha=245;

Wfit=Xfit.^(1/pr.a);
if intrinsic==1
    %Fit to ocupancy:
    %[simu,~,~]=heRunCovid19(pr,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,[ones(1,length(tvec)-1)],tvec(1:numInt+1),0,data);
    %Fit to admissions:
    %%BH
    %Fitting link function:
    [simu,simu2,~,rhohat]=beRunCovid19(pr,be,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,[ones(1,length(tvec)-1)],tvec(1:lx2+1),plotRun,data);
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
end