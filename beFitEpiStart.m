function [poptim,Ypred,delta,resnorm]=beFitEpiStart(ydata,X,data,thetaIn,Xfull,coeff)%,lb,ub)
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
tvec=[-68.7792,monthStart(2),93,monthStart(5:end)];
%
tvec=tvec(1:4);%[-40,365:368]; End=19
lt=length(tvec);
X=X(:,1:lt-1);
Xfull=Xfull(:,1:lt-1);

[lx1,lx2]=size(X);

xdata=85:tvec(end);
ydata=ydata(0+(1:length(xdata)));
%If data is just England:
%ydata=ydata*(sum(data.Npop)/56286961);%England, mid-2019 (ONS)

x0=thetaIn;
%
%R0, t0, t1, alpha, p1:
lb=[-210,60,.5,0,0];%zeros(1,lx-2)];
ub=[00,120,2,1,1];%zeros(1,lx-2)];
%}

%%
fun=@(params,xdata)sim2fit(params,data,xdata,X,intrinsic,Xfull,coeff,tvec,lx1,lx2);
plot(xdata,[fun(thetaIn,xdata);ydata'])
%{
tic
rng default;%for reproducibility
options=optimoptions(@lsqcurvefit,'MaxFunctionEvaluations',1e2,'MaxIterations',1e2);
problem=createOptimProblem('lsqcurvefit','x0',x0,'objective',fun,'xdata',xdata,'ydata',ydata','lb',lb,'ub',ub,'options',options);
ms=MultiStart;
[poptim,resnorm]=run(ms,problem,10);
toc
Ypred=1;%sim2fit(poptim,data,xdata,X,thetaIn,intrinsic,Xfull);
delta=1;
%}

options = optimoptions(@lsqcurvefit, ...
    'MaxFunctionEvaluations', 5e3, ...
    'MaxIterations',  5e3, ...
    'StepTolerance',  1e-8, ...
    'FunctionTolerance', 1e-8, ...
    'FiniteDifferenceType','forward', ...   % less fragile
    'FiniteDifferenceStepSize', 1e-3);      % don’t poke too hard
best = inf; best_out = []; poptim=nan(1,5); resnorm=inf;
for t0i = -90:-74
  for t1i = 85:97
    fun2 = @(z,xdata) sim2fit([t0i,t1i,z(1),z(2),z(3)], data, xdata, X, intrinsic, Xfull, coeff, tvec, lx1, lx2);
    z0  = thetaIn(3:5);%[1.0, 0.5, 0.2];   % [rH, alpha, p1] starts
    lb  = [0.4, 0.0, 0.0];   ub = [2.5, 1.0, 1.0];
    [zhat,~,res,~] = lsqcurvefit(fun2, z0, xdata, ydata', lb, ub, options);
    if res < best, best = res; poptim=[t0i,t1i,zhat]; resnorm=res; end %best_out = struct('t0',t0i,'t1',t1i,'z',zhat,'res',res); end
  end
end
Ypred=1;
delta=1;

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
R0=2.8;%params(1);
t0=params(1);
t1=params(2);
tvec(1)=t0;
tvec(3)=t1;
propIn=params(3);
alpha=params(4)*ones(1,3);%0.3536;%params(1);

%try
    %%BH [pr,be,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=bePrepCovid19(data,lx2,R0,[ones(1,lx2-2)],zeros(1,5),coeff,ones(5,lx2),alpha);%5=numPCA+2
    [pr,be,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=bePrepCovid19(data,R0,ones(1,lx2-2),zeros(1,5),ones(1,3)',zeros(5,lx2),alpha,propIn);
    pr.leak=0;
    pr.xfull=Xfull;
    Wfit=Xfit.^(1/pr.a);
    if intrinsic==1
        %Fit to ocupancy:
        %[simu,~,~]=heRunCovid19(pr,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,[ones(1,length(tvec)-1)],tvec(1:numInt+1),0,data);
        %Fit to admissions:
        %%BH [simu,simu2,~,rhohat]=beRunCovid19(pr,be,vx,n,ntot,NN,NNbar,beta,[ones(1,length(tvec)-1)],tvec(1:numInt+1),0,data);
        be.BiFirstFit=params(5);
        [simu,simu2,~,rhohat]=beRunCovid19(pr,be,vx,n,NN,NNbar,beta,[ones(1,length(tvec)-1)],tvec,0,data);
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
%catch
    % ANY failure in ODE/eigs/etc. → big finite penalty
%    f = zeros(size(xdata)); f(:) = 1e9; rhohat = NaN;
%end

%f(isinf(f))=-1e6;
%f(isnan(f))=-1e6;

end