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
tvec=[-68.7792,monthStart(2),93,110];%monthStart(4)];
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

%{
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
    %[zhat,~,res,~] = lsqcurvefit(fun2, z0, xdata, ydata', lb, ub, options);
    if res < best, best = res; poptim=[t0i,t1i,zhat]; resnorm=res; end %best_out = struct('t0',t0i,'t1',t1i,'z',zhat,'res',res); end
  end
end
Ypred=1;
delta=1;
%}

options_nl = optimoptions(@lsqnonlin, ...
    'MaxFunctionEvaluations',5e3,'MaxIterations',5e3, ...
    'StepTolerance',1e-8,'FunctionTolerance',1e-8, ...
    'FiniteDifferenceType','forward','FiniteDifferenceStepSize',1e-3, ...
    'Display','off');

best_res = inf; poptim = nan(1,5); best_z = []; best_t0 = NaN; best_t1 = NaN;

for t0i = -50:-50
  for t1i = 80:100
    % objective over z = [rH, alpha, p1] with t0,t1 held fixed
    obj = @(z) resid_peakaware_full([t0i, t1i, z(1), z(2)], ...%, z(3)
                                    data, xdata, X, intrinsic, Xfull, coeff, tvec, lx1, lx2, ydata');

    z0 = thetaIn(3:4);              % start for [rH, alpha, p1]
    lb = [0, 0.0];           % bounds for z
    ub = [0.4, 1.0];

    % solve (single-start is usually fine here)
    [zhat,~,resnorm] = lsqnonlin(obj, z0, lb, ub, options_nl);

    if resnorm < best_res
        best_res = resnorm;
        best_z   = zhat;
        best_t0  = t0i; best_t1 = t1i;
        poptim   = [best_t0, best_t1, best_z(:)'];
    end
  end
end

resnorm = best_res;

Ypred=1;
delta=1;

% Plot best fit
ymod = sim2fit(poptim, data, xdata, X, intrinsic, Xfull, coeff, tvec, lx1, lx2);
figure('Units','centimeters','Position',[0 0 20 20]); hold on;
bar(xdata, ydata);
plot(xdata, ymod, 'r-', 'LineWidth', 2.5);
xlabel('Time'); ylabel('Hospital Occupancy'); box on; grid on; axis square;
title(sprintf('Best t0=%d, t1=%d, rH=%.2f, alpha=%.2f', ... %p1=%.2f', ...
      best_t0, best_t1, best_z(1), best_z(2)));%, best_z(3)));


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
    R0 = 2.8;
    t0 = params(1); t1 = params(2);
    propIn = 1;%params(3);
    phi2=params(3);
    alpha  = params(4)*ones(1,3);
    p1     = 1;%params(5);

    % keep timeline sane
    tvec(1)=t0; tvec(3)=t1;
    if ~(t1 > tvec(2)+1 && t1 < tvec(4)-1)
        f = ones(size(xdata))*1e9; rhohat = NaN; return;
    end

    try
        [pr,be,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta] = ...
            bePrepCovid19(data,R0,ones(1,lx2-2),zeros(1,5),ones(1,3)',zeros(5,lx2),alpha,propIn);
        pr.leak=0; pr.xfull=Xfull; be.BiFirstFit=p1; pr.phi2=phi2;

        [simu,simu2,~,rhohat]=beRunCovid19(pr,be,vx,n,NN,NNbar,beta,[ones(1,length(tvec)-1)],tvec,0,data);
        t = simu(:,1)'; h = simu2';  % admissions
        f = interp1(t, h, xdata, 'linear');
        f(~isfinite(f)) = 1e9;
    catch
        f = ones(size(xdata))*1e9; rhohat = NaN;
    end
end

function r = resid_peakaware_full(params, data, xdata, Xfit, intrinsic, Xfull, coeff, tvec, lx1, lx2, y)
    [f,~] = sim2fit(params, data, xdata, Xfit, intrinsic, Xfull, coeff, tvec, lx1, lx2);

    % base residuals on log-scale (helps early/late)
    eps0 = 1e-6;
    r_ts = (f-y).^2;%log(max(f,eps0)) - log(max(y,eps0));

    % peak timing & height residuals
    [~, iyd] = max(y);  tpk_y = xdata(iyd);  hpk_y = y(iyd);
    [~, ifm] = max(f);  tpk_f = xdata(ifm);  hpk_f = f(ifm);

    wt = 0.0;   % timing weight
    wh = 0.0;   % height weight
    r_peak_t = sqrt(wt) * (tpk_f - tpk_y);
    r_peak_h = sqrt(wh) * (log(max(hpk_f,eps0)) - log(max(hpk_y,eps0)));

    % extra weight near observed peak (optional)
    sig = 7; wp = 0.0;
    wts = 1 + wp*exp(-0.5*((xdata - tpk_y)/sig).^2);
    r_ts = sqrt(wts(:)) .* r_ts(:);

    r = [r_ts; r_peak_t; r_peak_h];
    r(~isfinite(r)) = 0;
end