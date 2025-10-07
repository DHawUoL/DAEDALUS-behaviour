function [poptim,Ypred,delta,resnorm] = beFitEpiStartGA(ydata,X,data,thetaIn,Xfull,coeff)

intrinsic = 1;
numPCA    = size(coeff,1);

% ---- timeline/windows
monthDur   = [1,31,29,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31];
monthStart = cumsum(monthDur);
tvec       = [-68.7792, monthStart(2), 93, 110];
tvec       = tvec(1:4);

X     = X(:,1:numel(tvec)-1);
Xfull = Xfull(:,1:numel(tvec)-1);
[lx1,lx2] = size(X);

% ---- data slice
xdata = 85:tvec(end);
ydata = ydata(0+(1:numel(xdata)));        % admissions (column later)
ydata = ydata*(sum(data.Npop)/56286961);

% ---- helper (model forward)
fun = @(params,xd) sim2fit(params,data,xd,X,intrinsic,Xfull,coeff,tvec,lx1,lx2);

% ---- Objective for GA (SSE on admissions)
obj = @(z) obj_ga(z, data, xdata, X, Xfull, coeff, tvec, lx1, lx2, ydata(:));

% ---- Nonlinear constraint enforcing t1-window sanity
nonlc = @(z) timewin_con(z, tvec);

% ---- Decision vars: [t0, t1, phi2, alpha] with t0,t1 integer
nvars  = length(thetaIn);%3;
IntCon = [1 2];

% ---- Bounds (cover your manual point)
lb = [-100,  80, 0.00];%, 0.00];
ub = [-60, 110, 1.0];%, 1.00];

% ---- Manual seed (optional but helpful)
manual = thetaIn;%[-60, 92, 0.1637, 0.5282];

% ---- Initial population
popSize = 120;
initPop = repmat(lb, popSize, 1) + rand(popSize,nvars).*repmat(ub-lb, popSize, 1);
initPop(1,:) = manual;

% ---- GA options
opts = optimoptions('ga', ...
    'PopulationSize',       popSize, ...
    'MaxGenerations',       200, ...
    'EliteCount',           max(2, round(0.02*popSize)), ...
    'MutationFcn',          {@mutationadaptfeasible}, ...
    'CrossoverFraction',    0.8, ...
    'FunctionTolerance',    1e-6, ...
    'MaxStallGenerations',  60, ...
    'UseParallel',          false, ...   % set true if your sim is thread-safe
    'InitialPopulationMatrix', initPop, ...
    'Display',              'iter');

% =========================
% 1) Run GA
% =========================
rng default;   % for reproducibility if you want
[z_ga, sse_ga] = ga(obj, nvars, [],[],[],[], lb, ub, nonlc, IntCon, opts);

% SSE at manual (safety check / comparison)
sse_manual = obj(manual);

% =========================
% 2) Local polish of [phi2,alpha] with t0,t1 fixed to GA integers
% =========================
z0_local = z_ga;
obj_local = @(w) ...
    (fun([z0_local, w], xdata) - ydata(:));

lb_local = [lb(3:nvars)];%, lb(4)];
ub_local = [ub(3:nvars)];%), ub(4)];
opts_ls  = optimoptions('lsqnonlin','Display','off','MaxIterations',1e3, ...
                        'StepTolerance',1e-8,'FunctionTolerance',1e-8, ...
                        'FiniteDifferenceType','forward','FiniteDifferenceStepSize',1e-3);

[w_hat, ~, ~] = lsqnonlin(obj_local, z0_local(3:nvars), lb_local, ub_local, opts_ls);
z_hybrid   = [z0_local(1:2), w_hat(:)'];
sse_hybrid = obj(z_hybrid);

% =========================
% 3) Keep the best among manual / GA / hybrid
% =========================
cands   = [manual; z_ga; z_hybrid];
sses    = [sse_manual; sse_ga; sse_hybrid];
[resnorm, kbest] = min(sses);
z_best  = cands(kbest,:);   % [t0 t1 phi2 alpha]

% =========================
% 4) Prepare outputs & plot
% =========================
poptim = z_best;%[z_best(1), z_best(2), z_best(3), z_best(4)];  % no p1 here
Ypred  = 1; delta = 1;

% plot fitted admissions
ymod = fun([z_best, 1], xdata);
figure('Units','centimeters','Position',[0 0 20 20]); hold on;
bar(xdata, ydata, 'FaceAlpha',0.5, 'EdgeColor','none');
plot(xdata, ymod, 'r-', 'LineWidth', 2.5);
xlabel('Time'); ylabel('Hospital admissions'); box on; grid on; axis square;
%title(sprintf('Best: t0=%d, t1=%d, phi2=%.4f, alpha=%.4f | SSE=%.3g', ...
      %z_best(1), z_best(2), z_best(3), z_best(4), resnorm));
end

% ===== helpers =====

function [f,rhohat]=sim2fit(params,data,xdata,Xfit,intrinsic,Xfull,coeff,tvec,lx1,lx2)
    R0 = 2.8;
    t0 = params(1); t1 = params(2);
    propIn = 1;%params(3);
    phi2=0;%params(3);
    alpha  = params(3)*ones(1,3);
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

function sse = obj_ga(z, data, xdata, X, Xfull, coeff, tvec, lx1, lx2, y)
    % z = [t0, t1, phi2, alpha], p1 fixed to 1
    [f,~] = sim2fit([z, 1], data, xdata, X, 1, Xfull, coeff, tvec, lx1, lx2);
    if any(~isfinite(f))
        sse = 1e12;
    else
        e   = f(:) - y(:);
        sse = sum(e.^2);
    end
end

function [c,ceq] = timewin_con(z, tvec)
    % enforce same sanity you check inside sim2fit
    t1 = z(2);
    c   = [ (tvec(2)+1) - t1;     % <= 0
             t1 - (tvec(4)-1) ];  % <= 0
    ceq = [];
end
