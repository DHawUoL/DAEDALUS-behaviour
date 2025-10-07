function [poptim,Ypred,delta,resnorm] = beFitEpiFull_GA(ydata,X,data,thetaIn,Xfull,coeff)
hlag = -7;
plotRun = 0;
intrinsic = 1;
projection = 0;

% ===== timeline / windows (your original logic) =====
if projection==1
    tvec = [1    32    61    93   107   169   176   223   227   230   250   258   265   271   279   294   310   322   330   338   349   370   384    407   418   433   445   463   474   491   504   517   540   561   567   575   594   605   617   631   642   652   661   679   693   702    716   742   784];
    xdata = 85:763;
else
    tvec = [1,2,61,91,127,134,141,148,153,155,162,167,169,175,176,186,200,211,216,218,223,227,230,236,[250,258,265,266,271,279,288,294,305,310,322,330,337,338,349,354,356,361,370,372,384,397,407,418,433,445,454,463,468,474,491,503,504,517,540,561,566,567,575]+hlag];
    xdata = 85:tvec(end-15);
end

lt = numel(tvec);
X     = X(:,1:lt-1);
Xfull = Xfull(:,1:lt-1);
[lx1,lx2] = size(X);

% admissions slice and England scale (as in your code)
ydata = ydata(0+(1:numel(xdata)));
ydata = ydata*(sum(data.Npop)/56286961);
ymean=0;%mean(ydata);

% weights (same as your weighted SSE)
weights = 1./(1 + ydata(:));

% forward model handle EXACTLY as you call it
fun = @(params,xd) sim2fit(params, data, xd, X, intrinsic, Xfull, coeff, tvec, lx1, lx2, plotRun, ymean);

% ===== bounds (use exactly your choices) =====
% params meaning in your sim2fit:
%   alpha = params(1) replicated to [1,1,1]
%   reducedParams = [1, params(2:end)] -> passed into bePrepCovid19 as forFeedback vector
% Keep these bounds consistent with your modeling choices.
lb = [0,  -20, -20,  -20];
ub = [1,   20,  20,  20];
nvars = numel(lb);

% ===== objective: weighted SSE on admissions =====
obj = @(z) obj_full_ga(z, fun, xdata, ydata(:), weights);

% ===== seed with your current thetaIn =====
popSize = 120;
initPop = repmat(lb, popSize, 1) + rand(popSize, nvars).*repmat(ub-lb, popSize, 1);
initPop(1,:) = thetaIn(:)';  % manual seed first

% ===== GA options =====
opts = optimoptions('ga', ...
    'PopulationSize',       popSize, ...
    'MaxGenerations',       250, ...
    'EliteCount',           max(2, round(0.02*popSize)), ...
    'MutationFcn',          {@mutationadaptfeasible}, ...
    'CrossoverFraction',    0.8, ...
    'FunctionTolerance',    1e-6, ...
    'MaxStallGenerations',  80, ...
    'InitialPopulationMatrix', initPop, ...
    'UseParallel',          true, ...   % set true only if your sim is thread-safe
    'Display',              'iter');

% ===== 1) GA =====
rng default;
%[z_ga, sse_ga] = ga(obj, nvars, [],[],[],[], lb, ub, [], [], opts);
IntCon=[];
[z_ga, sse_ga] = ga(obj, nvars, [],[],[],[], lb, ub, [], IntCon, opts);

% also evaluate manual point (thetaIn) for apples-to-apples comparison
sse_manual = obj(thetaIn(:)');

% ===== 2) local polish of all continuous params with lsqnonlin =====
% (you can restrict to some subset if you prefer)
obj_vec = @(z) sqrt(weights).*( fun(z, xdata) - ydata(:) );  % vector residuals
opts_ls  = optimoptions('lsqnonlin','Display','off','MaxIterations',1e3, ...
                        'StepTolerance',1e-8,'FunctionTolerance',1e-8, ...
                        'FiniteDifferenceType','forward','FiniteDifferenceStepSize',1e-3);

[z_hyb, ~, ~] = lsqnonlin(obj_vec, z_ga, lb, ub, opts_ls);
sse_hyb = obj(z_hyb);

% ===== 3) keep the best among manual / GA / hybrid =====
cands = [thetaIn(:)'; z_ga; z_hyb];
sses  = [sse_manual;   sse_ga; sse_hyb];
[resnorm, kbest] = min(sses);
poptim = cands(kbest,:);

% ===== outputs & plot =====
Ypred = 1; delta = 1;

ymod = fun(poptim, xdata);
figure('Units','centimeters','Position',[0 0 20 20]); hold on;
bar(xdata, ydata, 'FaceAlpha',0.5, 'EdgeColor','none');
plot(xdata, ymod, 'r-', 'LineWidth', 2.2);
xlabel('Time'); ylabel('Hospital admissions'); box on; grid on; axis square;
title(sprintf('Full fit: alpha=%.4f | SSE=%.3g', poptim(1), resnorm));
legend({'Data','Model'}, 'Location','best');

end

function [f,rhohat]=sim2fit(params,data,xdata,Xfit,intrinsic,Xfull,coeff,tvec,lx1,lx2,plotRun,ymean)
R0=2.8;%2.2;
tvec(1)=-80;%-59;
alpha=params([1,1,1]);
propIn=1;
%reducedParams=[1,params(2:end)];
reducedParams=[1,params(2),0,params(3:4)];
%BH
%Fitting link function:
[pr,be,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=bePrepCovid19(data,R0,ones(1,lx2-2),reducedParams,coeff,zeros(5,lx2),alpha,propIn);
pr.leak=0; pr.xfull=Xfull; be.BiFirstFit=1; pr.phi2=0;%.186;
%[pr,be,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=bePrepCovid19(data,R0,ones(1,lx2-2),[params(2:end),0.8036*params(3)-0.3232],coeff,zeros(5,lx2),alpha);
%Interaction term:
%{
delta = params(4);                 % behaviour lead/lag in days
idx_start = 7;                             % first window to shift
Xfull_shift = shift_driver_by_days(Xfull, tvec, delta, idx_start, 'previous');
pr.xfull = Xfull_shift;
%}
pr.xfull=Xfull;
pr.ymean=ymean;

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

% ===== helper: scalar SSE for GA =====
function sse = obj_full_ga(z, fun, xdata, y, weights)
    f = fun(z, xdata);
    if any(~isfinite(f))
        sse = 1e12;
        return;
    end
    e   = f(:) - y(:);
    sse = sum( weights .* (e.^2) );
end

function Xfull_shift = shift_driver_by_days(Xfull, tvec, delta_days, idx_start, step_rule)
% Shift a piecewise-constant driver Xfull by delta_days (lead if >0, lag if <0)
% without changing tvec. Works with irregular tvec.
%
% Xfull:      [m x (lt-1)] driver, column j active on [tvec(j), tvec(j+1))
% tvec:       [1 x lt]     breakpoints (ascending)
% delta_days: scalar, shift in days (+ leads behaviour earlier)
% idx_start:  first window index to shift (e.g. 7), earlier windows untouched
% step_rule:  'previous' (left-continuous) or 'next' (right-continuous)
%
% Returns:
%   Xfull_shift with same size/shape as Xfull, re-averaged over the *original* windows.

    if nargin < 5 || isempty(step_rule), step_rule = 'previous'; end
    [m, L] = size(Xfull);  lt = numel(tvec);
    assert(L == lt-1, 'Xfull must have lt-1 columns.');

    % 1) Build original step change-times and values for each row
    t0 = tvec(1);  tN = tvec(end);
    Tdaily = (floor(t0) : ceil(tN)-1)';                 % daily support
    nd = numel(Tdaily);

    % “Announcement schedule” for behaviour (only shift from idx_start onward)
    tb = tvec(:);
    tb_shift = tb;
    tb_shift(idx_start:end) = tb_shift(idx_start:end) + delta_days;

    % Guard: keep strictly increasing (rare collisions with small deltas)
    if any(diff(tb_shift) <= 0)
        % Resolve by nudging by eps; behaviourally negligible, numerically safe
        epsb = 1e-6;
        d = diff(tb_shift);
        bad = find(d <= 0);
        for k = bad'
            tb_shift(k+1) = tb_shift(k) + max(epsb, tb_shift(k+1)-tb_shift(k) + epsb);
        end
    end

    % 2) Evaluate shifted step function daily, row-by-row
    Xfull_shift = zeros(m, L);
    for r = 1:m
        % row r values per original windows
        v = Xfull(r,:).';                          % Lx1
        % form left-continuous step knots/values for *shifted* schedule
        % Values apply on [tb_shift(j), tb_shift(j+1))
        tv = tb_shift(1:end-1);                    % Lx1 start times
        vv = v;                                    % Lx1 values

        % sample daily using step rule
        % map each day to its window index in shifted schedule
        % left-continuous: use previous value (bins with edges tv)
        switch lower(step_rule)
            case 'previous'
                % find last tv <= day
                idx = arrayfun(@(t) find(tv <= t, 1, 'last'), Tdaily, 'uni', 0);
            case 'next'
                % find first tv > day, then take previous window
                idx = arrayfun(@(t) find(tv > t, 1, 'first')-1, Tdaily, 'uni', 0);
            otherwise
                error('step_rule must be ''previous'' or ''next''.');
        end
        idx = cellfun(@(k) (isempty(k) * 0 + ~isempty(k) * k), idx);  % zeros for pre-first
        idx(idx < 1) = 1;                           % clamp before first change to first value
        idx(idx > L) = L;                           % clamp after last to last value
        v_daily = vv(idx);                          % nd x 1 daily values

        % 3) Re-aggregate back onto the *original* (unshifted) windows
        for j = 1:L
            dj = (Tdaily >= tvec(j)) & (Tdaily < tvec(j+1));
            if any(dj)
                Xfull_shift(r,j) = mean(v_daily(dj));  % average over days in that window
            else
                % If window narrower than a day, fall back to point sampling at its start
                Xfull_shift(r,j) = vv(max(1, find(tv <= tvec(j), 1, 'last')));
            end
        end
    end
end
