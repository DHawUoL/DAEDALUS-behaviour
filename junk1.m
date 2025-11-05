function [pbest, yhat, bestLagDays, resnorm, out] = beFitEpiFull_GA_with_lag(ydata, data, thetaIn, Xfull, tvec, b0, loss, lag_range)

%[tvec, Xfull_crude, Xfull_pca, pinfo, Xdaily_crude, Xdaily_pca, b0] = build_indices_from_changes(IndexTable, 487, struct('tfrom',91,'min_win',1,'negate_stringency',true));
%bePlotFitSimple(dataAdmissions, dataUK2, pbest, Xfull_crude, tvec, b0)
%plot(Xdaily_crude','linewidth',2); hold on; plot(85:tvec(end),dataAdmissions(1:tvec(end)-84)/4e3,'linewidth',2)

X=ones(size(Xfull,2));
coeff=ones(1,3)';
%Old version:
%[params5,Ypred,delta,resnorm]=beFitEpiFull_GA(dataAdmissions,ones(1,size(Xfull_crude,2)),dataUK2,[params5],Xfull_crude,ones(1,3)',tvec);
%This version:
%[popt, Ypred, delta_hat, res] = beFitEpiFull_GA(dataAdmissions, dataUK, thetaIn, Xfull_crude, tvec);

% Adds one extra integer decision variable: lag/lead in days (−7…+7).
% GA treats it as an integer index; internally we map to actual day-shift.
%
% Inputs:
%   loss: 'poiss' | 'wls1' | 'wls05' | 'ols'   (default 'poiss')
%   lag_range: e.g. -7:7
%
% Outputs:
%   pbest        : best continuous params (same length as thetaIn)
%   yhat         : model at best solution
%   bestLagDays  : chosen lag in days (positive = drivers lead)
%   resnorm      : objective value at best solution
%   out          : struct with fields (lag_idx, pop, etc.)

    if nargin < 9 || isempty(loss),      loss = 'poiss'; end
    if nargin < 10 || isempty(lag_range), lag_range = -7:7; end

    xdata = 85:tvec(end);
    y     = ydata(0+(1:numel(xdata))) * (sum(data.Npop)/56286961);
    y     = y(:);
    ymean = 0;

    [~,lx2] = size(X);
    first_real_col = 4;               % after the 3 forced-zero windows
    lags = lag_range(:)';             % row vector, e.g. [-7 -6 ... 6 7]

    % --- bounds for continuous params (as in your code) ---
    switch numel(thetaIn)
        case 5, lb = [0, 0, -40, 0,  1];  ub = [1, 40, 0, 40, 80];
        case 4, lb = [0, -40, -40,  1];       ub = [1, 40, 40, 120];
        case 3, lb = [0, -40, -40];           ub = [1, 40, 40];
        case 2, lb = [-40, -40];              ub = [40, 40];
        otherwise, error('Unexpected thetaIn length.');
    end

    % --- one extra integer var: lag index in 1..numel(lags) ---
    lbz = [lb, 1];
    ubz = [ub, numel(lags)];
    IntCon = numel(lbz);          % last variable is integer

    % --- objective closures (wrap simulator with lag) ---
    sim_core = @(p, Xf) sim2fit_global(p, data, xdata, 1, Xf, tvec, lx2, 0, ymean, b0);

    % scalar objective for GA (depends on loss)
    [obj_scalar, obj_resvec] = make_loss(sim_core, y, loss);

    % z = [params, lag_idx]
    fun_scalar = @(z) obj_scalar(@()sim_with_lag(z));
    fun_vec    = @(z) obj_resvec(@()sim_with_lag(z));

    % --- GA seeds/opts ---
    popSize = 120;
    initPop = repmat(lbz, popSize, 1) + rand(popSize, numel(lbz)).*repmat(ubz-lbz, popSize, 1);
    initPop(1,1:numel(thetaIn)) = thetaIn(:)';   % seed params
    initPop(:,end) = round(initPop(:,end));      % make integer column valid

    opts = optimoptions('ga', ...
        'PopulationSize', popSize, 'MaxGenerations', 250, ...
        'EliteCount', max(2, round(0.02*popSize)), ...
        'MutationFcn', {@mutationadaptfeasible}, ...
        'CrossoverFraction', 0.8, ...
        'FunctionTolerance', 1e-6, ...
        'MaxStallGenerations', 80, ...
        'InitialPopulationMatrix', initPop, ...
        'UseParallel', true, ...
        'Display','iter');

    % --- GA ---
    rng default
    [z_ga, f_ga] = ga(@(z) fun_scalar(z), numel(lbz), [], [], [], [], lbz, ubz, [], IntCon, opts);

    % --- local polish on continuous part ---
    opts_ls = optimoptions('lsqnonlin','Display','off','MaxIterations',1e3, ...
        'StepTolerance',1e-8,'FunctionTolerance',1e-8, ...
        'FiniteDifferenceType','forward','FiniteDifferenceStepSize',1e-3);

    lag_idx = max(1, min(numel(lags), round(z_ga(end))));
    lag_days = lags(lag_idx);
    Xfull_shift = shift_Xfull_by_days(Xfull, tvec, lag_days, first_real_col);

    objvec_cont = @(p) obj_resvec(@() sim_core(p, Xfull_shift));
    [p_loc, ~, ~] = lsqnonlin(objvec_cont, z_ga(1:end-1), lb, ub, opts_ls);

    % --- choose best of GA vs hybrid ---
    f_ga2   = fun_scalar(z_ga);
    f_hyb   = obj_scalar(@() sim_core(p_loc, Xfull_shift));
    if f_hyb <= f_ga2
        pbest = p_loc;
        bestLagDays = lag_days;
        [yhat, ~] = sim_core(pbest, Xfull_shift);
        resnorm = f_hyb;
    else
        pbest = z_ga(1:end-1);
        lag_idx = max(1, min(numel(lags), round(z_ga(end))));
        bestLagDays = lags(lag_idx);
        Xfull_shift = shift_Xfull_by_days(Xfull, tvec, bestLagDays, first_real_col);
        [yhat, ~] = sim_core(pbest, Xfull_shift);
        resnorm = f_ga2;
    end

    out = struct('lag_idx', find(lags==bestLagDays,1), 'lag_days', bestLagDays, ...
                 'ga_point', z_ga, 'ga_obj', f_ga);

    % ------- nested helpers -------
    function [f, rhohat] = sim_with_lag(z)
        lag_idx = max(1, min(numel(lags), round(z(end))));
        del     = lags(lag_idx);
        Xs      = shift_Xfull_by_days(Xfull, tvec, del, first_real_col);
        [f, rhohat] = sim_core(z(1:end-1), Xs);
    end
end

function [obj_scalar, obj_resvec] = make_loss(sim_core, y, scheme)
% sim_core: @(p, Xf) -> [yhat, rhohat]
    switch lower(scheme)
        case 'poiss'
            obj_scalar = @(simcall) pois_dev_scalar( vec(call_sim(simcall)), y );
            obj_resvec = @(simcall) pois_dev_components( vec(call_sim(simcall)), y );
        case 'wls1'
            w = 1./(1+y);
            obj_scalar = @(simcall) sse_scalar( vec(call_sim(simcall)), y, w );
            obj_resvec = @(simcall) sqrt(w).*( vec(call_sim(simcall)) - y );
        case 'wls05'
            w = 1./sqrt(1+y);
            obj_scalar = @(simcall) sse_scalar( vec(call_sim(simcall)), y, w );
            obj_resvec = @(simcall) sqrt(w).*( vec(call_sim(simcall)) - y );
        case 'ols'
            obj_scalar = @(simcall) sse_scalar( vec(call_sim(simcall)), y );
            obj_resvec = @(simcall) vec(call_sim(simcall)) - y;
        otherwise
            error('Unknown loss');
    end

    function yhat = call_sim(simcall)
        [yhat, ~] = simcall();
    end
end

function Xfull_shift = shift_Xfull_by_days(Xfull, tvec, delta_days, first_real_col)
% Re-aggregate a shifted step driver back onto the ORIGINAL tvec windows.
% Positive delta_days means behaviour leads (starts earlier).
    if delta_days == 0
        Xfull_shift = Xfull; return;
    end
    [m, L] = size(Xfull);
    Xfull_shift = zeros(m, L);

    tb = tvec(:);
    tb_shift = tb;
    tb_shift(first_real_col:end) = tb_shift(first_real_col:end) + delta_days;

    % guard monotonicity (rare collisions with small deltas)
    bad = find(diff(tb_shift) <= 0);
    for k = bad'
        tb_shift(k+1) = tb_shift(k) + 1e-6;
    end

    % daily sampling & re-aggregation
    Tdaily = (tvec(1) : tvec(end)-1)';         % 1-day grid
    for r = 1:m
        v = Xfull(r,:).';
        tv = tb_shift(1:end-1);
        vv = v;

        % left-continuous sampling
        idx = arrayfun(@(t) find(tv <= t, 1, 'last'), Tdaily, 'uni', 0);
        idx = cellfun(@(k) (isempty(k) * 1 + ~isempty(k) * k), idx);
        idx(idx < 1) = 1; idx(idx > L) = L;

        v_daily = vv(idx);
        for j = 1:L
            dj = (Tdaily >= tvec(j)) & (Tdaily < tvec(j+1));
            if any(dj)
                Xfull_shift(r,j) = mean(v_daily(dj));
            else
                Xfull_shift(r,j) = vv(max(1, find(tv <= tvec(j), 1, 'last')));
            end
        end
    end
end

% ---- small utilities used above ----
function v = vec(x), v = x(:); end
function val = sse_scalar(yhat, y, w)
    yhat = yhat(:);
    if nargin < 3, val = sum((yhat - y).^2);
    else,           val = sum(w(:).*(yhat - y).^2);
    end
end
function val = pois_dev_scalar(yhat, y)
    yhat = max(yhat(:), eps); y = max(y(:), eps);
    val  = 2*sum( yhat - y + y.*log(y./yhat) );
end
function r = pois_dev_components(yhat, y)
    yhat = max(yhat(:), eps); y = max(y(:), eps);
    d    = 2*( yhat - y + y.*log(y./yhat) );
    r    = sign(y - yhat).*sqrt(d);
end
