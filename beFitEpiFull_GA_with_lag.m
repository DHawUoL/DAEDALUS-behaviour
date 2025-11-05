function [pbest, yhat, bestLagDays, resnorm, out] = beFitEpiFull_GA_with_lag(ydata, data, thetaIn, Xfull, tvec, b0, loss, lag_range)

    if nargin < 9 || isempty(loss),       loss = 'poiss'; end
    if nargin < 10 || isempty(lag_range), lag_range = -7:7; end

    xdata = 85:tvec(end);
    y     = ydata(1:numel(xdata)) * (sum(data.Npop)/56286961);
    y     = y(:);
    ymean = 0;

    X = ones(size(Xfull,2));                 % (dummy; your sim ignores it)
    coeff = ones(1,3)';                      % (dummy; crude mode)
    [~,lx2] = size(X);

    first_real_col = 4;                      % after the 3 forced zeros
    anchor_day     = 250;                    % <-- start shifting from here
    lags           = lag_range(:)';

    % --- bounds (use your preferred ones) ---
    switch numel(thetaIn)
        case 5, lb = [0, -10, -10, -10,  1];  ub = [1, 10, 10, 10, 30];
        case 4, lb = [0, -40, -40,   1];      ub = [1, 80, 40, 120];
        case 3, lb = [0, -40, -40];           ub = [1, 40, 40];
        case 2, lb = [-40, -40];              ub = [40, 40];
        otherwise, error('Unexpected thetaIn length.');
    end

    % integer lag index
    lbz = [lb, 1];  ubz = [ub, numel(lags)];
    IntCon = numel(lbz);

    % ---- simulator core: takes tvec_in (NOT Xf) ----
    sim_core = @(p, tvec_in) sim2fit_global(p, data, xdata, 1, Xfull, tvec_in, lx2, 0, ymean, b0);

    % ---- choose loss ----
    [obj_scalar, obj_resvec] = make_loss(sim_core, y, loss);

    % z = [params, lag_idx]
    fun_scalar = @(z) obj_scalar(@() sim_with_lag(z));
    fun_vec    = @(z) obj_resvec(@() sim_with_lag(z));

    % ---- GA ----
    popSize = 120;
    initPop = repmat(lbz, popSize, 1) + rand(popSize, numel(lbz)).*repmat(ubz-lbz, popSize, 1);
    initPop(1,1:numel(thetaIn)) = thetaIn(:)';
    initPop(:,end) = round(initPop(:,end));

    opts = optimoptions('ga','PopulationSize',popSize,'MaxGenerations',250, ...
        'EliteCount',max(2,round(0.02*popSize)),'MutationFcn',{@mutationadaptfeasible}, ...
        'CrossoverFraction',0.8,'FunctionTolerance',1e-6,'MaxStallGenerations',80, ...
        'InitialPopulationMatrix',initPop,'UseParallel',true,'Display','iter');

    rng default
    [z_ga, f_ga] = ga(@(z) fun_scalar(z), numel(lbz), [], [], [], [], lbz, ubz, [], IntCon, opts);

    % ---- local polish on continuous vars at GA lag ----
    opts_ls = optimoptions('lsqnonlin','Display','off','MaxIterations',1e3, ...
        'StepTolerance',1e-8,'FunctionTolerance',1e-8, ...
        'FiniteDifferenceType','forward','FiniteDifferenceStepSize',1e-3);

    lag_idx   = max(1, min(numel(lags), round(z_ga(end))));
    lag_days  = lags(lag_idx);
    tvec_shift= shift_tvec_after(tvec, lag_days, anchor_day, first_real_col, true);

    objvec_cont = @(p) obj_resvec(@() sim_core(p, tvec_shift));
    p_loc = lsqnonlin(objvec_cont, z_ga(1:end-1), lb, ub, opts_ls);

    % ---- pick best ----
    f_ga2 = fun_scalar(z_ga);
    f_hyb = obj_scalar(@() sim_core(p_loc, tvec_shift));

    if f_hyb <= f_ga2
        pbest = p_loc; bestLagDays = lag_days; resnorm = f_hyb;
        [yhat, ~] = sim_core(pbest, tvec_shift);
    else
        pbest = z_ga(1:end-1); bestLagDays = lag_days; resnorm = f_ga2;
        [yhat, ~] = sim_core(pbest, tvec_shift);
    end

    out = struct('lag_days',bestLagDays, 'ga_point',z_ga, 'ga_obj',f_ga);

    % ---- nested ----
    function [f, rhohat] = sim_with_lag(z)
        li   = max(1, min(numel(lags), round(z(end))));
        d    = lags(li);
        tvs  = shift_tvec_after(tvec, d, anchor_day, first_real_col, true);
        [f, rhohat] = sim_core(z(1:end-1), tvs);
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

function tvec2 = shift_tvec_after(tvec, delta_days, anchor_day, first_real_col, keep_end)
% Shift breakpoints >= anchor_day (but not the first 3 windows) by +delta_days.
% keep_end=true keeps tvec(end) fixed.
    if nargin < 5, keep_end = true; end
    tvec2 = tvec(:)';
    j0 = max(first_real_col, find(tvec2 >= anchor_day, 1, 'first'));
    if isempty(j0), return; end
    if keep_end
        tvec2(j0:end-1) = tvec2(j0:end-1) + delta_days;
    else
        tvec2(j0:end)   = tvec2(j0:end)   + delta_days;
    end
    % enforce strictly increasing
    epsb = 1e-6;
    for k = j0:(numel(tvec2)-1)
        if tvec2(k+1) <= tvec2(k)
            tvec2(k+1) = tvec2(k) + epsb;
        end
    end
end

