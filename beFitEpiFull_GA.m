function [poptim,resnorm] = beFitEpiFull_GA(ydata,X,data,thetaIn,Xfull,coeff,tvec,b0)
%[tvec, Xfull_crude, Xfull_pca, pca_info, Xdaily_crude, Xdaily_pca] = build_indices_with_pca(finalData3, 487);
hlag = 0;
plotRun = 0;
intrinsic = 1;
%projection = 0;
%centre=0;
%oldParamsIn=0;

if numel(tvec) >= 19 && hlag ~= 0
    tvec(19:end) = tvec(19:end) + hlag;
end
xdata = 85:tvec(end);
[lx1,lx2] = size(X);

% admissions slice and England scale (as in your code)
ydata = ydata(0+(1:numel(xdata)));
ydata = ydata*(sum(data.Npop)/56286961);
ymean=0;%mean(ydata);

% weights (same as your weighted SSE)
weights = 1./(1 + ydata(:));

% forward model handle EXACTLY as you call it
fun = @(params,xd) sim2fit_global(params, data, xd, X, intrinsic, Xfull, coeff, tvec, lx1, lx2, plotRun, ymean, b0);

% ===== bounds (use exactly your choices) =====
% params meaning in your sim2fit:
%   alpha = params(1) replicated to [1,1,1]
%   reducedParams = [1, params(2:end)] -> passed into bePrepCovid19 as forFeedback vector
% Keep these bounds consistent with your modeling choices.
if length(thetaIn)==5
    lb = [0, -40, -40,  -40,  1];
    ub = [1,  40,  40,  200,  200];
elseif length(thetaIn)==4
    lb = [0, -40, -40,  1];
    ub = [1,  40,  40,  120];
elseif length(thetaIn)==3
    lb = [0, -40,  -40];
    ub = [1,  40,  40];
elseif length(thetaIn)==2
    lb = [-40,  -40];
    ub = [40,  40];
end
nvars = numel(lb);

% ===== objective: weighted SSE on admissions =====
% choose the loss you want:
[obj, obj_vec] = make_obj(fun, xdata, ydata(:), 'poiss');   % or 'ols','wls1','wls05' % scheme: 'wls1' (1/(1+y)), 'wls05' (1/sqrt(1+y)), 'ols', 'poiss'
%@(z) obj_full_ga(z, fun, xdata, ydata(:), weights);

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
% choose the loss you want:
[obj, obj_vec] = make_obj(fun, xdata, ydata(:), 'poiss');   % or 'ols','wls1','wls05'

% --- GA ---
[z_ga, sse_ga] = ga(obj, nvars, [],[],[],[], lb, ub, [], [], opts);

% --- manual point for comparison ---
sse_manual = obj(thetaIn(:)');

% --- local polish (use residuals that match the chosen loss) ---
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

ymod = fun(poptim, xdata);
figure('Units','centimeters','Position',[0 0 20 20]); hold on;
bar(xdata, ydata, 'FaceAlpha',0.5, 'EdgeColor','none');
plot(xdata, ymod, 'r-', 'LineWidth', 2.2);
xlabel('Time'); ylabel('Hospital admissions'); box on; grid on; axis square;
title(sprintf('Full fit: alpha=%.4f | SSE=%.3g', poptim(1), resnorm));
legend({'Data','Model'}, 'Location','best');

end

function [f,rhohat]=sim2fit(params,data,xdata,Xfit,intrinsic,Xfull,coeff,tvec,lx1,lx2,plotRun,ymean,arg)
R0=2.8;
tvec(1)=-80;
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
    %
    pc=2;%PC1 OR PC2
    kstar=params(2);
    k3=params(3);
    delta=params(4);%delta>eps_safe
    v0pca=-kstar*arg(pc)-delta;
    if pc==1
        reducedParams=[1,kstar,0,k3,v0pca];   
    else
        reducedParams=[1,0,kstar,k3,v0pca];                      
    end
    %}
elseif length(params)==3
    reducedParams=[1,0,0,params(2:3)];
    %{
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

function [objfun, resvec_fun] = make_obj(fun, xdata, y, scheme)
% Returns:
%   objfun(z)     -> scalar objective for GA/fmin*
%   resvec_fun(z) -> residual vector for lsqnonlin (matching the loss)
    y = y(:);  % column

    switch lower(scheme)
        case 'ols'
            objfun     = @(z) sse_scalar( vec(fun(z,xdata)), y );
            resvec_fun = @(z) vec(fun(z,xdata)) - y;

        case 'wls1'   % weights = 1/(1+y)
            w          = 1./(1+y);
            objfun     = @(z) sse_scalar( vec(fun(z,xdata)), y, w );
            resvec_fun = @(z) sqrt(w).*( vec(fun(z,xdata)) - y );

        case 'wls05'  % weights = 1/sqrt(1+y)
            w          = 1./sqrt(1+y);
            objfun     = @(z) sse_scalar( vec(fun(z,xdata)), y, w );
            resvec_fun = @(z) sqrt(w).*( vec(fun(z,xdata)) - y );

        case 'poiss'  % Poisson deviance
            objfun     = @(z) pois_dev_scalar( vec(fun(z,xdata)), y );
            resvec_fun = @(z) pois_dev_components( vec(fun(z,xdata)), y );

        otherwise
            error('Unknown scheme');
    end
end

% ---- helpers ----
function v = vec(x), v = x(:); end

function val = sse_scalar(yhat, y, w)
    yhat = yhat(:);
    if nargin < 3
        val = sum( (yhat - y).^2 );
    else
        val = sum( w(:) .* (yhat - y).^2 );
    end
end

function val = pois_dev_scalar(yhat, y)
    yhat = max(yhat(:), eps);
    y    = max(y(:),    eps);
    % 2 * sum( yhat - y + y*log(y/yhat) )
    val  = 2*sum( yhat - y + y .* log( y ./ yhat ) );
end

function r = pois_dev_components(yhat, y)
    % component-wise residuals for lsqnonlin: sign * sqrt(deviance_i)
    yhat = max(yhat(:), eps);
    y    = max(y(:),    eps);
    d    = 2*( yhat - y + y .* log( y ./ yhat ) );
    r    = sign(y - yhat) .* sqrt(d);
end


%{
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

if centre==1
    %means=mean(Xfull(:,4:end),2);%[0.5344    0.5247    0.1490];
    %Xfull(:,4:end)=Xfull(:,4:end)-repmat(means,1,size(Xfull,2)-3);
    
    % --- center + rotate v1,v2 to reduce VIF ---
    % Xfull is 2 x T (rows: v1, v2). Only center/rotate from the first “real” window.
    mu12   = mean(Xfull(1:2,4:end), 2);        % 2x1 means across time
    W      = Xfull(1:2,4:end) - mu12;          % center
    [U,~,~]= svd(W,'econ');                    % 2x2 orthonormal rotation (PCA/QR both OK)
    % Replace v1,v2 in Xfull by *centered+rotated* components (orthogonal)
    Xfull(1:2,4:end) = U' * (Xfull(1:2,4:end) - mu12);
    mu12_rot = U' * mu12;                      % needed for intercept correction
    
    if oldParamsIn==1
        alpha_old=thetaIn(1); k1_old=thetaIn(2); k2_old=thetaIn(3); k3_old=thetaIn(4); v0_old=thetaIn(5);
        % rotate old k's
        b = U' * [k1_old; k2_old];          % -> b1,b2 seed
        % pick eta0 so that softplus matches desired x(0) = v0_old - [k1_old,k2_old]*mu12
        x0_old = v0_old - [k1_old k2_old] * mu12;
        % solve -softplus(eta0) = x0_old  -> softplus(eta0) = -x0_old
        softplus_inv = @(y) log(exp(y) - 1);              % for y>0; guard numerically
        eta0_seed = softplus_inv(max(1e-6, -x0_old));
        thetaIn = [alpha_old, b(1), b(2), k3_old, eta0_seed];
    end
%else
    %mu12_rot=0;
end

% ===== timeline / windows (your original logic) =====
if projection==1
    tvec = [1    32    61    93   107   169   176   223   227   230   250   258   265   271   279   294   310   322   330   338   349   370   384    407   418   433   445   463   474   491   504   517   540   561   567   575   594   605   617   631   642   652   661   679   693   702    716   742   784];
    xdata = 85:763;
else
    tvec = [1,2,61,91,127,134,141,148,153,155,162,167,169,175,176,186,200,211,216,218,223,227,230,236,[250,258,265,266,271,279,288,294,305,310,322,330,337,338,349,354,356,361,370,372,384,397,407,418,433,445,454,463,468,474,491,503,504,517,540,561,566,567,575]+hlag];
    xdata = 85:tvec(end-7);
    tvec=tvec(1:end-7);
    %tvec=tvec(1:16);
    %tvec(end)=183;%1st July
end

lt = numel(tvec);
X     = X(:,1:lt-1);
Xfull = Xfull(:,1:lt-1);
%}
