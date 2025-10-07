function [xsto, outsto, history, accept_rate, covmat] = MCMC_adaptive_noz(F, x0, n, sigma, fixinds, blockind, displ)
% Adaptive Metropolis (Haario et al., 2001) — simple, robust implementation.
% F(x) must return a log posterior (use -Inf outside bounds).
% x0: 1xd or dx1; sigma: scalar step scale OR dxd covariance matrix.
% fixinds: [] or 2×k matrix [idx; values] to hold coordinates fixed.
% blockind is ignored (kept for API compatibility).
% displ: true/false for very light console output.

if nargin < 7 || isempty(displ), displ = false; end
d = numel(x0);
x = x0(:);                          % d×1
% Set up fixed coordinates
if ~isempty(fixinds)
    fixed_idx  = fixinds(1,:); fixed_idx = fixed_idx(:);
    fixed_vals = fixinds(2,:); fixed_vals = fixed_vals(:);
else
    fixed_idx  = [];
    fixed_vals = [];
end

% Initial covariance:
if isscalar(sigma)
    C0 = eye(d);
    step_scale = (2.38^2/d) * sigma^2;     % Haario scaling with user factor
else
    C0 = (sigma + sigma')/2;               % ensure symmetric
    C0 = C0 + 1e-12*eye(d);                % SPD jitter
    step_scale = 1.0;                      % already scaled by user
end

% If we have fixed coords, zero out their rows/cols in C0:
if ~isempty(fixed_idx)
    C0(fixed_idx,:) = 0; C0(:,fixed_idx) = 0;
end

eps_jit   = 1e-9;               % jitter for SPD
adapt_burn = max(3*d, 500);     % start adapting after some burn-in
adapt_every = 10;               % refresh adaptive covariance this often
shrink     = 0.05;              % shrinkage toward spherical
target_acc = 0.234;             % Robbins–Monro step-size adaptation target
gamma_rm   = 0.05;              % small learning rate for step scale

% Storage
xsto     = zeros(d, n);
outsto   = -inf(1, n);
history  = zeros(n, d+1);       % proposals + accept flag
xsto(:,1) = x;
FX = F(x');                     % F expects row or column; allow either
if ~isfinite(FX), error('Initial state has -Inf log posterior.'); end
outsto(1) = FX;
acc_count = 0;

% Running mean/covariance trackers (Welford)
mu_run = x;
S_run  = zeros(d);              % sum of outer products for covariance
% Keep also an identity scale to guard degeneracy
I_d    = eye(d);
C_adapt = C0;

for t = 2:n
    % Current proposal covariance (Gaussian random walk)
    if t <= adapt_burn
        Cprop = step_scale * (C0 + eps_jit*I_d);
    else
        Cprop = step_scale * (C_adapt + eps_jit*I_d);
    end

    % Propose
    L = chol((Cprop + Cprop')/2 + 1e-12*I_d, 'lower');
    Y = x + L*randn(d,1);

    % Enforce fixed coords
    if ~isempty(fixed_idx)
        Y(fixed_idx) = fixed_vals;
    end

    % Evaluate F at proposal
    FY = F(Y'); if ~isfinite(FY), FY = -Inf; end

    % MH accept/reject
    dlog = FY - FX;             % symmetric proposal → no Hastings term
    if dlog >= 0 || log(rand) < dlog
        x  = Y;
        FX = FY;
        acc = 1;
        acc_count = acc_count + 1;
    else
        acc = 0;
    end

    % Store
    xsto(:,t) = x;
    outsto(t) = FX;
    history(t,1:d) = x(:)';
    history(t,end) = acc;

    % Update running mean/cov (Welford) and build adaptive covariance
    % Only update on *accepted* OR *every step*? Standard: every step.
    dt    = x - mu_run;
    mu_run = mu_run + dt / t;
    S_run  = S_run + (dt * (x - mu_run)');   % = sum (x - mean)*(x - mean)'

    % Refresh adaptive covariance occasionally after burn-in
    if t > adapt_burn && mod(t, adapt_every)==0
        % Empirical covariance
        Cemp = S_run / (t - 1);
        % Zero fixed coords
        if ~isempty(fixed_idx)
            Cemp(fixed_idx,:) = 0; Cemp(:,fixed_idx) = 0;
        end
        % Shrinkage → guard against degeneracy and help mixing early
        mvar = mean(diag(Cemp));
        C_shrunk = (1 - shrink)*Cemp + shrink*mvar*I_d;
        C_adapt  = (C_shrunk + C_shrunk')/2;

        % Robbins–Monro adapt of global step scale (very gentle)
        if t <= 10*adapt_burn   % stop adapting step size eventually
            acc_recent = mean(history(max(2,t-199):t,end));  % last 200
            step_scale = step_scale * exp(gamma_rm*(acc_recent - target_acc));
            % Keep step_scale within sane bounds
            step_scale = min(max(step_scale, 1e-6), 1e3);
        end
    end

    if displ && mod(t, round(n/20))==0
        fprintf('Iter %6d / %6d | acc=%.3f | step=%.3g\n', t, n, acc_count/t, sqrt(step_scale));
    end
end

accept_rate = acc_count / n;
% Return last adaptive covariance (unscaled by step_scale for inspection)
covmat = (t > adapt_burn) * C_adapt + (t <= adapt_burn) * C0;

% Match your expected shapes
xsto    = xsto.';     % n × d
history = history;    % already n × (d+1)
end
