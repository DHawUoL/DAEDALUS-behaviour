% Adaptive MCMC, using Haario et al:
% https://link.springer.com/article/10.1007/s11222-008-9110-y
% and
% http://probability.ca/jeff/ftpdir/adaptex.pdf

% F:          Function giving log-posterior density for a parameter set x
% x0:         Initial value of parameter set x
% n:          Number of iterations
% cov0:       Initial covariance matrix
% fac:        Scaling factor for covariance matrix. Set fac = 1 for defaultar
% fixinds:    Elements of x that should be held fixed. Set fixinds = [] for full MCMC
% blockinds:  Number of 'epi parameters' (e.g. beta, X2) in x, if we want to vary epi and non-epi parameters as independent 'blocks'. Set blockinds = [] if we want a full covariance matrix
% displ:      Structure with display options. Set displ = true to show progress

function [xsto, outsto, history, accept_rate,covmat] = MCMC_adaptive(F, x0, n, sigma, fixinds, blockind, displ, ind_opts)%, bounds)

if nargin < 8 || isempty(ind_opts), ind_opts = struct; end
if ~isfield(ind_opts,'mu'), ind_opts.mu = x0(:)'; end  % safe default
if ~isfield(ind_opts,'C'),  ind_opts.C  = eye(length(x0)); end
if ~isfield(ind_opts,'p_ind'), ind_opts.p_ind = 0.00; end
if ~isfield(ind_opts,'p_t'),   ind_opts.p_t   = 0.00; end
if ~isfield(ind_opts,'nu'),    ind_opts.nu    = 5;    end

global_scale=0.75;

%d = length(x0); b = 0.05; sd = sigma*2.4^2/d; %DH commented out
%d = length(x0); b = 0.05; sd = sigma*2.4^2/d; %DH added to change b etc.
d = length(x0); b = .15; %sd = sigma*2.4^2/d; %DH added to change b etc. %xx
target_acc   = 0.24;          % aim for ~30%
log_scale    = -.05;           % proposal multiplier = exp(log_scale)
scale_mult = exp(log_scale);
%adapt_start  = max(600, 6*d); % start adapting after some burn-in
nfree      = d - size(fixinds,2);                 % if you ever fix bits
adapt_start = max(300, 4*max(1,nfree));           % was max(600,6*d)
window       = 1000;           % covariance window (last W samples)
shrink       = 0.03;          % shrinkage toward spherical (10%)
recompute_every = 10;
jitter       = 1e-12;         % SPD jitter
% ----- adaptation freeze controls -----
auto_freeze        = false;                % also allow auto-freeze
freeze_band        = [0.18 0.32];         % target acceptance band (acc200)
freeze_check_every = 200;                 % how often to check stability
stable_hits_needed = 3;                   % require K consecutive hits
stable_hits        = 0;
adaptation_frozen  = false;               % state flag

vridge = [-1; 0.5; 0.5; 0.5];  % same shape you used in the penalty
vridge = vridge / norm(vridge);
p_ridge = 0.05;                      % 5% of iterations
ridge_scale = 0.25;                  % tune: 0.2–0.5 works

if isscalar(sigma)
    sigmode = "scalar";
elseif isvector(sigma) && numel(sigma)==d
    sigmode = "vector";                 % per-dimension scales
    sigma = sigma(:)';                  % row
elseif ismatrix(sigma) && all(size(sigma)==[d d])
    sigmode = "matrix";                 % full covariance
    sigma = (sigma + sigma.')/2 + 1e-12*eye(d);  % ensure SPD
else
    error('sigma must be scalar, length-%d vector, or %d×%d matrix', d, d, d);
end

if ~isempty(fixinds)
    inds = fixinds(1,:); vals = fixinds(2,:);
else
    inds = []; vals = [];
end


% Checks on the initial covariance matrix
cov0 = eye(d);
if ~isempty(fixinds)
    inds = fixinds(1,:); 
    cov0(inds,:) = 0;
    cov0(:,inds) = 0;
end
covmat=cov0;

% Initiate the output matrices
xsto = zeros(d,n); outsto = zeros(1,n);
history = zeros(d+1,n);                                                    % Rows: 1:d Proposed values 4. Accept or reject

xsto(:,1) = x0(:); xbar = xsto;
FX = F(x0); outsto(1) = FX;
acc = 0;

if displ; figure; end

% --- Start the MCMC loop -------------------------------------------------

switch sigmode
  case "scalar",  baseCov = (global_scale^2/d) * cov0;
  case "vector",  S = diag(sigma); baseCov = (global_scale^2/d) * (S*cov0*S);
  case "matrix",  baseCov = (global_scale^2/d) * sigma;
end

% initialize throttled adaptive cache ONCE per call
baseCov_scaled   = baseCov * exp(log_scale);
sdmatrix_cached  = baseCov_scaled;% instead of (2.4^2/d)*cov0

% base probs
p_de_base     = 0.40;
p_1d_base     = 0.15;
p_pc2_base    = 0.15;
p_global_base = 0.07;

for t = 2:n
    X = xsto(:,t-1);



    FY = [];    % clear any carry-over; we'll fill it in below
    hastings = 0;   % proposal log-correction (only nonzero for independence)
    if t < adapt_start
        % warm-up: fixed kernel
        Y = mvnrnd(X', baseCov);
    else
        % --- refresh the independence (Laplace-like) proposal from recent z states ---
        if isfield(ind_opts,'refresh_every') && mod(t, ind_opts.refresh_every)==0
            % use a recent window of chain states in *z*-space
            W   = 1000;                          % window length
            i1  = max(2, t - W);                 % avoid the very first column
            Zwin = xsto(:, i1:t-1);              % d × (#window)
    
            % Option A: use only ACCEPTED states (recommended if available)
            acc_idx = find(history(end,1:t-1) == 1);          % indices with accept=1
            acc_idx = acc_idx(acc_idx >= i1);                 % restrict to window
            if numel(acc_idx) >= 10
                Zwin = xsto(:, acc_idx);
            end
    
            % centre and shape for independence proposal
            mu_ind = mean(Zwin, 2)';                           % 1×d
            C_ind  = cov(Zwin'); C_ind = (C_ind + C_ind')/2;   % d×d
            % regularise + inflate slightly
            C_ind  = 2.0^2 * (C_ind + 1e-9*eye(size(C_ind,1)));
    
            % update the options used by the mixture "independence" kernel
            ind_opts.mu = mu_ind;
            ind_opts.C  = C_ind;
        end
        % --- refresh adaptive covariance if due ---
        if ~adaptation_frozen && mod(t, recompute_every) == 0
            ind1 = t-1; ind0 = max(2, ind1 - window + 1);
            Cb   = cov(xsto(:,ind0:ind1)');                     % recent window
            if ~isempty(inds), Cb(inds,:) = 0; Cb(:,inds) = 0; end
            mtrace = mean(diag(Cb));
            Cb = (1 - shrink)*Cb + shrink*mtrace*eye(d);        % shrinkage
            Cb = (Cb + Cb')/2 + jitter*eye(d);                  % SPD

            % NEW: absolute variance floors in z-space (prevents α dimension collapsing)
            min_var_all   = 2^2;     % general floor for every dim in z
            min_var_alpha = 2^2;    % stronger floor for z_alpha specifically
            Cb = Cb + diag(max(0, min_var_all - diag(Cb)));
            Cb(1,1) = max(Cb(1,1), min_var_alpha);
    
            sdmatrix = (2.4^2/d) * exp(log_scale) * Cb;         % Haario scaling
            [V,D] = eig((sdmatrix + sdmatrix.')/2);
            evals   = max(diag(D), 1e-12);
            medEval = median(evals);
            low_abs   = 1e-6;               % absolute floor
            low_rel   = medEval/100;
            evals     = max(evals, max(low_abs, low_rel));
            evals   = min(max(evals, medEval/100), 10*medEval); % cap extremes
            sdmatrix_new = V*diag(evals)*V.' + 1e-12*eye(d);
    
            % smooth update to the cached proposal
            alpha_smooth      = 0.2;
            sdmatrix_cached   = (1-alpha_smooth)*sdmatrix_cached + alpha_smooth*sdmatrix_new;
            covmat            = Cb;                              % expose last C
        end
    
        % ---------- PROPOSAL MIXTURE (DE / 1D / 2D / GLOBAL / GAUSSIAN) ----------
        % probabilities for each kernel (must sum to <= 1)
        p_de     = p_de_base;   % differential-evolution jump
        p_1d     = p_1d_base;   % 1-D oriented along top PC
        p_pc2    = p_pc2_base;   % 2-D oriented in top PC plane          << new
        p_global = p_global_base;   % global "reset" from broad distribution << new
        %p_ind = 0.05;
        %p_t = 0.05;  nu = 5;
        lap_mu  = ind_opts.mu;            % centre at LSQ/MAP
        lap_C   = covmat * 2.0 + 1e-6*eye(d);  % inflate a bit; refresh every 500 iters
        
        % small helpers
        step_mult   = exp(log_scale);            % your adaptive step-size multiplier

        % --- decide which kernels are eligible this iteration ---
        acc_idx = find(history(end,1:t-1) == 1);
        acc_idx = acc_idx(:)';
        % require a minimum number of *distinct* accepted states
        uniq_ok = numel(unique(acc_idx)) >= max(30, 6*d);

        % recent acceptance (use same window you like elsewhere)
        w0  = min(200, t-1);
        accW = mean(history(end, t-w0:t-1));
        

        
        p_ind = ind_opts.p_ind;
        p_t   = ind_opts.p_t;



        
        % ... compute p_de, p_1d, p_pc2, p_global as you already do ...
        psum  = p_de + p_1d + p_pc2 + p_global + p_ind + p_t + p_ridge;
        scale = min(0.85 / max(psum, eps), 1);
        p_de   = p_de   * scale;  p_1d  = p_1d  * scale;
        p_pc2  = p_pc2  * scale;  p_global = p_global * scale;
        p_ind  = p_ind  * scale;  p_t   = p_t   * scale;
        p_ridge= p_ridge* scale;
        p_joint = 1 - (p_de + p_1d + p_pc2 + p_global + p_ind + p_t + p_ridge);

        %{
        % eligibility (as you already had)
        uniq_ok    = numel(unique(acc_idx)) >= max(30, 6*d);
        use_de     = (t >= adapt_start + 150)  && uniq_ok;
        use_global = (t >= adapt_start + 400)  && uniq_ok;

        % tilt weights: when acc is low, push toward DE/global; when high, favor Gaussian/subspace
        low  = 0.18; high = 0.28;
        tilt = max(0,min(1,(high - accW)/(high - low)));    % 0..1, larger when acc is low
        
        p_de     = p_de_base     * use_de     * (1 + 0.7*tilt);
        p_global = p_global_base * use_global * (1 + 1.0*tilt);
        p_1d     = p_1d_base     * (1 - 0.3*tilt);
        p_pc2    = p_pc2_base    * (1 - 0.3*tilt);
        
        % renormalize to leave room for Gaussian
        psum   = p_de + p_1d + p_pc2 + p_global;
        scale  = min(0.85/ max(psum, eps), 1);              % keep <= 0.85 total mass
        p_de   = p_de   * scale;
        p_1d   = p_1d   * scale;
        p_pc2  = p_pc2  * scale;
        p_global = p_global * scale;
        p_joint  = 1 - (p_de + p_1d + p_pc2 + p_global + p_ind + p_t);    % Gaussian gets the rest
        %}

        u = rand;
        
        if numel(acc_idx) < 2
            % Not enough accepted states yet → fallback
            Y = mvnrnd(X', sdmatrix_cached);   % or baseCov_scaled
        elseif u < p_de
            % ----- Differential Evolution (tempered + crossover) -----
            r1 = acc_idx(randi(numel(acc_idx)));
            r2 = acc_idx(randi(numel(acc_idx)));
            while r2 == r1, r2 = acc_idx(randi(numel(acc_idx))); end
            
            gamma0   = 2.38 / sqrt(2*d);
            warmfac  = min(1.0, (t - adapt_start)/1200);             % warm faster
            gamma_de = (0.85 + 0.25*warmfac) * gamma0;               % a tad larger
            
            diffv = xsto(:,r1) - xsto(:,r2);
            
            % binomial crossover: include each dim with prob Cr, but force at least one
            Cr  = 0.9;                                              % high crossover
            mask = rand(d,1) < Cr;
            if ~any(mask), mask(randi(d)) = true; end
            diffv = diffv .* mask;
            
            jit  = mvnrnd(zeros(1,d), 1e-6 * (sdmatrix_cached + 1e-12*eye(d)))';
            Y    = ( X + gamma_de * diffv + jit )';
            
        elseif u < p_de + p_1d
            % ----- 1-D oriented step along top PC of recent covariance
            % use accepted times for local PCs when possible
            win_pc  = 400;
            acc_win = acc_idx(acc_idx >= max(2, t-win_pc) & acc_idx <= t-1);
            if numel(acc_win) >= max(10, 2*d)
                ids = acc_win;
            else
                ids = max(2, t-win_pc):(t-1);
            end
            Cpc = cov(xsto(:, ids)');
            Cpc = (Cpc + Cpc')/2 + 1e-12*eye(d);
            [V1,D1] = eig(Cpc);
            [~,ix]  = max(diag(D1));
            v1 = V1(:,ix);
            s1 = sqrt(max(D1(ix,ix), 1e-12));
            step = step_mult * s1 * randn;
            Y = (X + step * v1)';
        
        elseif u < p_de + p_1d + p_pc2
            % ----- 2-D oriented step in the plane of the top two PCs (new)
            % use accepted times for local PCs when possible
            win_pc  = 400;
            acc_win = acc_idx(acc_idx >= max(2, t-win_pc) & acc_idx <= t-1);
            if numel(acc_win) >= max(10, 2*d)
                ids = acc_win;
            else
                ids = max(2, t-win_pc):(t-1);
            end
            Cpc = cov(xsto(:, ids)');
            Cpc = (Cpc + Cpc')/2 + 1e-12*eye(d);
            [V2,D2] = eig(Cpc);
            [eigs_sorted, order] = sort(diag(D2), 'descend');
            k1 = order(1); k2 = order(min(2,numel(order)));
            v1 = V2(:,k1); v2 = V2(:,k2);
            s1 = sqrt(max(eigs_sorted(1), 1e-12));
            s2 = sqrt(max(eigs_sorted(min(2,numel(eigs_sorted))), 1e-12));
            c1 = step_mult * s1 * randn;
            c2 = step_mult * s2 * randn;
            Y = (X + c1*v1 + c2*v2)';
        
        elseif u < p_de + p_1d + p_pc2 + p_global
            % ----- Global reset (new)
            r0 = acc_idx(randi(numel(acc_idx)));
            inflate = 3.0;
            Y = mvnrnd(xsto(:,r0)', inflate^2 * baseCov);
        elseif u < p_de + p_1d + p_pc2 + p_global + p_ind
            mu = ind_opts.mu(:);           % d×1
            C  = ind_opts.C;
            L  = chol(C + 1e-12*eye(d), 'lower');
        
            % pick a random scale for this proposal
            s = ind_opts.scales(randi(numel(ind_opts.scales)));
            ycol = mu + sqrt(s) * (L * randn(d,1));  % d×1
            Y    = ycol.';                           % 1×d
        
            % Hastings correction for mixture-of-Gaussians independence kernel
            % We approximate by the chosen component (good in practice).
            vX = (L \ (X - mu));   mX = (vX' * vX) / s;
            vY = (L \ (ycol - mu)); mY = (vY' * vY) / s;
            hastings = -0.5*(mX - mY);

        elseif u < p_de + p_1d + p_pc2 + p_global + p_ind + p_t
            nu = ind_opts.nu;
            L  = chol(sdmatrix_cached + 1e-12*eye(d), 'lower');
        
            g   = randn(d,1);
            chi = chi2rnd(nu);
            tstep = L * (g / sqrt(chi/nu)); % d×1
        
            ycol = X + tstep;               % d×1
            Y    = ycol.';                  % 1×d
            % symmetric in t-metric -> hastings stays 0
        elseif u < p_de + p_1d + p_pc2 + p_global + p_ind + p_t + p_ridge
            step = ridge_scale * randn;     % 1-D Gaussian along ridge
            ycol = X + step * vridge;       % d×1
            Y    = ycol.';                  % 1×d
            % symmetric -> hastings = 0 (already default)
        else
            % ----- Joint Gaussian: mostly adaptive, sometimes global to keep mobility
            if rand < b
                Y = mvnrnd(X', baseCov_scaled);           % global/spherical-ish
            else
                Y = mvnrnd(X', sdmatrix_cached);          % local adaptive
            end
        end
    end

    % enforce fixed coords & compute FY if not already done
    if ~isempty(inds), Y(inds) = vals; end
    if isempty(FY)
        FY = F(Y);
        if ~isfinite(FY) || ~isreal(FY), FY = -Inf; end
    end
    history(1:d,t) = Y;
    
    dlog = (FY - FX) + hastings;
    if dlog>=0 || log(rand)<dlog %rand < exp(FY-FX)
        % Accept
        xsel = Y(:);
        FX = FY;
        acc = acc+1;
        history(end,t) = 1;
    else
        % Reject
        xsel = xsto(:,t-1);
    end
    xsto(:,t) = xsel;
    outsto(t) = FX;
    xbar(:,t) = (xbar(:,t-1)*(t-1) + xsel)/t;
    % per-step step-size adaptation, gated by freeze
    % --- step-size adaptation using a moving window ---
    if ~adaptation_frozen && t >= adapt_start
        W   = 150;                                    % window length
        w0  = min(W, t-1);
        aW  = mean(history(end, (t-w0):t-1));         % recent acceptance
        old = log_scale;
        eta = 0.08;                                   % learning rate
        log_scale = log_scale + eta * (aW - target_acc);
        log_scale = min(max(log_scale, -1.2), 08);
        % rescale proposals in-place
        sdmatrix_cached = sdmatrix_cached * exp(log_scale - old);
        baseCov_scaled  = baseCov         * exp(log_scale - old);
    end
    % auto-freeze adaptation
    if auto_freeze && ~adaptation_frozen && t >= adapt_start + 400 && mod(t, freeze_check_every) == 0
        acc_recent = mean(history(end, max(2,t-199):t));
        if acc_recent >= freeze_band(1) && acc_recent <= freeze_band(2)
            stable_hits = stable_hits + 1;
        else
            stable_hits = 0;
        end
        if stable_hits >= stable_hits_needed
            adaptation_frozen = true;
            fprintf('*** Freezing step-size at t=%d: acc200=%.2f, scale=%.3f ***\n', ...
                    t, acc_recent, exp(log_scale));
        end
    end

    if mod(t, ind_opts.refresh_every)==0
        ind0 = max(2, t-800); ind1 = t-1;
        ind_opts.mu = mean(xsto(:,ind0:ind1), 2).';
        ind_opts.C  = cov(xsto(:,ind0:ind1)') + 1e-6*eye(d);
        ind_opts.C  = 2.0^2 * ind_opts.C;   % inflate
    end

        % Display options
    if displ && (mod(t,round(n/25))==0); fprintf('%0.5g ', t/n*25); end
    if displ && (mod(t,200)==0)
        plot(xsto(:,1:t-1)'); xlim([0 n]); drawnow;
    end

    if mod(t,200)==0
      acc200 = mean(history(end, max(2,t-199):t));
      fprintf('t=%d  acc200=%.2f  scale=%.3f\n', t, acc200, exp(log_scale));
    end

end

accept_rate = acc/n;
xsto = xsto';
history = history';
end