function [chains] = beFitEpiBayesian_final_2(ydata,X,data,Xfull,coeff,Diag)
% 4-parameter Bayesian fit:
% theta = [alpha, kstar, k3, v0]
% - Tight Normal prior on v0 (effectively “fixed but with uncertainty”)
% - Negative Binomial likelihood for admissions
% - Blocked, adaptive MH in z-space with logistic box transform

filename = "chains_blocked_ak3v0";
hlag     = 0;
plotRun  = 0;

%% ---- Timeline and data slice ----
tvec  = [1,2,61,91,[127,134,141,148,153,155,162,167,169,175,176,186,200,211,216,218,223,227,230,236,...
                    250,258,265,266,271,279,288,294,305,310,322,330,337,338,349,354,356,361,370,372,...
                    384,397,407,418,433,445,454,463,468,474,491,503,504,517,540,561,566,567,575]+hlag];
xdata = 85:tvec(end-15);

lt    = numel(tvec);
X     = X(:,1:lt-1);
Xfull = Xfull(:,1:lt-1);
[~,lx2] = size(X);

ydata = ydata(1:numel(xdata));
ydata = ydata*(sum(data.Npop)/56286961);  % scale England-only → UK

%% ---- Parameter box (plim = [ub; lb]) ----
% alpha in [0,1], k* and k3 fairly wide, v0 in a *tight* box around v0_fix
if isfield(Diag,'v0_hat'), v0_fix = Diag.v0_hat; else, v0_fix = 4.3792; end
v0_box = 0.75;   % half-width of the *box* for v0 (still enforced), SD set below

lb   = [0,   -20, -20, v0_fix - v0_box];
ub   = [1,    20,  20, v0_fix + v0_box];
plim = [ub; lb];

%% ---- LSQ anchors (used for weak normals on k*, k3) ----
% If not provided, default to zero-mean, broad SEs.
if ~isfield(Diag,'p_hat') || numel(Diag.p_hat)<4
    lsq_hat = [0.4, 0, 0, v0_fix];
else
    lsq_hat = Diag.p_hat(:).';
end
if ~isfield(Diag,'se') || numel(Diag.se)<4
    lsq_se  = [0.1, 8, 8, 0.3];
else
    lsq_se  = max(Diag.se(:).', 1e-6);
end

%% ---- Priors: alpha ~ Beta, k*,k3 ~ Normal(LSQ, SE), v0 ~ Normal(v0_fix, tau^2 tight) ----
% Beta for alpha: center near LSQ alpha with reasonable concentration
to_beta = @(m,K) deal(max(m*K,1e-6), max((1-m)*K,1e-6));
a_alpha = 2; b_alpha = 2; % default symmetric if LSQ alpha is silly
if lsq_hat(1)>0 && lsq_hat(1)<1
    [a_alpha,b_alpha] = to_beta(lsq_hat(1), 80); % quite informative around LSQ alpha
end
tau_v0 = 0.35;   % << tighten/loosen here; 0.25–0.5 works well

%% ---- Chain starts (in x-space) ----
nChains = 4;
startPoints = zeros(nChains,4);
jitter = [0.02, 0.8, 0.8, 0.05];
for c=1:nChains
    x0 = lsq_hat + jitter.*randn(1,4);
    x0(1) = min(max(x0(1), lb(1)+1e-6), ub(1)-1e-6);
    x0(2) = min(max(x0(2), lb(2)+1e-3), ub(2)-1e-3);
    x0(3) = min(max(x0(3), lb(3)+1e-3), ub(3)-1e-3);
    x0(4) = min(max(x0(4), lb(4)+1e-6), ub(4)-1e-6);
    startPoints(c,:) = x0;
end

%% ---- Independence/Mixture proposal options (in z-space) ----
mu_z = x2z(lsq_hat, plim);
ind_opts = struct;
ind_opts.mu    = mu_z(:).';
ind_opts.C     = (3.5^2)*eye(numel(mu_z));
ind_opts.p_ind = 0.20;
ind_opts.p_t   = 0.10;
ind_opts.nu    = 4;
ind_opts.refresh_every = 500;
ind_opts.scales = [1, 4, 16];

%% ---- Run blocked adaptive MCMC in z-space ----
chains = cell(nChains,1);
for c=1:nChains
    z0 = x2z(startPoints(c,:), plim);
    blockList = { [1,2,3,4] };   % single block over all 4 params
    nPerBlock = 3000;
    nCycles   = 4;
    sigma0    = [1.2, 1.4, 1.4, 0.6]; % per-dim z-scales to start
    flats     = [0.70, 0.90, 1.00];

    z_all = []; acc_rates = []; cov_last = eye(4);
    for i = 1:numel(flats)
        Fz = @(z) posteriorz(z, flats(i), data, xdata, ydata, X, Xfull, coeff, ...
                             tvec, lx2, plim, plotRun, lsq_hat, lsq_se, ...
                             a_alpha,b_alpha, v0_fix,tau_v0);
        [z_block, z_final, acc_blk, cov_last] = blockMCMC_wrapper(Fz, z0, sigma0, nPerBlock, nCycles, blockList, false, ind_opts);
        z_all     = [z_all; z_block];
        acc_rates = [acc_rates; acc_blk];
        z0        = z_final + 0.05*randn(size(z_final));
        sigma0    = cov_last*1.4 + 1e-6*eye(4);
    end

    chainData = struct;
    chainData.xsto       = mapSample(z_all, plim);
    chainData.x_final    = startPoints(c,:);
    chainData.acc_rates  = acc_rates;
    chains{c} = chainData;
end

save(filename, 'chains', 'startPoints');

end % beFitEpiBayesian_final


%% ====================== Likelihood & Prior =========================
function ll = log_like(theta, data, xdata, y, Xfit, Xfull, coeff, tvec, lx2, plotRun)
    mu = sim2fit(theta, data, xdata, Xfit, 1, Xfull, coeff, tvec, lx2, plotRun, 0);
    mu = mu(:); y = y(:);
    if any(~isfinite(mu)), ll = -Inf; return; end
    mu = max(mu, 1e-8);
    % Negative Binomial for overdispersed counts:
    r  = 50;  % 30–200 are typical; tune if needed
    ll = sum( gammaln(y+r) - gammaln(r) - gammaln(y+1) ...
            + r.*log(r) + y.*log(mu) - (y+r).*log(r+mu) );
end

function lp = log_prior_ak3v0(theta, lsq_hat, lsq_se, a_alpha,b_alpha, v0_fix,tau_v0, plim)
    % theta = [alpha, kstar, k3, v0]
    a = theta(1); ks = theta(2); k3 = theta(3); v0 = theta(4);
    % box
    if ~(plim(2,1) < a  && a  < plim(1,1)), lp = -Inf; return; end
    if ~(plim(2,2) < ks && ks < plim(1,2)), lp = -Inf; return; end
    if ~(plim(2,3) < k3 && k3 < plim(1,3)), lp = -Inf; return; end
    if ~(plim(2,4) < v0 && v0 < plim(1,4)), lp = -Inf; return; end

    % alpha ~ Beta(a_alpha,b_alpha) over [0,1]
    if a<=0 || a>=1, lp = -Inf; return; end
    lp = (a_alpha-1)*log(a) + (b_alpha-1)*log(1-a) - betaln(a_alpha,b_alpha);

    % k*, k3 ~ Normal(LSQ, SE^2) (weakly informative)
    lp = lp - 0.5*((ks - lsq_hat(2))/lsq_se(2))^2 - log(lsq_se(2)*sqrt(2*pi));
    lp = lp - 0.5*((k3 - lsq_hat(3))/lsq_se(3))^2 - log(lsq_se(3)*sqrt(2*pi));

    % v0 ~ Normal(v0_fix, tau_v0^2)  (tight)
    lp = lp - 0.5*((v0 - v0_fix)/tau_v0)^2 - log(tau_v0*sqrt(2*pi));
end

function f = posteriorx(x, flat, data, xdata, ydata, Xfit, Xfull, coeff, ...
                        tvec, lx2, plim, plotRun, lsq_hat, lsq_se, ...
                        a_alpha,b_alpha, v0_fix,tau_v0)
    ll = log_like(x, data, xdata, ydata, Xfit, Xfull, coeff, tvec, lx2, plotRun);
    if ~isfinite(ll), f = -Inf; return; end
    lp = log_prior_ak3v0(x, lsq_hat, lsq_se, a_alpha,b_alpha, v0_fix,tau_v0, plim);
    f  = flat*ll + lp;
end

function f = posteriorz(z, flat, data, xdata, ydata, Xfit, Xfull, coeff, ...
                        tvec, lx2, plim, plotRun, lsq_hat, lsq_se, ...
                        a_alpha,b_alpha, v0_fix,tau_v0)
    [x, logJ] = z2x(z, plim);
    ll = log_like(x, data, xdata, ydata, Xfit, Xfull, coeff, tvec, lx2, plotRun);
    if ~isfinite(ll), f = -Inf; return; end
    lp = log_prior_ak3v0(x, lsq_hat, lsq_se, a_alpha,b_alpha, v0_fix,tau_v0, plim);
    f  = flat*ll + lp + logJ;   % Jacobian for the box→R transform
end


%% ====================== Simulation wrapper ========================
function [f,rhohat]=sim2fit(params,data,xdata,Xfit,intrinsic,Xfull,coeff,tvec,lx2,plotRun,ymean)
    R0=2.8;
    tvec(1)=-80;

    alpha=params([1,1,1]);
    propIn=1;
    reducedParams=[1,params(2),0,params(3:4)];

    % Prep (this gives you NNbar)
    [pr,be,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta] = ...
        bePrepCovid19(data,R0,ones(1,lx2-2),reducedParams,coeff,zeros(5,lx2),alpha,propIn);

    pr.leak=0; 
    pr.xfull = Xfull;     % <-- your PCs / indices go here (as before)
    be.BiFirstFit=1; 
    pr.phi2=0;
    pr.ymean=ymean;

    % === IMPORTANT: build Xit with the correct number of sector rows ===
    % In beRunCovid19, they do: NNvec = repmat(NNbar(1:lx),1,lt-1) .* Xit;
    % NNbar has length (lx + lc). lc=4 in beRunCovid19. So sectors = length(NNbar)-4.
    lc = 4;
    lx_sect = numel(NNbar) - lc;      % inferred #sectors used inside beRunCovid19
    Xit_mat = ones(lx_sect, lx2);     % “no closures” but with correct shape

    % (Optional sanity: if your Xfull was accidentally transposed somewhere)
    % if size(Xfull,2) ~= lx2 && size(Xfull,1) == lx2
    %     pr.xfull = Xfull.';  % keep PCs as [nDrivers × (lt-1)]
    % end

    if intrinsic==1
        [simu,simu2,~,rhohat] = beRunCovid19( ...
            pr,be,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta, ...
            Xit_mat, tvec(1:lx2+1), plotRun, data);
    else
        Wfit = Xfit.^(1/pr.a);
        [simu,simu2,~,rhohat] = beRunCovid19( ...
            pr,be,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta, ...
            Wfit, tvec(1:lx2+1), plotRun, data);
    end

    t = simu(:,1)'; 
    h = simu2';
    f = interp1(t,h,xdata,'linear'); 
    f(~isfinite(f)) = 1e9;
end



%% ====================== MCMC helpers ==============================
function [all_xsto, final_x, accept_rates, last_cov] = blockMCMC_wrapper(F, x0, sigma, nPerBlock, nCycles, blockList, displ, ind_opts)
    x = x0; d = length(x0);
    all_xsto = []; accept_rates = zeros(nCycles, numel(blockList));
    last_cov = eye(d);
    for cycle = 1:nCycles
        fprintf('\nCycle %d / %d\n', cycle, nCycles);
        for b = 1:numel(blockList)
            block  = blockList{b};
            fixed  = setdiff(1:d, block);
            fixinds = [fixed; x(fixed)];
            sigma_b = promote_sigma(sigma, d, block);
            fprintf('  Block %d: [%s]\n', b, num2str(block));
            [xsto,~,~,acc_b,cov_b] = MCMC_adaptive(F, x, nPerBlock, sigma_b, fixinds, numel(block), displ, ind_opts);
            x = xsto(end,:);
            all_xsto = [all_xsto; xsto];
            accept_rates(cycle,b) = acc_b;
            last_cov = cov_b;
        end
    end
    final_x = x;
end

function sigma_b = promote_sigma(sigma, d, block)
    if isscalar(sigma),             sigma_b = sigma;
    elseif isvector(sigma)
        if numel(sigma)==d
            sigma_b = sigma(:).';
        elseif numel(sigma)==numel(block)
            sigma_b = zeros(1,d); sigma_b(block) = sigma(:).';
        else, error('sigma length mismatch');
        end
    elseif ismatrix(sigma)
        if all(size(sigma)==[d d])
            sigma_b = sigma;
        elseif all(size(sigma)==[numel(block) numel(block)])
            sigma_b = zeros(d); sigma_b(block,block) = sigma;
        else, error('sigma size mismatch');
        end
    else, error('Unsupported sigma type');
    end
end

% Adaptive MH (Haario) — unchanged logic from your working version
function [xsto, outsto, history, accept_rate, covmat] = MCMC_adaptive(F, x0, n, sigma, fixinds, blockind, displ, ind_opts)
    if nargin < 8 || isempty(ind_opts), ind_opts = struct; end
    if ~isfield(ind_opts,'mu'), ind_opts.mu = x0(:)'; end
    if ~isfield(ind_opts,'C'),  ind_opts.C  = eye(length(x0)); end
    if ~isfield(ind_opts,'p_ind'), ind_opts.p_ind = 0.00; end
    if ~isfield(ind_opts,'p_t'),   ind_opts.p_t   = 0.00; end
    if ~isfield(ind_opts,'nu'),    ind_opts.nu    = 5;    end
    if ~isfield(ind_opts,'refresh_every'), ind_opts.refresh_every = 500; end
    if ~isfield(ind_opts,'scales'), ind_opts.scales = [1,4,16]; end

    global_scale = 0.75;
    d = length(x0); b = 0.15;
    target_acc = 0.24; log_scale = -0.05;

    nfree = d - size(fixinds,2);
    adapt_start = max(300, 4*max(1,nfree));
    window = 1000; shrink = 0.03; jitter = 1e-12;

    if isscalar(sigma), sigmode="scalar";
    elseif isvector(sigma) && numel(sigma)==d, sigmode="vector"; sigma = sigma(:).';
    elseif ismatrix(sigma) && all(size(sigma)==[d d]), sigmode="matrix"; sigma = (sigma+sigma.')/2 + 1e-12*eye(d);
    else, error('sigma must be scalar, length-%d vector, or %d×%d matrix', d,d,d);
    end

    if ~isempty(fixinds), inds = fixinds(1,:); vals = fixinds(2,:);
    else, inds = []; vals = [];
    end

    cov0 = eye(d);
    if ~isempty(inds), cov0(inds,:) = 0; cov0(:,inds)=0; end
    covmat = cov0;

    xsto = zeros(n,d); outsto = zeros(1,n); history = zeros(n,d+1);
    xsto(1,:) = x0(:).'; xbar = xsto(1,:);
    FX = F(x0); outsto(1) = FX; acc = 0;

    switch sigmode
        case "scalar", baseCov = (global_scale^2/d) * cov0;
        case "vector", S = diag(sigma); baseCov = (global_scale^2/d) * (S*cov0*S);
        case "matrix", baseCov = (global_scale^2/d) * sigma;
    end
    baseCov_scaled  = baseCov * exp(log_scale);
    sdmatrix_cached = baseCov_scaled;

    p_de_base=0.40; p_1d_base=0.15; p_pc2_base=0.15; p_global_base=0.07;
    p_ridge=0.05; ridge_scale=0.25; vridge = ones(d,1); vridge = vridge/norm(vridge);

    for t=2:n
        X = xsto(t-1,:);
        FY = []; hastings=0;

        if t < adapt_start
            Y = mvnrnd(X, baseCov);
        else
            if isfield(ind_opts,'refresh_every') && mod(t, ind_opts.refresh_every)==0
                W = 1000; i1 = max(2, t-W);
                Zwin = xsto(i1:t-1,:).';
                mu_ind = mean(Zwin,2)'; C_ind = cov(Zwin'); C_ind=(C_ind+C_ind')/2;
                C_ind = 2.0^2*(C_ind + 1e-9*eye(d));
                ind_opts.mu = mu_ind; ind_opts.C = C_ind;
            end
            if mod(t,10)==0
                ind1 = t-1; ind0 = max(2, ind1 - window + 1);
                Cb = cov(xsto(ind0:ind1,:)); if ~isempty(inds), Cb(inds,:)=0; Cb(:,inds)=0; end
                mtrace = mean(diag(Cb));
                Cb = (1-shrink)*Cb + shrink*mtrace*eye(d);
                Cb = (Cb+Cb')/2 + jitter*eye(d);
                sdmatrix_new = (2.4^2/d)*exp(log_scale)*Cb + 1e-12*eye(d);
                sdmatrix_cached = 0.8*sdmatrix_cached + 0.2*sdmatrix_new;
                covmat = Cb;
            end
            p_ind = ind_opts.p_ind; p_t = ind_opts.p_t; nu = ind_opts.nu;
            psum = p_de_base+p_1d_base+p_pc2_base+p_global_base+p_ind+p_t+p_ridge;
            scale = min(0.85/max(psum,eps),1);
            p_de=p_de_base*scale; p_1d=p_1d_base*scale; p_pc2=p_pc2_base*scale;
            p_global=p_global_base*scale; p_ind=p_ind*scale; p_t=p_t*scale; p_ridge=p_ridge*scale;
            p_joint = 1 - (p_de+p_1d+p_pc2+p_global+p_ind+p_t+p_ridge);

            u=rand;
            if u < p_de
                acc_idx = find(history(1:t-1,end)==1); if numel(acc_idx)<2, Y = mvnrnd(X, sdmatrix_cached);
                else
                    r1 = acc_idx(randi(numel(acc_idx))); r2 = acc_idx(randi(numel(acc_idx)));
                    while r2==r1, r2 = acc_idx(randi(numel(acc_idx))); end
                    gamma0 = 2.38/sqrt(2*d); gamma_de = 0.85*gamma0;
                    diffv = (xsto(r1,:)-xsto(r2,:));
                    Cr=0.9; mask = rand(1,d)<Cr; if ~any(mask), mask(randi(d))=true; end
                    diffv = diffv.*mask;
                    jit = mvnrnd(zeros(1,d), 1e-6*(sdmatrix_cached + 1e-12*eye(d)));
                    Y = X + gamma_de*diffv + jit;
                end
            elseif u < p_de+p_1d
                ids = max(2,t-400):(t-1); Cpc = cov(xsto(ids,:)); Cpc=(Cpc+Cpc')/2 + 1e-12*eye(d);
                [V,D] = eig(Cpc); [~,ix]=max(diag(D)); v1=V(:,ix); s1=sqrt(max(D(ix,ix),1e-12));
                step = exp(log_scale)*s1*randn; Y = X + step*v1';
            elseif u < p_de+p_1d+p_pc2
                ids = max(2,t-400):(t-1); Cpc = cov(xsto(ids,:)); Cpc=(Cpc+Cpc')/2 + 1e-12*eye(d);
                [V,D] = eig(Cpc); [evals,ord]=sort(diag(D),'descend');
                k1=ord(1); k2=ord(min(2,numel(ord)));
                v1=V(:,k1); v2=V(:,k2); s1=sqrt(max(evals(1),1e-12)); s2=sqrt(max(evals(min(2,numel(evals))),1e-12));
                c1=exp(log_scale)*s1*randn; c2=exp(log_scale)*s2*randn; Y = X + c1*v1' + c2*v2';
            elseif u < p_de+p_1d+p_pc2+p_global
                inflate = 3.0; Y = mvnrnd(X, inflate^2*baseCov);
            elseif u < p_de+p_1d+p_pc2+p_global+p_ind
                mu = ind_opts.mu(:); C = ind_opts.C; L = chol(C+1e-12*eye(d),'lower');
                s = ind_opts.scales(randi(numel(ind_opts.scales)));
                ycol = mu + sqrt(s)*(L*randn(d,1)); Y = ycol';
                vX = L\(X-mu)'; vY = L\(ycol-mu); mX=(vX'*vX)/s; mY=(vY'*vY)/s; hastings = -0.5*(mX - mY);
            elseif u < p_de+p_1d+p_pc2+p_global+p_ind+p_t
                L = chol(sdmatrix_cached+1e-12*eye(d),'lower'); g = randn(d,1); chi = chi2rnd(nu);
                ycol = X + (L*(g/sqrt(chi/nu)))'; Y = ycol;
            elseif u < p_de+p_1d+p_pc2+p_global+p_ind+p_t+p_ridge
                step = ridge_scale*randn; Y = X + step*vridge';
            else
                if rand<b, Y = mvnrnd(X, baseCov_scaled); else, Y = mvnrnd(X, sdmatrix_cached); end
            end
        end

        if ~isempty(fixinds), Y(fixinds(1,:)) = fixinds(2,:); end
        if isempty(FY), FY = F(Y); if ~isfinite(FY), FY = -Inf; end, end
        history(t,1:d) = Y;

        dlog = (FY - FX) + hastings;
        if dlog>=0 || log(rand)<dlog
            xsel = Y; FX = FY; acc = acc+1; history(t,end)=1;
        else
            xsel = xsto(t-1,:);
        end
        xsto(t,:) = xsel; outsto(t) = FX; xbar = (xbar*(t-1)+xsel)/t;

        if t>=adapt_start
            W=150; w0=min(W,t-1); aW = mean(history((t-w0):(t-1),end));
            eta=0.08; old=log_scale; log_scale = min(max(log_scale + eta*(aW-0.24), -1.2), 0.8);
            sdmatrix_cached = sdmatrix_cached * exp(log_scale-old);
            baseCov_scaled  = baseCov         * exp(log_scale-old);
        end

        if mod(t,200)==0
            acc200 = mean(history(max(2,t-199):t,end));
            fprintf('t=%d  acc200=%.2f  scale=%.3f\n', t, acc200, exp(log_scale));
        end
    end

    accept_rate = acc/n;
end

%% ====================== Box transform & mapping ===================
function [z, logJ] = x2z(x, plim)
    x = x(:)'; d = numel(x); z = zeros(1,d); logJ=0;
    for j=1:d
        lb = plim(2,j); ub = plim(1,j); w = ub-lb;
        xj = min(max(x(j), lb+1e-12), ub-1e-12);
        u  = (xj - lb)/w; z(j) = log(u/(1-u));
        logJ = logJ + (log(w) + log(u) + log(1-u));
    end
end

function [x, logJ] = z2x(z, plim)
    z = z(:)'; d = numel(z); x = zeros(1,d); logJ=0;
    for j=1:d
        lb = plim(2,j); ub = plim(1,j); w = ub-lb;
        s  = 1./(1+exp(-z(j))); s = min(max(s,1e-12),1-1e-12);
        x(j) = lb + w*s;
        logJ = logJ + (log(w) + log(s) + log(1-s));
    end
end

function xsto = mapSample(zsto, plim)
    N = size(zsto,1); xsto = zeros(N, size(plim,2));
    for i=1:N, xsto(i,:) = z2x(zsto(i,:), plim); end
end
