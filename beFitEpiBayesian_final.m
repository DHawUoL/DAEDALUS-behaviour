function [chains]=beFitEpiBayesian_final(ydata,X,data,Xfull,coeff,Diag)
filename="chains_blocked1";
hlag=0;
plotRun=0;%Plot simulation at x0 - causes an error so a fit doesn't go ahead

%% RUN %%

%Stringency:
tvec=[1,2,61,94,[127,134,141,148,153,155,162,167,169,175,176,186,200,211,216,218,223,227,230,236,250,258,265,266,271,279,288,294,305,310,322,330,337,338,349,354,356,361,370,372,384,397,407,418,433,445,454,463,468,474,491,503,504,517,540,561,566,567,575]+hlag];
xdata=85:tvec(end-15);%-2 -7
lt=length(tvec);
X=X(:,1:lt-1);
Xfull=Xfull(:,1:lt-1);
[lx1,lx2]=size(X);
ydata=ydata(0+(1:length(xdata)));
%If data is just England:
ydata=ydata*(sum(data.Npop)/56286961);%England, mid-2019 (ONS)

% From diagnostics / last frequentist fit:
% ---- 4-parameter model: theta = [alpha, phi, kstar, k3] ----
lb   = [0,  -20, -20];
ub   = [1,   20,  20];
plim = [ub; lb];   % 2 x 4

% From diagnostics / last frequentist fit (length must be 4)
lsq_hat = Diag.p_hat(:).';     % [alpha_hat, phi_hat, kstar_hat, k3_hat]
lsq_se  = Diag.se(:).';        % [se_alpha, se_phi,  se_kstar,  se_k3]
% Beta priors for alpha, phi (map mean->(a,b))
to_beta = @(m,K) deal(max(m*K,1e-6), max((1-m)*K,1e-6));  % safe clamp
[a_alpha,b_alpha] = to_beta(min(max(lsq_hat(1),1e-6),1-1e-6), 80);
[a_phi,  b_phi]   = to_beta(min(max(lsq_hat(2),1e-6),1-1e-6), 150);
% Start points (nChains x 4)
nChains = 4;
jitter  = [0.02, 0.3, 0.2];           % per-dimension jitter
startPoints = repmat(lsq_hat, nChains, 1) + randn(nChains,3).*jitter;
for j=1:3
    startPoints(:,j) = min(max(startPoints(:,j), plim(2,j)+1e-6), plim(1,j)-1e-6);
end

%% RUN MCMC IN PARALLEL %%
%Make sure Parallel Computing Toolbox is available
%
if isempty(gcp('nocreate'))
    parpool('local', nChains);
end
%}
chains = cell(nChains,1);
mu_z = x2z(lsq_hat, plim);             % centre in z-space
ind_opts.mu    = mu_z(:).';            % 1x4
ind_opts.C     = (3.5^2) * eye(numel(mu_z));   % fairly broad
ind_opts.p_ind = 0.20;
ind_opts.p_t   = 0.10;
ind_opts.nu    = 4;
ind_opts.refresh_every = 500;
ind_opts.scales = [1, 4, 16];

for c = 1:nChains
    chainData=struct;
    %% Blocks:
    blockList = { [1,2,3]};%[2,4] };         % 4 parameters
    base   = 1;
    sigmas = [1, 1, 1];              % z-space scales (length 4)
    nPerBlock = 3000;
    nCycles   = 4;
    
    theta = startPoints(c,:);
    z0    = x2z(theta, plim);
    sigma0 = [1.2, 1.4, 1.4];              % 4-dim

    z_all = []; acc_rates = [];
    flats = [0.6 0.85 1.0];
    for i = 1:length(flats)
        Fz = @(z) posteriorz(z, flats(i), data, xdata, ydata, X, Xfull, coeff, ...
                     tvec, lx1, lx2, plim, plotRun, lsq_hat, lsq_se, ...
                     a_alpha,b_alpha, a_phi,b_phi);
    
        [z_alli, z_final, acc_ratesi, cov_last] = blockMCMC_wrapper(Fz, z0, sigma0, nPerBlock, nCycles, blockList, false, ind_opts);
    
        z_all     = [z_all; z_alli];
        acc_rates = [acc_rates; acc_ratesi];
    
        % update starting point and carry forward covariance (as a matrix)
        %z0     = z_final;
        %sigma0 = cov_last;  % this is d×d; MCMC_adaptive will treat it as "matrix" mode
        z0     = z_final + 0.05*randn(size(z_final));     % small re-centering jitter
        sigma0 = cov_last * 1.5 + 1e-6*eye(numel(z0));    % modest inflation
    end
    chainData.xsto=mapSample(z_all,plim);
    chainData.x_final=theta;
    chainData.acc_rates=acc_rates;
    %}
    chains{c}=chainData;
end
%Save results
save(filename, 'chains', 'startPoints');

end

%% SIMULATION %%

function [f,rhohat]=sim2fit(params,data,xdata,Xfit,Xfull,coeff,tvec,lx1,lx2,plotRun)
    % params = [alpha, phi, kstar, k3]
    R0 = 2.2;
    tvec(1) = -145;

    alpha = params([1,1,1]);   % broadcast into 3-group alpha vector if needed
    phi   = 1;%params(2);
    ks    = params(2);
    k3    = params(3);

    % --- linear collapse for k1,k2 from k* ---
    a1=-.814;   b1= 8.0161;
    a2=-.8067;  b2=-8.0887;
    k1 = a1*ks + b1;
    k2 = a2*ks + b2;

    % reducedParams = [1, k1, k2, k3, v0]; here v0=0 after removal
    reducedParams = [1, k1, k2, k3, 0];

    [pr,be,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta] = ...
        bePrepCovid19(data, R0, ones(1,lx2-2), reducedParams, coeff, zeros(5,lx2), alpha);

    pr.leak  = phi;      % implements (1 - phi*p)^2 in your transmission code
    pr.xfull = Xfull;

    % "intrinsic" fit: all-ones intervention input
    [simu,simu2,~,rhohat] = beRunCovid19(pr,be,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta, ...
                                         ones(lx1, length(tvec)-1), tvec(1:lx2+1), plotRun, data);

    t = simu(:,1)';
    h = simu2';                         % admissions
    f = interp1(t,h,xdata);
end


%% BLOCK WRAPPER %%

function [all_xsto, final_x, accept_rates, last_cov] = blockMCMC_wrapper(F, x0, sigma, nPerBlock, nCycles, blockList, displ, ind_opts)
    x = x0; d = length(x0);
    all_xsto = []; accept_rates = zeros(nCycles, length(blockList));
    last_cov = [];  % pass last learned Cb back to caller

    for cycle = 1:nCycles
        fprintf('\nCycle %d / %d\n', cycle, nCycles);
        for b = 1:length(blockList)
            block = blockList{b};
            fixed = setdiff(1:d, block);
            fixinds = [fixed; x(fixed)];

            % --- promote sigma to full-dim form ---
            if isscalar(sigma)
                sigma_b = sigma;  % OK
            elseif isvector(sigma)
                if numel(sigma)==d
                    sigma_b = sigma(:)';                        % length-d
                elseif numel(sigma)==numel(block)
                    sigma_b = zeros(1,d);                       % map into full length
                    sigma_b(block) = sigma(:)';
                    sigma_b(fixed) = 0;
                else
                    error('sigma vector must have length d or length(block)');
                end
            elseif ismatrix(sigma)
                if all(size(sigma)==[d d])
                    sigma_b = sigma;                            % full matrix
                    sigma_b(fixed,:) = 0; sigma_b(:,fixed) = 0;
                elseif all(size(sigma)==[numel(block) numel(block)])
                    sigma_b = zeros(d);                         % embed block matrix
                    sigma_b(block,block) = sigma;
                else
                    error('sigma matrix must be d×d or block×block');
                end
            else
                error('Unsupported sigma type');
            end

            fprintf('  Block %d: sampling params [%s]\n', b, num2str(block));
            [xsto,~,~,accept_rate,cov_b] = MCMC_adaptive(F, x, nPerBlock, sigma_b, fixinds, numel(block), displ, ind_opts);

            x = xsto(end, :);
            all_xsto = [all_xsto; xsto];
            accept_rates(cycle, b) = accept_rate;
            last_cov = cov_b;   % keep the most recent adaptive covariance
        end
    end
    final_x = x;
end


%% PRIOR %%

function lp = log_prior_4(theta, lsq_hat, lsq_se, a_alpha,b_alpha, a_phi,b_phi, plim)
    % theta = [alpha, phi, kstar, k3]
    alpha = theta(1); phi = theta(2);
    kstar = theta(3); k3  = theta(4);

    % Guard bounds (alpha & phi in (0,1); others inside plim)
    if ~(plim(2,1) < alpha && alpha < plim(1,1)), lp = -Inf; return; end
    if ~(0 < phi && phi < 1),                      lp = -Inf; return; end
    if ~(plim(2,3) < kstar && kstar < plim(1,3)),  lp = -Inf; return; end
    if ~(plim(2,4) < k3    && k3    < plim(1,4)),  lp = -Inf; return; end

    % Beta(alpha on [0,1]), Beta(phi), Gaussians for kstar,k3
    % alpha has bounds [lb,ub] but here lb=0, ub=1 so mapping is identity
    u = (alpha - plim(2,1)) / (plim(1,1)-plim(2,1));   % in (0,1)
    lp = (a_alpha-1)*log(u) + (b_alpha-1)*log(1-u) - betaln(a_alpha,b_alpha) ...
       - log(plim(1,1)-plim(2,1));

    lp = lp + (a_phi-1)*log(phi) + (b_phi-1)*log(1-phi) - betaln(a_phi,b_phi);

    % weakly-informative normals centred at LSQ
    lp = lp - 0.5*((kstar - lsq_hat(3))/lsq_se(3))^2 - log(lsq_se(3)*sqrt(2*pi));
    lp = lp - 0.5*((k3    - lsq_hat(4))/lsq_se(4))^2 - log(lsq_se(4)*sqrt(2*pi));
end

function lp = log_prior_3(theta, lsq_hat, lsq_se, a_alpha,b_alpha, a_phi,b_phi, plim)
    % theta = [alpha, phi, kstar, k3]
    alpha = theta(1); 
    kstar = theta(2); k3  = theta(3);

    % Guard bounds (alpha & phi in (0,1); others inside plim)
    if ~(plim(2,1) < alpha && alpha < plim(1,1)), lp = -Inf; return; end
    if ~(plim(2,2) < kstar && kstar < plim(1,2)),  lp = -Inf; return; end
    if ~(plim(2,3) < k3    && k3    < plim(1,3)),  lp = -Inf; return; end

    % Beta(alpha on [0,1]), Beta(phi), Gaussians for kstar,k3
    % alpha has bounds [lb,ub] but here lb=0, ub=1 so mapping is identity
    u = (alpha - plim(2,1)) / (plim(1,1)-plim(2,1));   % in (0,1)
    lp = (a_alpha-1)*log(u) + (b_alpha-1)*log(1-u) - betaln(a_alpha,b_alpha) ...
       - log(plim(1,1)-plim(2,1));

    % weakly-informative normals centred at LSQ
    lp = lp - 0.5*((kstar - lsq_hat(2))/lsq_se(2))^2 - log(lsq_se(2)*sqrt(2*pi));
    lp = lp - 0.5*((k3    - lsq_hat(3))/lsq_se(3))^2 - log(lsq_se(3)*sqrt(2*pi));
end

%% LIKELIHOOD %%
%{
function ll = log_like_4(theta, data, xdata, ydata, Xfit, Xfull, coeff, tvec, lx1, lx2, plotRun)
    ymodel = sim2fit(theta, data, xdata, Xfit, Xfull, coeff, tvec, lx1, lx2, plotRun);
    if any(~isfinite(ymodel)), ll = -Inf; return; end
    ymodel = max(ymodel, 1e-9);  % strictly positive Poisson mean
    % Poisson log-likelihood (up to constants)
    ll = sum( -ymodel' + ydata .* log(ymodel') - gammaln(ydata+1) );
end
%}
function ll = log_like_4(theta, data, xdata, y, Xfit, Xfull, coeff, tvec, lx1, lx2, plotRun)
    mu = sim2fit(theta, data, xdata, Xfit, Xfull, coeff, tvec, lx1, lx2, plotRun);
    mu = mu(:); y = y(:);
    if any(~isfinite(mu)), ll = -Inf; return; end
    mu = max(mu, 1e-8);
    r  = 50;  % try 30–200; lower = more overdispersion
    % y ~ NB(r, p=r/(r+mu))
    ll = sum( gammaln(y+r) - gammaln(r) - gammaln(y+1) ...
            + r.*log(r) + y.*log(mu) - (y+r).*log(r+mu) );
end

%% Variable transform x/z

function [z, logJ] = x2z(x, plim)
    x = x(:)'; d = numel(x);
    z = zeros(1,d);
    logJ = 0;
    for j = 1:d
        lb = plim(2,j); ub = plim(1,j); w = ub - lb;
        % keep strictly inside to avoid log(0)
        xj = min(max(x(j), lb+1e-12), ub-1e-12);
        u  = (xj - lb)/w;           % in (0,1)
        z(j) = log(u/(1-u));        % logit
        logJ = logJ + (log(w) + log(u) + log(1-u));  % |dz/dx| inverse = w*u*(1-u)
    end
end

function [x, logJ] = z2x(z, plim)
    z = z(:)'; d = numel(z);
    x = zeros(1,d);
    logJ = 0;
    for j = 1:d
        lb = plim(2,j); ub = plim(1,j); w = ub - lb;
        s = 1./(1+exp(-z(j)));             % (0,1)
        s = min(max(s,1e-12),1-1e-12);     % clamp for safety
        x(j) = lb + w*s;
        logJ = logJ + (log(w) + log(s) + log(1-s)); % sum of per-dim Jacobians
    end
end

%% Total posterior on z:

function f = posteriorx(x, flat, data, xdata, ydata, Xfit, Xfull, coeff, ...
                        tvec, lx1, lx2, plim, plotRun, lsq_hat, lsq_se, ...
                        a_alpha,b_alpha, a_phi,b_phi)
    ll = log_like_4(x, data, xdata, ydata, Xfit, Xfull, coeff, tvec, lx1, lx2, plotRun);
    if ~isfinite(ll), f = -Inf; return; end
    lp = log_prior_4(x, lsq_hat, lsq_se, a_alpha,b_alpha, a_phi,b_phi, plim);
    f  = flat*ll + lp;
end

function f = posteriorz(z, flat, data, xdata, ydata, Xfit, Xfull, coeff, ...
                        tvec, lx1, lx2, plim, plotRun, lsq_hat, lsq_se, ...
                        a_alpha,b_alpha, a_phi,b_phi)
    [x, logJ] = z2x(z, plim);
    ll = log_like_4(x, data, xdata, ydata, Xfit, Xfull, coeff, tvec, lx1, lx2, plotRun);
    if ~isfinite(ll), f = -Inf; return; end
    lp = log_prior_3(x, lsq_hat, lsq_se, a_alpha,b_alpha, a_phi,b_phi, plim);
    f  = flat*ll + lp + logJ;   % add Jacobian for transform
end

function xsto = mapSample(zsto, plim)
    N = size(zsto,1);
    xsto = zeros(N, size(plim,2));
    for i = 1:N
        [xi, ~] = z2x(zsto(i,:), plim);
        xsto(i,:) = xi;
    end
end
