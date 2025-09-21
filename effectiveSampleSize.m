function [ESS, tau]=effectiveSampleSize(chains)

burn=2000;

% Suppose you saved 10 chains in cell array `chains`, each with .xsto
% Make a 3-D array: [iters x params x chains]
T = size(chains{1}.xsto,1)-burn;
P = size(chains{1}.xsto,2);
C = numel(chains);

X3 = zeros(T, P, C);
for c = 1:C
    X3(:,:,c) = chains{c}.xsto(burn+1:end,:);   % ensure burn-in already removed
end

[ESS, tau] = ess_multi(X3);  % or ess_multi({chains{1}.xsto, chains{2}.xsto, ...})
end

function [ESS, tau_hat] = ess_multi(chains, maxLag)
% ESS and integrated autocorrelation time (tau) for multiple MCMC chains
% chains:  cell array {C} of [T x P] matrices, or 3D array [T x P x C]
% maxLag: optional, default = min(1000, floor(T/2))
%
% Returns:
%   ESS     : 1 x P vector of effective sample sizes across all chains
%   tau_hat : 1 x P vector of integrated autocorrelation times

    % ---- normalize input to cell array of [T x P] ----
    if isnumeric(chains) && ndims(chains) == 3
        C = size(chains,3);
        tmp = cell(C,1);
        for c = 1:C
            tmp{c} = chains(:,:,c);
        end
        chains = tmp;
    elseif ~iscell(chains)
        error('chains must be a cell array of [T x P] or a 3D numeric array [T x P x C].');
    end
    C = numel(chains);
    T = size(chains{1},1);
    P = size(chains{1},2);
    if nargin < 2 || isempty(maxLag)
        maxLag = min(1000, floor(T/2));
    end

    % basic checks
    for c = 1:C
        if ~isequal(size(chains{c},1), T) || ~isequal(size(chains{c},2), P)
            error('All chains must have the same [T x P] size.');
        end
    end

    % ---- compute tau per chain, per param; then combine ----
    tau_chain = nan(C, P);
    for c = 1:C
        X = chains{c};
        % center each parameter within the chain
        X = X - mean(X,1);
        for p = 1:P
            x = X(:,p);
            % autocorrelation up to maxLag (unbiased, using FFT-free formula)
            acf = acf_unbiased(x, maxLag);
            % Geyer's initial positive sequence on gamma_k = rho_{2k-1}+rho_{2k}
            gamma = acf(2:end); % drop lag 0
            % pair them: (1+2),(3+4),...
            K = floor(numel(gamma)/2);
            pairSums = gamma(2*(1:K)-1) + gamma(2*(1:K));  % rho1+rho2, rho3+rho4, ...
            % truncate at first negative pair sum
            idx = find(pairSums <= 0, 1, 'first');
            if ~isempty(idx)
                pairSums = pairSums(1:idx-1);
            end
            tau = 1 + 2*sum(pairSums);  % integrated autocorrelation time
            % guardrails
            if ~isfinite(tau) || tau < 1
                tau = 1;
            end
            tau_chain(c,p) = tau;
        end
    end

    % weight taus by chain length (all equal here)
    tau_hat = mean(tau_chain, 1, 'omitnan');
    Ntot = T * C;
    ESS = Ntot ./ tau_hat;

    % optional: tiny printout
    fprintf('ESS (per param):\n'); disp(ESS);
    fprintf('Tau  (per param):\n'); disp(tau_hat);
    fprintf('ESS fraction of total (ESS/N):\n'); disp(ESS / Ntot);
end

function acf = acf_unbiased(x, maxLag)
% Unbiased sample autocorrelation up to maxLag (lag 0..maxLag)
    x = x(:);
    n = numel(x);
    v = var(x, 1);             % population variance (1/n)
    if v == 0
        acf = [1; zeros(maxLag,1)];
        return
    end
    acf = zeros(maxLag+1,1);
    acf(1) = 1;
    for k = 1:maxLag
        % unbiased covariance estimator uses (n-k) in denominator
        cov_k = (x(1:n-k)' * x(1+k:n)) / (n - k);
        acf(k+1) = cov_k / v;
    end
end
