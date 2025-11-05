function S = diag_pca_equivalence(poptim5, V12_raw, coeff, mu, varargin)
% DIAG_PCA_EQUIVALENCE
% Verifies that replacing (k1,k2,v0) with (k_pca, v0_pca or eta0) reproduces x(t).
% 
% Inputs
%   poptim5 : [alpha,k1,k2,k3,v0] from the 5-param *raw indices* fit
%   V12_raw : 2 x T raw index matrix [v1; v2] for all windows used in the fit
%   coeff   : 2 x 2 PCA loadings C  (from your pcaReduction)
%   mu      : 1 x 2 PCA mean of raw features (same as from pcaReduction, row vector)
%   varargin{1} (optional): b0 = -mu*coeff, 1x2, for reporting only
%
% Outputs (struct S)
%   k_raw, k_pca, v0_raw, v0_pca, eta0
%   max_abs_diff_direct : max |x_raw - x_pca_direct|
%   max_abs_diff_eta    : max |x_raw - x_pca_using_eta|
%   also prints a compact report

    softplus     = @(z) log1p(exp(-abs(z))) + max(z,0);
    softplus_inv = @(y) log(exp(max(y,1e-12)) - 1);

    if size(mu,1) == 1, mu = mu(:); end     % columnize
    if size(V12_raw,1) ~= 2
        error('V12_raw must be 2 x T (rows: v1, v2).');
    end
    if ~isequal(size(coeff), [2 2])
        error('coeff must be 2x2.');
    end

    % unpack raw
    k_raw  = poptim5(2:3).';           % 2x1
    v0_raw = poptim5(5);               % scalar
    T      = size(V12_raw,2);

    % build Z = C'*(x - mu)
    Z = coeff' * (V12_raw - mu*ones(1,T));     % 2xT

    % textbook mapping:
    % x_raw(t)  = k_raw' * x(t) + v0_raw
    % x_pca(t)  = (C k_raw)' * z(t) + (v0_raw + k_raw' mu)
    k_pca  = coeff' * k_raw;                    % 2x1
    v0_pca = v0_raw + dot(k_raw, mu);          % scalar

    x_raw = (k_raw.' * V12_raw) + v0_raw;      % 1xT
    x_pca_direct = (k_pca.' * Z) + v0_pca;     % 1xT

    % eta0 path (what sim uses internally): v0 = -softplus(eta0)
    eta0  = softplus_inv(-v0_pca);             % requires v0_pca < 0
    x_pca_eta = (k_pca.' * Z) - softplus(eta0);

    % basic diffs
    d1 = x_raw - x_pca_direct;
    d2 = x_raw - x_pca_eta;

    % package + print
    S.k_raw  = k_raw;    S.k_pca = k_pca;
    S.v0_raw = v0_raw;   S.v0_pca = v0_pca;   S.eta0 = eta0;
    S.max_abs_diff_direct = max(abs(d1));
    S.max_abs_diff_eta    = max(abs(d2));
    if ~isempty(varargin), S.b0 = varargin{1}; end

    fprintf('\n=== PCA mapping diagnostics ===\n');
    fprintf('k_raw  = [% .4f % .4f]^T\n', k_raw);
    fprintf('k_pca  = C*k_raw = [% .4f % .4f]^T\n', k_pca);
    fprintf('v0_raw = % .4f\n', v0_raw);
    fprintf('v0_pca = v0_raw + k_raw''*mu = % .4f\n', v0_pca);
    fprintf('eta0   = softplus^{-1}(-v0_pca) = % .4f\n', eta0);
    fprintf('max |x_raw - x_pca_direct| = %.3g\n', S.max_abs_diff_direct);
    fprintf('max |x_raw - x_pca_eta|    = %.3g\n', S.max_abs_diff_eta);
    if ~isempty(varargin)
        fprintf('b0 (=-mu*C) = [% .4f % .4f]\n', S.b0);
    end
    fprintf('================================\n');
end
