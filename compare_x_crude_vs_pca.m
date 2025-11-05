function out = compare_x_crude_vs_pca(Vraw, Vpc, qoptim5, varargin)
%compare_x_crude_vs_pca(Vraw, Vpc, qoptim5, 'COEFF', C, 'mu', mu, 'Plot', true)

%COMPARE_X_CRUDE_VS_PCA  Compare x(t) in raw vs 2-PC bases on the SAME windows.
%
% Usage (mapping known):
%   out = compare_x_crude_vs_pca(Vraw, Vpc, qoptim5, 'COEFF', C, 'mu', mu, 'Plot', true);
%
% Usage (PCA params already mapped):
%   out = compare_x_crude_vs_pca(Vraw, Vpc, qoptim5, 'k_pca', kpc, 'v0_pca', v0pc, 'Plot', true);
%
% Inputs:
%   Vraw     2×N  raw drivers (rows: v1, v2) for EXACT windows used by the fit
%   Vpc      2×N  PC drivers (rows: PC1, PC2) on the SAME windows
%   qoptim5  1×5  crude params [alpha, k1_raw, k2_raw, k3, v0_raw]
%
% Name-Value options (choose one path):
%   'COEFF'  2×2  PCA loadings from Vraw'  (columns = PCs)
%   'mu'     1×2  PCA means used with COEFF
%     (If COEFF & mu are given, k_pca and v0_pca are derived automatically.)
%
%   'k_pca'  2×1  PCA-space weights (if you already have them)
%   'v0_pca' 1×1  PCA-space intercept
%
%   'Plot'   logical (default=false) show comparison plot
%
% Output struct:
%   out.x_crude, out.x_pca, out.diff, out.max_abs_diff
%   out.k_raw, out.v0_raw, out.k_pca, out.v0_pca
%   out.has_mapping (true if COEFF/mu provided)
%
% Notes:
%   - Vraw and Vpc must correspond to the EXACT same N windows the simulator uses.
%   - If COEFF/mu are supplied, this function also checks that
%     Vpc ≈ ((Vraw' - mu) * COEFF)' and reports the max deviation.

% ---------- parse & validate ----------
p = inputParser;
p.addParameter('COEFF', [], @(x)isnumeric(x)&&isequal(size(x),[2,2]));
p.addParameter('mu',    [], @(x)isnumeric(x)&&isequal(size(x),[1,2]));
p.addParameter('k_pca', [], @(x)isnumeric(x)&&isequal(size(x),[2,1]));
p.addParameter('v0_pca',[], @(x)isnumeric(x)&&isscalar(x));
p.addParameter('Plot',  false, @(x)islogical(x)&&isscalar(x));
p.parse(varargin{:});
C      = p.Results.COEFF;
mu     = p.Results.mu;
k_pca  = p.Results.k_pca;
v0_pca = p.Results.v0_pca;
doPlot = p.Results.Plot;

assert(isnumeric(Vraw)&&size(Vraw,1)==2, 'Vraw must be 2×N.');
assert(isnumeric(Vpc) &&size(Vpc,1)==2,  'Vpc must be 2×N.');
assert(size(Vraw,2)==size(Vpc,2),        'Vraw and Vpc must have the same N columns.');
assert(isnumeric(qoptim5)&&numel(qoptim5)==5, 'qoptim5 must be [alpha k1 k2 k3 v0].');

% ---------- crude-basis x(t) ----------
k_raw  = qoptim5(2:3).';   % 2×1
v0_raw = qoptim5(5);       % scalar
x_crude = (k_raw.' * Vraw) + v0_raw;  % 1×N

% ---------- determine PCA parameters ----------
hasMapping = ~isempty(C) && ~isempty(mu);
if hasMapping
    % derive k_pca, v0_pca from mapping + sanity check Vpc construction
    k_pca_calc  = C.' * k_raw;           % 2×1
    v0_pca_calc = v0_raw + k_raw.' * mu.'; % scalar

    % if user also passed k_pca/v0_pca, check consistency; else use derived
    if ~isempty(k_pca)
        assert(isequal(size(k_pca),[2,1]), 'k_pca must be 2×1.');
        if max(abs(k_pca - k_pca_calc)) > 1e-8
            warning('k_pca differs from C''*k_raw by %.2e (using derived value).', ...
                    max(abs(k_pca - k_pca_calc)));
        end
    end
    if ~isempty(v0_pca)
        if abs(v0_pca - v0_pca_calc) > 1e-8
            warning('v0_pca differs from v0_raw + k_raw''*mu by %.2e (using derived value).', ...
                    abs(v0_pca - v0_pca_calc));
        end
    end
    k_pca  = k_pca_calc;
    v0_pca = v0_pca_calc;

    % optional: verify that Vpc matches PCA transform of Vraw
    Vpc_from_raw = ((Vraw.' - mu) * C).'; % 2×N
    pc_recon_err = max(abs(Vpc(:) - Vpc_from_raw(:)));
else
    % mapping not provided; require k_pca & v0_pca
    assert(~isempty(k_pca) && ~isempty(v0_pca), ...
        'Provide either COEFF & mu, or k_pca & v0_pca.');
    pc_recon_err = NaN; % not checked
end

% ---------- PCA-basis x(t) ----------
x_pca = (k_pca.' * Vpc) + v0_pca;  % 1×N

% ---------- stats & optional plot ----------
diff_vec = x_crude - x_pca;
max_abs_diff = max(abs(diff_vec));

if doPlot
    figure('Name','x(t): crude vs PCA');
    plot(x_crude, '-', 'LineWidth', 1.6); hold on;
    plot(x_pca,   '--','LineWidth', 1.6);
    grid on; box on; xlabel('window index'); ylabel('x(t)');
    legend('crude','PCA','Location','best');
    title(sprintf('max |x_{crude}-x_{pca}| = %.2e', max_abs_diff));
end

% ---------- package output ----------
out = struct();
out.x_crude      = x_crude;
out.x_pca        = x_pca;
out.diff         = diff_vec;
out.max_abs_diff = max_abs_diff;
out.k_raw        = k_raw;
out.v0_raw       = v0_raw;
out.k_pca        = k_pca;
out.v0_pca       = v0_pca;
out.has_mapping  = hasMapping;
out.pc_recon_err = pc_recon_err;  % only meaningful if COEFF/mu provided

% helpful console prints
fprintf('max |x_crude - x_pca| = %.3e\n', max_abs_diff);
if hasMapping
    fprintf('max |Vpc - ((Vraw''-mu)*COEFF)''| = %.3e\n', pc_recon_err);
end
end
