function Xfull2 = orthogonalise_pc2_to_adm_data(Xfull, ydata, xdata, tvec, idx_pc2, first_real_col, hlag)
% Residualise PC2 against observed admissions over the calibration window.
tend=503;
tvec(tvec>tend)=[];
xdata(xdata>tend)=[];
ydata=ydata(1:length(xdata));
    if nargin < 7, hlag = 0; end
    if nargin < 6 || isempty(first_real_col), first_real_col = 4; end

    % ---- Trim to a finite horizon (your edit) ----
    tend = 503;
    tvec = tvec(tvec <= tend);
    if numel(tvec) < 2
        Xfull2 = Xfull; return;  % nothing to do
    end
    xdata = xdata(xdata <= tend);
    ydata = ydata(1:numel(xdata));
    Xfull = Xfull(:, 1:(numel(tvec)-1));

    Xfull2 = Xfull;
    L  = numel(tvec) - 1;
    j0 = max(1, min(first_real_col, L));  % safety clamp

    % ---- Build window-level admissions (respecting hlag) ----
    y_win    = zeros(1, L);
    y_smooth = movmean(ydata, 7, 'Endpoints','shrink');
    for j = 1:L
        t0 = tvec(j)     + hlag;
        t1 = tvec(j+1)-1 + hlag;
        m  = (xdata >= t0) & (xdata <= t1);
        if any(m)
            y_win(j) = mean(y_smooth(m));
        else
            y_win(j) = 0;
        end
    end

    % ---- Slices used for residualisation ----
    p = Xfull(idx_pc2, j0:end);   % 1 x M (intended)
    a = y_win(       j0:end);     % 1 x M

    % Guard: if nothing to process, bail out cleanly
    if isempty(p) || isempty(a) || numel(p) ~= numel(a)
        % fall back to no-orthogonalisation
        Xfull2(idx_pc2, :) = Xfull(idx_pc2, :);
        return;
    end

    % Force to row vectors with matching length
    p = reshape(p, 1, []);
    a = reshape(a, 1, []);

    % ---- Standardise & residualise p against a ----
    a_mu = mean(a); a_sd = std(a);
    if a_sd < 1e-12
        % Admissions window series is (near) constant → nothing to project out
        Xfull2(idx_pc2, j0:end) = p;
        return;
    end
    a0 = (a - a_mu) / a_sd;

    p_mu = mean(p); p_sd = std(p);
    if p_sd < 1e-12
        % PC2 already (near) constant, keep as-is
        Xfull2(idx_pc2, j0:end) = p;
        return;
    end
    p0 = (p - p_mu) / p_sd;

    % Scalar projection using dot products (orientation-safe)
    denom  = max(dot(a0, a0), 1e-12);
    coeff  = dot(p0, a0) / denom;
    p0_orth = p0 - coeff * a0;

    % Rescale back so k_* remains comparable to pre-orth values
    p_orth = p_mu + p_sd * p0_orth;

    % Write back only the usable columns
    Xfull2(idx_pc2, j0:end) = p_orth;
end

