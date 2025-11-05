function [tvec, Xfull_crude, Xfull_pca, pca_info, Xdaily_crude, Xdaily_pca,b0] = ...
    build_indices_with_pca(tab, t_end, opts)

% BUILD_INDICES_WITH_PCA_FROM_FLAGS
% Table layout:
%   col1 = day number, col3 = trust, col4 = stringency, col5 = selected flag (1 = breakpoint)
%
% Inputs
%   tab   : table with columns as above
%   t_end : last day included (so tvec(end) = t_end)
%   opts.negate_stringency (false) -> if true, uses -stringency
%
% Outputs
%   tvec          : breakpoints, first three windows are forced zeros; tvec(end)=t_end
%   Xfull_crude   : [2 x (lt-1)] (trust, stringency window means), first 3 cols = 0
%   Xfull_pca     : [2 x (lt-1)] (PC1, PC2 window means), first 3 cols = 0
%   pca_info      : struct with fields .coeff (2x2), .mu (1x2), .explained (1x2)
%   Xdaily_crude  : [2 x t_end] daily (trust, stringency)
%   Xdaily_pca    : [2 x t_end] daily (PC1, PC2), zeroed for days < tfrom
%   (hard coded below)

    arguments
        tab table
        t_end (1,1) double
        opts.negate_stringency (1,1) logical = false
    end
    
    tfrom=92;

    % ----- 0) Pull and pre-process daily series -----
    day   = tab{:,1};
    trust = tab{:,3};
    strg  = tab{:,4};
    flag  = tab{:,5};

    %if opts.negate_stringency
        strg = -strg;
    %end

    if max(day) < t_end
        error('t_end (%d) exceeds max day in table (%d).', t_end, max(day));
    end

    use_mask = (day >= day(1)) & (day <= t_end);
    day   = day(use_mask);
    trust = trust(use_mask);
    strg  = strg(use_mask);
    flag  = flag(use_mask);

    % ----- 1) Build tvec from flags, plus the 3 zero windows and the end day -----
    % First 3 windows: [1,2), [2,61), [61,tfrom) -> breakpoints at [1,2,61,tfrom]
    forced = [1; 2; 61; tfrom];

    % Use flagged rows as additional breakpoints (within range and from tfrom onward)
    flag_days = day(flag == 1);
    flag_days = flag_days(flag_days >= tfrom & flag_days <= t_end);

    % Merge and finalize
    tvec = unique([forced; flag_days; t_end]);
    tvec = tvec(tvec >= day(1));     % ensure not before the start
    tvec = sort(tvec(:));
    if tvec(end) ~= t_end
        tvec(end+1) = t_end;
    end

    L  = numel(tvec) - 1;     % windows
    lt = numel(tvec);

    % ----- 2) Daily matrices -----
    last_day = t_end;
    Xdaily_crude = zeros(2, last_day);
    % Assign daily (assumes 'day' are consecutive or at least unique)
    Xdaily_crude(1, day) = trust;
    Xdaily_crude(2, day) = strg;

    % ----- 3) Aggregate CRUDE window means; zero first 3 columns -----
    Xfull_crude = zeros(2, L);
    for j = 1:L
        d0 = tvec(j);
        d1 = tvec(j+1) - 1;
        if d1 < d0, d1 = d0; end
        seg = d0:d1;
        Xfull_crude(1,j) = mean(Xdaily_crude(1, seg), 'omitnan');
        Xfull_crude(2,j) = mean(Xdaily_crude(2, seg), 'omitnan');
    end
    k0 = min(3, L);
    Xfull_crude(:,1:k0) = 0;

    % ----- 4) PCA fit on days >= tfrom, project all days, zero < tfrom -----
    fit_mask = (day >= tfrom) & (day <= t_end);
    Xfit = [trust(fit_mask), strg(fit_mask)];   % N x 2

    if size(Xfit,1) < 5
        error('Too few days >= tfrom to run PCA (got %d).', size(Xfit,1));
    end

    [coeff, ~, ~, ~, explained, mu] = pca(Xfit, 'Centered', true);

    % Project *all* days 1..t_end with same mu/coeff
    Xall = [Xdaily_crude(1,:).', Xdaily_crude(2,:).'];    % t_end x 2
    PCs  = (Xall - mu) * coeff;                            % t_end x 2

    b0 = -mu  * coeff;

    Xdaily_pca = zeros(2, last_day);
    Xdaily_pca(1,:) = PCs(:,1)';   % PC1
    Xdaily_pca(2,:) = PCs(:,2)';   % PC2
    % Zero before day tfrom to match model convention
    Xdaily_pca(:, 1:tfrom) = 0;

    % Aggregate PCA window means; zero first 3 columns
    Xfull_pca = zeros(2, L);
    for j = 1:L
        d0 = tvec(j);
        d1 = tvec(j+1) - 1;
        if d1 < d0, d1 = d0; end
        seg = d0:d1;
        Xfull_pca(1,j) = mean(Xdaily_pca(1, seg), 'omitnan'); % PC1
        Xfull_pca(2,j) = mean(Xdaily_pca(2, seg), 'omitnan'); % PC2
    end
    Xfull_pca(:,1:k0) = 0;

    % ----- 5) Pack PCA info -----
    pca_info = struct('coeff', coeff, 'mu', mu, 'explained', explained);

    tvec=tvec';
end
