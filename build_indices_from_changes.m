function [tvec, Xfull_crude, Xfull_pca, pca_info, Xdaily_crude, Xdaily_pca, b0] = build_indices_from_changes(tab, t_end, varargin)
% TABLE COLUMNS (as per your example):
%   Day, Date, Trust, Stringrncy, Select
% Outputs:
%   tvec             : breakpoints, starting [1,2,61,91,...,t_end]
%   Xfull_crude      : [2 x (lt-1)] window means (Trust, -Stringrncy), first 3 cols = 0
%   Xfull_pca        : [2 x (lt-1)] window means (PC1, PC2), first 3 cols = 0
%   pca_info         : .coeff, .mu, .explained
%   Xdaily_crude     : [2 x t_end] daily (Trust, -Stringrncy) after fill
%   Xdaily_pca       : [2 x t_end] daily (PC1, PC2), zero for days < tfrom
%   b0               : -mu*coeff (your “extra negative bit” base for PCA)

% ---- options
opts.tfrom = 91;               % day where drivers first become nonzero
opts.min_win = 1;               % min window length to keep (>= tfrom)
opts.negate_stringency = true;  % use -Stringrncy
if ~isempty(varargin)
    if isstruct(varargin{1})
        s = varargin{1}; fn = fieldnames(s);
        for k=1:numel(fn), opts.(fn{k}) = s.(fn{k}); end
    else
        for k=1:2:numel(varargin), opts.(varargin{k}) = varargin{k+1}; end
    end
end
tfrom = opts.tfrom;

% ---- read columns robustly
day   = tab{:,'Day'};
trust = tab{:,'Trust'};
if ismember('Stringrncy', tab.Properties.VariableNames)
    strg = tab{:,'Stringrncy'};
else
    strg = tab{:,'Stringency'}; % fallback if spelled differently
end

if opts.negate_stringency, strg = -strg; end
if max(day) < t_end, error('t_end (%d) > max Day in table (%d).', t_end, max(day)); end

% ---- build dense daily series 1..t_end
Xdaily_crude = nan(2, t_end);
valid = (day >= 1) & (day <= t_end);
Xdaily_crude(1, day(valid)) = trust(valid).';
Xdaily_crude(2, day(valid)) = strg(valid).';

% forward fill from first non-NaN; set anything before that to 0; then back-fill
for r = 1:2
    first_non_nan = find(~isnan(Xdaily_crude(r,:)), 1, 'first');
    if isempty(first_non_nan)
        Xdaily_crude(r,:) = 0;
    else
        % set leading NaNs to 0
        if first_non_nan > 1
            Xdaily_crude(r,1:first_non_nan-1) = 0;
        end
        % forward-fill and back-fill remaining gaps
        Xdaily_crude(r,:) = fillmissing(Xdaily_crude(r,:), 'previous');
        Xdaily_crude(r,:) = fillmissing(Xdaily_crude(r,:), 'next');
    end
end

% ---- forced early breakpoints and change points from 92+
forced = [1; 2; 61; 91];
cp = forced;
tol = 1e-12;
for d = max(92,1):t_end
    if any(abs(Xdaily_crude(:,d) - Xdaily_crude(:,d-1)) > tol)
        cp(end+1,1) = d; %#ok<AGROW>
    end
end
tvec = unique([cp; t_end]);
tvec = sort(tvec(:));
L = numel(tvec) - 1;

% ---- optionally merge tiny windows (>= tfrom)
if opts.min_win > 1
    j = 1;
    while j < numel(tvec)
        wlen = tvec(j+1) - tvec(j);
        if tvec(j) >= tfrom && wlen < opts.min_win
            tvec(j+1) = [];       % merge forward
        else
            j = j + 1;
        end
    end
    L = numel(tvec) - 1;
end

% ---- aggregate CRUDE window means; zero first 3 columns
Xfull_crude = zeros(2, L);
for j = 1:L
    seg = tvec(j):max(tvec(j+1)-1, tvec(j));
    Xfull_crude(1,j) = mean(Xdaily_crude(1, seg), 'omitnan');
    Xfull_crude(2,j) = mean(Xdaily_crude(2, seg), 'omitnan');
end
k0 = min(3, L);
Xfull_crude(:,1:k0) = 0;

%Flip trust:
Xfull_crude(1,:)=1-Xfull_crude(1,:);
Xdaily_crude(1,:)=1-Xdaily_crude(1,:);

% ---- PCA on days >= tfrom, project all days, zero < tfrom
fit_mask = (1:t_end) >= tfrom;
Xfit = [Xdaily_crude(1,fit_mask).', Xdaily_crude(2,fit_mask).'];
[coeff,~,~,~,explained,mu] = pca(Xfit, 'Centered', true);

PCs = ([Xdaily_crude(1,:).', Xdaily_crude(2,:).'] - mu) * coeff;
b0  = -mu * coeff;

Xdaily_pca = zeros(2, t_end);
Xdaily_pca(1,:) = PCs(:,1)';   % PC1
Xdaily_pca(2,:) = PCs(:,2)';   % PC2
Xdaily_pca(:, 1:tfrom-1) = 0;

Xfull_pca = zeros(2, L);
for j = 1:L
    seg = tvec(j):max(tvec(j+1)-1, tvec(j));
    Xfull_pca(1,j) = mean(Xdaily_pca(1, seg), 'omitnan');
    Xfull_pca(2,j) = mean(Xdaily_pca(2, seg), 'omitnan');
end
Xfull_pca(:,1:k0) = 0;

pca_info = struct('coeff', coeff, 'mu', mu, 'explained', explained);
end
