function plotmat = bePlotFrequentist_final(ydata, data, Xfull, pointEst, Diag, tvec, b0, bestLagDays)
tvec=tvec(:).';
tvec(tvec>=250)=tvec(tvec>=250)+bestLagDays;

% Frequentist fit visualization:
%   - bounded parametric fan using Diag.pcov (+bounds clipping)
%   - linearized delta-method band around MLE
%   - in-sample vs projection split

% -------------------- switches / styling --------------------
%hlag        = 0;
intrinsic   = 1;
seed_draws  = 42;
n_draws     = 100;
z_delta     = 1.96;%0.674;%1.96;
bar_alpha   = 0.35;
in_alpha    = 0.10;
out_alpha   = 0.25;
line_col    = [0 0.447 0.741];  % lines(1)

xdata = 85:tvec(end);
[~,lx2] = size(Xfull);
lt = numel(tvec);                      % <-- FIX: was missing

% -------------------- scale data (England) --------------------
ydata = ydata(0+(1:numel(xdata)));
ydata = ydata * (sum(data.Npop)/56286961);
ymean = 0;

% -------------------- simulator handle --------------------
% Use the function name you actually have available:
fun_full  = @(params) sim2fit_global(params, data, xdata, intrinsic, Xfull, tvec, lx2, 0, ymean, b0);
fun_model = @(p) fun_full(p);   % single-output context -> returns f only

% -------------------- parametric fan --------------------
fan_opts = struct('truncate','clip','seed',seed_draws);
sample   = drawParamSamples(Diag, n_draws, fan_opts);

plotmat  = zeros(n_draws, numel(xdata));
for i = 1:n_draws
    [fi, ~]       = fun_full(sample(i,:));
    plotmat(i,:)  = fi;
end

% -------------------- MLE trajectory & p(t) -----------------
[f_hat, rhohat] = fun_full(Diag.p_hat);

% piecewise p(t) per day from rhohat window-means
tvecPlus = [1, tvec(2:end)];
value    = zeros(1, tvec(end));
for j = 1:lt-1
    tj = round(tvecPlus(j)) : round(tvecPlus(j+1));
    tj(tj<1 | tj>numel(value)) = [];   % just in case
    if ~isempty(tj), value(tj) = rhohat(j); end
end
value(1:round(tvec(3))-1) = 0;

% -------------------- delta-method band --------------------
[mu_lin, lo_lin, hi_lin] = delta_band(fun_model, Diag.p_hat, Diag.pcov, xdata, z_delta);

% -------------------- split in-sample vs projection --------
t_cut  = tvec(max(1, lt-7));          % <-- guard for short tvec
ix_in  = xdata <= t_cut;
ix_out = xdata >  t_cut;

% prctile can be 3D -> squeeze to [3 x T]
pr_in  = squeeze(prctile(plotmat(:,ix_in),  [25 50 75], 1));   % [3 x Tin]
pr_out = squeeze(prctile(plotmat(:,ix_out), [25 50 75], 1));   % [3 x Tout]

% -------------------- scaling ------------------------------
ycap   = max([ydata(:); prctile(plotmat(:),99); hi_lin(:)]);
factor = 5e3;

ydata_s  = ydata / factor;
plotmat_s= plotmat / factor;
mu_s     = mu_lin / factor;  lo_s = lo_lin / factor;  hi_s = hi_lin / factor;
pr_in    = pr_in / factor;   pr_out = pr_out / factor;
ycap_s   = ycap / factor;

% -------------------- plotting ------------------------------
figure('Color','w'); hold on;

% Data
bar(xdata, ydata_s, 'FaceColor', 0.5*[1 1 1], ...
    'EdgeColor', 0.5*[1 1 1], 'LineWidth', .5, 'FaceAlpha', bar_alpha);

% Linearized band
fill([xdata fliplr(xdata)], [lo_s' fliplr(hi_s')], line_col, ...
     'FaceAlpha',0.12,'EdgeColor','none');
plot(xdata, mu_s, 'Color', line_col, 'LineWidth', 2.2);

% Parametric fan (in-sample)
if any(ix_in)
    fill([xdata(ix_in) fliplr(xdata(ix_in))], ...
         [pr_in(1,:) fliplr(pr_in(3,:))], line_col, ...
         'FaceAlpha', in_alpha, 'EdgeColor','none');
    plot(xdata(ix_in), pr_in(2,:), 'Color', line_col, 'LineWidth', 1.2);
end

% Parametric fan (projection)
if any(ix_out)
    fill([xdata(ix_out) fliplr(xdata(ix_out))], ...
         [pr_out(1,:) fliplr(pr_out(3,:))], line_col, ...
         'FaceAlpha', out_alpha, 'EdgeColor','none');
    plot(xdata(ix_out), pr_out(2,:), 'Color', line_col, 'LineWidth', 2.0);
end

% Boundary
plot(t_cut*[1 1], [0 1], 'k:', 'LineWidth', 1.5);

% p(t) on right axis
yyaxis right
vx = 1:numel(value);  % ensure domain covers xdata
plot(vx, value, 'k--', 'LineWidth', 1.8);
ylim([0 1]);
ylabel('p(t)')
yyaxis left

% Horizontal effectiveness line (paper’s \delta in pointEst(1))
plot(xdata, pointEst(1)*ones(size(xdata)), '--', 'Color', [0.5 0 0], 'LineWidth', 1.6);

% Axes & labels
xlim([xdata(1), xdata(end)]);
ylim([0, 1.05*max(1.0, ycap_s)]);
xlabel('Time'); ylabel('Hospital Admissions / 5k');

% Month ticks
monthDur   = [1,31,29,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31,31];
monthStart = cumsum(monthDur);
xnames = {'Jan 2020','Feb 2020','Mar 2020','Apr 2020','May 2020','Jun 2020','Jul 2020','Aug 2020','Sep 2020','Oct 2020','Nov 2020','Dec 2020', ...
          'Jan 2021','Feb 2021','March 2021','Apr 2021'};
xvec   = monthStart(1:numel(xnames));
xticks(xvec); xticklabels(xnames); xtickangle(45);

grid on; box on;

% Legend
hData  = plot(NaN,NaN,'s','MarkerSize',8,'MarkerFaceColor',0.5*[1 1 1], 'MarkerEdgeColor',0.5*[1 1 1]);
hMLE   = plot(NaN,NaN,'-','Color',line_col,'LineWidth',2.2);
hPT    = plot(NaN,NaN,'k--','LineWidth',1.8);
hDelta = plot(NaN,NaN,'--','Color',[0.5 0 0],'LineWidth',1.6);
legend([hData,hMLE,hPT,hDelta], {'Data','MLE (delta-band)','p(t)','Effectiveness \delta'}, 'Location','southeast');

end


% ======================================================================
% =======  DELTA-METHOD LINEARIZED CONFIDENCE BAND  =====================
% ======================================================================
function [mu, lo, hi] = delta_band(fun, p_hat, C, xdata, z)
% Linearize model at p_hat: f(p) ≈ f(p_hat) + J (p - p_hat)
% Var(f_t) = J_t * C * J_t'   -> band mu ± z*se
    f0  = fun(p_hat); f0 = f0(:);
    K   = numel(p_hat);  T = numel(f0);
    J   = zeros(T, K);
    eps = 1e-5*(1+abs(p_hat(:)'));
    for k = 1:K
        dk = zeros(size(p_hat)); dk(k) = eps(k);
        fp = fun(p_hat + dk);
        fm = fun(p_hat - dk);
        J(:,k) = (fp(:) - fm(:)) / (2*eps(k));
    end
    Vt = sum((J*C).*J, 2);
    se = sqrt(max(Vt, 0));
    mu = f0;
    lo = mu - z*se;
    hi = mu + z*se;
end

% ======================================================================
% =======  DRAW PARAMETER SAMPLES WITH BOUNDS & COVARIANCE  =============
% ======================================================================
function sample = drawParamSamples(Diag, n, opts)
% drawParamSamples  Draw θ ~ N(θ̂, pcov), with optional clipping to bounds.
    if nargin < 3, opts = struct; end
    if ~isfield(opts,'truncate'), opts.truncate = 'clip'; end
    if ~isfield(opts,'seed'),     opts.seed     = [];    end
    if ~isfield(Diag,'p_hat'),    error('Diag.p_hat missing'); end

    mu = Diag.p_hat(:);
    if isfield(Diag,'pcov') && ~isempty(Diag.pcov)
        C = (Diag.pcov + Diag.pcov.')/2;
    elseif isfield(Diag,'se') && ~isempty(Diag.se)
        C = diag(Diag.se(:).^2);
    else
        error('No covariance or SEs in Diag.');
    end

    % make PD
    [L,p] = chol(C,'lower');
    if p~=0
        [V,D] = eig(C);
        D     = diag(max(diag(D),0));
        Cfix  = V*D*V' + 1e-10*eye(numel(mu));
        L     = chol((Cfix+Cfix.')/2, 'lower');
    end

    if ~isempty(opts.seed), rng(opts.seed); end
    Z      = randn(numel(mu), n);
    draws  = mu + L*Z;            % K x n
    sample = draws.';             % n x K

    % bounds from Diag if present
    if ~isfield(opts,'lb') || isempty(opts.lb)
        opts.lb = [];
        if isfield(Diag,'lb') && ~isempty(Diag.lb), opts.lb = Diag.lb(:)'; end
    end
    if ~isfield(opts,'ub') || isempty(opts.ub)
        opts.ub = [];
        if isfield(Diag,'ub') && ~isempty(Diag.ub), opts.ub = Diag.ub(:)'; end
    end

    if strcmpi(opts.truncate,'clip')
        if ~isempty(opts.lb), sample = max(sample, opts.lb); end
        if ~isempty(opts.ub), sample = min(sample, opts.ub); end
    end
end
