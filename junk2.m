function h = bePlotFitSimple(ydata, data, params, Xfull, tvec, b0, varargin)
X=ones(size(Xfull,2));
coeff=ones(1,3)';

% bePlotFitSimple  Plot model fit (single line) vs admissions data.
% Usage:
%   bePlotFitSimple(hosp, X, dataUK, popt, Xfull, coeff, tvec, b0);
%
% Optional name/value:
%   'Factor'     (default 5e3)   scale divisor for readability (Admissions/5k)
%   'Color'      (default [0 0.447 0.741])
%   'BarAlpha'   (default 0.35)
%   'LineWidth'  (default 2.2)
%   'ShowDelta'  (default false)  if true, draws horizontal line at params(1)
%
% Assumes your simulator: sim2fit_global(params, data, xdata, X, intrinsic, Xfull, coeff, tvec, lx1, lx2, plotRun, ymean, b0)

% ---- options ----
p = inputParser;
addParameter(p,'Factor',5e3);
addParameter(p,'Color',[0 0.447 0.741]);
addParameter(p,'BarAlpha',0.35);
addParameter(p,'LineWidth',2.2);
addParameter(p,'ShowDelta',false);
parse(p,varargin{:});
Factor    = p.Results.Factor;
col       = p.Results.Color;
barAlpha  = p.Results.BarAlpha;
lw        = p.Results.LineWidth;
showDelta = p.Results.ShowDelta;

% ---- prep ----
xdata = 85:tvec(end);
[lx1,lx2] = size(X);

% England scaling (as in your code)
y = ydata(0+(1:numel(xdata)));
y = y * (sum(data.Npop)/56286961);
ymean = 0;

% Simulator handle
intrinsic = 1; plotRun = 0;
fun = @(prm) sim2fit_global(prm, data, xdata, intrinsic, Xfull, tvec, lx2, plotRun, ymean, b0);

% Evaluate model
[yhat, ~] = fun(params);

% ---- plot ----
figure('Color','w'); hold on;
bar(xdata, y/Factor, 'FaceColor', 0.5*[1 1 1], 'EdgeColor', 0.5*[1 1 1], ...
    'LineWidth', .5, 'FaceAlpha', barAlpha);
plot(xdata, yhat/Factor, '-', 'Color', col, 'LineWidth', lw);

% Optional horizontal effectiveness line (your paper's delta)
if showDelta
    plot(xdata, params(1)*ones(size(xdata)), '--', 'Color', [0.5 0 0], 'LineWidth', 1.6);
end

% Axes & labels
xlabel('Time'); ylabel(sprintf('Hospital Admissions / %.0f', Factor));
grid on; box on;
xlim([xdata(1), xdata(end)]);

% Month ticks (Jan 2020 start)
monthDur   = [1,31,29,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31,31];
monthStart = cumsum(monthDur);
names = {'Jan 2020','Feb 2020','Mar 2020','Apr 2020','May 2020','Jun 2020','Jul 2020','Aug 2020', ...
         'Sep 2020','Oct 2020','Nov 2020','Dec 2020','Jan 2021','Feb 2021','Mar 2021','Apr 2021', ...
         'May 2021','Jun 2021','Jul 2021','Aug 2021','Sep 2021','Oct 2021','Nov 2021','Dec 2021','Jan 2022'};
mx = monthStart(1:min(numel(names), find(monthStart<=xdata(end),1,'last')));
xticks(mx); xticklabels(names(1:numel(mx))); xtickangle(45);

% Legend
h = struct;
h.data = plot(nan,nan,'s','MarkerSize',8,'MarkerFaceColor',0.5*[1 1 1],'MarkerEdgeColor',0.5*[1 1 1]);
h.fit  = plot(nan,nan,'-','Color',col,'LineWidth',lw);
if showDelta
    h.delta = plot(nan,nan,'--','Color',[0.5 0 0],'LineWidth',1.6);
    legend([h.data h.fit h.delta], {'Data','Model','Effectiveness \delta'}, 'Location','southeast');
else
    legend([h.data h.fit], {'Data','Model'}, 'Location','southeast');
end
end
