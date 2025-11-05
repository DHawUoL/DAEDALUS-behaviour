function h = bePlotFitSimple(ydata, data, params, Xfull, tvec, b0, hosp, bestLagDays, varargin)
tvec=tvec(:).';
tvec(tvec>=250)=tvec(tvec>=250)+bestLagDays;

% bePlotFitSimple  Plot model fit (single line) vs admissions data,
% and (optionally) p(t)=sigma(x(t)) on a right axis.
%
% Required:
%   ydata, data, params, Xfull, tvec, b0
%
% Name/Value options:
%   'Factor'     (default 5e3)          scale for admissions (Admissions/5k)
%   'Color'      (default [0 0.447 0.741])
%   'BarAlpha'   (default 0.35)
%   'LineWidth'  (default 2.2)
%   'ShowDelta'  (default false)        dashed horizontal line at params(1)
%   'ShowPT'     (default false)        plot p(t)=sigmoid(x(t)) on right axis
%   'Xdaily12'   (default [])           2xT daily inputs [v1; v2]
%   'V3'         (default [])           v3; daily 1xT or windowed 1x(L) (L=numel(tvec)-1)
%   'V3Mode'     (default 'window')     'daily' or 'window'
%
% Assumes sim2fit_global signature:
%   sim2fit_global(params,data,xdata,X,intrinsic,Xfull,coeff,tvec,lx1,lx2,plotRun,ymean,b0)

% ---- options ----
p = inputParser;
addParameter(p,'Factor',5e3);
addParameter(p,'Color',[0 0.447 0.741]);
addParameter(p,'BarAlpha',0.35);
addParameter(p,'LineWidth',2.2);
addParameter(p,'ShowDelta',false);
addParameter(p,'ShowPT',false);
addParameter(p,'Xdaily12',[]);

%addParameter(p,'V3',[]);
tvech=tvec(4:end-1)-84;
hvec=[0,0,0,hosp(tvech)']/5e3;
addParameter(p,'V3',hvec);

addParameter(p,'V3Mode','window');  % 'daily' or 'window'
parse(p,varargin{:});
Factor    = p.Results.Factor;
col       = p.Results.Color;
barAlpha  = p.Results.BarAlpha;
lw        = p.Results.LineWidth;
showDelta = p.Results.ShowDelta;
showPT    = p.Results.ShowPT;
Xdaily12  = p.Results.Xdaily12;
V3        = p.Results.V3;
V3Mode    = lower(p.Results.V3Mode);

% ---- prep for simulator ----
xdata = 85:tvec(end);
%X      = ones(1, size(Xfull,2));
%coeff  = ones(1,3)';            % crude model
[~,lx2] = size(Xfull);

% England scaling (as in your code)
y = ydata(0+(1:numel(xdata)));
y = y * (sum(data.Npop)/56286961);
ymean = 0;

% Simulator handle (correct full signature)
intrinsic = 1; plotRun = 0;
fun = @(prm) sim2fit_global(prm, data, xdata, intrinsic, Xfull, tvec, lx2, plotRun, ymean, b0);

% Evaluate model
[yhat, ~] = fun(params);

% ---- main plot ----
figure('Color','w'); hold on;
bar(xdata, y/Factor, 'FaceColor', 0.5*[1 1 1], 'EdgeColor', 0.5*[1 1 1], ...
    'LineWidth', .5, 'FaceAlpha', barAlpha);
plot(xdata, yhat/Factor, '-', 'Color', col, 'LineWidth', lw);

% Optional horizontal delta line (paper's effectiveness baseline if you want it)
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

% ---- optional p(t) from WINDOW values (Xfull & V3win) ----
if showPT
    if numel(params) < 4
        warning('ShowPT requested but params length < 5; skipping p(t).');
    else
        if numel(params)==5
            k1 = params(2); k2 = params(3); k3 = params(4); delta = params(5);
        elseif numel(params)==4
            k1 = params(2); 
            k2 = params(3); 
            k3 = 0;%params(3); 
            delta = params(4);
        else
            error('Too many parameters')

        L = numel(tvec)-1;                % number of windows
        % v1, v2 from Xfull windows:
        if size(Xfull,1) < 2 || size(Xfull,2) < L
            warning('Xfull must be at least 2 x (numel(tvec)-1). Skipping p(t).');
        else
            v1w = Xfull(1,1:L);
            v2w = Xfull(2,1:L);

            % v3 windowed: prefer external V3win if present, else Xfull(3,:)
            %if exist('V3win','var') && numel(V3win) == L
                v3w = p.Results.V3(:).';
            %elseif size(Xfull,1) >= 3
                %v3w = Xfull(3,1:L);
            %else
                %warning('No V3win provided and Xfull has <3 rows; assuming v3=0.');
                %v3w = zeros(1,L);
            %end

            % x(window) and p(window)
            x_win = k1*v1w + k2*v2w + k3*v3w - delta;
            p_win = 1 ./ (1 + exp(-x_win));

            % Force p(t)=0 in first 3 windows
            zN = min(3, L);
            p_win(1:zN) = 0;

            % Build step vectors that align with tvec edges
            tt = [tvec(1:end-1), tvec(end)];       % edges, length L+1
            pp = [p_win, p_win(end)];              % hold last value to the end

            % Plot on right axis (only show in visible x-range)
            yyaxis right
            stairs(tt, pp, 'k--', 'LineWidth', 1.6);
            ylim([0 1]); ylabel('p(t)');
            yyaxis left

            % Update legend entries gracefully
            if showDelta
                legend({'Data','Model','Effectiveness \delta','p(t)'}, 'Location','southwest');
            else
                legend({'Data','Model','p(t)'}, 'Location','southwest');
            end
        end
    end
end
end