function plot_mcmc_correlation_poster_spearman(X)
    % X: matrix of MCMC samples [N_samples x 7]
    % Columns: alpha, m, k1, k2, k3, v0, phi

    varNames = {'\alpha','m','k_1','k_2','k_3','v_0','\phi'};
    nVars = length(varNames);

    figure('Position',[100 100 1000 800]);
    [h,ax,bigax] = gplotmatrix(X,[],[],[],[],[],[],[]);
    
    for i = 1:nVars
        for j = 1:nVars
            if i ~= j
                xdata = X(:,j);
                ydata = X(:,i);

                % Spearman correlation
                R = corr(xdata, ydata, 'type','Spearman');
                axes(ax(i,j));
                text(0.05, 0.85, sprintf('%.2f',R), 'Units','normalized', ...
                    'FontSize',10,'FontWeight','bold','Color','k');

                % Regression line (still linear for visual guidance)
                hold on
                p = polyfit(xdata, ydata, 1);
                xfit = linspace(min(xdata), max(xdata), 100);
                yfit = polyval(p, xfit);
                plot(xfit, yfit, 'r-', 'LineWidth', 1);
                hold off
            end
        end
    end

    % Axis labels
    for i = 1:nVars
        xlabel(ax(end,i), varNames{i}, 'Interpreter','tex', 'FontSize', 12)
        ylabel(ax(i,1), varNames{i}, 'Interpreter','tex', 'FontSize', 12)
    end

    % Tighten layout
    set(ax,'FontSize',10,'LooseInset',get(gca,'TightInset'))
    set(bigax, 'Visible','off')
    set(gcf, 'Color','w');
    tightfig;
end

% --- helper function to minimize white space ---
function tightfig()
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
end
