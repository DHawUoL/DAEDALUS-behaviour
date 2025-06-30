function plot_mcmc_correlation_full(X)
    % X: matrix of MCMC samples [N_samples x 7]
    % Columns: alpha, m, k1, k2, k3, v0, phi

    % Variable names
    varNames = {'\alpha','m','k_1','k_2','k_3','v_0','\phi'};
    nVars = length(varNames);

    % Create figure
    figure;
    [h,ax,bigax] = gplotmatrix(X,[],[],[],[],[],[],[]);

    % Loop through each subplot to add correlation coefficients and regression lines
    for i = 1:nVars
        for j = 1:nVars
            if i ~= j
                % Get data for this subplot
                xdata = X(:,j);
                ydata = X(:,i);

                % Calculate Pearson correlation
                R = corr(xdata, ydata, 'type','Pearson');

                % Add correlation text
                axes(ax(i,j));
                text(0.05, 0.85, sprintf('%.2f',R), 'Units','normalized', 'FontSize',9, 'FontWeight','bold');

                % Add regression line
                hold on
                p = polyfit(xdata, ydata, 1);
                xfit = linspace(min(xdata), max(xdata), 100);
                yfit = polyval(p, xfit);
                plot(xfit, yfit, 'r-', 'LineWidth', 1);
                hold off
            end
        end
    end

    % Label axes
    for i = 1:nVars
        xlabel(ax(end,i), varNames{i}, 'Interpreter','tex', 'FontSize', 10)
        ylabel(ax(i,1), varNames{i}, 'Interpreter','tex', 'FontSize', 10)
    end

    % General formatting
    set(ax, 'FontSize', 9)
    title(bigax, 'Posterior Correlation Matrix', 'FontSize', 12)
end
