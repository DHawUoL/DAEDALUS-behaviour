function [a,b]=plotRidgeLine(xsto)
burn=1e3;
xpost=xsto(burn+1:end,:);
param_index_alpha=1;
param_index_m=2;

% === INPUT: posterior samples ===
% Replace these with your actual parameter indices
alpha_samples = xpost(:, param_index_alpha);  % e.g. xpost(:,1)
m_samples     = xpost(:, param_index_m);      % e.g. xpost(:,2)

% === Combine and run PCA ===
X = [alpha_samples, m_samples];
[coeff, score, ~] = pca(X);

% === Extract ridge direction: m ≈ a*alpha + b ===
v = coeff(:,1);            % Principal component direction
a = v(2) / v(1);           % Slope of the line
b = mean(m_samples) - a * mean(alpha_samples);  % Intercept

% === Plot posterior samples and ridge line ===
figure;
scatter(alpha_samples, m_samples, 10, 'filled', 'MarkerFaceAlpha', 0.1);
hold on;

% Plot the ridge line
alpha_range = linspace(min(alpha_samples), max(alpha_samples), 100);
m_fit = a * alpha_range + b;
plot(alpha_range, m_fit, 'r-', 'LineWidth', 2);

xlabel('\alpha');
ylabel('m');
title('Posterior Ridge: m ≈ a·\alpha + b');
legend('Posterior samples', sprintf('m = %.3f·\\alpha + %.3f', a, b), 'Location', 'best');
grid on;
