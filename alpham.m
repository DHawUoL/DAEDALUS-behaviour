function [a,b]=alpham(xsto)
burn=500;
alpha_samples=xsto(burn+1:end,1);
m_samples=xsto(burn+1:end,2);
[coeff, score, ~] = pca([alpha_samples, m_samples]);
% Direction vector of PC1
v = coeff(:,1);  % [v_alpha; v_m]
% Express m = a*alpha + b
a = v(2) / v(1);  % slope
b = mean(m_samples) - a * mean(alpha_samples);  % intercept
