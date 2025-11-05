function f=make4paramSeed(theta5,Xfull)

alpha5  = theta5(1);
k1_5    = theta5(2);
k2_5    = theta5(3);
k3_5    = theta5(4);
delta5  = theta5(5);

PC1series = Xfull(1,4:end);   % 1 x T
PC2series = Xfull(2,4:end);   % 1 x T

%{
%Option 1:
y = k1_5 .* PC1series + k2_5 .* PC2series;
z = PC1series;
k1_4_seed = sum(z .* y) / sum(z .* z);   % <-- LS collapse of PC1+PC2 onto PC1
k3_4_seed = k3_5;
delta4_seed = delta5;
alpha4_seed = alpha5;
%}
%
%Other option:
y = k1_5 * PC1series + k2_5 * PC2series;
z = PC2series;
k1_4_seed = sum(z .* y) / sum(z .* z);
%}

theta4_seed = [alpha4_seed, k1_4_seed, k3_4_seed, delta4_seed];
f=theta4_seed;
end
