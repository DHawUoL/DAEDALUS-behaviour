function Rhat = gelmanRubin(chains)
% chains: cell array {m x 1}, each m is an n x p matrix
% Returns: Rhat (1 x p) vector of PSRF values

burn=2000;

m = numel(chains);        % number of chains
n = size(chains{1}.xsto, 1);   % iterations per chain
p = size(chains{1}.xsto, 2);   % parameters

% Check consistent size
for i = 2:m
    if size(chains{i}.xsto,1) ~= n || size(chains{i}.xsto,2) ~= p
        error('All chains must have same length and parameter count');
    end
end

% Compute mean per chain
chainMeans = zeros(m,p);
for i = 1:m
    chainMeans(i,:) = mean(chains{i}.xsto(burn+1:end,:), 1);
end

% Overall mean
overallMean = mean(chainMeans, 1);

% Between-chain variance B
B = n * var(chainMeans, 0, 1); % 0 → normalise by m-1

% Within-chain variance W
W = zeros(1,p);
for i = 1:m
    W = W + var(chains{i}.xsto(burn+1:end,:), 0, 1);
end
W = W / m;

% Estimate of marginal posterior variance
varPlus = ((n-1)/n) .* W + (1/n) .* B;

% Potential scale reduction factor Rhat
Rhat = sqrt(varPlus ./ W);

end
