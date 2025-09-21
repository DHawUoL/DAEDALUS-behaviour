function f=MCMCdiagnostics(chains)

figure;
for c = 1:numel(chains)
    subplot(numel(chains),1,c);
    plot(chains{c}.xsto(:,1));  % first column = alpha
    ylabel(sprintf('Chain %d',c));
end
xlabel('Iteration');
sgtitle('Trace plots for \alpha');



paramIdx = 3; % which parameter to check
figure; hold on;
for c = 1:numel(chains)
    histogram(chains{c}.xsto(:,paramIdx), 'Normalization','pdf', 'DisplayStyle','stairs');
end
xlabel(sprintf('Parameter %d',paramIdx));
ylabel('Posterior density');
legend(arrayfun(@(i) sprintf('Chain %d',i),1:numel(chains),'UniformOutput',false));



% Usage:
RhatVals = gelman_rubin(chains);
disp(RhatVals);

end



function Rhat = gelman_rubin(chains)
    m = numel(chains);
    n = size(chains{1}.xsto,1);
    p = size(chains{1}.xsto,2);
    Rhat = zeros(1,p);

    for j = 1:p
        % Extract parameter j from all chains
        X = cell2mat(cellfun(@(c) c.xsto(:,j), chains, 'UniformOutput',false));
        X = reshape(X, n, m);

        % Between-chain variance
        chainMeans = mean(X);
        B = n * var(chainMeans, 1);

        % Within-chain variance
        W = mean(var(X, 0, 1));

        % Estimate of marginal posterior variance
        varHat = ((n-1)/n) * W + (1/n) * B;

        % Potential scale reduction factor
        Rhat(j) = sqrt(varHat / W);
    end
end



