function stats = info_criteria(y, yhat, k, weights)
% INFO_CRITERIA  Compute AIC, AICc, and BIC for (possibly) weighted LS fits.
%
% Inputs
%   y       : n×1 data vector
%   yhat    : n×1 model predictions at the optimum
%   k       : number of free parameters estimated in the model
%   weights : (optional) n×1 nonnegative weights; default = ones(n,1)
%
% Interpretation
%   Uses the Gaussian GLS/Laplace form:  -2 log L  ~  n*log(RSS_w/n) + const
%   where RSS_w = sum_i w_i * (y_i - yhat_i)^2.
%   Constants cancel when comparing models fit to the same data/weights.
%
% Output (struct)
%   .n, .k, .RSSw, .sigma2hat, .AIC, .AICc, .BIC

    if nargin < 4 || isempty(weights)
        weights = ones(size(y));
    end
    y     = y(:);
    yhat  = yhat(:);
    w     = weights(:);
    assert(numel(y)==numel(yhat) && numel(y)==numel(w), 'Length mismatch');

    n        = numel(y);
    resid    = y - yhat;
    RSSw     = sum(w .* (resid.^2));        % your weighted SSE
    sigma2   = RSSw / n;                    % MLE for weighted Gaussian scale

    % Information criteria (constants dropped; fine for comparisons)
    AIC  = n*log(sigma2) + 2*k;
    % small-sample correction (only valid if n > k+1)
    if n > (k + 1)
        AICc = AIC + (2*k*(k+1)) / (n - k - 1);
    else
        AICc = NaN;
    end
    BIC  = n*log(sigma2) + k*log(n);

    stats = struct('n',n,'k',k,'RSSw',RSSw,'sigma2hat',sigma2, ...
                   'AIC',AIC,'AICc',AICc,'BIC',BIC);
end
