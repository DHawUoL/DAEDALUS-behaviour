function [xsto, outsto, history, accept_rate, covmat] = MCMC_adaptive_row(F, x0, n, sigma, fixinds, blockind, displ)

symmetrize = 0;
x0 = x0(:)';                       % ensure row
d  = size(x0,2);
b  = 0.2;

long_every   = 50;
long_factor1 = 4.0;
long_factor2 = 4.0;

% per-parameter scaling support
if isscalar(sigma)
    scaleLeft = 1;
else
    if isvector(sigma) && numel(sigma)==d
        scaleLeft = diag(sigma(:));
    else
        error('sigma must be scalar or a length-d vector');
    end
end
sd_scalar = (2.4^2 / d);

if ~isscalar(sigma)
    sigma = (sigma + sigma.')/2;
    sigma = sigma + 1e-12*eye(size(sigma));
end

if nargin < 6 || isempty(blockind), blockind = 0; end  % handle "no blocks"

inds = []; vals = [];
if ~isempty(fixinds)
    inds = fixinds(1,:); 
    vals = fixinds(2,:); 
    vals = vals(:)';                  % ensure row
end

% Initial covariance mask
cov0 = eye(d);
if blockind > 0
    cov0(1:blockind, (blockind+1):end) = 0;
    cov0((blockind+1):end, 1:blockind) = 0;
end
if ~isempty(inds)
    cov0(inds,:) = 0; cov0(:,inds) = 0;
end
if symmetrize==1
    cov0 = (cov0 + cov0.')/2; 
    cov0 = cov0 + 1e-12*eye(size(cov0));
end

% Outputs (row orientation: n x d)
xsto    = zeros(n,d);
outsto  = zeros(n,1);
history = zeros(n,d+1);  % [proposal , accepted_flag]

xsto(1,:) = x0;
FX = F(x0); outsto(1) = FX;
if ~isscalar(FX) || ~isfinite(FX), error('F(x0) must be finite scalar.'); end
acc = 0;

B1 = 1:blockind;
B2 = (blockind+1):d;

for t = 2:n
    X = xsto(t-1,:);                 % row

    % --- proposal
    Y0 = mvnrnd(X, cov0*(0.4^2)* (isscalar(sigma)*sigma + ~isscalar(sigma)) / d);  % row

    if t < 2*d
        Y = Y0;
        if ~isempty(inds), Y(inds) = vals; end
    else
        covmat = cov(xsto(1:t-1,:));                % d x d
        if ~isempty(inds), covmat(inds,:) = 0; covmat(:,inds) = 0; end
        if blockind > 0
            covmat(1:blockind, (blockind+1):end) = 0;
            covmat((blockind+1):end, 1:blockind) = 0;
        end
        covmat = (covmat + covmat.')/2 + 1e-12*eye(d);

        if isscalar(scaleLeft)
            sdmatrix = sd_scalar * covmat;
        else
            sdmatrix = sd_scalar * (scaleLeft * covmat * scaleLeft);
        end
        if ~isscalar(sigma)
            sdmatrix = (sdmatrix + sdmatrix.')/2 + 1e-12*eye(d);
        end
        if symmetrize==1
            [V,D] = eig(sdmatrix);
            D(D<1e-12) = 1e-12;
            sdmatrix = V*D*V.';
            sdmatrix = (sdmatrix + sdmatrix.')/2;
        end

        do_long = (t > 2*d) && (mod(t,long_every)==0) && (blockind>0) && (blockind<d);
        if do_long
            if rand < 0.5 && ~isempty(B1)
                Cb = cov(xsto(1:t-1, B1));
                Cb = (Cb+Cb')/2 + 1e-12*eye(numel(B1));
                [V,D] = eig(Cb);
                [~,ix] = max(diag(D)); v1 = V(:,ix);
                step = long_factor1 * sqrt(max(D(ix,ix),1e-12)) * randn;
                Y = X; 
                Y(B1) = X(B1) + (step*v1)';        % keep row
            else
                if ~isempty(B2)
                    Cb = cov(xsto(1:t-1, B2));
                    Cb = (Cb+Cb')/2 + 1e-12*eye(numel(B2));
                    [V,D] = eig(Cb);
                    [~,ix] = max(diag(D)); v1 = V(:,ix);
                    step = long_factor2 * sqrt(max(D(ix,ix),1e-12)) * randn;
                    Y = X;
                    Y(B2) = X(B2) + (step*v1)';    % keep row
                else
                    Y = (1-b)*mvnrnd(X, sdmatrix) + b*mvnrnd(X, cov0*(0.5^2)* (isscalar(sigma)*sigma + ~isscalar(sigma))/d);
                end
            end
        else
            Yg = mvnrnd(X, cov0*(0.5^2)* (isscalar(sigma)*sigma + ~isscalar(sigma))/d);
            Y  = (1-b)*mvnrnd(X, sdmatrix) + b*Yg;   % row
        end

        % optional DE-style kick (enable if desired)
        % if (t > 200) && (mod(t,50) == 0)
        %     r1 = randi([max(2,t-200), t-100]);
        %     r2 = randi([max(2,t-200), t-100]);
        %     gamma = 2.38/sqrt(2*d);
        %     Y = X + gamma*(xsto(r1,:) - xsto(r2,:)) + 0.01*randn(1,d);
        % end

        if ~isempty(inds), Y(inds) = vals; end
    end

    history(t,1:d) = Y;

    FY = F(Y);
    if ~isscalar(FY) || ~isfinite(FY)
        error('Target F returned non-scalar or non-finite at t=%d (size=%s).', t, mat2str(size(FY)));
    end

    dlog = FY - FX;
    if dlog >= 0 || log(rand) < dlog
        xsel = Y; FX = FY; acc = acc + 1; history(t,end) = 1;
    else
        xsel = xsto(t-1,:);
    end
    xsto(t,:)  = xsel;
    outsto(t)  = FX;

    if displ && (mod(t,round(max(10,n/25)))==0), fprintf('%0.5g ', t/n*25); end
end

accept_rate = acc/n;
covmat = cov(xsto);   % final empirical cov (n x d input)
end
