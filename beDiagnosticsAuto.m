function Diag = beDiagnosticsAuto(ydata, X, data, pointEst, Xfull, coeff, paramNames)
% Parameter-agnostic diagnostics around a provided point estimate.
% Returns H = J'*J, pcov (Gauss-Newton), pcov_robust (sandwich), etc.

%% -------------------- User-editable bounds (ONE place) --------------------
K  = numel(pointEst);
% Example for current 4-param model: [alpha, phi, kstar, k3]
lb = [0,  -25, -25];   % length must be K
ub = [1,   25,  25];
assert(numel(lb)==K && numel(ub)==K, 'lb/ub must have length equal to numel(pointEst).');
%% -------------------------------------------------------------------------

if nargin < 7 || isempty(paramNames)
    paramNames = arrayfun(@(i)sprintf('p%d',i), 1:K, 'uni',0);
end

% Housekeeping
hlag = 0; projection = 0; plotRun = 0;

% Time grid + slicing (matches your convention)
if projection==1
    tvec  = [1,2,61,93,134,141,148,155,162,169,176,186,200,211,218,223,227,236,250,258,266,271,279,294,310,322,330,338,349,370,384,397,407,418,433,445,463,468,474,491,504,517,540,561,567,575];
    xdata = 85:tvec(end);
else
    tvec  = [1,2,61,94,[127,134,141,148,153,155,162,167,169,175,176,186,200,211,216,218,223,227,230,236,250,258,265,266,271,279,288,294,305,310,322,330,337,338,349,354,356,361,370,372,384,397,407,418,433,445,454,463,468,474,491,503,504,517,540,561,566,567,575]+hlag];
    xdata = 85:tvec(end-7);
end
lt    = length(tvec);
X     = X(:,1:lt-1);
Xfull = Xfull(:,1:lt-1);
[~,lx2] = size(X);

% Data (England scaling once)
ydata = ydata(1:length(xdata));
ydata = ydata * (sum(data.Npop)/56286961);

% Model handle for lsqcurvefit (column output guaranteed)
fun = @(params, x) sim2fit(params, data, x, X, 1, Xfull, coeff, tvec, size(X,1), lx2, plotRun);

% Local refit at the provided point estimate to get a clean Jacobian
p0   = pointEst(:).';
opts = optimoptions('lsqcurvefit', ...
    'FiniteDifferenceType','central', ...
    'OptimalityTolerance',1e-9, ...
    'StepTolerance',1e-12, ...
    'FunctionTolerance',1e-12, ...
    'MaxFunctionEvaluations',2e5, ...
    'Display','off', ...
    'SpecifyObjectiveGradient',false);

[p_hat,~,res,~,~,~,J] = lsqcurvefit(fun, p0, xdata, ydata', lb, ub, opts); % note ydata'
J = full(J);

% Core diagnostics
N   = numel(res);
k   = numel(p_hat);
dof = max(N - k, 1);
rss = sum(res.^2);
s2  = rss / dof;

% SVD for conditioning / ridge direction
[U,S,V]  = svd(J,'econ'); %#ok<ASGLU>
s        = diag(S);
condJ    = s(1) / s(end);
v_ridge  = V(:,end);
frac_info = s(end)^2 / sum(s.^2);

% Column correlation & VIFs
Z  = zscore(J,0,1);
R  = corrcoef(Z);
if rcond(R) < 1e-10
    VIF = diag(pinv(R));
else
    VIF = diag(inv(R));
end
maxAbsCorr = max(max(abs(R - diag(diag(R)))));


% Gauss–Newton Hessian and covariance
H     = J.'*J;                 % Hessian approx of 1/2 * ||res||^2
pcov  = s2 * pinv(H);          % covariance of params under LS assumptions
se    = sqrt(max(diag(pcov), 0));
tcrit = tinv(0.975, dof);
CI    = [p_hat(:) - tcrit.*se(:),  p_hat(:) + tcrit.*se(:)];

% Robust (sandwich) covariance
Smeat = zeros(k);
for i=1:N
    Ji = J(i,:).';
    Smeat = Smeat + (res(i)^2) * (Ji*Ji.');
end
pcov_robust = pinv(H) * Smeat * pinv(H);

% Print concise tables
T_ridge = table(paramNames(:), v_ridge, 'VariableNames', {'param','ridge_loading'});
disp(T_ridge);
fprintf('cond(J) = %.2e,   worst-dir info fraction = %.3g\n', condJ, frac_info);

T_vif = array2table([diag(R), VIF], 'RowNames', paramNames, ...
                    'VariableNames', {'selfcorr','VIF'});
disp('Max |corr| among columns:'), disp(maxAbsCorr);
disp(T_vif);

T_ci = table((1:k).', p_hat(:), se(:), CI(:,1), CI(:,2), ...
    'VariableNames', {'index','estimate','se','ci_lo','ci_hi'});
disp(T_ci);

% Pack outputs
Diag = struct();
Diag.p_hat       = p_hat(:).';
Diag.residuals   = res(:);
Diag.rss         = rss;
Diag.sigma2      = s2;
Diag.J           = J;
Diag.S           = s;
Diag.V           = V;
Diag.H           = H;              % <— Gauss–Newton Hessian
Diag.pcov        = pcov;           % <— covariance (LS assumptions)
Diag.pcov_robust = pcov_robust;    % <— robust covariance
Diag.se          = se(:);
Diag.CI          = CI;
Diag.paramNames  = paramNames(:);
Diag.condJ       = condJ;
Diag.v_ridge     = v_ridge(:);
Diag.frac_info   = frac_info;
Diag.R           = R;
Diag.maxAbsCorr  = maxAbsCorr;
Diag.VIF         = VIF(:);
Diag.dof         = dof;
Diag.tcrit       = tcrit;
Diag.lb          = lb(:).';
Diag.ub          = ub(:).';
Diag.xdata       = xdata(:);
Diag.ydata       = ydata(:);
end


function [f,rhohat]=sim2fit(params,data,xdata,Xfit,intrinsic,Xfull,coeff,tvec,lx1,lx2,plotRun)
R0=2.2;
tvec(1)=-145;
alpha=params([1,1,1]);
a1=-.814;
b1=8.0161;
a2=-.8067;
b2=-8.0887;
ks=params(2);%ksrat0=0.1613
reducedParams=[1,a1*ks+b1,a2*ks+b2,params(3),0];
[pr,be,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=bePrepCovid19(data,R0,ones(1,lx2-2),reducedParams,coeff,zeros(5,lx2),alpha);
pr.leak=1;%params(2);
pr.xfull=Xfull;
Wfit=Xfit.^(1/pr.a);
if intrinsic==1
    [simu,simu2,~,rhohat]=beRunCovid19(pr,be,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,[ones(1,length(tvec)-1)],tvec(1:lx2+1),plotRun,data);
else

    [simu,simu2,~,rhohat]=beRunCovid19(pr,be,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,Wfit,tvec(1:lx2+1),0,data);
end
t=simu(:,1)';
h=simu2';
f=interp1(t,h,xdata); 
end

function [a,b,kstar_hat,u,ab_info] = collapse_k1k2(J, idx_k1, idx_k2, v1, v2, p_hat)
% COLLAPSE_K1K2  Find the most identifiable linear combo k* = a k1 + b k2
% using the Jacobian columns for k1,k2, and build the composite index
% u_t = a v1_t + b v2_t (RMS-normalized so its overall scale is 1).
%
% Inputs:
%   J        : N x K Jacobian from lsqcurvefit at the refit optimum
%   idx_k1   : column index of k1 inside the parameter vector
%   idx_k2   : column index of k2 inside the parameter vector
%   v1, v2   : vectors (or same-length columns) of the two indices over time
%   p_hat    : 1 x K vector of the refit parameter estimates
%
% Outputs:
%   a,b         : coefficients for the composite index (RMS-normalized u)
%   kstar_hat   : suggested value a*k1_hat + b*k2_hat
%   u           : composite index, u = a v1 + b v2, RMS(u)=1
%   ab_info     : diagnostics (G, eigenvalues, ratio, and the dropped dir)

    % --- 1) Build the 2-col sensitivity block and its Gram matrix
    J12 = J(:, [idx_k1, idx_k2]);                 % N x 2
    G   = J12.' * J12;                             % 2 x 2, Fisher in (k1,k2)

    % --- 2) Top eigenvector of G gives the best-informed linear combo
    [V,D] = eig((G+G')/2);
    [lam,ix] = max(diag(D));
    ab = V(:,ix);                                  % 2x1, up to scale
    ab = ab / norm(ab);                            % unit Euclidean norm

    % --- 3) Build the composite index and RMS-normalize its scale
    u_raw = ab(1)*v1(:) + ab(2)*v2(:);             % T x 1
    s = sqrt(mean(u_raw.^2));                      % RMS
    if s > 0
        ab = ab / s;                               % fold scale into (a,b)
        u  = u_raw / s;                            % so RMS(u)=1
    else
        u  = u_raw;                                % degenerate fallback
    end
    a = ab(1); b = ab(2);

    % --- 4) Suggested single-parameter value at the optimum
    kstar_hat = a * p_hat(idx_k1) + b * p_hat(idx_k2);

    % --- 5) Diagnostics (strength of kept vs dropped directions)
    ab_perp = [-ab(2); ab(1)];                     % orthogonal in R^2
    info_keep = ab.' * G * ab;
    info_drop = ab_perp.' * G * ab_perp;
    ab_info = struct('G',G, 'lambda_max',lam, ...
                     'info_keep',info_keep, 'info_drop',info_drop, ...
                     'info_ratio', info_keep / max(info_drop, eps), ...
                     'ab_perp', ab_perp);
end

function [a, b, kstar_hat] = collapse_k12_from_hat(k1_hat, k2_hat, a_guess)
% a_guess = [a1; a2] (can be your ridge direction; need not be unit length)
a = a_guess(:);
khat = [k1_hat; k2_hat];
kstar_hat = (a' * khat) / (a' * a);
b = khat - a * kstar_hat;  % b = [b1; b2]
end

%{
[a1,a2,kstar_hat,u,ab_info] = collapse_k1k2(J, idx_k1, idx_k2, v1, v2, p_hat);
%disp([a1,a2,kstar_hat])
%disp(ab_info.info_ratio)   % >>  how much more informed the kept combo is

[k1_hat, k2_hat] = deal(pointEst(3), pointEst(4));   % from your last good fit
[a, b, kstar_hat] = collapse_k12_from_hat(k1_hat, k2_hat, [a1,a2]);
a1 = a(1); a2 = a(2);  b1 = b(1); b2 = b(2);

disp([a1,b1,a2,b2,kstar_hat])
disp(ab_info.info_ratio)   % >>  how much more informed the kept combo is
%}