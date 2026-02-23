function [full_draw,bhat,qfit,eps,full_draw_SAVS,bhat_SAVS,qfit_SAVS,eps_SAVS,Z_fit] = BQVAR_Estim_plag(Y,taus,eq_i,iter,burn,Prior_choice,p)
% BQVAR_Estim_plag  Estimate a single BQVAR equation with p lags.
%
%   Generalises BQVAR_Estim to an arbitrary lag order p.
%   The regressor matrix Z contains lags 1..p of all K variables,
%   plus contemporaneous values of equations 1..(eq_i-1) for
%   triangular identification.
%
%   Inputs:
%       Y            (T x K)   data matrix
%       taus         (1 x Q)   quantile grid
%       eq_i         scalar    equation index (1..K)
%       iter         scalar    total MCMC iterations
%       burn         scalar    burn-in iterations
%       Prior_choice string    prior label passed to QVP_wrapper
%       p            scalar    number of lags (default = 1)
%
%   Outputs:  same as BQVAR_Estim

if nargin < 7 || isempty(p)
    p = 1;
end

[T, K] = size(Y);
Teff   = T - p;                        % effective sample size after losing p obs

% --- Build lag matrix (Teff x K*p) ---
Z = [];
for lag = 1:p
    Z = [Z, Y(p+1-lag:end-lag, :)];    % lag 1, lag 2, ..., lag p
end

% --- Contemporaneous regressors for triangular identification ---
if eq_i > 1
    for j = 1:(eq_i - 1)
        Z = [Z, Y(p+1:end, j)];
    end
end

% Design matrix including intercept
Z_fit = [ones(Teff, 1), Z];

% Dependent variable
y = Y(p+1:end, eq_i);

% Estimate
[bhat, full_draw, bhat_SAVS, full_draw_SAVS, time] = ...
    QVP_wrapper(y, Z, taus, iter, burn, Prior_choice);

% In-sample fit and residuals
qfit      = Z_fit * bhat;
qfit_SAVS = Z_fit * bhat_SAVS;
eps       = y - qfit;
eps_SAVS  = y - qfit_SAVS;

end
