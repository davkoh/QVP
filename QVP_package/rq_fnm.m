function b = rq_fnm(X, y, p)
% Construct the dual problem of quantile regression
% Solve it with lp_fnm
%
%
[m n] = size(X);
u = ones(m, 1);
a = (1 - p) .* u;
b = -lp_fnm_qr(X', -y', X' * a, u, a)';
end