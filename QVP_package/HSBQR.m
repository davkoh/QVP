function [bhat,full_draw,time]=HSBQR(y,X,quant,iter,nburn,chain_i)
% HSBQR is a function that applies the Horseshoe prior to the Bayesian
% Quantile Regression. Please cite the authors (David Kohns and Tibor
% Szendrei) paper (found in zipped file) in your work.The function
% generates a vector of betas which is an average for all draws.

% This code is free to use for academic purposes only, provided that the
% paper is cited as:
%
% Kohns, D.E. and Szendrei, T. (2024). Horseshoe Prior Bayesian Quantile
% Regression, JRSSC
%
% This code comes without technical support of any kind. It is expected to
% reproduce the results reported in the paper. Under no circumstances will
% the authors be held responsible for any use (or misuse) of this code in
% any way.

%
% The input arguments are the following:
% 1. y=LHS variable
% 2. X=Matrix of explanatory variables
% 3. quantile=set the quantile (between 0 and 1)
% 4. inter=Total Number of simulations
% 5. nburn=Number of burn-in runs
% 6: chain_i=chain we are currently on


nsave = iter - nburn;   % Number of total draws
T=size(X,1);
y = y;
x = X;
[n,p] = size(x);
n_q=length(quant);

% ==============| Define priors
V = 9*eye(p);
Vinv = inv(V);

%% output %%
betaout=zeros(p,nsave);
lambdaout=zeros(p,nsave);
tauout=zeros(nsave,1);
sigmaSqout=zeros(nsave,1);

%% matrices %%
I_n=eye(n);
l=ones(n,1);

% prior for sigma2 ~ IG(a0,b0)
a0 = 0.1;
b0 = 0.1;

% ==============| Initialize vectors
if chain_i==1
    beta = zeros(p,n_q);
    z = ones(n,n_q);
    sigma = ones(1,n_q);
    theta = zeros(1,n_q);
    tau_sq = zeros(1,n_q);

    plambda=ones(p,n_q);
    ptau=ones(1,n_q);
else
    beta = -2 + 4*rand(p,n_q);
    z = exp(-2 + 4*rand(n,n_q));
    sigma = exp(-2 + 4*rand(n_q,1))';
    theta = zeros(1,n_q);
    tau_sq = zeros(1,n_q);

    plambda= exp(-2 + 4*rand(p,n_q));
    ptau= exp(-2 + 4*rand(1,n_q));
end

% ==============| Storage matrices
beta_draws = zeros(p,n_q,nsave);

tic
for irep = 1:iter
    % Print every "iter" iterations on the screen
    if rem(irep,1000)==0
        if irep~=iter
            fprintf("Iteration %i out of %i \n", irep,iter);
        end
    end
    for q= 1:n_q
        tau_sq(:,q) = 2/(quant(q)*(1-quant(q)));
        theta(:,q) = (1-2*quant(q))/(quant(q)*(1-quant(q)));
        lambda=plambda(:,q);
        tau=ptau(:,q);


        %% HS prior Implementation

        % Sample regression coefficients beta
        U = diag( 1./(sigma(1,q).*tau_sq(:,q).*z(:,q)) );
        try
            y_tilde = chol(U)*(y-theta(:,q)*z(:,q));
            X_tilde= chol(U)*x;
        catch
            y_tilde = chol(nearestSPD(U))*(y-theta(:,q)*z(:,q));
            X_tilde = chol(nearestSPD(U))*x;
        end

        lambda_star=tau*lambda;
        U_bar=bsxfun(@times,(lambda_star.^2),X_tilde');

        u=normrnd(0,lambda_star);
        v=X_tilde*u+normrnd(0,l);

        v_star=(X_tilde*U_bar+I_n)\(y_tilde-v);
        beta(:,q)=(u+U_bar*v_star);

        %% update lambda_j's in a block using slice sampling %%
        eta = 1./(lambda.^2);
        upsi = unifrnd(0,1./(1+eta));
        tempps = beta(:,q).^2/(2*tau^2);
        ub = (1-upsi)./upsi;

        % now sample eta from exp(tempv) truncated between 0 & upsi/(1-upsi)
        Fub = 1 - exp(-tempps.*ub); % exp cdf at ub
        Fub(Fub < (1e-4)) = 1e-4;  % for numerical stability
        up = unifrnd(0,Fub);
        eta = -log(1-up)./tempps;
        lambda = 1./sqrt(eta);

        %% update tau %%
        tempt = sum((beta(:,q)./lambda).^2)/(2);
        et = 1/tau^2;
        utau = unifrnd(0,1/(1+et));
        ubt = (1-utau)/utau;
        Fubt = gamcdf(ubt,(p+1)/2,1/tempt);
        Fubt = max(Fubt,1e-8); % for numerical stability
        ut = unifrnd(0,Fubt);
        et = gaminv(ut,(p+1)/2,1/tempt);
        tau = 1/sqrt(et);

        % Sample regression variance sigma2
        a1 = a0 + 3*n/2;
        sse = (y-x*beta(:,q) - theta(:,q)*z(:,q)).^2;
        a2 = b0 + sum(sse./(2*z(:,q)*tau_sq(:,q))) + sum(z(:,q));
        sigma(1,q) = 1./gamrnd(a1,1./a2);

        % Sample latent variables z_{t}
        for t = 1:n
            k1 = sqrt(theta(:,q).^2 + 2*tau_sq(:,q))/abs(y(t,:)-x(t,:)*beta(:,q));
            k2 = (theta(:,q).^2 + 2*tau_sq(:,q))/(sigma(1,q)*tau_sq(:,q));
            z(t,q) = max(1./Draw_IG(k1,k2),1e-4);
        end
        plambda(:,q)=lambda;
        ptau(:,q)=tau;
    end
    %end

    if irep > nburn
        beta_draws(:,:,irep-nburn) = beta;
    end
end
time=toc

full_draw=beta_draws;
bhat=squeeze(mean(beta_draws,3));
end
