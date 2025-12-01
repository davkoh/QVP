function [bhat,full_draw,pomega,time] = HS(ysmall,Xsmall,tau,iter,burn,chain_i)
if nargin<6
    chain_i=1;
end

%% Collect terms needed
mcmc = iter-burn;
tau=unique(sort(tau));
Q = length(tau);
[T, K]=size(Xsmall); %MAKE SURE NO CONSTANT IN X!
dim_mat=[Q,K,T];

%% Bring to System form
% Large Form Data
Y = sparse(repmat(ysmall,Q,1));
bigX = sparse(kron(eye(Q),Xsmall));

% Difference matrix
H = speye(K*Q,K*Q) - sparse(K+1:K*Q,1:(Q-1)*K,...
    ones(1,(Q-1)*K),K*Q,K*Q); %% This would be the difference matrix that I need!

%% Storage Matrices & Initialisation
% Initialisation
if chain_i==1
    %beta = ones(K*Q,1);
    mu = ones(T*Q,1);
    omega=ones(T,Q);
    bigalpha = sparse(zeros(Q*T,1));

    % Local level Paremeters
    lam_local = ones(K*Q,1);
    eta1_beta = ones(K*Q,1);
    % Global level Parameters
    nu_global = ones(Q,1);
    eta2_beta = ones(Q,1);
else
    %beta = ones(K*Q,1);
    mu = (4*(rand(T*Q,1)-0.5));
    omega=exp(4*(rand(T,Q)-0.5));
    bigalpha = sparse(exp(4*(rand(Q*T,1)-0.5)));

    % Local level Paremeters
    lam_local = exp(4*(rand(K*Q,1)-0.5));
    eta1_beta = exp(4*(rand(K*Q,1)-0.5));
    % Global level Parameters
    nu_global = exp(4*(rand(Q,1)-0.5));
    eta2_beta = exp(4*(rand(Q,1)-0.5));
end

%Omega = diag(vec(omega));
mu_tau = zeros(T*Q,1); % storage for deterministic part of mixture
% Pre-compute H'*H
for j = 1:Q
    lam(j*K-K + 1 : j*K) = lam_local(j*K-K + 1 : j*K) * nu_global(j);
end
iSigm_beta=sparse(1:K*(Q),1:K*(Q),1./(lam));
HTH = iSigm_beta;

for i = 1:length(tau)
    mu_tau((i*T-T +1)  :i*T,:) = (1-2*tau(i));
end
mu = sparse(mu.*mu_tau);

iOm= sparse(1:T*Q,1:T*Q,1./vec(omega));

% Storage
beta_store = zeros(K*Q,mcmc);
alpha_store = zeros(Q*T,mcmc);
omega_store = zeros(T*Q,mcmc);

tic
for i = 1:iter
    if rem(i,1000)==0
        if i~=iter
            fprintf("Iteration %i out of %i \n", i,iter);
        end
    end

    % Compute y* = y - X*beta - alpha
    ystar = Y - bigalpha - mu;
    % Sample beta
    Kbeta = (HTH + bigX'*iOm*bigX);

    try
        Cbeta = chol(Kbeta,'lower');
    catch
        Cbeta = chol(nearestSPD(full(Kbeta)),'lower'); %SVD does not support sparse matrix so convert Kbeta to normal first
        Cbeta = sparse(Cbeta); %Convert back to sparse matrix
    end
    try
        beta_hat = Kbeta\(bigX'*iOm*ystar);
    catch
        beta_hat = pinv(full(Kbeta))*(bigX'*iOm*ystar); %Incase Kbeta is badly scaled.
    end
    beta = beta_hat + Cbeta'\randn(K*Q,1);

    % Sample Hyperparameters
    for j = 1:Q
        % Level Parameters
        % sample Local
        lam_local(j*K-K + 1 : j*K) = 1./gamrnd(1, 1./( 1./eta1_beta(j*K-K + 1 : j*K) + 0.5*beta(j*K-K + 1 : j*K).^2/nu_global(j)));
        % sample Global
        nu_global(j) = 1/gamrnd( 0.5*(K), 1/( 1/eta2_beta(j) + 0.5*sum(sum(beta(j*K-K + 1 : j*K).^2./lam_local(j*K-K + 1 : j*K)))  ) );
        % sample mixture local
        eta1_beta(j*K-K + 1 : j*K) = 1./gamrnd(1, 1./(1 + 1./lam_local(j*K-K + 1 : j*K)));
        % samplel mixture
        eta2_beta(j) = 1/gamrnd(1, 1/( 1/sqrt(T) + 1/nu_global(j) )); %Orig: 1/gamrnd(1, 1/( 1 + 1/nu_global(j) ));

        lam(j*K-K + 1 : j*K) = lam_local(j*K-K + 1 : j*K) * nu_global(j);
    end
    iSigm_beta = sparse(1:K*(Q),1:K*(Q),1./(lam));
    HTH = iSigm_beta;

    % Sample alpha
    ystar = Y - bigX*beta - mu;
    for j = 1:Q
        Kalpha = (sum(1./omega(:,j)));
        alpha_hat = Kalpha\(sum(sparse(1:T,1:T,1./omega(:,j))*ystar(j*T-T+1:j*T)));
        bigalpha(j*T-T+1 : j*T) = alpha_hat + sqrt(Kalpha)\randn;
    end

    % Sample Omega
    a = Y - bigX*beta - bigalpha;
    for j = 1:Q
        for t = 1:T
            omega(t,j) =  gigrnd(1,1,abs(a(t + (j-1)*T)),1);
        end
    end

    iOm= sparse(1:T*Q,1:T*Q,1./vec(omega));
    mu = mu_tau.*vec(omega);

    if i> burn
        isave = i - burn;
        beta_store(:,isave) = beta;
        alpha_store(:,isave) = bigalpha;
        omega_store(:,isave) = vec(omega);
    end

end
time=toc

%Collect terms
pbeta = mean(beta_store,2);
palpha = mean(alpha_store,2);
pomega = mean(omega_store,2);

palpha = palpha(1:T:Q*T);
pbeta = reshape(pbeta,[K,Q]);
pomega = reshape(pomega,[T,Q]);

bhat=[palpha';pbeta];

full_draw=nan([size(bhat),size(beta_store,2)]);
for i=1:size(beta_store,2)
    temp_beta=beta_store(:,i);
    temp_beta=reshape(temp_beta,[K,Q]);
    temp_alpha=alpha_store(:,i);
    temp_alpha=temp_alpha(1:T:Q*T);

    full_draw(:,:,i)=[temp_alpha';temp_beta];
end

end