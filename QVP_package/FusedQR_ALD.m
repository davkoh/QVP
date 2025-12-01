function [bhat,full_draw,pomega,time]=FusedQR_ALD(ysmall,Xsmall,tau,iter,burn,chain_i)
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

% Pre-compute H'*H
HTH = H'*H;

% Pre-allocate the new diagonal of Sigma_beta
sigma = zeros(Q*K,1);

a0 = 0.1;
b0 = 0.1;

%% Priors

%sig_alpha = 2;
%mean_alpha = 0;


%% Storage Matrices & Initialisation
% Initialisation
if chain_i==1
    omega=ones(T,Q);
    beta0 = zeros(K,1);
    bigalpha = sparse(zeros(Q*T,1));
    mu = ones(T*Q,1);

        % Scale of ALD
    sigma_ald = ones(Q,1);

    % Locals for level
    lam_local = ones(K,1);
    eta1_beta = ones(K,1);
    % Global for level
    nu_global = 1;
    eta2_beta = 1;
    % Locals for difference
    lam_diff_local = ones(K*(Q),1);
    eta1_diff_beta = ones(K*(Q),1);
    % Global for differences
     nu_diff_global  = ones(Q,1); %nu_diff_global = 1;
     eta2_diff_beta = ones(Q,1); %eta2_diff_beta =1;
    %Sig_beta = sparse(1:K*(Q),1:K*(Q),[lam_diff_local*nu_diff_global]);
else
    mu = (4*(rand(T*Q,1)-0.5));
    omega=exp(4*(rand(T,Q)-0.5));
    beta0 = zeros(K,1);
    bigalpha = sparse(exp(4*(rand(Q*T,1)-0.5)));

            % Scale of ALD
    sigma_ald = exp(4*(rand(Q,1)-0.5));

    % Locals for level
    lam_local = exp(4*(rand(K,1)-0.5));
    eta1_beta = exp(4*(rand(K,1)-0.5));
    % Global for level
    nu_global = exp(4*(rand(1,1)-0.5));
    eta2_beta = exp(4*(rand(1,1)-0.5));
    % Locals for difference
    lam_diff_local = exp(4*(rand(K*(Q),1)-0.5));
    eta1_diff_beta = exp((4*rand(K*(Q),1)-0.5));
    % Global for differences
    nu_diff_global = exp(4*(rand(Q,1)-0.5));%nu_diff_global = exp(4*(rand(1,1)-0.5));
    eta2_diff_beta = exp(4*(rand(Q,1)-0.5));%eta2_diff_beta = exp(4*(rand(1,1)-0.5));
    %Sig_beta = sparse(1:K*(Q),1:K*(Q),[lam_diff_local*nu_diff_global]);
end

mu_tau = zeros(T*Q,1); % storage for deterministic part of mixture
tau2 = zeros(Q,1); 
theta = zeros(Q,1);
for i = 1:length(tau)
    mu_tau((i*T-T +1)  :i*T,:) = (1-2*tau(i))/(tau(i)*(1-tau(i)));
    tau2(i) = 2/(tau(i)*(1-tau(i)));
    theta(i) = (1-2*tau(i))/(tau(i)*(1-tau(i)));
end
mu = sparse(mu.*mu_tau);

iOm= sparse(1:T*Q,1:T*Q,1./vec(omega));
%beta = ones(K*Q,1);
%Omega = diag(vec(omega));
%iSigm_beta = sparse(1:K*Q,1:K*Q,ones(K*Q,1));

% Storage
beta_store = zeros(K*Q,mcmc);
beta0_store = zeros(K,mcmc);
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
        beta_hat = Kbeta\(HTH*kron(ones(Q,1),beta0) + bigX'*iOm*ystar);
    catch
        beta_hat = pinv(full(Kbeta))*(HTH*kron(ones(Q,1),beta0) + bigX'*iOm*ystar); %Incase Kbeta is badly scaled.
    end
    beta = beta_hat + Cbeta'\randn(K*Q,1);

    del_beta = [beta(1:K) - beta0;beta(1+K : end) - beta(1: end - K)];
    %del_beta = [beta(1:K) - beta0;beta(1+K : end-K) - beta(1: end - K*2);beta(end-K+1:end)];

    % Sample Hyperparameters
    % Difference Parameters
    for qq = 1:Q
    % sample Local
    lam_diff_local(qq*K-K+1:qq*K) = 1./gamrnd(1, 1./( 1./eta1_diff_beta(qq*K-K+1:qq*K) + 0.5*del_beta(qq*K-K+1:qq*K).^2/nu_diff_global(qq)));
    % sample Global
    nu_diff_global(qq) = 1/gamrnd( 0.5*K, 1/( 1/eta2_diff_beta(qq) + 0.5*sum(sum(del_beta(qq*K-K+1:qq*K).^2./lam_diff_local(qq*K-K+1:qq*K))))); %nu_diff_global = 1/gamrnd( 0.5*(K*(Q)), 1/( 1/eta2_diff_beta + 0.5*sum(sum(del_beta.^2./lam_diff_local))  ) );
    % samplel mixture
    eta2_diff_beta(qq) = 1/gamrnd(1, 1/( 1/sqrt(T) + 1/(nu_diff_global(qq)) ));%eta2_diff_beta = 1/gamrnd(1, 1/( 1/sqrt(T*Q) + 1/nu_diff_global )); %eta2_diff_beta = 1/gamrnd(1, 1/( 1/(1 + nu_diff_global) ));%eta2_diff_beta = 1/gamrnd(1, 1/( 1/sqrt(T*Q) + 1/nu_diff_global ));
    % Expand to new diagonal
    sigma(qq*K-K+1 : qq*K) = lam_diff_local(qq*K-K+1 : qq*K)*nu_diff_global(qq);
    end


    % sample mixture local
    eta1_diff_beta = 1./gamrnd(1, 1./(1 + 1./lam_diff_local));
    

    iSigm_beta = sparse(1:K*(Q),1:K*(Q),[1./(sigma)]);
    %iSigm_beta = sparse(1:K*(Q),1:K*(Q),[1./(1000*ones(K*Q,1))]);
    HTH = H'*iSigm_beta*H;

    % Sample beta0
    beta0 = zeros(K,1);

    % Sample alpha
    ystar = Y - bigX*beta - mu;
    for j = 1:Q
        Kalpha = ( sum(1./(omega(:,j)*sigma_ald(j)^1/1*tau2(j)).^(1/1)));
        alpha_hat = Kalpha\(sum(sparse(1:T,1:T,1./(omega(:,j)*sigma_ald(j)^1/1*tau2(j)).^(1/1))*ystar(j*T-T+1:j*T)));
        bigalpha(j*T-T+1 : j*T) = alpha_hat + sqrt(Kalpha)\randn;
    end

        % Sample Omega
    a = Y - bigX*beta - bigalpha;
    for j = 1:Q
    for t = 1:T
       omega(t,j) = max(gigrnd(1/2,2/sigma_ald(j)+theta(j)^2/(tau2(j)*sigma_ald(j)),a(t + (j-1)*T)^2/(tau2(j)*sigma_ald(j)),1),1e-4);
    end
    end

    mu = mu_tau.*vec(omega);

        % Sample scale of the ALD 
    a = Y - bigX*beta - bigalpha - mu;
    for j = 1 : Q
    a1 = a0 + 3*T/2;
    sse = (a(j*T-T+1:j*T)).^2;
    a2 = b0 + sum(sse./(2*omega(:,j)*tau2(j))) + sum(omega(:,j));
    sigma_ald(j) = 1./gamrnd(a1,1./a2);
    end

    iOm= sparse(1:T*Q,1:T*Q,1./vec(((omega'.*sigma_ald.*tau2)')));

    if i> burn
        isave = i - burn;
        beta_store(:,isave) = beta;
        beta0_store(:,isave) = beta0;
        alpha_store(:,isave) = bigalpha;
        omega_store(:,isave) = vec(omega);
    end

end
time=toc

pbeta = mean(beta_store,2);
palpha = mean(alpha_store,2);
pomega = mean(omega_store,2);
pbeta0 = mean(beta0_store,2);

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