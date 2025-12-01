function [bhat,bhat_unspars,full_draw_sparse,full_draw,pomega,time] = LQVP_SAVS_ALD(ysmall,Xsmall,tau,iter,burn,chain_i)
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

% Priors
a0 = 0.1;
b0 = 0.1;

%% Storage Matrices & Initialisation
% Initialisation
if chain_i==1
    omega=ones(T,Q);
    %beta = zeros(K,1);
    bigalpha = sparse(zeros(Q*T,1));
    %Omega = diag(vec(omega));
    mu = ones(T*Q,1);

            % Scale of ALD
    sigma_ald = ones(Q,1);

    % Prior for Quantile invariant coefficients
    % Locals for level
    lam_level = ones(K,1);
    eta1_level= ones(K,1);
    % Global for level
    nu_level = 1;
    eta2_level = 1;

    % Prior for Non-centred differences
    % Locals for initial conditions
    lam_init_diff = ones(K,1);
    eta1_init_diff = ones(K,1);
    % Global for level
    nu_init_diff = 1;
    eta2_init_diff = 1;
    % Locals for difference
    lam_diff= ones(K*(Q),1);
    eta1_diff = ones(K*(Q),1);
    % Global for differences
    nu_diff = ones(Q,1);
    eta2_diff =ones(Q,1);
    sigma = ones(K*(Q),1);

    V = zeros(K*Q,1);
else
    omega=exp(4*(rand(T,Q)-0.5));
    %beta = zeros(K,1);
    bigalpha = sparse(4*(rand(Q*T,1)-0.5));
    %Omega = diag(vec(omega));
    mu = (4*(rand(T*Q,1)-0.5));

                % Scale of ALD
    sigma_ald = exp(4*(rand(Q,1)-0.5));

    % Prior for Quantile invariant coefficients
    % Locals for level
    lam_level = exp(4*(rand(K,1)-0.5));
    eta1_level= exp(4*(rand(K,1)-0.5));
    % Global for level
    nu_level = exp(4*(rand(1,1)-0.5));
    eta2_level = exp(4*(rand(1,1)-0.5));

    % Prior for Non-centred differences
    % Locals for initial conditions
    lam_init_diff = exp(4*(rand(K,1)-0.5));
    eta1_init_diff = exp(4*(rand(K,1)-0.5));
    % Global for level
    nu_init_diff = exp(4*(rand(1,1)-0.5));
    eta2_init_diff = exp(4*(rand(1,1)-0.5));
    % Locals for difference
    lam_diff= exp(4*(rand(K*(Q),1)-0.5));
    eta1_diff = exp(4*(rand(K*(Q),1)-0.5));
    % Global for differences
    nu_diff = exp(4*(rand(Q,1)-0.5));
    eta2_diff =exp(4*(rand(Q,1)-0.5));
    sigma = exp(4*(rand(K*Q,1)-0.5));

    V = 4*(rand(K*Q,1)-0.5);
end

%iSig_beta_level = sparse(1:K,1:K,ones(K,1));
iSig_beta_level = sparse(1:K,1:K,1./(lam_level*nu_level));
iSig_beta_diff = sparse(1:K*(Q),1:K*(Q),1./(sigma));
bigX_tilde = bigX*sparse(diag(sigma).^2);

beta = zeros(K*Q,1);

bigX_star = bigX;

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

% Storage
beta_level_store = zeros(K,mcmc);
beta_store = zeros(K*Q,mcmc);
beta_savs_store = zeros(K,mcmc);
beta_fullsavs_store = zeros(K*Q,mcmc);
V_store = zeros(K*Q,mcmc);
alpha_store = zeros(Q*T,mcmc);
omega_store = zeros(T*Q,mcmc);

tic
for i = 1:iter
    if rem(i,1000)==0
        if i~=iter
            fprintf("Iteration %i out of %i \n", i,iter);
        end
    end
    %% Update Level beta
    ystar = Y - bigalpha - mu - bigX_star*V;
    Mom2 = iSig_beta_level;
    Mom1 = zeros(K,1);
    for j = 1:Q
        ybar = ystar(j*T-T+1:j*T);
        % Update Second Moment
        Mom2 = Xsmall'*diag(1./(omega(:,j)*sigma_ald(j)*tau2(j)))*Xsmall + Mom2;
        % Update First Moment
        Mom1 = Xsmall'*diag(1./(omega(:,j)*sigma_ald(j)*tau2(j)))*ybar + Mom1;
    end
    Kbeta = Mom2;
    try
        Cbeta = chol(Kbeta,'lower');
    catch
        Cbeta = chol(nearestSPD(full(Kbeta)),'lower'); %SVD does not support sparse matrix so convert Kbeta to normal first
        Cbeta = sparse(Cbeta); %Convert back to sparse matrix
    end
    beta_hat = Kbeta\Mom1;
    beta_level = beta_hat + Cbeta'\randn(K,1);

    % Sample Hyperparameters for Level beta
    % Level Parameters
    % sample Local
    lam_level = 1./gamrnd(1, 1./( 1./eta1_level + 0.5*beta_level.^2/nu_level));
    % sample Global
    nu_level = 1/gamrnd( 0.5*(K), 1/( 1/eta2_level + 0.5*sum(sum(beta_level.^2./lam_level))  ) );
    % sample mixture local
    eta1_level = 1./gamrnd(1, 1./(1 + 1./lam_level));
    % samplel mixture
    eta2_level = 1/gamrnd(1, 1/( 1/sqrt(T*Q) + 1/nu_level )); %Orig: 1/gamrnd(1, 1/( 1 + 1/nu_level )); 1/gamrnd(1, 1/( 1/sqrt(T*Q) + 1/nu_level ));

    % Assemble inverse Sigma Level
    iSig_beta_level = sparse(1:K,1:K,1./(lam_level*nu_level));


    %% Update Difference Beta

    ystar = Y - bigalpha - mu - bigX*repmat(beta_level,Q,1);

    % Sample beta
    Kbeta = (HTH + bigX_tilde'*iOm*bigX_tilde);
    try
        Cbeta = chol(Kbeta,'lower');
    catch
        Cbeta = chol(nearestSPD(full(Kbeta)),'lower'); %SVD does not support sparse matrix so convert Kbeta to normal first
        Cbeta = sparse(Cbeta); %Convert back to sparse matrix
    end
    try
        beta_hat = Kbeta\(bigX_tilde'*iOm*ystar);
    catch
        beta_hat = pinv(full(Kbeta))*(bigX_tilde'*iOm*ystar); %Incase Kbeta is badly scaled.
    end
    beta = beta_hat + Cbeta'\randn(K*Q,1);

    %% Update Difference Hyperparameters
    bigX_star = bigX.*repmat(beta',T*Q,1);
    % Sample V
    KV = (iSig_beta_diff + bigX_star'*iOm*bigX_star);
    %CV = chol(KV,'lower');
    try
        CV = chol(KV,'lower');
    catch
        CV = chol(nearestSPD(full(KV)),'lower'); %SVD does not support sparse matrix so convert Kbeta to normal first
        CV = sparse(CV); %Convert back to sparse matrix
    end
    %V_hat = KV\(bigX_star'*iOm*ystar);
    %V_hat = pinv(full(KV))*(bigX_star'*iOm*ystar); %Incase Kbeta is badly scaled.
    try
       V_hat = KV\(bigX_star'*iOm*ystar);
    catch
       V_hat = pinv(full(KV))*(bigX_star'*iOm*ystar); %Incase Kbeta is badly scaled.
    end
    V = V_hat + CV'\randn(K*Q,1);

    % Sample Hyperparameters

    % Difference Parameters
    for qq = 1:Q
    % sample Local
    lam_diff(qq*K-K+1:qq*K) = 1./gamrnd(1, 1./( 1./eta1_diff(qq*K-K+1:qq*K) + 0.5*V(qq*K-K+1:qq*K).^2/nu_diff(qq)));
    % sample Global
    nu_diff(qq) = 1/gamrnd( 0.5*(K), 1/( 1/eta2_diff(qq) + 0.5*sum(sum(V(qq*K-K+1:qq*K).^2./lam_diff(qq*K-K+1:qq*K)))  ) );
    % samplel mixture
    eta2_diff(qq) =  1/gamrnd(1, 1/( 1/sqrt(T) + 1/nu_diff(qq) )); % 1/gamrnd(1, 1/( 1/sqrt(T*Q) + 1/nu_diff )); %Orig:
    sigma(qq*K-K+1 : qq*K) = lam_diff(qq*K-K+1 : qq*K)*nu_diff(qq);
    end
    
    % sample mixture local
    eta1_diff = 1./gamrnd(1, 1./(1 + 1./lam_diff));
    
    iSig_beta_diff = sparse(1:K*(Q),1:K*(Q),1./(sigma));

    % Update bigX_tilde
    %Sig_beta_diff = sparse(diag([lam_init_diff*nu_init_diff;lam_diff*nu_diff].^(0.5)));
    bigX_tilde = bigX*diag(V);


    %% Sample alpha
    ystar = Y - bigX*repmat(beta_level,Q,1) - bigX_star*V - mu;
    for j = 1:Q
        Kalpha = ( sum(1./(omega(:,j)*sigma_ald(j)^1/1*tau2(j)).^(1/1)));
        alpha_hat = Kalpha\(sum(sparse(1:T,1:T,1./(omega(:,j)*sigma_ald(j)^1/1*tau2(j)).^(1/1))*ystar(j*T-T+1:j*T)));
        bigalpha(j*T-T+1 : j*T) = alpha_hat + sqrt(Kalpha)\randn;
    end

    %% Sample Omega
    a = Y - bigX*repmat(beta_level,Q,1) - bigX_star*V - bigalpha;
    for j = 1:Q
    for t = 1:T
       omega(t,j) = max(gigrnd(1/2,2/sigma_ald(j)+theta(j)^2/(tau2(j)*sigma_ald(j)),a(t + (j-1)*T)^2/(tau2(j)*sigma_ald(j)),1),1e-4);
    end
    end

    mu = mu_tau.*vec(omega);

        % Sample scale of the ALD 
    a = Y - bigX*repmat(beta_level,Q,1) - bigX_star*V - mu -bigalpha;
    for j = 1 : Q
    a1 = a0 + 3*T/2;
    sse = (a(j*T-T+1:j*T)).^2;
    a2 = b0 + sum(sse./(2*omega(:,j)*tau2(j))) + sum(omega(:,j));
    sigma_ald(j) = 1./gamrnd(a1,1./a2);
    end

    iOm= sparse(1:T*Q,1:T*Q,1./vec(((omega'.*sigma_ald.*tau2)')));

    %% End Loop

    if i> burn
        isave = i - burn;
        beta_level_store(:,isave) = beta_level;
        beta_store(:,isave) = repmat(beta_level,Q,1) + diag(V)*beta;
        beta_savs_store(:,isave) = savs(Xsmall,beta_level,-2,1);
        %beta_fullsavs_store(:,isave) =
        %repmat(savs(Xsmall,beta_level,-2,1),Q,1) +
        %diag(savs(bigX_star,V,-2,1))*beta; This was the original one
        beta_fullsavs_store(:,isave) = repmat(beta_level,Q,1) + diag(savs(bigX_star,V,-2,1))*beta;
        %beta_fullsavs_store(:,isave) = repmat(savs(Xsmall,beta_level,-2,1),Q,1) + diag(savs(bigX,V,-2,1))*beta;
        V_store(:,isave) = V;
        alpha_store(:,isave) = bigalpha;
        omega_store(:,isave) = vec(omega);
    end

end
time=toc

%Collect terms
pomega = mean(omega_store,2);
palpha = mean(alpha_store,2);
pbeta_level = mean(beta_level_store,2);
pbeta_full = mean(beta_store,2);
pbeta_savs = mean(beta_savs_store,2);
pbeta_savs_diff = mean(beta_fullsavs_store,2);

palpha = palpha(1:T:Q*T);
pbeta = reshape(mean(beta_fullsavs_store,2),[K,Q]);
pbeta_dense = reshape(mean(beta_store,2),[K,Q]);
pomega = reshape(pomega,[T,Q]);

bhat=[palpha';pbeta];
bhat_unspars=[palpha';pbeta_dense];

full_draw=nan([size(bhat),size(beta_store,2)]);
full_draw_sparse=nan([size(bhat),size(beta_store,2)]);
for i=1:size(beta_store,2)
    temp_beta=beta_store(:,i); %use beta_store instead of beta_fullsavs_store for full draw: beta_fullsavs_store can have chains of 0 which leads to Neff of NaN.
    temp_beta=reshape(temp_beta,[K,Q]);
    temp_beta_sparse=beta_fullsavs_store(:,i);
    temp_beta_sparse=reshape(temp_beta_sparse,[K,Q]);
    temp_alpha=alpha_store(:,i);
    temp_alpha=temp_alpha(1:T:Q*T);

    full_draw(:,:,i)=[temp_alpha';temp_beta];
    full_draw_sparse(:,:,i)=[temp_alpha';temp_beta_sparse];
end


end


