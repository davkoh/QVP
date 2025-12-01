function [bhat,full_draw,bhat_SAVS,full_draw_SAVS,time] = QVP_wrapper(ysmall,Xsmall,tau,iter,burn,Prior_choice,chain_i,stdize_flg,datwht_flg)
% Prior choice:
% 1     FreqQR              is normal QR
% 2     BRW                 is Bondell, Reich, Wang 2010 esitmator
% 3     BQR                 is Bayesian QR with ALD setup
% 4     HSBQR               is HS-BQR of Kohns, Szendrei (2024)
% 5     QVP_HS              is HS diff shrinkage within QVP setup
% 6     LQVP_HS             is HS priors with QVP setup and shared coefficients
% 7     LQVP_HS_SAVS        is HS priors with QVP setup and shared coefficients and SAVS
% 8     HS                  is HS with QVP setup but no fused shrinkage
% 9     QVP_alpha           is QVP with difference in constant influencing global

warning('error', 'MATLAB:nearlySingularMatrix') %Turn nearly singular warning into error. This activates try-catch loop on warnings too!
%% Collect terms needed
if nargin<9
    datwht_flg=1;
end
if nargin<8
    stdize_flg=1;
end
if nargin<7
    chain_i=1;
end

tau=unique(sort(tau));
Q = length(tau);
[T, K]=size(Xsmall); %MAKE SURE NO CONSTANT IN X!

if datwht_flg==1
    gramm_m = 1/T * (Xsmall'*Xsmall);
    [U,S,D] = svd(gramm_m);
    L= U;
    Lam_inv_sqr = diag(diag(S).^(-1/2));
    Xuse = Xsmall*L*Lam_inv_sqr/sqrt(T);
else
    Xuse=Xsmall;
end

if stdize_flg==1
    [Xuse,Centre,Scale]=normalize(Xuse);
    transform=[1,-Centre./Scale];
    transform=[transform;zeros(K,1),diag(1./Scale)];
end

bhat=nan(K+1,Q);
bhat_SAVS=nan(K+1,Q);
full_draw=nan([size(bhat),iter-burn]);
full_draw_SAVS=nan([size(bhat),iter-burn]);
if strcmp(Prior_choice,"FreqQR")==1 || strcmp(Prior_choice,"BRW")==1
    tic;
    if strcmp(Prior_choice,"FreqQR")==1
        for i=1:Q
            bhat(:,i)=rq_fnm([ones(T,1),Xuse],ysmall,tau(i));
        end
    elseif strcmp(Prior_choice,"BRW")==1
        bhat=BRW(ysmall,Xuse,tau);
    end
    time=toc
elseif strcmp(Prior_choice,"BQR")==1
    [bhat,full_draw,time]=BQR(ysmall,[ones(T,1),Xuse],tau,iter,burn,chain_i);
elseif strcmp(Prior_choice,"HSBQR")==1
    [bhat,full_draw,time]=HSBQR(ysmall,[ones(T,1),Xuse],tau,iter,burn,chain_i);
elseif strcmp(Prior_choice,"QVP_HS")==1
    [bhat,full_draw,~,time]=QVP_ALD(ysmall,Xuse,tau,iter,burn,chain_i);
elseif strcmp(Prior_choice,"QVP_alpha")==1
    [bhat,full_draw,~,time]=QVP_alpha(ysmall,Xuse,tau,iter,burn,chain_i);
elseif strcmp(Prior_choice,"LQVP_HS")==1 || strcmp(Prior_choice,"LQVP_HS_SAVS")==1
    [bhat_SAVS,bhat,full_draw_SAVS,full_draw,~,time]=LQVP_SAVS_ALD(ysmall,Xuse,tau,iter,burn,chain_i); %If you want to use sparsified full_draws: [bhat_SAVS,bhat,full_draw,~,pomega,time]=
elseif strcmp(Prior_choice,"HS")==1
    [bhat,full_draw,~,time]=HS(ysmall,Xuse,tau,iter,burn,chain_i);
elseif strcmp(Prior_choice,"CQR_HS")==1
    [bhat_SAVS,bhat,full_draw_SAVS,full_draw,~,time] = CQR_SAVS_ALD(ysmall,Xuse,tau,iter,burn,chain_i);
elseif strcmp(Prior_choice,"FusedQR_ALD")==1
    [bhat,full_draw,~,time]=FusedQR_ALD(ysmall,Xuse,tau,iter,burn,chain_i);
end


if stdize_flg==1
    bhat=transform*bhat;
    bhat_SAVS=transform*bhat_SAVS;
    for saved_iter=1:size(full_draw,3)
        full_draw(:,:,saved_iter)=transform*squeeze(full_draw(:,:,saved_iter));
        full_draw_SAVS(:,:,saved_iter)=transform*squeeze(full_draw_SAVS(:,:,saved_iter));
    end
end

if datwht_flg==1
    for l = 1:Q
        bhat(2:end,l) = L*Lam_inv_sqr*bhat(2:end,l)/sqrt(T);
        bhat_SAVS(2:end,l) = L*Lam_inv_sqr*bhat_SAVS(2:end,l)/sqrt(T);
        for saved_iter=1:size(full_draw,3)
            full_draw(2:end,l,saved_iter)=L*Lam_inv_sqr*squeeze(full_draw(2:end,l,saved_iter))/sqrt(T);
            full_draw_SAVS(2:end,l,saved_iter)=L*Lam_inv_sqr*squeeze(full_draw_SAVS(2:end,l,saved_iter))/sqrt(T);
        end
    end
end

end
