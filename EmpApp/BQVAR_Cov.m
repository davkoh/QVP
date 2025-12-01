function [Sigma,vect_mat] = BQVAR_Cov(All_Draws)

    [K,Q,iter]=size(All_Draws);
    vect_mat=nan(iter,K*Q);

    it=0;
    for k_it=1:K
        for q_it=1:Q
            it=it+1;
            vect_mat(:,it)=squeeze(All_Draws(k_it,q_it,:))';
        end
    end
    Sigma=cov(vect_mat);
end