function [full_draw,bhat,qfit,eps,full_draw_SAVS,bhat_SAVS,qfit_SAVS,eps_SAVS,Z_fit] = BQVAR_Estim(Y,taus,eq_i,iter,burn,Prior_choice)
[T,~]= size(Y);
Z=[Y(1:end-1,:)];
if eq_i>1
    for j=1:(eq_i-1)
        Z=[Z,Y(2:end,j)];
    end
end
Z_fit= [ones(T-1,1) Z];
[bhat,full_draw,bhat_SAVS,full_draw_SAVS,time]=QVP_wrapper(Y(2:end,eq_i),Z, taus,iter,burn,Prior_choice);
qfit=Z_fit*bhat;
qfit_SAVS=Z_fit*bhat_SAVS;
eps=Y(2:end,eq_i)-qfit;
eps_SAVS=Y(2:end,eq_i)-qfit_SAVS;
end


%[bhat,full_draw,bhat_SAVS,time] = QVP_wrapper(ysmall,Xsmall,tau,iter,burn,Prior_choice,chain_i,stdize_flg,datwht_flg)