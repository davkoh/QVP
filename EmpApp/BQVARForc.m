function [Y_mean, CI, Y_N, M]=BQVARForc(Varib_num,Y,theta_full,taus,sample,All_Draws,h,N,k_choice,CI_quant,given_theta_flg)

if nargin<11
    given_theta_flg=0;
end

K=size(Y,2);
Q=length(taus);

Y_N=zeros(h+1,N);
M=zeros(K^2,K^2);
for i=1:N
    W=Y(sample,:)';
    Y_N(1,i)=W(k_choice);
    if given_theta_flg==0
        theta_full=zeros(K,h);
        for l=1:h
            theta_full(:,l)=randsample(taus,2,true)';
        end
    end
    for l=1:h
        theta=theta_full(:,l);
        bhat_use=squeeze(All_Draws(:,:,randi(size(All_Draws,3)))); %Bring in model uncert as well given that we have posterior draws

        omega=zeros(K,1,Q);
        A_0=zeros(K,K,Q);
        A_1=zeros(K,K,Q);
        for eq_i=1:K
            if eq_i==1
                bhat_temp=bhat_use(1:Varib_num(eq_i),:);
            else
                bhat_temp=bhat_use(Varib_num(eq_i-1)+1:Varib_num(eq_i-1)+Varib_num(eq_i),:);
            end
            omega(eq_i,1,:)=bhat_temp(1,:);
            for quant_i=1:Q
                A_1(eq_i,:,quant_i)=bhat_temp(2:(K+1),quant_i);
                if eq_i>1
                    A_0(eq_i,1:eq_i-1,quant_i)=bhat_temp(K+2:(K+1+(eq_i-1)),quant_i);
                end
            end
        end

        A_0_use=zeros(K,K);
        A_1_use=zeros(K,K);
        omega_use=zeros(K,1);
        for eq_i=1:K
            [~,loc]=min(abs(theta(eq_i)-taus));


            A_1_use(eq_i,:)=squeeze(A_1(eq_i,:,loc));
            A_0_use(eq_i,:)=squeeze(A_0(eq_i,:,loc));
            omega_use(eq_i,:)=squeeze(omega(eq_i,:,loc));
        end
        %[omega_use, A_0_use, A_1_use,~,~]=QuantileVAR_System(Y,randsample(taus,2,true));
        Gamma=pinv(eye(K)-A_0_use);
        Phi=Gamma*A_1_use;
        phi=Gamma*omega_use;
        W=phi+Phi*W;
        Y_N(1+l,i)=W(k_choice);
    end
    M=M+kron(Phi,Phi);
end
M=M./N;

Y_mean=mean(Y_N,2);
CI=[quantile(Y_N',CI_quant/2)',quantile(Y_N',1-(CI_quant/2))'];
end