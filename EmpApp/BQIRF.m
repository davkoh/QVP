function [QIRF]=BQIRF(Theta,taus,bhat_use,h,eps_full,Varib_num)
K=size(eps_full,3);
Q=size(bhat_use,2);
IRF=zeros(K,h);
%[~, A_0, A_1,eps,b]=QuantileVAR_System(Y,theta);

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

eps=[];
for eq_i=1:K
    [~,loc]=min(abs(Theta(eq_i)-taus));
    
    A_1_use(eq_i,:)=squeeze(A_1(eq_i,:,loc));
    A_0_use(eq_i,:)=squeeze(A_0(eq_i,:,loc));
    eps=[eps,squeeze(eps_full(:,loc,eq_i))];
end

Shock=[0, std(eps(:,2))]';
IRF(:,1)=(eye(K)-A_0_use)\Shock;
Phi=(eye(K)-A_0_use)\A_1_use;
for l=2:h
    IRF(:,l)=Phi*IRF(:,l-1);
end
QIRF=IRF';
end