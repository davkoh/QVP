function[savs_out] = savs(X,pbeta,pen,iter)
K = size(X,2);
pbeta_savs=zeros(K,iter);
for s = 1:iter
for j = 1:K
    xtx = X(:,j)'*X(:,j);
    mu = abs(pbeta(j,s))^(pen);
    if mu >= abs(pbeta(j,s))*xtx
      pbeta_savs(j,s)=0;
    else
      pbeta_savs(j,s) = sign(pbeta(j,s))*(abs(pbeta(j,s))*xtx - mu)/xtx;
    end
end
end
 savs_out = pbeta_savs;