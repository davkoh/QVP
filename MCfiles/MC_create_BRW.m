function constr_montecarlo = MC_create_BRW(T, const, tau, quantspec_cutoff)
%case=1 (same as BRW)
%case=2 (same as BRW)
%case=3 (same as BRW except quant variation set to 0.1)
%case=4 (case 2 with quant spec sparsity)
%case=5 (case 1 with half covariates 0.1 quant variation, and other half set to 1)


if nargin < 4
    quantspec_cutoff = 0.1;
end

if quantspec_cutoff>0.5
    quantspec_cutoff=1-quantspec_cutoff;
end

T = T + 100;
ynum=5;
max_beta=10;
beta_mat=nan(max_beta+1,length(tau),ynum);
Data_y=nan(T,ynum);
Data_x=nan(T,max_beta,ynum);

for ydes=1:ynum
    beta_full=nan(10+1,length(tau));
    beta_full(1,:)=norminv(tau)+1;
    eps=normrnd(0,1,[T,1]);
    if ydes==1
        for row=2:4+1
            beta_full(row,:)=norminv(tau)*0.1+1;
        end
        X=unifrnd(0,1,[T,4]);
        beta_use=[1,1,1,1];
        theta_use=[0.1,0.1,0.1,0.1];
    elseif ydes==2
        for row=2:size(beta_full,1)
            if row<=5
                beta_full(row,:)=norminv(tau)*0.1+1;
            else
                beta_full(row,:)=norminv(tau)*0;
            end
        end
        X=unifrnd(0,1,[T,10]);
        beta_use=[1,1,1,1,0,0,0,0,0,0];
        theta_use=[0.1,0.1,0.1,0.1,0,0,0,0,0,0];
    elseif ydes==3
        for row=2:7+1
            if row<=4
                beta_full(row,:)=norminv(tau)*0.1+1;
            else
                beta_full(row,:)=norminv(tau)*0+1;
            end
        end
        X=unifrnd(0,1,[T,7]);
        beta_use=[1,1,1,1,1,1,1];
        theta_use=[0.1,0.1,0.1,0,0,0,0];
    elseif ydes==4
        quant_spec_mat=(round(tau,3)<=repmat(round(quantspec_cutoff,3),1,length(tau)))+(round(tau,3)>=round(1-quantspec_cutoff,3));
        for row=2:size(beta_full,1)
            if row<=5
                beta_full(row,:)=norminv(tau)*0.1+1;
            else
                beta_full(row,:)=norminv(tau).*quant_spec_mat;
            end
        end
        beta_full(end-1:end,:)=zeros(2,length(tau));
        X=unifrnd(0,1,[T,10]);
        beta_use=[1,1,1,1,0,0,0,0,0,0];
        theta_use=[0.1,0.1,0.1,0.1,1,1,1,1,0,0];
        qcheck=(eps<=norminv(quantspec_cutoff))+(eps>=norminv(1-quantspec_cutoff));
    elseif ydes==5
        for row=2:4+1
            if row<=3
                beta_full(row,:)=norminv(tau)*0.1+1;
            else
                beta_full(row,:)=norminv(tau)*1+1;
            end
        end
        X=unifrnd(0,1,[T,4]);
        beta_use=[1,1,1,1];
        theta_use=[0.1,0.1,1,1];
    end
    beta_mat(:,:,ydes)=beta_full;

    if ydes==4
        y=const*ones(T,1)+X*beta_use'+(ones(T,1)+X(:,(theta_use==0.1))*theta_use(theta_use==0.1)'+qcheck.*X(:,(theta_use==1))*theta_use(theta_use==1)').*eps;
    else
        y=const*ones(T,1)+X*beta_use'+(ones(T,1)+X*theta_use').*eps;
    end

    Data_y(:,ydes)=y;
    Data_x(:,1:size(X,2),ydes)=X;
end

constr_montecarlo = struct();
constr_montecarlo.Y = Data_y;
constr_montecarlo.X = Data_x;
constr_montecarlo.coeff_full = beta_mat;
end