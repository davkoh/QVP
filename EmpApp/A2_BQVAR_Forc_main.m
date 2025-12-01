clear all
close all
clc

%Parameters
taus=0.05:0.05:0.95;
iter = 15000;
burn = 10000;
Prior_choice=["BRW","BQR","HSBQR","QVP_HS","LQVP_HS","LQVP_HS_SAVS"];
fig_row_num=2;
fig_col_num=3;
setfontsize_main=18;
T_Mat=96:224; %the expanding window for forecast (Original: Start in T=1:96 and expand window until T=1:224->T_Mat=96:224)

%Choose the horizon h and the number of repetitions N
h=10;
N=10^3;
AllDraws_Flg=1; %If 1, for forecast all posterior draws will be used. If it is 0, only posterior means are used

%Paths
mainpath=cd;
dataloc=mainpath+"/data/";
funcs=fileparts(mainpath)+"/QVP_package/";

addpath(funcs,dataloc);

%Read data
xlsfilename = 'DataQVAR.xlsx';
[X, TEXT]   = xlsread(xlsfilename);
dateorig    = (TEXT(3:end,1));

%Define variables
Headers = {'IP growth' 'CISS'};
Y_Full=[100*(log(X(2:end,1))-log(X(1:end-1,1))) X(2:end,2)];
date=datenum(dateorig(1:end,1));

[T, K]=size(Y_Full);
Q=length(taus);



for T_use=T_Mat
    Y=Y_Full(1:T_use,:);
    Y_datename=dateorig(2:T_use,:);
    Y_forc=Y_Full(T_use:T_use+h,:);
    Y_forc_datename=dateorig(T_use:T_use+h,:);

    sample=T_use;

    Forecast_Save=nan(h+1,Q,K,length(Prior_choice));
    for Prior_i=1:length(Prior_choice) %Swap to parfor loop if parallel toolbox available
        %% Estimation
        bhat_full=[];
        Varib_num=[];

        omega=zeros(K,1,Q);
        A_0=zeros(K,K,Q);
        A_1=zeros(K,K,Q);

        for eq_i=1:K
            [full_draw,bhat,qfit,eps,full_draw_SAVS,bhat_SAVS,qfit_SAVS,eps_SAVS,Z_fit] = BQVAR_Estim(Y,taus,eq_i,iter,burn,Prior_choice(Prior_i));

            if strcmp(Prior_choice(Prior_i),"LQVP_HS_SAVS")==0 && strcmp(Prior_choice(Prior_i),"LQVP_HS_ASIS_SAVS")==0
                bhat_use=bhat;
                eps_use=eps;
                full_draw_use=full_draw;
            else
                bhat_use=bhat_SAVS;
                eps_use=eps_SAVS;
                full_draw_use=full_draw_SAVS;
            end
            omega(eq_i,1,:)=bhat_use(1,:);
            for quant_i=1:Q
                A_1(eq_i,:,quant_i)=bhat_use(2:(K+1),quant_i);
                if eq_i>1
                    A_0(eq_i,1:eq_i-1,quant_i)=bhat_use(K+2:(K+1+(eq_i-1)),quant_i);
                end
            end
            bhat_full=[bhat_full;bhat_use];
            Varib_num=[Varib_num;size(bhat_use,1)]; %tracking the # of variables in an equation
            if eq_i==1
                All_Draws=full_draw;
                eps_full=eps_use;
            else
                All_Draws=cat(1,All_Draws,full_draw);
                eps_full=cat(3,eps_full,eps_use);
            end

        end
        [V_T,~] = BQVAR_Cov(All_Draws);

        %% Forecast
        All_Draws_Mean=bhat_full; %Means for SAVS variants for the "multiforecast" setups we use unsparsified draws
        if strcmp(Prior_choice(Prior_i),"BRW")==0%strcmp(Prior_choice(Prior_i),"LQVP_HS_SAVS")==0 && strcmp(Prior_choice(Prior_i),"LQVP_HS_ASIS_SAVS")==0  && strcmp(Prior_choice(Prior_i),"BRW")==0
            All_Draws_Post=All_Draws;
        else
            All_Draws_Post=bhat_full; %For SAVS Alldraws is without sparsification. So using the mean bhat for these. No conf band
        end
        if AllDraws_Flg==1
            All_Draws_Use=All_Draws_Post;
        else
            All_Draws_Use=All_Draws_Mean;
        end

        %Forecast each eq
        for varib_i=1:K
            theta_full=zeros(K,h);
            [~, ~, Y_N, M]=BQVARForc(Varib_num,Y,theta_full,taus,sample,All_Draws_Use,h,N,varib_i,0.1,0);
            % Evaluate covariance stationarity of the system
            Cov_stationarity=max(abs(eig(M)))

            Forecast_Save(:,:,varib_i,Prior_i)=prctile(Y_N , taus*100, 2 ); %h-ahead multistep forecasts for the quantiles of interest
        end
    end
    file=mainpath+"/ZZZ_Forc_T="+T_use+".mat";
    save(file,"Forecast_Save","Y_forc","taus","Y","Y_forc_datename","Y_datename","Prior_choice","sample");
end


