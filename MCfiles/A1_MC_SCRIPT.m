clear
clc
close all

tau = 0.05:0.05:0.95;
iter = 15000;
burn = 10000;
MC_T=100;
MC_num=2;
rho = 0.0;
Chain_num=4;
Prior_mat=["FreqQR","BRW","BQR","HSBQR","CQR_HS","FusedQR_ALD","QVP_HS","LQVP_HS","LQVP_HS_SAVS"];

%% Set seed (as a function of the Job ID)
rng(1998,'twister');  % Pour one out to the undertaker

mainpath=cd;
funcs=fileparts(mainpath)+"/QVP_package/";

addpath(funcs);

%% Rearrange Prior_mat to make sure SAVS is always saved for LQVP_HS and LQVP_HS_ASIS
%Order also matters! Since we run the LQVP anyways, the savs is run
%automatically too. To ensure that the location of runs and the prior_mat
%is correct we remove all entries with "xxx_SAVS" in them

loc=(ones(1,length(Prior_mat))-contains(Prior_mat,"_SAVS"))==1;
Prior_temp=Prior_mat(loc);
ASIS_flg=sum(contains(Prior_temp,"LQVP_HS_ASIS"));
LQVP_flg=sum(contains(Prior_temp,"LQVP_HS"))-ASIS_flg;

%Keep only non-SAVS
loc=(ones(1,length(Prior_temp))-contains(Prior_temp,"LQVP_HS"))==1;
Prior_temp=Prior_temp(loc);

if LQVP_flg==1 && ASIS_flg==1
    Prior_mat=[Prior_temp,"LQVP_HS","LQVP_HS_SAVS","LQVP_HS_ASIS","LQVP_HS_ASIS_SAVS"];
elseif LQVP_flg==1 && ASIS_flg==0
    Prior_mat=[Prior_temp,"LQVP_HS","LQVP_HS_SAVS"];
elseif LQVP_flg==0 && ASIS_flg==1
    Prior_mat=[Prior_temp,"LQVP_HS_ASIS","LQVP_HS_ASIS_SAVS"];
elseif LQVP_flg==0 && ASIS_flg==0
    Prior_mat=Prior_temp;
end

for MC_it=1:MC_num
    %% Prelocate parfor loop variables
    if rho==0
        datause=MC_create_BRW(MC_T, 1, tau);
    else
        datause=MC_create_BRW_normal(MC_T, 1, tau,0.1,rho);
    end
    True_beta=datause.coeff_full; %Call true_beta here. Easier to preallocate space
    case_tot=size(datause.X,3); %Save case totals here.
    Prio_mat_length=length(Prior_mat);

    %% Run Monte Carlo experiments
    % Pre-alocation
    bhat_ALL=nan([size(True_beta),length(Prior_mat)]);
    size(bhat_ALL)
    bhat_deviation=nan([size(True_beta),length(Prior_mat)]);
    Rhat=nan([size(True_beta),length(Prior_mat)]);
    Neff=nan([size(True_beta),length(Prior_mat)]);
    wQS_ALL=nan([4,size(True_beta,3),length(Prior_mat)]);
    Cross_ALL=nan([size(True_beta,2),size(True_beta,3),length(Prior_mat)]);
    Time_mat=nan([case_tot,Prio_mat_length,Chain_num]);

    % Run the MC instance
    Prelocated_y=nan(MC_T,size(True_beta,3));
    Prelocated_x=nan(MC_T,size(True_beta,1)-1,size(True_beta,3));
    Forecast_y=nan(100,size(True_beta,3));
    Forecast_x=nan(100,size(True_beta,1)-1,size(True_beta,3));

    if rho==0
        datause=MC_create_BRW(MC_T, 1, tau);
    else
        datause=MC_create_BRW_normal(MC_T, 1, tau,0.1,rho);
    end
    %Save data that will be used for estimation
    Prelocated_y(:,:)=datause.Y(1:MC_T,:);
    Prelocated_x(:,:,:)=datause.X(1:MC_T,:,:);

    %Save data that will be used for OOS evaluation
    Forecast_y(:,:)=datause.Y(MC_T+1:end,:);
    Forecast_x(:,:,:)=datause.X(MC_T+1:end,:,:);

    for Prior_i=1:Prio_mat_length
        Prior_choice=Prior_mat(Prior_i);

        if strcmp(Prior_choice,"LQVP_HS_SAVS")==0 && strcmp(Prior_choice,"LQVP_HS_ASIS_SAVS")==0
            Rhat_temp=nan([size(True_beta)]);
            Neff_temp=nan([size(True_beta)]);

            %Extra Preallocation for faster ruunning
            bhat_ALL_temp_UNSPARS=nan([size(True_beta)]);
            bhat_ALL_temp_SAVS=nan([size(True_beta)]);
            bhat_deviation_temp_UNSPARS=nan([size(True_beta)]);
            bhat_deviation_temp_SAVS=nan([size(True_beta)]);
            wQS_temp_UNSPARS=nan([4,3]);
            wQS_temp_SAVS=nan([4,3]);
            Cross_temp_UNSPARS=nan([size(True_beta,2),size(True_beta,3)]);
            Cross_temp_SAVS=nan([size(True_beta,2),size(True_beta,3)]);

            True_beta_temp=True_beta;
            for case_num=1:case_tot
                ysmall=squeeze(Prelocated_y(:,case_num));
                Xsmall=squeeze(Prelocated_x(:,:,case_num));
                Xsmall(:, all(isnan(Xsmall),1)) = []; %Remove columns with only NaN
                bhat_chains=nan([size(Xsmall,2)+1,length(tau),Chain_num]);
                bhat_chains_SAVS=nan([size(Xsmall,2)+1,length(tau),Chain_num]);
                posterior_ALL_temp=nan([size(Xsmall,2)+1,length(tau),iter-burn,Chain_num]);
                time_temp=nan(Chain_num,1);
                for chain_i=1:Chain_num %Change this to parfor loop if parallel toolbox available
                    fprintf("Draws for model: %s, Sim-Case: %i, Chain: %i \n", Prior_mat(Prior_i),case_num,chain_i);
                    try
                        [bhat,full_draw,bhat_SAVS,~,time] = QVP_wrapper(ysmall,Xsmall,tau,iter,burn,Prior_choice,chain_i,1,1);
                    catch
                        bhat=nan([size(Xsmall,2)+1,length(tau)]);
                        bhat_SAVS=nan([size(Xsmall,2)+1,length(tau)]);
                        full_draw=nan([size(Xsmall,2)+1,length(tau),iter-burn]);
                        time=nan(1,1);
                    end

                    %Store chains
                    bhat_chains(:,:,chain_i)=bhat; %Save seperately to average later
                    bhat_chains_SAVS(:,:,chain_i)=bhat_SAVS; %Save seperately to average later
                    posterior_ALL_temp(:,:,:,chain_i)=full_draw;
                    time_temp(chain_i)=time;
                end
                Time_mat(case_num,Prior_i,:)=time_temp;
                if strcmp(Prior_choice,"LQVP_HS")==1 || strcmp(Prior_choice,"LQVP_HS_ASIS")==1
                    Time_mat(case_num,Prior_i+1,:)=time_temp;
                end

                %Calculate converge diag
                Rhat_mat=nan(size(bhat_chains,1:2));
                Neff_mat=nan(size(bhat_chains,1:2));

                for quant_i = 1:length(tau)
                    beta_temp = squeeze(posterior_ALL_temp(:,quant_i,:,:));
                    beta_temp = permute(beta_temp,[2,1,3]);
                    [Rhat_calc,Neff_calc] = psrf(beta_temp);
                    Rhat_mat(:,quant_i) = Rhat_calc;
                    Neff_mat(:,quant_i) = Neff_calc;
                end

                Rhat_temp(1:size(bhat_chains,1),:,case_num)=Rhat_mat;
                Neff_temp(1:size(bhat_chains,1),:,case_num)=Neff_mat;

                %Calculations for SAVS and non-SAVS variants
                for bhat_it=1:2
                    if bhat_it==1
                        bhat_chain_use=bhat_chains;
                    else
                        bhat_chain_use=bhat_chains_SAVS;
                    end

                    %We take average of all chains and store those too
                    bhat=squeeze(mean(bhat_chain_use,3));
                    bhat_save=nan(size(True_beta_temp,1),size(True_beta_temp,2));
                    bhat_save(1:size(bhat,1),:)=bhat;

                    %Calculate deviation from true beta
                    True_beta_use=squeeze(True_beta_temp(:,:,case_num));
                    True_beta_use(all(isnan(True_beta_use),2),:) = []; %Remove rows with only NaN
                    bhat_dev_save=nan(size(True_beta_temp,1),size(True_beta_temp,2));
                    bhat_dev_save(1:size(bhat,1),:)=(bhat-True_beta_use).^2;

                    %Cross calc
                    Crossfit=[ones(size(Xsmall,1),1),Xsmall]*bhat;

                    %OOS eval
                    y_OOS=squeeze(Forecast_y(:,case_num));
                    xfit=squeeze(Forecast_x(:,:,case_num));
                    xfit(:, all(isnan(xfit),1)) = []; %Remove columns with only NaN
                    yfit=[ones(size(y_OOS,1),1),xfit]*bhat;
                    wQS_temp=nan(4,1); %New line
                    for w=1:4
                        wQS_temp(w,1)=mean(wQS(y_OOS, yfit, tau, w));
                    end

                    if bhat_it==1
                        bhat_ALL_temp_UNSPARS(:,:,case_num)=bhat_save;
                        bhat_deviation_temp_UNSPARS(:,:,case_num)=bhat_dev_save;
                        wQS_temp_UNSPARS(:,case_num)=wQS_temp;
                        Cross_temp_UNSPARS(:,case_num)=mean(sort(Crossfit,2)~=Crossfit);
                    else
                        bhat_ALL_temp_SAVS(:,:,case_num)=bhat_save;
                        bhat_deviation_temp_SAVS(:,:,case_num)=bhat_dev_save;
                        wQS_temp_SAVS(:,case_num)=wQS_temp;
                        Cross_temp_SAVS(:,case_num)=mean(sort(Crossfit,2)~=Crossfit);
                    end
                end
            end
            bhat_ALL(:,:,:,Prior_i)=bhat_ALL_temp_UNSPARS;
            bhat_deviation(:,:,:,Prior_i)=bhat_deviation_temp_UNSPARS;
            Rhat(:,:,:,Prior_i)=Rhat_temp;
            Neff(:,:,:,Prior_i)=Neff_temp;
            wQS_ALL(:,:,Prior_i)=wQS_temp_UNSPARS;
            Cross_ALL(:,:,Prior_i)=Cross_temp_UNSPARS;

            if strcmp(Prior_choice,"LQVP_HS")==1 || strcmp(Prior_choice,"LQVP_HS_ASIS")==1
                bhat_ALL(:,:,:,Prior_i+1)=bhat_ALL_temp_SAVS;
                bhat_deviation(:,:,:,Prior_i+1)=bhat_deviation_temp_SAVS;
                Rhat(:,:,:,Prior_i+1)=Rhat_temp;
                Neff(:,:,:,Prior_i+1)=Neff_temp;
                wQS_ALL(:,:,Prior_i+1)=wQS_temp_SAVS;
                Cross_ALL(:,:,Prior_i+1)=Cross_temp_SAVS;
            end
        end
    end

    % Save output
    file="ZZZ_Results_it="+MC_it+"_T="+MC_T+".mat";
    save(file,'wQS_ALL','Cross_ALL','bhat_deviation','bhat_ALL','Rhat','Neff','Prior_mat','tau','Time_mat','MC_T','rho');

    close all
end
