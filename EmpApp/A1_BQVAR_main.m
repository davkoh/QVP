clear all
close all
clc

%Parameters
taus=0.05:0.05:0.95;
iter = 15000;
burn = 10000;
Prior_choice=["BQR","HSBQR","BRW","QVP_HS","LQVP_HS","LQVP_HS_SAVS"];
fig_row_num=2;
fig_col_num=3;
setfontsize_main=18;

%Paths
mainpath=cd;
dataloc=mainpath+"/data/";
funcs=fileparts(mainpath)+"/QVP_package/";

addpath(funcs,dataloc);

%% CODE STARTS HERE

Prior_mat=Prior_choice;
for mod_i=1:length(Prior_mat)
    if strcmp(Prior_mat(mod_i),"QVP_HS")==1
        Prior_mat(mod_i)="QVP";
    elseif strcmp(Prior_mat(mod_i),"LQVP_HS")==1
        Prior_mat(mod_i)="NC-QVP";
    elseif strcmp(Prior_mat(mod_i),"LQVP_HS_SAVS")==1
        Prior_mat(mod_i)="NC-QVP_{SAVS}";
    elseif strcmp(Prior_mat(mod_i),"LQVP_HS_ASIS")==1
        Prior_mat(mod_i)="NC-QVP_{ASIS}";
    elseif strcmp(Prior_mat(mod_i),"LQVP_HS_ASIS_SAVS")==1
        Prior_mat(mod_i)="NC-QVP_{ASIS,SAVS}";
    end
end

%Read data
xlsfilename = 'DataQVAR.xlsx';
[X, TEXT]   = xlsread(xlsfilename);
date        = datenum(TEXT(2:end,1));

%Define variables
Headers = {'IP growth' 'CISS'};
Y=[100*(log(X(2:end,1))-log(X(1:end-1,1))) X(2:end,2)];
date=date(2:end,1);

[T, K]=size(Y);
Q=length(taus);

%For Stress test forecast
sample = find(year(date) == 2008 & month(date) ==  9);

%DataPlot
fig=figure();
set(gcf, 'Position', [1,1,1920/2,1080/2],'DefaultAxesFontSize',setfontsize_main);
yyaxis left
plot(date, Y(:, 1),'k','LineWidth',2)
ax.YColor = 'k';
ylim([-5, 5])
ylabel('Industrial production growth');
datetick('x', 'yyyy')
yyaxis right
plot(date, Y(:, 2),'r--','LineWidth',2)
ax.YColor = 'r';
ylabel('Financial Stress Index');
legend('Industrial production', 'CISS', 'Location', 'NorthWest')
filename=mainpath+"/FigureData.jpg";
saveas(fig,filename);
close all

%% Estimation
All_Pr_collect{1,2}="bhat_full";
All_Pr_collect{1,3}="All_draws";
All_Pr_collect{1,4}="V_T";
All_Pr_collect{1,5}="eps";
All_Pr_collect{1,6}="omega";
All_Pr_collect{1,7}="A1";
All_Pr_collect{1,8}="A0";
All_Pr_collect{1,9}="qfit";
All_Pr_collect{1,10}="Varib_num";


for Prior_i=1:length(Prior_mat) %Swap to parfor loop if parallel toolbox available
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
            All_Draws=full_draw_use;
            eps_full=eps_use;
            qfit_full=qfit;
        else
            All_Draws=cat(1,All_Draws,full_draw_use);
            eps_full=cat(3,eps_full,eps_use);
            qfit_full=cat(3,qfit_full,qfit);
        end

    end
    [V_T,~] = BQVAR_Cov(All_Draws);

    temp={Prior_mat(Prior_i),bhat_full,All_Draws,V_T,eps_full,omega,A_1,A_0,qfit_full,Varib_num};
    [All_Pr_collect{Prior_i+1,1:10}]=temp{:};
end

Varib_num=All_Pr_collect{2,10};

filename=mainpath+"/EUruns.mat";
save(filename,"All_Pr_collect","Varib_num","taus","Prior_mat");

%% Crossnum
CrossNum=NaN(length(Prior_mat),3);
for Prior_i=1:length(Prior_mat)
    temp=All_Pr_collect{Prior_i+1,9};
    for eq_i=1:K
        qfit_use=squeeze(temp(:,:,eq_i));
        qfit_use_sort=sort(qfit_use,2);

        CrossNum(Prior_i,eq_i)=mean(mean((qfit_use~=qfit_use_sort)));

        if eq_i==1
            qfit_overall=qfit_use;
            qfit_overall_sort=qfit_use_sort;
        else
            qfit_overall=[qfit_overall;qfit_use];
            qfit_overall_sort=[qfit_overall_sort;qfit_use_sort];
        end
    end
    CrossNum(Prior_i,end)=mean(mean((qfit_overall~=qfit_overall_sort)));
end

[Prior_mat',round(CrossNum*100,2)]

%% Beta Profile
Parameters = {'b_1(\tau_q)','a_{11}(\tau_q)','a_{12}(\tau_q)','b_2(\tau_q)','a_{21}(\tau_q)','a_{22}(\tau_q)','a_{021}(\tau_q)'};
for l=1:size(All_Pr_collect{2,2},1)
    fig=figure();
    set(gcf, 'Position', [1,1,1920,1080/2],'DefaultAxesFontSize',setfontsize_main);
    tiledlayout(1,fig_row_num*fig_col_num-1);
    for Prior_i=2:length(Prior_mat)
        nexttile;
        bhat_use=All_Pr_collect{1+1,2};
        All_draw_use=All_Pr_collect{1+1,3};

        b=squeeze(bhat_use(l,:));
        temp=squeeze(All_draw_use(l,:,:))';
        bands=nan(2,length(taus));
        for tau_i=1:length(taus)
            temp_tau=squeeze(temp(:,tau_i));
            bands(1,tau_i)=prctile(temp_tau,5);
            bands(2,tau_i)=prctile(temp_tau,95);
        end
        %sigma=std(temp,'omitnan');

        plot(taus,b','*-b','LineWidth',2);
        hold on
        %plot(taus,[b-2*sigma;b+2*sigma]','--c','LineWidth',2);
        plot(taus,bands','--c','LineWidth',2);

        bhat_use=All_Pr_collect{Prior_i+1,2};
        All_draw_use=All_Pr_collect{Prior_i+1,3};

        b=squeeze(bhat_use(l,:));
        temp=squeeze(All_draw_use(l,:,:))';
        bands=nan(2,length(taus));
        for tau_i=1:length(taus)
            temp_tau=squeeze(temp(:,tau_i));
            bands(1,tau_i)=prctile(temp_tau,5);
            bands(2,tau_i)=prctile(temp_tau,95);
        end
        %sigma=std(temp,'omitnan');
        plot(taus,b','*-r','LineWidth',2);
        %plot(taus,[b-2*sigma;b+2*sigma]','--k','LineWidth',2);
        plot(taus,bands','--k','LineWidth',2);

        set(gcf,'Color','w')
        axis tight; box on;
        hold off;
        xlim([taus(1) taus(end)])
        xlabel('\tau')
        ylabel('Coefficient')
        if l==1
            ylim([-1 2]);
        elseif l==2
            ylim([-0.5 0.5]);
        elseif l==3
            ylim([-6 2]);
        elseif l==4
            ylim([-0.05 0.15]);
        elseif l==5
            ylim([-0.06 0.04]);
        elseif l==6
            ylim([0.5 1.5]);
        elseif l==7
            ylim([-0.02 0.05]);
        end
        title(Prior_mat(Prior_i));
        subtitle(Parameters(l));
    end
    filename=mainpath+"/FigureCoeff"+l+".jpg";
    saveas(fig,filename);
    close all
end

%% QIRF
h=30;
Theta=[taus;0.5*ones(size(taus))];

fig=figure();
set(gcf, 'Position', [1,1,1920,1080],'DefaultAxesFontSize',setfontsize_main);
tiledlayout(fig_row_num,fig_col_num);
for Prior_i=1:length(Prior_mat)
    bhat_use=All_Pr_collect{Prior_i+1,2};
    eps_full=All_Pr_collect{Prior_i+1,5};
    Varib_num=All_Pr_collect{Prior_i+1,10};

    QIRF_IP=zeros(h,size(Theta,2));
    for i=1:length(taus)
        %[QIRF, CI]=BQIRF(theta,taus,b_T,h,V_T,A_0_full,A_1_full,eps_full,Varib_num)
        IR=BQIRF(Theta(:,i),taus,bhat_use,h,eps_full,Varib_num);
        QIRF_IP(:,i)=IR(:,1);
    end

    nexttile;
    H=(2:h);
    surf(H', taus', QIRF_IP(2:end,:)', QIRF_IP(2:end,:)')
    colormap summer;
    xlabel('h')
    zlabel('QIRF')
    zlim([-0.3 0.05]);
    ylabel('\tau')
    title(Prior_mat(Prior_i));
end
filename=mainpath+"/QIRF.jpg";
saveas(fig,filename);
close all

%% Quantile Specific QIRF
N=10^3;
Prior_use=["BQR","QVP_{HS}","QVP"];
tau_choice=[0.05,0.5,0.95];

fig=figure();
set(gcf, 'Position', [1,1,1920,1080],'DefaultAxesFontSize',setfontsize_main);
tiledlayout(2,length(tau_choice));
for Prior_i=1:size(All_Pr_collect,1)
    Prior_iter=All_Pr_collect{Prior_i,1};
    if sum(strcmp(Prior_use,Prior_iter))>0
        if sum(strcmp("QVP_{HS}",Prior_iter))>0
            Prior_iter="QVP";
        end
        Draws_use=All_Pr_collect{Prior_i,3};
        eps_full=All_Pr_collect{Prior_i,5};
        Varib_num=All_Pr_collect{Prior_i,10};
        for tau_i=1:length(tau_choice)
            Theta=[tau_choice(tau_i);0.5];
            QIRF_IP=nan(N,h);
            for iter=1:N
                bhat_use=squeeze(Draws_use(:,:,randi(size(Draws_use,3))));
                IR=BQIRF(Theta,taus,bhat_use,h,eps_full,Varib_num);
                QIRF_IP(iter,:)=IR(:,1)';
            end
            mean_path=mean(QIRF_IP);
            bands=nan(2,h);
            for h_i=1:h
                temp_h=squeeze(QIRF_IP(:,h_i));
                bands(1,h_i)=prctile(temp_h,5);
                bands(2,h_i)=prctile(temp_h,95);
            end
            nexttile;
            subtitle(Prior_iter)
            plot([1:h],mean_path(1:end)','-r','LineWidth',2);
            hold on
            plot([1:h],bands(:,1:end)','--k','LineWidth',2);
            hold off
            yline(0);
            title(Prior_iter)
            subtitle(sprintf("\\tau=%.2f",round(tau_choice(tau_i),2)));
            if tau_i==1
                ylim([-0.3 0]);
            elseif tau_i==2
                ylim([-0.1 0.05]);
            else
                ylim([-0.03 0.06]);
            end
        end
    end
end
filename=mainpath+"/QIRF_QuantSpec.jpg";
saveas(fig,filename);
close all

%% Stress test

%Choose the horizon h and the number of repetitions N
h=10;
N=10^3;

fig = figure;
set(gcf, 'Position', [1,1,1920,1080],'DefaultAxesFontSize',setfontsize_main);
tiledlayout(fig_row_num,fig_col_num);
for Prior_i=1:length(Prior_mat)
    if strcmp(Prior_mat(Prior_i),"BRW")==0
        All_Draws1=All_Pr_collect{Prior_i+1,3}; %BRW is not Bayesian
    else
        All_Draws1=All_Pr_collect{Prior_i+1,2};
    end
    if strcmp(Prior_mat(Prior_i),"BRW")==0
        All_Draws2=All_Pr_collect{Prior_i+1,3};
    else
        All_Draws2=All_Pr_collect{Prior_i+1,2}; %For SAVS Alldraws is without sparsification. So using the mean bhat for these. No conf band
    end
    theta_full=zeros(K,h);
    [~, ~, Y_N, M]=BQVARForc(Varib_num,Y,theta_full,taus,sample,All_Draws1,h,N,1,0.1,0);
    % Evaluate covariance stationarity of the system
    Cov_stationarity=max(abs(eig(M)))

    %Calculate forecasts using repetitive median and stress scenarios
    theta_full=0.5*ones(K,h);
    [Y_median, CI_median, ~, ~]=BQVARForc(Varib_num,Y,theta_full,taus,sample,All_Draws2,h,N,1,0.05,1);

    theta_full=[[0.1;0.9].*ones(K,6),0.5*ones(K,4)];
    [Y_stress, CI_stress, ~, ~]=BQVARForc(Varib_num,Y,theta_full,taus,sample,All_Draws2,h,N,1,0.05,1);

    nexttile;
    f=zeros(2,1);
    plot(0:h, prctile(Y_N , sort([taus*100,1,99]), 2 ), 'Color', 1/300*[200,200,200],'linestyle','--');
    hold on
    f(1)=plot(0:h, Y_median(:,1),'r-','LineWidth',1);
    hold on
    plot(0:h, CI_median(:,1:2),'k--','LineWidth',1);
    hold on
    f(2)=plot(0:h, Y_stress(:,1),'*-b','LineWidth',1);
    hold on
    plot(0:h, CI_stress(:,1:2),'k--','LineWidth',1);
    hold off
    axis tight;
    box on;
    ylim([-8, 4])
    xlim([0,10])
    title(Prior_mat(Prior_i));
end
lgd=legend([f(1),f(2)],'Median','Stress');
lgd.Location='southoutside';
lgd.Orientation='horizontal';
lgd.Layout.Tile = 'south';

filename=mainpath+"/StressTest.jpg";
saveas(fig,filename);
close all

