clear all
close all
clc

%Parameters
h_choice=1;
eq_choice=1;
T_Mat=96:224; %the expanding window for forecast (Original: Start in T=1:96 and expand window until T=1:224->T_Mat=96:224)
fig_row_num=2;
fig_col_num=3;
setfontsize_main=18;
plot_prior=["BRW","BQR","HSBQR","QVP_HS","LQVP_HS"];

%Paths
mainpath=cd;
dataloc=mainpath+"/data/";
funcs=fileparts(mainpath)+"/QVP_package/";

addpath(funcs,dataloc);

%% CODE STARTS HERE

% Read data
xlsfilename = 'DataQVAR.xlsx';
[X, TEXT]   = xlsread(xlsfilename);


%% Read forecast
ForcFile=mainpath+"/ZZZ_Forc_T="+T_Mat(1)+".mat";
load(ForcFile); %Just load for priors and taus

Prior_mat=Prior_choice;
for mod_i=1:length(Prior_mat)
    if strcmp(Prior_mat(mod_i),"QVP_HS")==1
        Prior_mat(mod_i)="QVP_{HS}";
    elseif strcmp(Prior_mat(mod_i),"LQVP_HS")==1
        Prior_mat(mod_i)="NCQVP_{HS}";
    elseif strcmp(Prior_mat(mod_i),"LQVP_HS_SAVS")==1
        Prior_mat(mod_i)="NCQVP_{HS,SAVS}";
    elseif strcmp(Prior_mat(mod_i),"LQVP_HS_ASIS")==1
        Prior_mat(mod_i)="NCQVP_{HS,ASIS}";
    elseif strcmp(Prior_mat(mod_i),"LQVP_HS_ASIS_SAVS")==1
        Prior_mat(mod_i)="NCQVP_{HS,ASIS,SAVS}";
    end
end


Y_forc_use=nan(length(T_Mat),length(taus)+1,length(Prior_choice)); %For eq_choice only
Y_forc_fullmat=nan(length(T_Mat),length(taus)+1,2,length(Prior_choice)); %For all
for mod_i=1:length(Prior_choice)
    it_loc=0;
    for t=T_Mat
        it_loc=it_loc+1;
        ForcFile=mainpath+"/ZZZ_Forc_T="+t+".mat";
        load(ForcFile)

        for eq_i=1:2
            Y_temp=squeeze(Forecast_Save(:,:,eq_i,mod_i));
            Y_use=[Y_forc(h_choice+1,eq_i),Y_temp(h_choice+1,:)];
            Y_forc_fullmat(it_loc,:,eq_i,mod_i)=Y_use;
            if eq_i==eq_choice
                Y_forc_use(it_loc,:,mod_i)=Y_use;
            end
        end
    end
end

%Figure
date=datetime(TEXT(2:end,1));
date=date(2:end,1);
date_use=2007:1/12:2020;
date_use=date_use(1:length(T_Mat));
medloc=find(0.5==round(taus,2));
fig=figure();
set(gcf, 'Position', [1,1,1920,1080],'DefaultAxesFontSize',setfontsize_main);
tiledlayout(fig_row_num,fig_col_num)
for i=1:length(Prior_choice)
    curr_prior=Prior_choice(i);

    if sum(strcmp(curr_prior,plot_prior))==1

        nexttile;
        temp=squeeze(Y_forc_use(:,:,i));
        temp_y=temp(:,1);
        temp_mat=temp(:,2:end);

        plot(date_use,temp(:,1),Color="#A2142F") %Maroon
        hold on
        plot(date_use,temp_mat(:,medloc),LineStyle="--",Color="#000000") %Black
        for j=1:(length(taus)-1)/2
            inBetween = [temp_mat(:,j); flipud(temp_mat(:,end-(j-1)))];
            patch([date_use'; flipud(date_use')], inBetween, 'b','FaceAlpha',.1,'EdgeAlpha',0);
        end
        xlim([min(date_use), max(date_use)]);
        if eq_choice==1
            ylim([-6, 4]);
        else
            ylim([0, 1]);
        end
        %xlim([min(date_use), 2014]);
        hold off

        title(Prior_mat(i))
    end
end
if h_choice==1
    file=mainpath+"/ForcFit_h=1.jpg";
    saveas(fig,file);
end
close all

%wqs
%EqSpec
wqs_mat=[];
for eq_i=1:2
    wqs_mat=[wqs_mat;["Equation "+eq_i,"CRPS","Centre","Left","Right"]];
    for i=1:length(Prior_choice)
        temp=squeeze(Y_forc_fullmat(:,:,eq_i,i));
        w_mat=Prior_mat(i);
        for w=1:4
            w_mat=[w_mat,mean(wQS(temp(:,1),temp(:,2:end),taus,w))];
        end
        wqs_mat=[wqs_mat;w_mat];
    end
end
%Overall
wqs_mat=[wqs_mat;["Overall","CRPS","Centre","Left","Right"]];
for i=1:length(Prior_choice)
    w_mat=Prior_mat(i);
    for w=1:4
        w_mat=[w_mat,mean([double(wqs_mat(i+1,w+1)),double(wqs_mat(i+length(Prior_mat)+eq_i,w+1))])];
    end
    wqs_mat=[wqs_mat;w_mat];
end
Sheet_Name="h="+h_choice;
file=mainpath+"/ForcRes.xlsx";
writematrix(wqs_mat,file,Sheet=Sheet_Name);