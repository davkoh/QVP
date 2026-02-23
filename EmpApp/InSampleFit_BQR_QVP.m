%% In-Sample Fit: BQR vs QVP_HS
% Estimates BQR and QVP_HS on the EU data and plots the in-sample quantile
% fit (X*beta_q) for three quantiles against the observed outcome for both
% IP growth and CISS.

clear all; close all; clc;

%% Settings
taus_full  = 0.05:0.05:0.95;          % full quantile grid (needed by estimator)
tau_plot   = [0.10, 0.50, 0.90];      % quantiles to highlight in plots
iter       = 15000;
burn       = 10000;
Prior_list = ["BQR","QVP_HS"];

%% Paths
mainpath = cd;
dataloc  = mainpath + "/data/";
funcs    = fileparts(mainpath) + "/QVP_package/";
addpath(funcs, dataloc);

%% Load data (same transformation as A1_BQVAR_main)
xlsfilename = 'DataQVAR.xlsx';
[X, TEXT]    = xlsread(xlsfilename);
date         = datenum(TEXT(2:end,1));

Headers = {'IP growth', 'CISS'};
Y       = [100*(log(X(2:end,1)) - log(X(1:end-1,1))),  X(2:end,2)];
date    = date(2:end,1);

[T, K] = size(Y);
Q      = length(taus_full);

% Indices of the three plot-quantiles inside the full grid
[~, tau_idx] = ismember(tau_plot, taus_full);

%% Estimation ----------------------------------------------------------
% Store results: rows = priors, columns = equations
qfit_store = cell(length(Prior_list), K);

for pr = 1:length(Prior_list)
    fprintf('\n===  Estimating %s  ===\n', Prior_list(pr));
    for eq_i = 1:K
        fprintf('  Equation %d / %d  (%s) ...\n', eq_i, K, Headers{eq_i});
        [~, ~, qfit, ~, ~, ~, ~, ~, ~] = ...
            BQVAR_Estim(Y, taus_full, eq_i, iter, burn, Prior_list(pr));
        qfit_store{pr, eq_i} = qfit;          % (T-1) x Q
    end
end

%% Plotting ------------------------------------------------------------
date_fit = date(2:end);   % qfit is based on Y(2:end,:), one obs lost for lag
Y_fit    = Y(2:end,:);    % observed outcomes aligned with qfit

colours  = lines(length(tau_plot));   % distinct colours per quantile
lw_q     = 1.4;                       % line width for quantile fits

for eq_i = 1:K
    fig = figure('Position', [50, 50, 1400, 800], 'Color', 'w');
    sgtitle(sprintf('In-sample quantile fit  -  %s', Headers{eq_i}), ...
            'FontSize', 16, 'FontWeight', 'bold');

    for pr = 1:length(Prior_list)
        subplot(length(Prior_list), 1, pr);
        hold on; box on;

        % Observed series
        h_obs = plot(date_fit, Y_fit(:, eq_i), 'k', 'LineWidth', 1.6);

        % Quantile fits
        qf = qfit_store{pr, eq_i};
        h_q = gobjects(length(tau_plot), 1);
        for j = 1:length(tau_plot)
            h_q(j) = plot(date_fit, qf(:, tau_idx(j)), ...
                          'Color', colours(j,:), 'LineWidth', lw_q, ...
                          'LineStyle', '--');
        end
        hold off;

        datetick('x', 'yyyy');
        ylabel(Headers{eq_i});
        title(Prior_list(pr), 'FontSize', 14);
        axis tight;

        % Legend in the first subplot only
        if pr == 1
            leg_labels = arrayfun(@(t) sprintf('\\tau = %.2f', t), ...
                                  tau_plot, 'UniformOutput', false);
            legend([h_obs; h_q], ['Observed', leg_labels], ...
                   'Location', 'best', 'FontSize', 10);
        end
    end

    % Save figure
    filename = fullfile(mainpath, sprintf('InSampleFit_%s.jpg', ...
               strrep(Headers{eq_i}, ' ', '_')));
    saveas(fig, filename);
    fprintf('Saved %s\n', filename);
end

fprintf('\nDone.\n');
