%% In-Sample Fit: BQR vs QVP_HS  (p-lag version)
% Estimates BQR and QVP_HS on the EU data with a user-specified number of
% lags p and plots the in-sample fit (X*beta_q) for three quantiles
% against the observed outcome for IP growth and CISS.

clear all; close all; clc;

%% Settings
p          = 1;                        % <--- number of lags
taus_full  = 0.1:0.1:0.9;          % full quantile grid
tau_plot   = [0.10, 0.50, 0.90];      % quantiles to highlight in plots
iter       = 1500;
burn       = 1000;
Prior_list   = ["BQR","QVP_HS"];
Prior_labels = ["BQR","QVP_{HS}"];

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
tau_idx = zeros(1, length(tau_plot));
for ti = 1:length(tau_plot)
    [~, tau_idx(ti)] = min(abs(taus_full - tau_plot(ti)));
end

%% Estimation ----------------------------------------------------------
nPr   = length(Prior_list);
nJobs = nPr * K;
qfit_store = cell(nPr, K);

% Linearise the (prior, equation) grid for parfor
pr_vec  = repelem(1:nPr, K);
eq_vec  = repmat(1:K, 1, nPr);

parfor idx = 1:nJobs
    pr   = pr_vec(idx);
    eq_i = eq_vec(idx);
    fprintf('  [job %d/%d]  %s  -  eq %d  (p=%d) ...\n', ...
            idx, nJobs, Prior_list(pr), eq_i, p);
    [~, ~, qfit, ~, ~, ~, ~, ~, ~] = ...
        BQVAR_Estim_plag(Y, taus_full, eq_i, iter, burn, Prior_list(pr), p);
    qfit_store{idx} = qfit;
end

% Reshape the linear cell back to (nPr x K)
qfit_store = reshape(qfit_store, K, nPr)';

%% Plotting ------------------------------------------------------------
date_fit = date(p+1:end);             % p obs lost to lags
Y_fit    = Y(p+1:end, :);

colours  = lines(length(tau_plot));
lw_q     = 1.4;

for eq_i = 1:K
    fig = figure('Position', [50, 50, 1400, 800], 'Color', 'w');
    sgtitle(sprintf('In-sample quantile fit  -  %s  (p = %d)', ...
            Headers{eq_i}, p), 'FontSize', 20, 'FontWeight', 'bold');

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

        % Highlight quantile crossings
        qf_plot = qf(:, tau_idx);
        cross_mask = false(size(qf_plot, 1), 1);
        for j = 1:length(tau_plot)-1
            cross_mask = cross_mask | (qf_plot(:, j) > qf_plot(:, j+1));
        end
        if any(cross_mask)
            yl = get(gca, 'YLim');
            cross_idx = find(cross_mask);
            for ci = 1:length(cross_idx)
                h_cr = patch([date_fit(cross_idx(ci))-15, date_fit(cross_idx(ci))+15, ...
                              date_fit(cross_idx(ci))+15, date_fit(cross_idx(ci))-15], ...
                             [yl(1), yl(1), yl(2), yl(2)], ...
                             'r', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
            end
        end

        hold off;

        datetick('x', 'yyyy');
        ylabel(Headers{eq_i}, 'FontSize', 14);
        xlabel('', 'FontSize', 14);
        title(Prior_labels(pr), 'FontSize', 16, 'Interpreter', 'tex');
        set(gca, 'FontSize', 16);
        axis tight;

        % Legend in the first subplot only
        if pr == 1
            leg_labels = arrayfun(@(t) sprintf('\\tau = %.2f', t), ...
                                  tau_plot, 'UniformOutput', false);
            legend([h_obs; h_q], [{'Observed'}, leg_labels], ...
                   'Location', 'best', 'FontSize', 16);
        end
    end

    % Save figure
    filename = fullfile(mainpath, sprintf('InSampleFit_%s_p%d.jpg', ...
               strrep(Headers{eq_i}, ' ', '_'), p));
    print(fig, filename, '-djpeg', '-r300');
    fprintf('Saved %s\n', filename);
end

fprintf('\nDone.\n');
