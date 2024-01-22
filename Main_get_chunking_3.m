clear; clc; close all

% This version of RRT is 1) removed the first RT, the cue period, and 2)
% combine same sequence and different days to make longer dataset (480
% trials per seq)


%% Inputs and conditions
addpath('src/');
load('example_data/dsp_example.mat');

% patient_id_simple = 'gRCS02';
% patient_id_simple = 'gRCS04';
% patient_id_simple = 'gRCS05';
% patient_id_simple = 'RCS17';
patient_id_simple = 'RCS20';

smooth_data_window = 10;

% Change from frame to ms
unit_multiplier = 1/500*1000; % 500 frames

% save_figs = 0;
save_figs = 1;
save_mat_data = 0;
% save_mat_data = 1;

condition = 'DBS-off High-Med';
parsing_condition = 'actual_key_presses_max time'; % Folder path

% trial_begin_condition = 'cue'; subfolder_name = 'Behavior_Starting from cue_begin';
trial_begin_condition = '1st key press'; subfolder_name = 'Behavior_Starting from 1st key press';

plot_mode = 'key_on_and_off'; % key begin and end will be strictly key presses
% plot_mode = 'overall_on_and_off'; % key begin and end will be overall key presses (consider capacitance and force)

% Load the target data
master_path = ['/Users/aprilxing/Library/CloudStorage/Box-Box/Motor Chunking Analysis/', patient_id_simple, '/', subfolder_name, '/', parsing_condition];
file_name = [patient_id_simple, '_seq_data_log, ', condition, ' Day 2-4', ', begin with ', trial_begin_condition, ', ', plot_mode, ', ', '07-Jan-2024']; %char(datetime('today'))
mat_name = [master_path, '/', file_name, '.mat'];
load(mat_name, 'seq_data_log')

%% Plot the IPIs over 3 days for each seq
h = figure;
h.Position = [1 1 1920 1029];

long_IPI = cell(1);

for i = 1:2 % 2 seqs
    max_Y = 1000; % ms
    num_trials = length(seq_data_log.IPIs_2el{1});

    long_IPI{i} = reshape([seq_data_log.IPIs_2el{i*(1:3)}], 3*num_trials, size(seq_data_log.IPIs_2el{1}, 2))*unit_multiplier;
    long_IPI{i} = [zeros(480, 1), long_IPI{i}];

    subplot(2, 1, i)
    plot(smoothdata(long_IPI{i}, "movmean", smooth_data_window), 'LineWidth', 2)
    ylim([0 max_Y])

    line([num_trials, num_trials], [0, max_Y], 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
    line([num_trials*2, num_trials*2], [0, max_Y], 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');

    xlabel('Trial over 3 days');
    ylabel('Reaction time (ms)');
    fig_name = [patient_id_simple, ', Reaction time over 3 days (480 trials) for sequence ', num2str(i), ', Smooth window = ', num2str(smooth_data_window)];
    title(fig_name);
    legend({'cue', 'Transition 1', 'Transition 2', 'Transition 3', 'Transition 4'}, 'Location', 'eastoutside')

    fontsize(h, 24, "points")
end

% Save the figure
if save_figs == 1
    saveas(gcf, [patient_id_simple, ', Reaction time over 3 days (480 trials), Smooth window = ', num2str(smooth_data_window), '.png']);
    saveas(gcf, [patient_id_simple, ', Reaction time over 3 days (480 trials), Smooth window = ', num2str(smooth_data_window), '.fig']);
end

%% For each sequence

for sequence = 1:2
    er_seq_all = cell(1);

    %% Build er_seq
    er_seq = cell(1);

    for number = (1:3) + (sequence-1)*3 % first 3 are seq 1, 4-6 are seq 3 (3 days each)
        % Choose target day and trial
        target = seq_data_log.Trials{number};

        er_seq{number} = zeros(length(target), 5);

        % RT = onset of n+1 key - offset of n key; er_seq = errored elements
        for i = 1:length(target) % Total trials
            trial_data = target{i};

            for j = 1:5 % 5 elements

                % Generate ER matrix
                checker = seq_data_log.Seqs{number}(j);

                if trial_data.response(j+1) ~= checker
                    er_seq{number}(i, j) = 1;
                end
            end
        end
    end

    er_seq_all{sequence} = er_seq;

    %% Fitting exponential decay model for RTs
    exponential_model = ['movement_time ~ a0' ...
        '+ a1*exp((b1/100)*(sequence_trial-1) + ' ...
        'b2/10*(sequence_trial-1)*(sequence_press-1))'];

    % initial_values = [0.18 0.38 -0.17, -0.1];
    % opts = statset('Display','off','TolFun',1e-5, ...
    %     'MaxIter', 100);

    initial_values = [0 0 0 0];
    opts = statset('Display','off','TolFun',1e-5, ...
        'MaxIter', 1000);

    rt_er_data = rt_er_data(1:2400, :);
    rt_er_data.movement_time = long_IPI{sequence}(:, 1);
    rt_er_data.day_count = [ones(800, 1)*1; ones(800, 1)*2; ones(800, 1)*3];
    rt_er_data.within_day_sequence_trial = repmat(repelem(1:160, 5), 1, 3)';
    rt_er_data.movement_time = reshape(long_IPI{sequence}.', [], 1);
    rt_er_data.sequence_trial = repelem(1:480, 5)';
    rt_er_data.sequence_press = repmat([1 2 3 4 5], 1, 480)';

    % Reshape error seq
    temp_er = cell(1);
    for i = 1:size(er_seq, 2)
        temp_er{i} = reshape(er_seq{i}.', [], 1);
    end
    rt_er_data.error = vertcat(temp_er{:});

    rt_er_data = rt_er_data(1:2400, :); % make sure size is correct

    % Fit
    nlmf = NonLinearModel.fit(rt_er_data, ...
        exponential_model, initial_values, 'Options', opts);

    % Preliminary visualizations
    h1 = figure;
    h1.Position = [1 1 1920 1029];
    clf;
    subplot(2, 1, 1);
    plot(rt_er_data.sequence_trial, rt_er_data.movement_time, '.');
    hold on;
    h = plot(rt_er_data.sequence_trial, smooth(rt_er_data.movement_time, 500), 'r-');
    xlabel('Trial over 3 days');
    ylabel('Reaction time (ms)');
    ylim([0 max_Y])

    line([num_trials, num_trials], [0, max_Y], 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
    line([num_trials*2, num_trials*2], [0, max_Y], 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');

    fig_name = [patient_id_simple, ', Original Reponse time, Sequence ', num2str(sequence)];
    title(fig_name);
    legend(h, 'smoothed reponse time');

    fontsize(gcf, 24, "points")


    subplot(2, 1, 2);
    plot(rt_er_data.sequence_trial, rt_er_data.error, '.');
    hold on;
    h = plot(rt_er_data.sequence_trial, smooth(rt_er_data.error, 500), 'r-');
    xlabel('Trial over 3 days');
    ylabel('Reaction time (ms)');
    ylim([0 1])

    line([num_trials, num_trials], [0, max_Y], 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
    line([num_trials*2, num_trials*2], [0, max_Y], 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');

    fig_name = [patient_id_simple, ', Original Errors, Sequence ', num2str(sequence)];
    title(fig_name);
    legend(h, 'smoothed error rate');

    fontsize(gcf, 24, "points")

    % Save the figure
    if save_figs == 1
        saveas(gcf, [patient_id_simple, ', Original Reponse time, Sequence ', num2str(sequence), '.png']);
        saveas(gcf, [patient_id_simple, ', Original Reponse time, Sequence ', num2str(sequence), '.fig']);
    end

    % Plot detrended response times
    h2 = figure;
    h2.Position = [1 1 1920 1029];
    clf;
    plot(rt_er_data.sequence_trial, table2array(nlmf.Residuals(:, 'Raw')), '.');
    hold on;
    plot(rt_er_data.sequence_trial, ...
        smooth(table2array(nlmf.Residuals(:, 'Raw')), 500), 'r-');
    xlabel('Trial over 3 days');
    ylabel('Residual response time (ms)');
    ylim([0 max_Y])

    line([num_trials, num_trials], [0, max_Y], 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
    line([num_trials*2, num_trials*2], [0, max_Y], 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');

    fig_name = [patient_id_simple, ', Residual response time, Sequence ', num2str(sequence)];
    title(fig_name);

    fontsize(gcf, 24, "points")

    % Save the figure
    if save_figs == 1
        saveas(gcf, [fig_name, '.png']);
        saveas(gcf, [fig_name, '.fig']);
    end

    %% Response time and error matrices

    % Create a response time and error matrix of trials vs element
    [rt_seq, er_seq] = mt_to_seq(rt_er_data, ...
        table2array(nlmf.Residuals(:, 'Raw')), ...
        rt_er_data.error);

    % Plot results
    h3 = figure;
    h3.Position = [1 1 1920 1029];
    clf;
    subplot(2, 1, 1);
    imagesc(rt_seq', [-0.25 0.25]);
    colormap('cool');
    colorbar;
    xlabel('Trial over 3 days');
    ylabel('Elements');

    line([num_trials, num_trials], [0, max_Y], 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
    line([num_trials*2, num_trials*2], [0, max_Y], 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');

    fig_name = [patient_id_simple, ', Response times, Sequence ', num2str(sequence)];
    title(fig_name);
    fontsize(gcf, 24, "points")

    subplot(2, 1, 2);
    imagesc(er_seq', [0 1]);
    colorbar;
    xlabel('Trial over 3 days');
    ylabel('Elements');

    line([num_trials, num_trials], [0, max_Y], 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
    line([num_trials*2, num_trials*2], [0, max_Y], 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');

    fig_name = [patient_id_simple, ', Errors, Sequence ', num2str(sequence)];
    title(fig_name);
    fontsize(gcf, 24, "points")

    % Save the figure
    if save_figs == 1
        saveas(gcf, [fig_name, '.png']);
        saveas(gcf, [fig_name, '.fig']);
    end
    %% Create space of chunk structures
    chunk_structures = create_chunks_nospace('n_seqlen', size(rt_seq, 2));
    figure(4);
    clf;
    % Plot space of chunks
    imagesc(chunk_structures');
    xlabel('Chunk structure index');
    ylabel('Element');
    colormap('jet');
    fig_name = ['All possible chunk structures, Sequence ', num2str(sequence)];
    title(fig_name);
    fontsize(figure(4), 24, "points")

    % Save the figure
    if save_figs == 1
        saveas(gcf, ['All possible chunk structures, Sequence ', num2str(sequence), '.png']);
        saveas(gcf, ['All possible chunk structures, Sequence ', num2str(sequence), '.fig']);
    end
    
    %% Find which chunk structure is present at each trial
    [rho, self_t, log_like, fm, T, rho_er, v, v_er, ...
        initial_dist, mean_pause, mean_inchunk, ...
        mean_pause_er, mean_inchunk_er, ...
        chunks, cor_chunks, gamma] = ...
        chunk_hmm_learn_param(rt_seq, er_seq, 'verbose', true, ...
        'fit_rt', true, 'fit_rt_rt', true, 'fit_er', true, 'fit_er_er', true, ...
        'fit_T', true, 'fit_rho', true, 'fit_rho_er', true, ...
        'chunks', chunk_structures);

    % compute mean and covariance of each chunk structure
    [chunk_means_rt, rt_cov, chunk_means_er, er_cov] = ...
        create_chunk_means_covs(chunks, cor_chunks, ...
        mean_pause, mean_inchunk, v, rho, ...
        mean_pause_er, mean_inchunk_er, v_er, rho_er);

    %% Plot results of algorithm

    % Expected chunking structure
    h1 = figure(5);
    h1.Position = [1 192 1470 785];

    clf;
    subplot(2, 1, 1);

    % For each trial, get the most possible chunking structure
    final_chunk_structure = zeros(length(gamma), 5);

    [maxValues, rowIndices] = max(gamma,[], 2);

    for i = 1:length(gamma)
        final_chunk_structure(i, :) = chunks(rowIndices(i), :)';
    end

    imagesc(final_chunk_structure');
    colormap('jet');
    xlabel('Trial');
    ylabel('Element');
    fig_name = ['Expected chunking structure, Sequence ', num2str(sequence)];
    title(fig_name);
    colorbar;

    subplot(2, 1, 2);
    n_chunks = apply(@(x)(length(unique(x))), chunks);

    % Expected number of chunk per trial (with some smoothing)
    plot(smooth(gamma * n_chunks, 100), '-')
    xlabel('Trial');
    ylabel('Number of chunks');
    fig_name = ['Expected number of chunks per trial, Sequence ', num2str(sequence)];
    title(fig_name);
    axis tight;
    colorbar;
    fontsize(figure(5), 24, "points")

    % Save the figure
    if save_figs == 1
        saveas(gcf, [fig_name, '.png']);
        saveas(gcf, [fig_name, '.fig']);
    end


    % Meean response times and errors fitted by model
    expected_rt = gamma * chunk_means_rt;
    expected_er = gamma * chunk_means_er;
    h6 = figure;
    h6.Position = [1 1 1920 1029];
    clf;
    subplot(2, 1, 1);
    imagesc(expected_rt', [-0.25 0.25]);
    colormap('cool');
    colorbar;
    xlabel('Trial over 3 days');
    ylabel('Element');

    line([num_trials, num_trials], [0, max_Y], 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
    line([num_trials*2, num_trials*2], [0, max_Y], 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');

    fig_name = [patient_id_simple, ', Expected response times, Sequence ', num2str(sequence)];
    title(fig_name);
    fontsize(gcf, 24, "points")
    
    subplot(2, 1, 2);
    imagesc(expected_rt', [0 0.2]);
    colorbar;
    xlabel('Trial over 3 days');
    ylabel('Element');

    line([num_trials, num_trials], [0, max_Y], 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
    line([num_trials*2, num_trials*2], [0, max_Y], 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');

    fig_name = [patient_id_simple, ', Expected error rate, Sequence ', num2str(sequence)];
    title(fig_name);
    fontsize(gcf, 24, "points")

    % Save the figure
    if save_figs == 1
        saveas(gcf, [fig_name, '.png']);
        saveas(gcf, [fig_name, '.fig']);
    end
end

close all
