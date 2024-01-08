clear; clc; close all

%% Inputs and conditions
addpath('src/');

% patient_id_simple = 'gRCS02';
% patient_id_simple = 'gRCS04';
patient_id_simple = 'gRCS05';
% patient_id_simple = 'RCS17';
% patient_id_simple = 'RCS20';

save_figs = 1;

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


for number = 2
    %% Build data structure

    % % Preallocate sample_table (temporarily called rt_er_data table to correspond with the demo file)
    % subject_id = nan(1, 1);
    % day_count = nan(1, 1);
    % sequence_id = nan(1, 1);
    % within_day_sequence_trial = nan(1, 1);
    % movement_time = nan(1, 1);
    % sequence_trial = nan(1, 1);
    % sequence_press = nan(1, 1);
    % error = nan(1, 1);
    % rt_er_data = table(subject_id, day_count, sequence_id, within_day_sequence_trial, movement_time, sequence_trial, sequence_press, error);


    % Choose target day and trial
    target = seq_data_log.Trials{number};

    % Reaction time matrix = size(trials, elements)
    rt_seq = zeros(length(target), 5);
    er_seq = zeros(length(target), 5);

    % RT = onset of n+1 key - offset of n key; er_seq = errored elements
    for i = 1:length(target) % Total trials
        trial_data = target{i};

        for j = 1:5 % 5 elements
            rt_seq(i, j) = trial_data.fk_onset(j+1) - trial_data.fk_offset(j);

            % Generate ER matrix
            checker = seq_data_log.Seqs{number}(j);

            if trial_data.response(j+1) ~= checker
                er_seq(i, j) = 1;
            end
        end
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
    title('All possible chunk structures');
    fontsize(figure(4), 24, "points")

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
    title('Expected chunking structure');
    colorbar;
    
    subplot(2, 1, 2);
    n_chunks = apply(@(x)(length(unique(x))), chunks);

    % Expected number of chunk per trial (with some smoothing)
    plot(smooth(gamma * n_chunks, 100), '-')
    xlabel('Trial');
    ylabel('Number of chunks');
    title('Expected number of chunks per trial');
    axis tight;
    colorbar;
    fontsize(figure(5), 24, "points")

    % Save the figure
    if save_figs == 1
        saveas(gcf, [patient_id_simple, num2str(number), '.png']);
        saveas(gcf, [patient_id_simple, num2str(number), '.fig']);
    end

    % Meean response times and errors fitted by model
    expected_rt = gamma * chunk_means_rt;
    expected_er = gamma * chunk_means_er;
    figure(6);
    clf;
    subplot(2, 1, 1);
    imagesc(expected_rt', [-0.25 0.25]);
    colormap('cool');
    colorbar;
    xlabel('Trial');
    ylabel('Element');
    title('Expected response times');
    subplot(2, 1, 2);
    imagesc(expected_rt', [0 0.2]);
    colorbar;
    xlabel('Trial');
    ylabel('Element');
    title('Expected error rate');

end