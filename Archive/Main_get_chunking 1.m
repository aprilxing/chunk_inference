clear; clc; close all

%% Load Demo data
addpath('src/');
load('example_data/dsp_example.mat', 'rt_er_data');

%% Inputs and conditions

% patient_id_simple = 'gRCS02';
% patient_id_simple = 'gRCS04';
patient_id_simple = 'gRCS05';
% patient_id_simple = 'RCS17';
% patient_id_simple = 'RCS20';

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
target = seq_data_log.Trials{1};

% Reaction time matrix = size(trials, elements)
rt_seq = zeros(length(target), 5);

% RT = onset of n+1 key - offset of n key
for i = 1:length(target) % Total trials 
    trial_data = target{i};
    
    for j = 1:5 % 5 elements 
        rt_seq(i, j) = trial_data.fk_onset(j+1) - trial_data.fk_offset(j);
    end 
end 


% Create er_seq -- psudo data for now 
er_seq = zeros(160, 5);
er_seq(:, 1) = rt_er_data.error(1:160);
er_seq(:, 2) = rt_er_data.error(161:320);
er_seq(:, 3) = rt_er_data.error(321:480);
er_seq(:, 4) = rt_er_data.error(481:640);
er_seq(:, 5) = rt_er_data.error(641:800);

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
figure(5);
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
subplot(2, 1, 2);

n_chunks = apply(@(x)(length(unique(x))), chunks);

% Expected number of chunk per trial (with some smoothing)
plot(smooth(gamma * n_chunks, 100), '-')
xlabel('Trial');
ylabel('Number of chunks');
title('Expected number of chunks per trial');
axis tight;

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



% figure(5).Position = [1 200 2240 1009];
fontsize(figure(4), 24, "points")
fontsize(figure(5), 24, "points")
