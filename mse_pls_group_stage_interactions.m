function mse_pls_group_stage_interactions(base_path)
% Group comparison analysis for MSE data by sleep stage
% Tests Sham vs TBI differences separately for each sleep stage

% Add PLS toolbox to path
addpath(genpath('PLS'));

% Parameters
conditions = {'Sham', 'TBI'};
stages = ["Wake", "NREM", "REM"];
channels = [1, 2];

% Expected sample sizes
expected_n = struct('Sham', 10, 'TBI', 10);

%% Define RGB colormap
rgb = [ ...
    94    79   162
    72   104   175
    50   136   189
    76   160   177
    102   194   165
    136   210   165
    171   221   164
    200   233   158
    230   245   152
    243   250   199
    255   255   255
    255   230   195
    254   224   139
    253   200   120
    253   174    97
    249   141    82
    244   109    67
    230    87    88
    213    62    79
    182    33    72
    158     1    66  ] / 255;

%% Load MSE data
fprintf('Loading MSE data...\n');
all_data = struct();

for ch = 1:length(channels)
    channel = channels(ch);
    fprintf('Processing Channel %d...\n', channel);

    for c = 1:length(conditions)
        condition = conditions{c};
        condition_path = fullfile(base_path, condition, 'Results');

        if exist(condition_path, 'dir')
            dataset_files = dir(fullfile(condition_path, '*_dataset.mat'));

            if ~isempty(dataset_files)
                mouse_data = struct();
                mouse_data.mouse_averages = [];
                mouse_data.mouse_ids = {};
                mouse_data.stages = [];
                mouse_data.scales = [];
                mouse_data.n_epochs_per_mouse = [];

                for f = 1:length(dataset_files)
                    file_path = fullfile(condition_path, dataset_files(f).name);
                    [~, filename, ~] = fileparts(dataset_files(f).name);                 
                    filename_parts = strsplit(filename, '_');
                    mouse_id = filename_parts{1};

                    data = load(file_path);
                    dataset = data.dataset;

                    % Get valid epochs based on channel
                    if channel == 1
                        valid_idx = ~cellfun(@isempty, dataset.MSE_EEG1);
                        valid_mse = dataset.MSE_EEG1(valid_idx);
                    else
                        valid_idx = ~cellfun(@isempty, dataset.MSE_EEG2);
                        valid_mse = dataset.MSE_EEG2(valid_idx);
                    end

                    valid_stages = dataset.Sleep_Stage(valid_idx);

                    % Convert MSE cell array to matrix
                    if ~isempty(valid_mse)
                        n_valid = length(valid_mse);
                        mse_matrix = zeros(n_valid, length(dataset.mse_scales));
                        for i = 1:n_valid
                            if ~isempty(valid_mse{i})
                                mse_matrix(i, :) = valid_mse{i};
                            end
                        end

                        % Average across epochs for each stage
                        unique_stages = unique(valid_stages);

                        for stage_idx = 1:length(unique_stages)
                            stage = unique_stages(stage_idx);

                            combo_idx = (valid_stages == stage);

                            if sum(combo_idx) >= 3 % Minimum 3 epochs
                                mouse_avg = mean(mse_matrix(combo_idx, :), 1);
                                n_epochs = sum(combo_idx);

                                % Store mouse-level data
                                mouse_data.mouse_averages = [mouse_data.mouse_averages; mouse_avg];
                                mouse_data.mouse_ids{end+1} = sprintf('%s_%s', mouse_id, stage);
                                mouse_data.stages = [mouse_data.stages; stage];
                                mouse_data.n_epochs_per_mouse = [mouse_data.n_epochs_per_mouse; n_epochs];
                                mouse_data.scales = dataset.mse_scales;
                            end
                        end
                    end
                end

                % Store mouse-level data
                field_name = sprintf('ch%d', channel);
                all_data.(condition).(field_name) = mouse_data;
                fprintf('    %s Channel %d: %d mouse-level averages\n', ...
                    condition, channel, size(mouse_data.mouse_averages, 1));
            end
        end
    end
end

%% Create results directory
results_dir = fullfile(base_path, 'PLS_ByStage');
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

%% PLS Analysis: Separate for each sleep stage
fprintf('\nRunning PLS Analysis...\n');

for ch = 1:length(channels)
    channel = channels(ch);
    field_name = sprintf('ch%d', channel);
    
    fprintf('Processing Channel %d...\n', channel);
    
    % Loop through each sleep stage separately
    for st = 1:length(stages)
        stage = stages(st);
        fprintf('  Analyzing %s...\n', stage);
        
        % Check data completeness for this stage
        stage_complete = true;
        for c_check = 1:length(conditions)
            condition_check = conditions{c_check};
            
            actual_count = 0;
            if isfield(all_data, condition_check) && ...
                    isfield(all_data.(condition_check), field_name)
                
                data_struct_check = all_data.(condition_check).(field_name);
                if ~isempty(data_struct_check.mouse_averages)
                    stage_idx_check = data_struct_check.stages == stage;
                    actual_count = sum(stage_idx_check);
                end
            end
            
            expected_count = expected_n.(condition_check);
            if actual_count ~= expected_count
                stage_complete = false;
                fprintf('    Warning: %s %s has %d subjects (expected %d)\n', ...
                    condition_check, stage, actual_count, expected_count);
                break;
            end
        end
        
        if ~stage_complete
            fprintf('    Skipping %s - incomplete data\n', stage);
            continue;
        end
        
        % Collect data for this stage only
        stage_data = [];
        group_labels_pls = {};
        group_sizes = [];
        
        for c = 1:length(conditions)
            condition = conditions{c};
            
            if isfield(all_data, condition) && ...
                    isfield(all_data.(condition), field_name)
                
                data_struct = all_data.(condition).(field_name);
                
                if ~isempty(data_struct.mouse_averages)
                    stage_idx = data_struct.stages == stage;
                    
                    if sum(stage_idx) > 0
                        condition_data = data_struct.mouse_averages(stage_idx, :);
                        
                        stage_data = [stage_data; condition_data];
                        group_sizes(end+1) = size(condition_data, 1);
                        group_labels_pls{end+1} = condition;
                    end
                end
            end
        end
        
        % Run PLS analysis for this stage
        if ~isempty(stage_data) && length(group_sizes) == 2 && min(group_sizes) >= 3
            % Clean data
            stage_data(isnan(stage_data) | isinf(stage_data)) = 0;
            
            % PLS options
            clear option
            option.method = 1;     % Mean-centred PLS
            option.num_perm = 500;
            option.num_boot = 1000;
            
            ncond = length(group_sizes);
            pls_result = pls_analysis({stage_data}, {group_sizes}, ncond, option);
            
            if isfield(pls_result, 'perm_result') && isfield(pls_result.perm_result, 'sprob')
                pvals = pls_result.perm_result.sprob;
                
                fprintf('    P-values: ');
                for pv = 1:length(pvals)
                    if pvals(pv) < 0.001
                        fprintf('LV%d<0.001 ', pv);
                    elseif pvals(pv) < 0.05
                        fprintf('LV%d=%.4f* ', pv, pvals(pv));
                    else
                        fprintf('LV%d=%.4f ', pv, pvals(pv));
                    end
                end
                fprintf('\n');
                
                % Plot LV1 only
                if length(pvals) >= 1
                    LV = 1;
                    p = pvals(LV);
                    headline = sprintf('%s: Sham vs TBI, Ch%d, LV%d p=%.4f', ...
                        stage, channel, LV, p);
                    
                    figure('Position', [100 100 1000 500]);
                    
                    % Bar plot
                    subplot(1,2,1)
                    z = pls_result.boot_result.orig_usc;
                    bar_handle = bar(z(:,LV));
                    hold on
                    
                    % Color bars
                    bar_colors = [0.3 0.6 0.9; 0.9 0.4 0.3]; % Sham blue, TBI red
                    bar_handle.CData = bar_colors;
                    
                    % Error bars
                    yneg = pls_result.boot_result.llusc(:,LV);
                    ypos = pls_result.boot_result.ulusc(:,LV);
                    errorbar(1:length(z), z(:,LV), yneg - z(:,LV), ypos - z(:,LV), '.', 'Color', 'black')
                    
                    xticks(1:length(group_labels_pls))
                    xticklabels(group_labels_pls)
                    xtickangle(45)
                    grid on
                    title('Group Contrast Scores')
                    ylabel('Contrast Scores')
                    xlabel('Condition')
                    
                    % Heatmap
                    subplot(1,2,2)
                    x = pls_result.boot_result.compare_u(:,LV);
                    x(abs(x)<2.3) = 0;
                    data_struct = all_data.(conditions{1}).(field_name);
                    scales = data_struct.scales;
                    plotdata = reshape(x, length(scales), []);
                    
                    imagesc(plotdata)
                    colormap(rgb)
                    clim([-7 7])
                    colorbar
                    yticks(1:5:length(scales))
                    yticklabels(1:5:length(scales))
                    xticks([])
                    title(sprintf('Salience LV%d', LV))
                    ylabel('Temporal Scale')
                    
                    sgtitle(headline, 'FontSize', 16, 'Interpreter', 'none')
                    set(gcf, 'Color', 'w')
                    
                    % Save figure
                    fig_filename = fullfile(results_dir, sprintf('%s_PLS_ch%d_LV%d.png', ...
                        stage, channel, LV));
                    saveas(gcf, fig_filename)
                    
                    fig_filename_highres = fullfile(results_dir, sprintf('%s_PLS_ch%d_LV%d.fig', ...
                        stage, channel, LV));
                    saveas(gcf, fig_filename_highres)
                    
                    fprintf('    Saved: %s\n', fig_filename);
                end
                
                % Save results
                save_filename = fullfile(results_dir, sprintf('%s_PLS_ch%d.mat', ...
                    stage, channel));
                save(save_filename, 'pls_result', 'pvals', 'group_sizes', 'conditions', ...
                    'stage', 'channel', 'stage_data', 'group_labels_pls');
            end
        end
    end
end

fprintf('\nAnalysis complete!\n');
end
