function mse_pls_group_cohort_interactions(base_path)
% Group-by-session interaction analysis for MSE data
% Tests whether Sham vs TBI differences change across sessions

% Add PLS toolbox to path
addpath(genpath('PLS'));

% Parameters
sessions = {'MSE_NLF_12m_chronic', 'MSE_NLF_6m_chronic'};
conditions = {'Sham', 'TBI'};
stages = ["Wake", "NREM", "REM"];
channels = [1, 2];
channel_names = struct('ch1', 'Parietal', 'ch2', 'Frontal');

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
fprintf('Loading MSE data across all sessions...\n');
all_data = struct();

for ch = 1:length(channels)
    channel = channels(ch);
    fprintf('Processing Channel %d (%s)...\n', channel, channel_names.(sprintf('ch%d', channel)));

    for s = 1:length(sessions)
        session = sessions{s};
        fprintf('  Loading session %s...\n', session);

        for c = 1:length(conditions)
            condition = conditions{c};
            condition_path = fullfile(base_path, session, condition, 'Results');

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

                                stage_idx_mask = (valid_stages == stage);

                                if sum(stage_idx_mask) >= 3 % Minimum 3 epochs
                                    mouse_avg = mean(mse_matrix(stage_idx_mask, :), 1);
                                    n_epochs = sum(stage_idx_mask);

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
                    all_data.(session).(condition).(field_name) = mouse_data;
                    fprintf('    %s Channel %d (%s): %d mouse-level averages\n', ...
                        condition, channel, channel_names.(field_name), size(mouse_data.mouse_averages, 1));
                end
            end
        end
    end
end

%% Create results directory
results_dir = fullfile(base_path, 'PLS_GroupSession_Interactions');
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

%% Group-by-Session Interaction Analysis
fprintf('\nRunning PLS Analysis...\n');

for ch = 1:length(channels)
    channel = channels(ch);
    field_name = sprintf('ch%d', channel);
    channel_label = channel_names.(field_name);

    for st = 1:length(stages)
        stage = stages(st);
        
        fprintf('Processing %s - %s...\n', stage, channel_label);
        
        % Check data completeness
        combination_complete = true;
        for s_check = 1:length(sessions)
            session_check = sessions{s_check};
            for c_check = 1:length(conditions)
                condition_check = conditions{c_check};
                
                actual_count = 0;
                if isfield(all_data, session_check) && isfield(all_data.(session_check), condition_check) && ...
                        isfield(all_data.(session_check).(condition_check), field_name)
                    
                    data_struct_check = all_data.(session_check).(condition_check).(field_name);
                    if ~isempty(data_struct_check.mouse_averages)
                        stage_idx_check = data_struct_check.stages == stage;
                        actual_count = sum(stage_idx_check);
                    end
                end
                
                expected_count = expected_n.(condition_check);
                if actual_count ~= expected_count
                    combination_complete = false;
                    break;
                end
            end
            if ~combination_complete
                break;
            end
        end
        
        if ~combination_complete
            fprintf('  Skipping - incomplete data\n');
            continue;
        end
        
        % Collect data for PLS analysis
        interaction_data = [];
        group_labels = [];
        session_labels = [];
        
        for s = 1:length(sessions)
            session = sessions{s};
            
            for c = 1:length(conditions)
                condition = conditions{c};
                
                if isfield(all_data, session) && isfield(all_data.(session), condition) && ...
                        isfield(all_data.(session).(condition), field_name)
                    
                    data_struct = all_data.(session).(condition).(field_name);
                    
                    if ~isempty(data_struct.mouse_averages)
                        stage_idx = data_struct.stages == stage;
                        
                        if sum(stage_idx) > 0
                            condition_data = data_struct.mouse_averages(stage_idx, :);
                            
                            interaction_data = [interaction_data; condition_data];
                            group_code = strcmp(condition, 'TBI');
                            group_labels = [group_labels; repmat(group_code, size(condition_data, 1), 1)];
                            session_labels = [session_labels; repmat(s, size(condition_data, 1), 1)];
                        end
                    end
                end
            end
        end
        
        % Run PLS analysis
        if ~isempty(interaction_data)
            % Clean data
            interaction_data(isnan(interaction_data) | isinf(interaction_data)) = 0;
            
            % Organize data for PLS
            pls_data = [];
            group_sizes = [];
            group_labels_pls = {};
            
            unique_sessions = unique(session_labels);
            unique_sessions = sort(unique_sessions);
            
            % All Sham first, then all TBI
            for s = 1:length(unique_sessions)
                session_idx_num = unique_sessions(s);
                session_name = sessions{session_idx_num};
                
                % Sham
                session_idx = session_labels == session_idx_num;
                group_idx = group_labels == 0; % Sham = 0
                combo_idx = session_idx & group_idx;
                
                if sum(combo_idx) >= 3
                    combo_data = interaction_data(combo_idx, :);
                    pls_data = [pls_data; combo_data];
                    group_sizes(end+1) = size(combo_data, 1);
                    group_labels_pls{end+1} = sprintf('Sham %s', session_name);
                end
            end
            
            for s = 1:length(unique_sessions)
                session_idx_num = unique_sessions(s);
                session_name = sessions{session_idx_num};
                
                % TBI
                session_idx = session_labels == session_idx_num;
                group_idx = group_labels == 1; % TBI = 1
                combo_idx = session_idx & group_idx;
                
                if sum(combo_idx) >= 3
                    combo_data = interaction_data(combo_idx, :);
                    pls_data = [pls_data; combo_data];
                    group_sizes(end+1) = size(combo_data, 1);
                    group_labels_pls{end+1} = sprintf('TBI %s', session_name);
                end
            end
            
            % PLS options
            if length(group_sizes) >= 2 && min(group_sizes) >= 3
                clear option
                option.method = 1;     % Mean-centred PLS
                option.num_perm = 500;
                option.num_boot = 1000;
                
                ncond = length(group_sizes);
                pls_result = pls_analysis({pls_data}, {group_sizes}, ncond, option);
                
                if isfield(pls_result, 'perm_result') && isfield(pls_result.perm_result, 'sprob')
                    pvals = pls_result.perm_result.sprob;
                    
                    fprintf('  P-values: ');
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
                        headline = sprintf('Group×Session: %s %s, LV%d p=%.4f', ...
                            stage, channel_label, LV, p);
                        
                        figure('Position', [100 100 1000 500]);
                        
                        % Bar plot
                        subplot(1,2,1)
                        z = pls_result.boot_result.orig_usc;
                        bar_handle = bar(z(:,LV));
                        hold on
                        
                        % Color bars
                        bar_colors = zeros(length(group_labels_pls), 3);
                        for i = 1:length(group_labels_pls)
                            if contains(group_labels_pls{i}, 'Sham')
                                bar_colors(i, :) = [0.3 0.6 0.9];
                            else
                                bar_colors(i, :) = [0.9 0.4 0.3];
                            end
                        end
                        bar_handle.CData = bar_colors;
                        
                        % Error bars
                        yneg = pls_result.boot_result.llusc(:,LV);
                        ypos = pls_result.boot_result.ulusc(:,LV);
                        errorbar(1:length(z), z(:,LV), yneg - z(:,LV), ypos - z(:,LV), '.', 'Color', 'black')
                        
                        xticks(1:length(group_labels_pls))
                        xticklabels(group_labels_pls)
                        xtickangle(45)
                        grid on
                        title('Session×Group Contrast Scores')
                        ylabel('Contrast Scores')
                        xlabel('Condition by Session')
                        
                        % Separator line between Sham and TBI groups
                        num_sham = sum(contains(group_labels_pls, 'Sham'));
                        xline(num_sham + 0.5, '--k', 'LineWidth', 2, 'Alpha', 0.7);
                        
                        % Heatmap
                        subplot(1,2,2)
                        x = pls_result.boot_result.compare_u(:,LV);
                        x(abs(x)<2.3) = 0;
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
                        fig_filename = fullfile(results_dir, sprintf('GroupSession_PLS_%s_%s_LV%d.png', ...
                            stage, channel_label, LV));
                        saveas(gcf, fig_filename)
                        
                        fig_filename_highres = fullfile(results_dir, sprintf('GroupSession_PLS_%s_%s_LV%d.fig', ...
                            stage, channel_label, LV));
                        saveas(gcf, fig_filename_highres)
                        
                        fprintf('  Saved: %s\n', fig_filename);
                    end
                    
                    % Save results
                    save_filename = fullfile(results_dir, sprintf('GroupSession_PLS_%s_%s.mat', ...
                        stage, channel_label));
                    save(save_filename, 'pls_result', 'pvals', 'group_sizes', 'conditions', ...
                        'sessions', 'stage', 'channel', 'channel_label', 'interaction_data', ...
                        'group_labels_pls', 'unique_sessions');
                end
            end
        end
    end
end

fprintf('\nAnalysis complete!\n');
end
