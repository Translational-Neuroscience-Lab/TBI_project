function plot_sessions_mse_comparison_cohorts(base_path, sessions, conditions)
% Plot MSE curves for individual sessions comparing TBI vs Sham conditions across sleep stages
% Shows mean curves with shaded SEM calculated from individual subjects
sessions = {'MSE_NLF_12m_chronic', 'MSE_NLF_6m_chronic'};
conditions = {'TBI', 'Sham'};

fprintf('Using subject-based SEM with shaded error regions...\n');

% Initialize storage
all_data = struct();
stages = ["Wake", "NREM", "REM"];
% Colors: TBI (red shades), Sham (blue shades)
% 6m group uses darker shades, 12m group uses lighter shades
colors_6m = {[0.63, 0, 0], [0, 0, 0.5]};  % Dark red and Dark blue for 6m
colors_12m = {[1, 0.4, 0.4], [0.4, 0.4, 1]};  % Light red and Light blue for 12m
channels = [1, 2];
channel_names = {'Parietal', 'Frontal'};

fprintf('Loading data from multiple sessions for both channels...\n');

% Load all data for both channels
for ch = 1:length(channels)
    channel = channels(ch);
    channel_name = channel_names{ch};
    fprintf('Processing Channel %d (%s)...\n', channel, channel_name);

    for s = 1:length(sessions)
        session = sessions{s};
        fprintf('  Processing %s...\n', session);

        for c = 1:length(conditions)
            condition = conditions{c};

            % Build path to Results folder
            condition_path = fullfile(base_path, session, condition, 'Results');

            if exist(condition_path, 'dir')
                % Find all dataset files in this condition
                dataset_files = dir(fullfile(condition_path, '*_dataset.mat'));

                if ~isempty(dataset_files)
                    % Initialize arrays for this session/condition/channel
                    session_condition_data = struct();
                    session_condition_data.individual_means = [];  % Store mean MSE per subject
                    session_condition_data.scales = [];
                    session_condition_data.n_subjects = 0;

                    % Load and compute mean for each subject
                    for f = 1:length(dataset_files)
                        file_path = fullfile(condition_path, dataset_files(f).name);

                        try
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
                                
                                % Store individual subject data by sleep stage
                                subject_data = struct();
                                subject_data.mse_data = mse_matrix;
                                subject_data.stages = valid_stages;
                                subject_data.scales = dataset.mse_scales;
                                subject_data.filename = dataset_files(f).name;
                                
                                % Calculate mean MSE for each sleep stage for this subject
                                for st = 1:length(stages)
                                    stage = stages(st);
                                    stage_idx = valid_stages == stage;
                                    
                                    if sum(stage_idx) > 0
                                        stage_mse = mse_matrix(stage_idx, :);
                                        subject_mean = mean(stage_mse, 1);
                                        
                                        % Initialize storage if needed
                                        if ~isfield(session_condition_data, char(stage))
                                            session_condition_data.(char(stage)) = [];
                                        end
                                        
                                        % Store this subject's mean for this stage
                                        session_condition_data.(char(stage)) = [session_condition_data.(char(stage)); subject_mean];
                                    end
                                end
                                
                                session_condition_data.scales = dataset.mse_scales;
                                session_condition_data.n_subjects = session_condition_data.n_subjects + 1;
                            end

                        catch ME
                            fprintf('Warning: Could not load %s: %s\n', dataset_files(f).name, ME.message);
                        end
                    end

                    % Store combined data with channel info
                    field_name = sprintf('ch%d', channel);
                    all_data.(session).(condition).(field_name) = session_condition_data;
                    fprintf('    %s Channel %d (%s): %d subjects loaded\n', condition, channel, ...
                           channel_name, session_condition_data.n_subjects);
                    
                else
                    fprintf('    Warning: No dataset files found in %s\n', condition_path);
                end
            else
                fprintf('    Warning: Directory not found: %s\n', condition_path);
            end
        end
    end
end

% Create plots - combined across sessions
fprintf('\nCreating combined session plots...\n');

for ch = 1:length(channels)
    channel = channels(ch);
    channel_name = channel_names{ch};
    fprintf('Processing Channel %d (%s) plots...\n', channel, channel_name);

    % Create mse_plots directory in base path
    plot_dir = fullfile(base_path, 'mse_plots_combined');
    if ~exist(plot_dir, 'dir')
        mkdir(plot_dir);
    end

    % Create plots for each sleep stage (combining all sessions)
    for st = 1:length(stages)
        stage = stages(st);
        
        fprintf('  Creating %s plot for Channel %d (%s)...\n', stage, channel, channel_name);
        
        % Create individual plot for this stage
        figure('Position', [100 100 900 600], 'Visible', 'off');

        legend_entries = {};
        legend_handles = [];

        % Loop through sessions first, then conditions
        for s = 1:length(sessions)
            session = sessions{s};
            
            % Extract session label (6m or 12m)
            if contains(session, '6m')
                session_label = '6m';
            elseif contains(session, '12m')
                session_label = '12m';
            else
                session_label = session;
            end

            for c = 1:length(conditions)
                condition = conditions{c};
                field_name = sprintf('ch%d', channel);

                if isfield(all_data, session) && isfield(all_data.(session), condition) && ...
                   isfield(all_data.(session).(condition), field_name)

                    data_struct = all_data.(session).(condition).(field_name);
                    
                    % Check if we have data for this stage
                    if isfield(data_struct, char(stage)) && ~isempty(data_struct.(char(stage)))
                        
                        % Get subject means for this stage
                        subject_means = data_struct.(char(stage));  % Each row is one subject's mean
                        n_subjects = size(subject_means, 1);
                        
                        % Calculate mean and SEM across subjects
                        grand_mean = mean(subject_means, 1);
                        grand_std = std(subject_means, 0, 1);
                        grand_sem = grand_std / sqrt(n_subjects);
                        
                        % Debug output
                        fprintf('    %s %s: %d subjects, SEM range: %.6f to %.6f\n', ...
                               session_label, condition, n_subjects, min(grand_sem), max(grand_sem));

                        % Get color based on session and condition
                        if contains(session, '6m')
                            color_rgb = colors_6m{c};
                        else
                            color_rgb = colors_12m{c};
                        end
                        
                        % Create shaded error region
                        scales = data_struct.scales;
                        upper_bound = grand_mean + grand_sem;
                        lower_bound = grand_mean - grand_sem;
                        
                        % Plot shaded area
                        fill([scales, fliplr(scales)], [upper_bound, fliplr(lower_bound)], ...
                             color_rgb, 'FaceAlpha', 0.33, 'EdgeColor', 'none', ...
                             'HandleVisibility', 'off');
                        hold on;
                        
                        % Plot mean line
                        h = plot(scales, grand_mean, '-', 'Color', color_rgb, ...
                                'LineWidth', 3, 'MarkerSize', 8);
                        
                        % Store legend info
                        legend_handles(end+1) = h;
                        legend_entries{end+1} = sprintf('%s %s (n=%d)', session_label, condition, n_subjects);
                    end
                end
            end
        end

        % Only finalize plot if we have data
        if ~isempty(legend_handles)
            xlabel('Temporal Scale', 'FontSize', 12);
            ylabel('Sample Entropy', 'FontSize', 12);
            title(sprintf('%s Stage - %s (All Sessions)', stage, channel_name), 'FontSize', 14);
            legend(legend_handles, legend_entries, 'Location', 'best', 'FontSize', 10);
            grid on;

            % Save plot
            plot_filename_png = fullfile(plot_dir, sprintf('%s_%s_combined.png', ...
                                       stage, channel_name));
            plot_filename_fig = fullfile(plot_dir, sprintf('%s_%s_combined.fig', ...
                                       stage, channel_name));
            
            % Save as PNG
            print(gcf, plot_filename_png, '-dpng', '-r300');
            fprintf('    Saved PNG: %s\n', plot_filename_png);
            
            % Save as FIG
            savefig(gcf, plot_filename_fig, 'compact');
            fprintf('    Saved FIG: %s\n', plot_filename_fig);
        else
            fprintf('    No data for %s %s - skipping plot\n', stage, channel_name);
        end

        close(gcf);
    end
end

fprintf('\nAll plots completed!\n');
end
