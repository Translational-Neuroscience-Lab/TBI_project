function plot_mse_comparison_average(base_path, conditions)
% Plot MSE curves comparing conditions across sleep stages for both channels
% Shows mean with shaded error bars (SEM)
% Conditions: TBI vs Sham
% Analysis without session or lighting condition separation

conditions = {'TBI', 'Sham'};
fprintf('Using shaded error bars with SEM...\n');

% Initialize storage
all_data = struct();
stages = ["Wake", "NREM", "REM"];
colors_conditions = {[0.63, 0, 0], [0, 0, 0.5]};  % Dark red for TBI, Dark blue for Sham
channels = [1, 2];

fprintf('Loading data for both channels...\n');

% Load all data for both channels
for ch = 1:length(channels)
    channel = channels(ch);
    fprintf('Processing Channel %d...\n', channel);

    for c = 1:length(conditions)
        condition = conditions{c};

        % Build path to Results folder
        condition_path = fullfile(base_path, condition, 'Results');

        if exist(condition_path, 'dir')
            % Find all dataset files in this condition
            dataset_files = dir(fullfile(condition_path, '*_dataset.mat'));

            if ~isempty(dataset_files)
                % Initialize arrays for this condition/channel
                condition_data = struct();
                condition_data.mse_data = [];
                condition_data.stages = [];
                condition_data.scales = [];
                condition_data.subject_ids = [];

                % Load and combine all subjects for this condition
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
                            
                            % Concatenate data for overall statistics
                            condition_data.mse_data = [condition_data.mse_data; mse_matrix];
                            condition_data.stages = [condition_data.stages; valid_stages];
                            condition_data.scales = dataset.mse_scales;
                            
                            % Add subject IDs for each epoch
                            subject_ids = repmat(f, n_valid, 1);
                            condition_data.subject_ids = [condition_data.subject_ids; subject_ids];
                        end

                    catch ME
                        fprintf('Warning: Could not load %s: %s\n', dataset_files(f).name, ME.message);
                    end
                end

                % Store combined data with channel info
                field_name = sprintf('ch%d', channel);
                all_data.(condition).(field_name) = condition_data;
                fprintf('  %s Channel %d: %d epochs from %d subjects loaded\n', condition, channel, ...
                       size(condition_data.mse_data, 1), length(dataset_files));
                
            else
                fprintf('  Warning: No dataset files found in %s\n', condition_path);
            end
        else
            fprintf('  Warning: Directory not found: %s\n', condition_path);
        end
    end
end

% Create plots - separate plots for each sleep stage
fprintf('\nCreating plots...\n');

for ch = 1:length(channels)
    channel = channels(ch);
    fprintf('Processing Channel %d plots...\n', channel);

    % Create mse_plots directory
    plot_dir = fullfile(base_path, 'mse_plots');
    if ~exist(plot_dir, 'dir')
        mkdir(plot_dir);
    end

    % Create separate plots for each sleep stage
    for st = 1:length(stages)
        stage = stages(st);
        
        fprintf('  Creating %s plot for Channel %d...\n', stage, channel);
        
        % Create individual plot for this stage
        figure('Position', [100 100 1000 700], 'Visible', 'off');

        legend_entries = {};
        legend_handles = [];

        for c = 1:length(conditions)
            condition = conditions{c};
            field_name = sprintf('ch%d', channel);

            if isfield(all_data, condition) && isfield(all_data.(condition), field_name)

                data_struct = all_data.(condition).(field_name);
                if ~isempty(data_struct.mse_data)
                    
                    % Filter by stage only
                    stage_idx = data_struct.stages == stage;

                    if sum(stage_idx) > 0
                        stage_mse = data_struct.mse_data(stage_idx, :);
                        mean_mse = mean(stage_mse, 1);
                        std_mse = std(stage_mse, 0, 1);
                        sem_mse = std_mse / sqrt(size(stage_mse, 1));

                        % Debug output
                        fprintf('    Data points: %d, SEM range: %.6f to %.6f\n', ...
                               size(stage_mse, 1), min(sem_mse), max(sem_mse));

                        % Get color
                        color_rgb = colors_conditions{c};
                        
                        % Create shaded error bar
                        x = data_struct.scales;
                        y = mean_mse;
                        err = abs(sem_mse);
                        
                        % Upper and lower bounds
                        upper = y + err;
                        lower = y - err;
                        
                        % Create patch for shaded area with 33% saturation
                        x_patch = [x, fliplr(x)];
                        y_patch = [upper, fliplr(lower)];
                        
                        patch_color = color_rgb + (1 - color_rgb) * (1 - 0.33);
                        h_patch = fill(x_patch, y_patch, patch_color, 'EdgeColor', 'none', ...
                                      'FaceAlpha', 0.33, 'HandleVisibility', 'off');
                        hold on;
                        
                        % Plot mean line on top
                        h = plot(x, y, 'o-', 'Color', color_rgb, 'LineWidth', 3, ...
                                'MarkerSize', 8, 'MarkerFaceColor', color_rgb);
                        
                        % Store legend info
                        legend_handles(end+1) = h;
                        legend_entries{end+1} = sprintf('%s (Mean ± SEM)', condition);
                    end
                end
            end
        end

        % Only create plot if we have data
        if ~isempty(legend_handles)
            xlabel('Temporal Scale', 'FontSize', 12);
            ylabel('Sample Entropy', 'FontSize', 12);
            title(sprintf('%s Stage - Channel %d', stage, channel), 'FontSize', 14);
            legend(legend_handles, legend_entries, 'Location', 'best', 'FontSize', 10);
            grid on;

            % Save plot
            plot_filename_png = fullfile(plot_dir, sprintf('%s_ch%d_average.png', stage, channel));
            plot_filename_fig = fullfile(plot_dir, sprintf('%s_ch%d_average.fig', stage, channel));
            
            % Save as PNG
            print(gcf, plot_filename_png, '-dpng', '-r300');
            fprintf('    Saved PNG: %s\n', plot_filename_png);
            
            % Save as FIG
            savefig(gcf, plot_filename_fig, 'compact');
            fprintf('    Saved FIG: %s\n', plot_filename_fig);
        else
            fprintf('    No data for %s Channel %d - skipping plot\n', stage, channel);
        end

        close(gcf);
    end
end

fprintf('\nAll plots completed!\n');
end
