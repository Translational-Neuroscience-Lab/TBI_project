function mse_pls_lights_sex(base_path_Adult)
% 6m subacute ("adult") cohort PLS analysis for MSE data
% Tests whether Sham vs TBI differences in the Adult cohort
% Runs only Channel 2 and plots both LV1 and LV2
% Also runs sex-split PLS with 4 groups:
%   Sham Male, Sham Female, TBI Male, TBI Female

% Add PLS toolbox to path
addpath(genpath('PLS'));

% Parameters
sessions = {'Session1'};
conditions = {'TBI', 'Sham'};
stages = ["Wake", "NREM", "REM"];
lighting_conditions = ["Lights_On", "Lights_Off"];
channels = [2];

% Define cohorts
cohorts = {'Adult'};
cohort_sessions = struct();
cohort_sessions.Adult = sessions;
cohort_paths = struct();
cohort_paths.Adult = base_path_Adult;

% --- Sex lookup map (hardcoded from MTN_Sex.txt) ---
sex_map = containers.Map(...
    {'Mouse#5','Mouse#18','Mouse#19','Mouse#20', ...                         % Sham Male
     'Mouse#3','Mouse#4','Mouse#12','Mouse#13','Mouse#16','Mouse#17', ...    % TBI Male
     'Mouse#1','Mouse#2','Mouse#10','Mouse#11', ...                          % Sham Female
     'Mouse#6','Mouse#7','Mouse#14','Mouse#15'}, ...                         % TBI Female
    {'Male','Male','Male','Male', ...
     'Male','Male','Male','Male','Male','Male', ...
     'Female','Female','Female','Female', ...
     'Female','Female','Female','Female'});

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

%% Load MSE data and combine into cohorts
fprintf('Loading MSE data and combining into cohorts...\n');
all_data = struct();
all_data_sex = struct();

for ch = 1:length(channels)
    channel = channels(ch);
    fprintf('Processing Channel %d...\n', channel);

    for co = 1:length(cohorts)
        cohort = cohorts{co};
        session_list = cohort_sessions.(cohort);
        base_path = cohort_paths.(cohort);

        fprintf('  Processing %s cohort (combining sessions %s)...\n', cohort, strjoin(session_list, ', '));

        % ----------------------------------------------------------
        % ORIGINAL: load per-condition (Sham, TBI)
        % ----------------------------------------------------------
        for c = 1:length(conditions)
            condition = conditions{c};

            % Initialize combined data structure
            combined_mouse_data = struct();
            combined_mouse_data.mouse_averages = [];
            combined_mouse_data.mouse_ids = {};
            combined_mouse_data.stages = [];
            combined_mouse_data.lighting = [];
            combined_mouse_data.scales = [];
            combined_mouse_data.n_epochs_per_mouse = [];

            % Combine data across all sessions in this cohort
            for s = 1:length(session_list)
                session = session_list{s};
                condition_path = fullfile(base_path, session, condition, 'Results');

                if exist(condition_path, 'dir')
                    dataset_files = dir(fullfile(condition_path, '*_dataset.mat'));

                    if ~isempty(dataset_files)
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

                            valid_stages   = dataset.Sleep_Stage(valid_idx);
                            valid_lighting = dataset.Lighting_Period(valid_idx);

                            % Convert MSE cell array to matrix
                            if ~isempty(valid_mse)
                                n_valid    = length(valid_mse);
                                mse_matrix = zeros(n_valid, length(dataset.mse_scales));
                                for i = 1:n_valid
                                    if ~isempty(valid_mse{i})
                                        mse_matrix(i, :) = valid_mse{i};
                                    end
                                end

                                % Average across epochs for each stage/lighting combination
                                unique_stages   = unique(valid_stages);
                                unique_lighting = unique(valid_lighting);

                                for stage_idx = 1:length(unique_stages)
                                    stage = unique_stages(stage_idx);
                                    for light_idx = 1:length(unique_lighting)
                                        lighting = unique_lighting(light_idx);

                                        combo_idx = (valid_stages == stage) & (valid_lighting == lighting);

                                        if sum(combo_idx) >= 3
                                            mouse_avg = mean(mse_matrix(combo_idx, :), 1);
                                            n_epochs  = sum(combo_idx);

                                            combined_mouse_data.mouse_averages = [combined_mouse_data.mouse_averages; mouse_avg];
                                            combined_mouse_data.mouse_ids{end+1} = sprintf('%s_%s_%s_%s', mouse_id, session, stage, lighting);
                                            combined_mouse_data.stages   = [combined_mouse_data.stages; stage];
                                            combined_mouse_data.lighting = [combined_mouse_data.lighting; lighting];
                                            combined_mouse_data.n_epochs_per_mouse = [combined_mouse_data.n_epochs_per_mouse; n_epochs];
                                            combined_mouse_data.scales   = dataset.mse_scales;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end

            % Store combined cohort data
            field_name = sprintf('ch%d', channel);
            all_data.(cohort).(condition).(field_name) = combined_mouse_data;
            fprintf('    %s %s Channel %d: %d mouse-level averages\n', ...
                cohort, condition, channel, size(combined_mouse_data.mouse_averages, 1));
        end

        % ----------------------------------------------------------
        % SEX-SPLIT: load per sex × condition group
        % ----------------------------------------------------------
        field_name = sprintf('ch%d', channel);
        sex_groups = {'Sham_Male', 'Sham_Female', 'TBI_Male', 'TBI_Female'};
        sex_group_names = {'Sham Male', 'Sham Female', 'TBI Male', 'TBI Female'};

        % Initialise all 4 sex-split group structs
        for g = 1:length(sex_groups)
            grp = sex_groups{g};
            all_data_sex.(cohort).(grp).(field_name).mouse_averages      = [];
            all_data_sex.(cohort).(grp).(field_name).mouse_ids           = {};
            all_data_sex.(cohort).(grp).(field_name).stages              = [];
            all_data_sex.(cohort).(grp).(field_name).lighting            = [];
            all_data_sex.(cohort).(grp).(field_name).scales              = [];
            all_data_sex.(cohort).(grp).(field_name).n_epochs_per_mouse  = [];
        end

        for c = 1:length(conditions)
            condition = conditions{c};

            for s = 1:length(session_list)
                session = session_list{s};
                condition_path = fullfile(base_path, session, condition, 'Results');

                if exist(condition_path, 'dir')
                    dataset_files = dir(fullfile(condition_path, '*_dataset.mat'));

                    if ~isempty(dataset_files)
                        for f = 1:length(dataset_files)
                            file_path = fullfile(condition_path, dataset_files(f).name);
                            [~, filename, ~] = fileparts(dataset_files(f).name);
                            filename_parts = strsplit(filename, '_');
                            mouse_id = filename_parts{1};

                            % Determine sex from lookup map
                            if isKey(sex_map, mouse_id)
                                sex = sex_map(mouse_id);
                            else
                                fprintf('Warning: Mouse ID %s not found in sex_map, skipping.\n', mouse_id);
                                continue;
                            end

                            % Build group key e.g. 'Sham_Male'
                            group_key = sprintf('%s_%s', condition, sex);

                            data = load(file_path);
                            dataset = data.dataset;

                            if channel == 1
                                valid_idx = ~cellfun(@isempty, dataset.MSE_EEG1);
                                valid_mse = dataset.MSE_EEG1(valid_idx);
                            else
                                valid_idx = ~cellfun(@isempty, dataset.MSE_EEG2);
                                valid_mse = dataset.MSE_EEG2(valid_idx);
                            end

                            valid_stages   = dataset.Sleep_Stage(valid_idx);
                            valid_lighting = dataset.Lighting_Period(valid_idx);

                            if ~isempty(valid_mse)
                                n_valid    = length(valid_mse);
                                mse_matrix = zeros(n_valid, length(dataset.mse_scales));
                                for i = 1:n_valid
                                    if ~isempty(valid_mse{i})
                                        mse_matrix(i, :) = valid_mse{i};
                                    end
                                end

                                unique_stages   = unique(valid_stages);
                                unique_lighting = unique(valid_lighting);

                                for stage_idx = 1:length(unique_stages)
                                    stage = unique_stages(stage_idx);
                                    for light_idx = 1:length(unique_lighting)
                                        lighting = unique_lighting(light_idx);

                                        combo_idx = (valid_stages == stage) & (valid_lighting == lighting);

                                        if sum(combo_idx) >= 3
                                            mouse_avg = mean(mse_matrix(combo_idx, :), 1);
                                            n_epochs  = sum(combo_idx);

                                            all_data_sex.(cohort).(group_key).(field_name).mouse_averages = ...
                                                [all_data_sex.(cohort).(group_key).(field_name).mouse_averages; mouse_avg];
                                            all_data_sex.(cohort).(group_key).(field_name).mouse_ids{end+1} = ...
                                                sprintf('%s_%s_%s_%s', mouse_id, session, stage, lighting);
                                            all_data_sex.(cohort).(group_key).(field_name).stages = ...
                                                [all_data_sex.(cohort).(group_key).(field_name).stages; stage];
                                            all_data_sex.(cohort).(group_key).(field_name).lighting = ...
                                                [all_data_sex.(cohort).(group_key).(field_name).lighting; lighting];
                                            all_data_sex.(cohort).(group_key).(field_name).n_epochs_per_mouse = ...
                                                [all_data_sex.(cohort).(group_key).(field_name).n_epochs_per_mouse; n_epochs];
                                            all_data_sex.(cohort).(group_key).(field_name).scales = dataset.mse_scales;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end

        % Report sex-split counts
        for g = 1:length(sex_groups)
            grp = sex_groups{g};
            fprintf('    %s %s Channel %d: %d mouse-level averages\n', ...
                cohort, grp, channel, size(all_data_sex.(cohort).(grp).(field_name).mouse_averages, 1));
        end
    end
end

%% Create results directories
results_dir = fullfile(base_path_Adult, 'PLS_Adult_Cohort');
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

results_dir_sex = fullfile(base_path_Adult, 'PLS_Adult_Sex');
if ~exist(results_dir_sex, 'dir')
    mkdir(results_dir_sex);
end

%% =========================================================
%% ORIGINAL 2-GROUP PLS (Sham vs TBI) — UNCHANGED
%% =========================================================
fprintf('\nRunning PLS Analysis for Adult cohort...\n');

for ch = 1:length(channels)
    channel = channels(ch);
    field_name = sprintf('ch%d', channel);

    for st = 1:length(stages)
        stage = stages(st);

        for lc = 1:length(lighting_conditions)
            lighting = lighting_conditions(lc);

            fprintf('Processing %s - %s - Channel %d...\n', stage, lighting, channel);

            % Check data completeness
            combination_complete = true;
            for co_check = 1:length(cohorts)
                cohort_check = cohorts{co_check};
                for c_check = 1:length(conditions)
                    condition_check = conditions{c_check};

                    actual_count = 0;
                    if isfield(all_data, cohort_check) && isfield(all_data.(cohort_check), condition_check) && ...
                            isfield(all_data.(cohort_check).(condition_check), field_name)

                        data_struct_check = all_data.(cohort_check).(condition_check).(field_name);
                        if ~isempty(data_struct_check.mouse_averages)
                            stage_idx_check    = data_struct_check.stages   == stage;
                            lighting_idx_check = data_struct_check.lighting == lighting;
                            combined_idx_check = stage_idx_check & lighting_idx_check;
                            actual_count = sum(combined_idx_check);
                        end
                    end

                    if actual_count == 0
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
            group_labels     = [];
            cohort_labels    = [];

            for co = 1:length(cohorts)
                cohort = cohorts{co};

                for c = 1:length(conditions)
                    condition = conditions{c};

                    if isfield(all_data, cohort) && isfield(all_data.(cohort), condition) && ...
                            isfield(all_data.(cohort).(condition), field_name)

                        data_struct = all_data.(cohort).(condition).(field_name);

                        if ~isempty(data_struct.mouse_averages)
                            stage_idx    = data_struct.stages   == stage;
                            lighting_idx = data_struct.lighting == lighting;
                            combined_idx = stage_idx & lighting_idx;

                            if sum(combined_idx) > 0
                                condition_data = data_struct.mouse_averages(combined_idx, :);

                                interaction_data = [interaction_data; condition_data];
                                group_code  = strcmp(condition, 'TBI');
                                group_labels = [group_labels; repmat(group_code, size(condition_data, 1), 1)];
                                cohort_code  = strcmp(cohort, 'Younger');
                                cohort_labels = [cohort_labels; repmat(cohort_code, size(condition_data, 1), 1)];
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
                % Groups: Sham, TBI
                pls_data    = [];
                group_sizes = [];
                group_labels_pls = {};

                group_combinations = [0 0; 0 1]; % [cohort(0=Adult), group(0=Sham,1=TBI)]
                group_names = {'Sham', 'TBI'};

                for g = 1:size(group_combinations, 1)
                    cohort_code = group_combinations(g, 1);
                    group_code  = group_combinations(g, 2);

                    combo_idx = (cohort_labels == cohort_code) & (group_labels == group_code);

                    if sum(combo_idx) >= 3
                        combo_data  = interaction_data(combo_idx, :);
                        pls_data    = [pls_data; combo_data];
                        group_sizes(end+1) = size(combo_data, 1);
                        group_labels_pls{end+1} = group_names{g};
                    end
                end

                % PLS options
                if length(group_sizes) >= 2 && min(group_sizes) >= 3
                    clear option
                    option.method   = 1;    % Mean-centred PLS
                    option.num_perm = 1000;
                    option.num_boot = 500;

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

                        % Plot both LV1 and LV2
                        max_LV = min(2, length(pvals));

                        for LV = 1:max_LV
                            p = pvals(LV);
                            headline = sprintf('Adult Cohort: %s %s Ch%d, LV%d p=%.4f', ...
                                stage, lighting, channel, LV, p);

                            figure('Position', [100 100 1000 500]);

                            % Bar plot
                            subplot(1,2,1)
                            z = pls_result.boot_result.orig_usc;
                            bar_handle = bar(z(:,LV));
                            hold on

                            % Color bars: dark blue (Sham), dark red (TBI)
                            bar_colors = zeros(length(group_labels_pls), 3);
                            for i = 1:length(group_labels_pls)
                                if contains(group_labels_pls{i}, 'Sham')
                                    bar_colors(i, :) = [0.1, 0.3, 0.7]; % Dark blue
                                else
                                    bar_colors(i, :) = [0.7, 0.1, 0.1]; % Dark red
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
                            title('Group Contrast Scores')
                            ylabel('Contrast Scores')
                            xlabel('Condition')

                            % Heatmap
                            subplot(1,2,2)
                            x = pls_result.boot_result.compare_u(:,LV);
                            x(abs(x)<2.3) = 0;
                            scales   = combined_mouse_data.scales;
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
                            fig_filename = fullfile(results_dir, sprintf('Adult_PLS_%s_%s_ch%d_LV%d.png', ...
                                stage, lighting, channel, LV));
                            saveas(gcf, fig_filename)

                            fig_filename_fig = fullfile(results_dir, sprintf('Adult_PLS_%s_%s_ch%d_LV%d.fig', ...
                                stage, lighting, channel, LV));
                            saveas(gcf, fig_filename_fig)

                            fprintf('  Saved: %s\n', fig_filename);

                            close(gcf);
                        end
                    end
                end
            end
        end
    end
end

%% =========================================================
%% SEX-SPLIT 4-GROUP PLS
%% (Sham Male, Sham Female, TBI Male, TBI Female)
%% =========================================================
fprintf('\nRunning Sex-Split PLS Analysis for Adult cohort...\n');

sex_groups      = {'Sham_Male', 'Sham_Female', 'TBI_Male', 'TBI_Female'};
sex_group_names = {'Sham Male', 'Sham Female', 'TBI Male', 'TBI Female'};

% Colors matching the MSE sex-split plots
sex_bar_colors = [0.1, 0.3, 0.7;   % Sham Male    - dark blue
                  0.4, 0.7, 1.0;   % Sham Female  - light blue
                  0.7, 0.1, 0.1;   % TBI Male  - dark red
                  1.0, 0.5, 0.5];  % TBI Female - light pink

for ch = 1:length(channels)
    channel = channels(ch);
    field_name = sprintf('ch%d', channel);

    for st = 1:length(stages)
        stage = stages(st);

        for lc = 1:length(lighting_conditions)
            lighting = lighting_conditions(lc);

            fprintf('Sex-Split PLS: Processing %s - %s - Channel %d...\n', stage, lighting, channel);

            % Check that every group has at least some data for this
            % stage × lighting combination before proceeding
            combination_complete = true;
            for g = 1:length(sex_groups)
                grp = sex_groups{g};
                actual_count = 0;

                if isfield(all_data_sex, 'Adult') && isfield(all_data_sex.Adult, grp) && ...
                        isfield(all_data_sex.Adult.(grp), field_name)

                    ds = all_data_sex.Adult.(grp).(field_name);
                    if ~isempty(ds.mouse_averages)
                        stage_idx_check    = ds.stages   == stage;
                        lighting_idx_check = ds.lighting == lighting;
                        actual_count = sum(stage_idx_check & lighting_idx_check);
                    end
                end

                if actual_count == 0
                    combination_complete = false;
                    break;
                end
            end

            if ~combination_complete
                fprintf('  Skipping sex-split - incomplete data\n');
                continue;
            end

            % Collect data for the 4-group sex-split PLS
            pls_data_sex    = [];
            group_sizes_sex = [];
            group_labels_sex = {};
            scales_sex = [];

            for g = 1:length(sex_groups)
                grp = sex_groups{g};
                ds  = all_data_sex.Adult.(grp).(field_name);

                stage_idx    = ds.stages   == stage;
                lighting_idx = ds.lighting == lighting;
                combined_idx = stage_idx & lighting_idx;

                if sum(combined_idx) >= 3
                    grp_data = ds.mouse_averages(combined_idx, :);
                    % Clean data
                    grp_data(isnan(grp_data) | isinf(grp_data)) = 0;

                    pls_data_sex    = [pls_data_sex; grp_data];
                    group_sizes_sex(end+1) = size(grp_data, 1);
                    group_labels_sex{end+1} = sex_group_names{g};

                    if isempty(scales_sex)
                        scales_sex = ds.scales;
                    end

                    fprintf('    %s: %d subjects\n', sex_group_names{g}, size(grp_data, 1));
                else
                    fprintf('    %s: only %d subjects (< 3), skipping this combination\n', ...
                        sex_group_names{g}, sum(combined_idx));
                    combination_complete = false;
                    break;
                end
            end

            if ~combination_complete || length(group_sizes_sex) < 2 || min(group_sizes_sex) < 3
                fprintf('  Skipping sex-split PLS - insufficient group sizes\n');
                continue;
            end

            % Run sex-split PLS
            clear option
            option.method   = 1;    % Mean-centred PLS
            option.num_perm = 1000;
            option.num_boot = 500;

            ncond_sex = length(group_sizes_sex);
            pls_result_sex = pls_analysis({pls_data_sex}, {group_sizes_sex}, ncond_sex, option);

            if isfield(pls_result_sex, 'perm_result') && isfield(pls_result_sex.perm_result, 'sprob')
                pvals_sex = pls_result_sex.perm_result.sprob;

                fprintf('  Sex-Split P-values: ');
                for pv = 1:length(pvals_sex)
                    if pvals_sex(pv) < 0.001
                        fprintf('LV%d<0.001 ', pv);
                    elseif pvals_sex(pv) < 0.05
                        fprintf('LV%d=%.4f* ', pv, pvals_sex(pv));
                    else
                        fprintf('LV%d=%.4f ', pv, pvals_sex(pv));
                    end
                end
                fprintf('\n');

                % Plot both LV1 and LV2
                max_LV = min(2, length(pvals_sex));

                for LV = 1:max_LV
                    p = pvals_sex(LV);
                    headline = sprintf('Adult Cohort (Sex Split): %s %s Ch%d, LV%d p=%.4f', ...
                        stage, lighting, channel, LV, p);

                    figure('Position', [100 100 1000 500]);

                    % Bar plot
                    subplot(1,2,1)
                    z = pls_result_sex.boot_result.orig_usc;
                    bar_handle = bar(z(:,LV));
                    hold on

                    % Color bars using sex_bar_colors
                    bar_colors_plot = zeros(length(group_labels_sex), 3);
                    for i = 1:length(group_labels_sex)
                        grp_idx = find(strcmp(sex_group_names, group_labels_sex{i}), 1);
                        if ~isempty(grp_idx)
                            bar_colors_plot(i, :) = sex_bar_colors(grp_idx, :);
                        end
                    end
                    bar_handle.CData = bar_colors_plot;

                    % Error bars
                    yneg = pls_result_sex.boot_result.llusc(:,LV);
                    ypos = pls_result_sex.boot_result.ulusc(:,LV);
                    errorbar(1:length(z), z(:,LV), yneg - z(:,LV), ypos - z(:,LV), '.', 'Color', 'black')

                    xticks(1:length(group_labels_sex))
                    xticklabels(group_labels_sex)
                    xtickangle(45)
                    grid on
                    title('Group Contrast Scores')
                    ylabel('Contrast Scores')
                    xlabel('Condition')

                    % Heatmap
                    subplot(1,2,2)
                    x = pls_result_sex.boot_result.compare_u(:,LV);
                    x(abs(x)<2.3) = 0;
                    plotdata = reshape(x, length(scales_sex), []);

                    imagesc(plotdata)
                    colormap(rgb)
                    clim([-7 7])
                    colorbar
                    yticks(1:5:length(scales_sex))
                    yticklabels(1:5:length(scales_sex))
                    xticks([])
                    title(sprintf('Salience LV%d', LV))
                    ylabel('Temporal Scale')

                    sgtitle(headline, 'FontSize', 16, 'Interpreter', 'none')
                    set(gcf, 'Color', 'w')

                    % Save figure
                    fig_filename_png = fullfile(results_dir_sex, sprintf('Adult_Sex_PLS_%s_%s_ch%d_LV%d.png', ...
                        stage, lighting, channel, LV));
                    fig_filename_fig = fullfile(results_dir_sex, sprintf('Adult_Sex_PLS_%s_%s_ch%d_LV%d.fig', ...
                        stage, lighting, channel, LV));

                    saveas(gcf, fig_filename_png)
                    saveas(gcf, fig_filename_fig)

                    fprintf('  Saved: %s\n', fig_filename_png);

                    close(gcf);
                end
            end
        end
    end
end

fprintf('\nAdult cohort PLS analysis complete!\n');
end
