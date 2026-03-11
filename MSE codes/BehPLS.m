function mse_pls_behavioural(base_path_older, neuropath_path)
% Behavioural PLS analysis for MSE data correlated with neuropathology
% All mice treated as one group (ncond=1)
% Behaviour matrix: 18 variables
%   - Iba1 IHC (3): TCTX, DHPC, ENT
%   - Western Blot (3): DA9/Actin, 6E10/Actin, PHF1/Actin
%   - Abeta ELISA (4): Soluble Abeta42, Insoluble Abeta42, Soluble Abeta40, Insoluble Abeta40
%   - Abeta IHC (3): TCTX, DHPC, ENT
%   - ThioS (3): TCTX, DHPC, ENT
%   - Treatment dummy (1): Control=0, Trazodone=1
%   - Sex dummy (1): Male=0, Female=1
%
% MSE data is averaged across all sessions per mouse per stage x lighting
% combo before entering PLS (one row per mouse).
% Verification .mat files are saved before each PLS run.

% Add PLS toolbox to path
addpath(genpath('PLS'));

% Parameters
older_sessions = {'Session1', 'Session2', 'Session4', 'Session5', 'Session6', 'Session7'};
conditions     = {'Trazodone', 'Control'};
stages         = ["Wake", "NREM", "REM"];
lighting_conditions = ["Lights_On", "Lights_Off"];
channels       = [2];

% Behaviour variable names (18 total) — used for axis labels
behav_column_names = {...
    'Iba1 TCTX',    'Iba1 DHPC',    'Iba1 ENT', ...
    'DA9/Actin',    '6E10/Actin',   'PHF1/Actin', ...
    'Sol Abeta42',  'Insol Abeta42','Sol Abeta40','Insol Abeta40', ...
    'Abeta IHC TCTX','Abeta IHC DHPC','Abeta IHC ENT', ...
    'ThioS TCTX',   'ThioS DHPC',   'ThioS ENT', ...
    'Treatment',    'Sex'};

%% -----------------------------------------------------------------------
%% Sex and treatment lookup maps
%% -----------------------------------------------------------------------
% Sex map (Male=0, Female=1)
sex_map = containers.Map(...
    {'MTN05','MTN06','MTN08', ...
     'MTN04','MTN07','MTN10','MTN11','MTN12','MTN13', ...
     'MTN02','MTN16','MTN19', ...
     'MTN01','MTN03','MTN17','MTN18','MTN20'}, ...
    {0,0,0, ...
     0,0,0,0,0,0, ...
     1,1,1, ...
     1,1,1,1,1});

% Treatment map (Control=0, Trazodone=1)
% Derived from which condition folder the mouse lives in.
% Control mice:   MTN05, MTN06, MTN08, MTN09, MTN02, MTN16, MTN19
% Trazodone mice: MTN04, MTN07, MTN10, MTN11, MTN12, MTN13,
%                 MTN01, MTN03, MTN17, MTN18, MTN20
treatment_map = containers.Map(...
    {'MTN05','MTN06','MTN08','MTN09', ...
     'MTN02','MTN16','MTN19', ...
     'MTN04','MTN07','MTN10','MTN11','MTN12','MTN13', ...
     'MTN01','MTN03','MTN17','MTN18','MTN20'}, ...
    {0,0,0,0, ...
     0,0,0, ...
     1,1,1,1,1,1, ...
     1,1,1,1,1});

%% -----------------------------------------------------------------------
%% Load neuropathology data into lookup maps
%% -----------------------------------------------------------------------
fprintf('Loading neuropathology data...\n');


% --- Iba1 IHC ---
iba1_data = readtable(fullfile(neuropath_path, 'MTN_Iba1_IHC.txt'), 'Delimiter', '\t');
iba1_map  = containers.Map();
for r = 1:height(iba1_data)
    mid = strtrim(iba1_data{r,1}{1});
    iba1_map(mid) = [iba1_data{r,2}, iba1_data{r,3}, iba1_data{r,4}];
end

% --- Western Blot ---
wb_data  = readtable(fullfile(neuropath_path, 'MTN_Western_Blot.txt'), 'Delimiter', '\t');
wb_map   = containers.Map();
for r = 1:height(wb_data)
    mid = strtrim(wb_data{r,1}{1});
    wb_map(mid) = [wb_data{r,2}, wb_data{r,3}, wb_data{r,4}];
end

% --- Abeta ELISA ---
elisa_data = readtable(fullfile(neuropath_path, 'MTN_Abeta_ELISA.txt'), 'Delimiter', '\t');
elisa_map  = containers.Map();
for r = 1:height(elisa_data)
    mid = strtrim(elisa_data{r,1}{1});
    elisa_map(mid) = [elisa_data{r,2}, elisa_data{r,3}, elisa_data{r,4}, elisa_data{r,5}];
end

% --- Abeta IHC ---
abeta_ihc_data = readtable(fullfile(neuropath_path, 'MTN_Abeta_IHC.txt'), 'Delimiter', '\t');
abeta_ihc_map  = containers.Map();
for r = 1:height(abeta_ihc_data)
    mid = strtrim(abeta_ihc_data{r,1}{1});
    abeta_ihc_map(mid) = [abeta_ihc_data{r,2}, abeta_ihc_data{r,3}, abeta_ihc_data{r,4}];
end

% --- ThioS ---
thios_data = readtable(fullfile(neuropath_path, 'MTN_ThioS.txt'), 'Delimiter', '\t');
thios_map  = containers.Map();
for r = 1:height(thios_data)
    mid = strtrim(thios_data{r,1}{1});
    thios_map(mid) = [thios_data{r,2}, thios_data{r,3}, thios_data{r,4}];
end

fprintf('Neuropathology data loaded.\n');

%% -----------------------------------------------------------------------
%% RGB colormap (unchanged)
%% -----------------------------------------------------------------------
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

%% -----------------------------------------------------------------------
%% Load MSE data — average across sessions per mouse
%% -----------------------------------------------------------------------
% For each mouse x stage x lighting combination, collect all valid epoch
% averages across ALL sessions, then compute one mean (Option 1).
%
% Data structure per channel:
%   all_data_avg.(field_name).mouse_averages  — one row per mouse per combo
%   all_data_avg.(field_name).mouse_ids_clean — MTN ID for each row
%   all_data_avg.(field_name).stages          — stage label for each row
%   all_data_avg.(field_name).lighting        — lighting label for each row
%   all_data_avg.(field_name).scales          — MSE scales vector

fprintf('Loading MSE data and averaging across sessions per mouse...\n');

all_data_avg = struct();

for ch = 1:length(channels)
    channel    = channels(ch);
    field_name = sprintf('ch%d', channel);
    fprintf('  Processing Channel %d...\n', channel);

    % Accumulator: keyed on "MTNXX_stage_lighting"
    % stores a matrix of epoch-averaged MSE rows (one per session)
    accum = struct();

    for c = 1:length(conditions)
        condition = conditions{c};

        for s = 1:length(older_sessions)
            session        = older_sessions{s};
            condition_path = fullfile(base_path_older, session, condition, 'Results');

            if ~exist(condition_path, 'dir')
                continue;
            end

            dataset_files = dir(fullfile(condition_path, '*_dataset.mat'));
            if isempty(dataset_files)
                continue;
            end

            for f = 1:length(dataset_files)
                file_path      = fullfile(condition_path, dataset_files(f).name);
                [~, filename, ~] = fileparts(dataset_files(f).name);
                filename_parts = strsplit(filename, '_');
                mouse_id       = filename_parts{1};  % e.g. 'MTN05'

                % Skip mice not in sex_map or treatment_map
                if ~isKey(sex_map, mouse_id) || ~isKey(treatment_map, mouse_id)
                    fprintf('    Warning: %s not in sex_map or treatment_map, skipping.\n', mouse_id);
                    continue;
                end

                % Skip mice missing any neuropath data
                if ~isKey(iba1_map, mouse_id) || ~isKey(wb_map, mouse_id) || ...
                   ~isKey(elisa_map, mouse_id) || ~isKey(abeta_ihc_map, mouse_id) || ...
                   ~isKey(thios_map, mouse_id)
                    fprintf('    Warning: %s missing neuropath data, skipping.\n', mouse_id);
                    continue;
                end

                data    = load(file_path);
                dataset = data.dataset;

                % Get valid epochs for the chosen channel
                if channel == 1
                    valid_idx = ~cellfun(@isempty, dataset.MSE_EEG1);
                    valid_mse = dataset.MSE_EEG1(valid_idx);
                else
                    valid_idx = ~cellfun(@isempty, dataset.MSE_EEG2);
                    valid_mse = dataset.MSE_EEG2(valid_idx);
                end

                valid_stages   = dataset.Sleep_Stage(valid_idx);
                valid_lighting = dataset.Lighting_Period(valid_idx);

                if isempty(valid_mse)
                    continue;
                end

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
                        lighting   = unique_lighting(light_idx);
                        combo_idx  = (valid_stages == stage) & (valid_lighting == lighting);

                        if sum(combo_idx) >= 3
                            session_avg = mean(mse_matrix(combo_idx, :), 1);

                            % Build accumulator key
                            accum_key = sprintf('%s__%s__%s', mouse_id, stage, lighting);
                            accum_key = strrep(accum_key, ' ', '_');

                            if isfield(accum, accum_key)
                                accum.(accum_key).rows = [accum.(accum_key).rows; session_avg];
                            else
                                accum.(accum_key).rows     = session_avg;
                                accum.(accum_key).mouse_id = mouse_id;
                                accum.(accum_key).stage    = stage;
                                accum.(accum_key).lighting = lighting;
                                accum.(accum_key).scales   = dataset.mse_scales;
                            end
                        end
                    end
                end
            end
        end
    end

    % Now collapse accum: average all session rows per key -> one row per mouse
    mouse_averages  = [];
    mouse_ids_clean = {};
    stage_labels    = [];
    lighting_labels = [];
    scales_out      = [];

    accum_keys = fieldnames(accum);
    for k = 1:length(accum_keys)
        key   = accum_keys{k};
        entry = accum.(key);
        mouse_mean = mean(entry.rows, 1);  % average across sessions

        mouse_averages  = [mouse_averages; mouse_mean];
        mouse_ids_clean{end+1} = entry.mouse_id;
        stage_labels    = [stage_labels;   entry.stage];
        lighting_labels = [lighting_labels; entry.lighting];
        if isempty(scales_out)
            scales_out = entry.scales;
        end
    end

    all_data_avg.(field_name).mouse_averages  = mouse_averages;
    all_data_avg.(field_name).mouse_ids_clean = mouse_ids_clean;
    all_data_avg.(field_name).stages          = stage_labels;
    all_data_avg.(field_name).lighting        = lighting_labels;
    all_data_avg.(field_name).scales          = scales_out;

    fprintf('  Channel %d: %d total mouse x stage x lighting rows\n', ...
        channel, size(mouse_averages, 1));
end

%% -----------------------------------------------------------------------
%% Create output directories
%% -----------------------------------------------------------------------
results_dir = fullfile(base_path_older, 'PLS_Behavioural');
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

verify_dir = fullfile(results_dir, 'Verification');
if ~exist(verify_dir, 'dir')
    mkdir(verify_dir);
end

%% -----------------------------------------------------------------------
%% Run Behavioural PLS per channel x stage x lighting
%% -----------------------------------------------------------------------
fprintf('\nRunning Behavioural PLS...\n');

for ch = 1:length(channels)
    channel    = channels(ch);
    field_name = sprintf('ch%d', channel);

    data_struct = all_data_avg.(field_name);

    for st = 1:length(stages)
        stage = stages(st);

        for lc = 1:length(lighting_conditions)
            lighting = lighting_conditions(lc);

            fprintf('Processing %s - %s - Channel %d...\n', stage, lighting, channel);

            % Find rows matching this stage x lighting combo
            stage_idx    = data_struct.stages   == stage;
            lighting_idx = data_struct.lighting == lighting;
            combo_idx    = stage_idx & lighting_idx;

            if sum(combo_idx) < 3
                fprintf('  Skipping - fewer than 3 mice for this combination\n');
                continue;
            end

            % Get MSE rows for this combo
            mse_rows   = data_struct.mouse_averages(combo_idx, :);
            mouse_ids  = data_struct.mouse_ids_clean(combo_idx);
            n_mice     = size(mse_rows, 1);
            scales     = data_struct.scales;

            % Clean MSE data
            mse_rows(isnan(mse_rows) | isinf(mse_rows)) = 0;

            % ----------------------------------------------------------
            % Build behaviour matrix row by row in the SAME order
            % as mse_rows — guarantees exact correspondence
            % ----------------------------------------------------------
            stacked_behavdata = zeros(n_mice, 18);
            row_labels        = cell(n_mice, 1);
            mouse_id_per_row  = cell(n_mice, 1);

            skip_flag = false;
            for m = 1:n_mice
                mid = mouse_ids{m};

                % Final check: all maps must have this mouse
                if ~isKey(iba1_map,      mid) || ~isKey(wb_map,       mid) || ...
                   ~isKey(elisa_map,     mid) || ~isKey(abeta_ihc_map, mid) || ...
                   ~isKey(thios_map,     mid) || ~isKey(sex_map,       mid) || ...
                   ~isKey(treatment_map, mid)
                    fprintf('  Warning: %s missing data, skipping entire combination\n', mid);
                    skip_flag = true;
                    break;
                end

                iba1_vals     = iba1_map(mid);       % 1x3
                wb_vals       = wb_map(mid);          % 1x3
                elisa_vals    = elisa_map(mid);       % 1x4
                abeta_vals    = abeta_ihc_map(mid);   % 1x3
                thios_vals    = thios_map(mid);       % 1x3
                treat_val     = treatment_map(mid);   % scalar
                sex_val       = sex_map(mid);         % scalar

                stacked_behavdata(m, :) = [...
                    iba1_vals, wb_vals, elisa_vals, ...
                    abeta_vals, thios_vals, ...
                    treat_val, sex_val];

                row_labels{m}       = sprintf('%s | %s | %s', mid, stage, lighting);
                mouse_id_per_row{m} = mid;
            end

            if skip_flag
                continue;
            end

            % ----------------------------------------------------------
            % Save verification file for this combination
            % ----------------------------------------------------------
            pls_data            = mse_rows;   % rename for clarity in .mat
            verify_filename     = fullfile(verify_dir, ...
                sprintf('verify_%s_%s_ch%d.mat', stage, lighting, channel));

            save(verify_filename, ...
                'pls_data', ...
                'stacked_behavdata', ...
                'row_labels', ...
                'behav_column_names', ...
                'mouse_id_per_row');

            fprintf('  Saved verification file: %s\n', verify_filename);

            % ----------------------------------------------------------
            % PLS options — Behavioural PLS
            % ----------------------------------------------------------
            clear option
            option.method            = 3;       % Behavioural PLS
            option.num_perm          = 1000;
            option.num_boot          = 500;
            option.stacked_behavdata = stacked_behavdata;

            ncond  = 1;
            nparts = n_mice;

            pls_result = pls_analysis({pls_data}, {nparts}, ncond, option);

            if ~isfield(pls_result, 'perm_result') || ~isfield(pls_result.perm_result, 'sprob')
                fprintf('  PLS did not return p-values, skipping.\n');
                continue;
            end

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

            % ----------------------------------------------------------
            % Plot LV1 and LV2
            % ----------------------------------------------------------
            max_LV = min(2, length(pvals));

            for LV = 1:max_LV
                p        = pvals(LV);
                headline = sprintf('Behavioural PLS: %s %s Ch%d, LV%d p=%.4f', ...
                    stage, lighting, channel, LV, p);

                figure('Position', [100 100 1200 500]);

                % --- Subplot 1: Behaviour correlations (18 bars) ---
                subplot(1,2,1)
                z    = pls_result.boot_result.orig_corr(:, LV);
                yneg = pls_result.boot_result.llcorr(:,  LV);
                ypos = pls_result.boot_result.ulcorr(:,  LV);

                bar(z)
                hold on
                grid on
                errorbar(1:length(z), z, yneg - z, ypos - z, '.', 'Color', 'black')
                xticks(1:length(z))
                xticklabels(behav_column_names)
                xtickangle(45)
                ylabel('Correlation with LV')
                title(sprintf('Behaviour Correlations, LV%d', LV))
                set(gca, 'FontSize', 12)
                hold off

                % --- Subplot 2: MSE salience heatmap ---
                subplot(1,2,2)
                x = pls_result.boot_result.compare_u(:, LV);
                x(abs(x) < 2.3) = 0;
                plotdata = reshape(x, length(scales), []);

                imagesc(plotdata)
                colormap(rgb)
                clim([-7 7])
                colorbar
                yticks(1:5:length(scales))
                yticklabels(1:5:length(scales))
                xticks([])
                title(sprintf('MSE Salience LV%d', LV))
                ylabel('Temporal Scale')
                set(gca, 'FontSize', 12)

                sgtitle(headline, 'FontSize', 16, 'Interpreter', 'none')
                set(gcf, 'Color', 'w')

                % Save figure
                fig_png = fullfile(results_dir, sprintf('Beh_PLS_%s_%s_ch%d_LV%d.png', ...
                    stage, lighting, channel, LV));
                fig_fig = fullfile(results_dir, sprintf('Beh_PLS_%s_%s_ch%d_LV%d.fig', ...
                    stage, lighting, channel, LV));

                saveas(gcf, fig_png)
                saveas(gcf, fig_fig)
                fprintf('  Saved: %s\n', fig_png);

                close(gcf);
            end

        end % lighting
    end % stage
end % channel

fprintf('\nBehavioural PLS analysis complete!\n');
end
