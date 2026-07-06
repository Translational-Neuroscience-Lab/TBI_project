% === Main Script ===
%First part takes the EDF file and transforms it into power values and
%creates PSD table
%Needs the TBI and sham EDF and tsv files to be store in separate folders
%within the "processed_EDF" folder
cd('S:\LAB_Kent\Studies\TBI\Scored_Data\NLF_6m_subacute\NLF_6m_tsv-edf_files\Processed_EDF')
epoch_length = 10; % seconds
folder_name = 'Processed EDF';
if ~exist(folder_name, 'dir')
    mkdir(folder_name);
end

% Grab EDF files from both TBI and Sham subfolders
base_path = 'S:\LAB_Kent\Studies\TBI\Scored_Data\NLF_6m_subacute\NLF_6m_tsv-edf_files\Processed_EDF';
edf_files = [dir(fullfile(base_path, 'TBI',  '*.edf')); ...
             dir(fullfile(base_path, 'Sham', '*.edf'))];

for f = 1:length(edf_files)
    file_path = fullfile(edf_files(f).folder, edf_files(f).name);
    fprintf('Processing %s\n', edf_files(f).name)

    preproc = pre_process(file_path, epoch_length);
    cleaned_data = preprocessing_artifact_p2p_matrix(preproc);
    psd_struct = compute_full_psd(cleaned_data);
    %psd = compute_full_psd(cleaned_data);

    % Initialize matrices
    freq = psd_struct.freq;
    power_matrix_EEG1 = psd_struct.EEG1_psd;  % Matrix: epochs x frequencies, NaN for rejected epochs
    power_matrix_EEG2 = psd_struct.EEG2_psd;  % Matrix: epochs x frequencies, NaN for rejected epochs

    % Create variable names
    freq_str = strrep(string(freq), '.', '_');
    var_names1 = strcat('EEG1_f_', freq_str, 'Hz');
    var_names2 = strcat('EEG2_f_', freq_str, 'Hz');
    all_var_names = [var_names1, var_names2];

    % Combine into table
    power_matrix = [power_matrix_EEG1, power_matrix_EEG2];
    psd_table = array2table(power_matrix, 'VariableNames', all_var_names);
    
    % Add time column
    psd_table.StartTime = cleaned_data.ts;
    psd_table = movevars(psd_table, 'StartTime', 'Before', 1);

    % Save
    [~, base_name, ~] = fileparts(edf_files(f).name);        % Get base filename without extension
    base_name = strrep(base_name, '_export', '');            % Remove '_export' if present
    out_name = fullfile(folder_name, [base_name '_psd.tsv']); % Construct output filename
    writetable(psd_table, out_name, 'FileType', 'text', 'Delimiter', '\t');

% --- CLEAR TO FREE MEMORY ---
    clear preproc cleaned_data psd_struct power_matrix power_matrix_EEG1 power_matrix_EEG2 psd_table
end

%%
%Second part: takes PSD table and scores, adds the treatment, mouse ID,
%vigilant state, light period, aggregates, normalizes and creates power spectra file

% Set the directory containing your data files
base_dir = 'S:\LAB_Kent\Studies\TBI\Scored_Data\NLF_6m_subacute\NLF_6m_tsv-edf_files\Processed_EDF';
data_dirs = {fullfile(base_dir, 'TBI'), fullfile(base_dir, 'Sham')};

fprintf('Starting PSD data processing...\n');

% Load and process all files
aggregated_data = process_all_psd_files(data_dirs);

% Save results to Excel
save_to_excel(aggregated_data, 'Power_Spectra_Clean_subacute6_lights.xlsx');

% New separate EEG files
create_separate_eeg_files(aggregated_data);

fprintf('\n=== PROCESSING COMPLETE ===\n');

%%
% === Function: Preprocessing ===
function output = pre_process(input_file, epoch_len)
    info = edfinfo(input_file);
    data = edfread(input_file);

    % Sampling rate
    fs = info.NumSamples(1) / seconds(info.DataRecordDuration); 
    epoch_dur = round(fs * epoch_len);  % samples per epoch

    % Flatten signals
    raw1 = data{:, 1};  
    raw2 = data{:, 2};  

    if iscell(raw1)
        EEG1 = vertcat(raw1{:});
        EEG2 = vertcat(raw2{:});
    else
        EEG1 = raw1;
        EEG2 = raw2;
    end

    % Epoch count
    total_samples = info.NumDataRecords * info.NumSamples(1);
    total_epochs = floor(total_samples / epoch_dur);

    % Cut signal
    EEG1 = EEG1(1 : total_epochs * epoch_dur);
    EEG2 = EEG2(1 : total_epochs * epoch_dur);

    % Reshape
    EEG1_epochs = reshape(EEG1, epoch_dur, [])';
    EEG2_epochs = reshape(EEG2, epoch_dur, [])';


    % Start time
    try
        start_time = datetime(strcat(info.StartDate, {' '}, info.StartTime), ...
                              'InputFormat', 'dd.MM.yy HH.mm.ss');
    catch
        % fallback format
        start_time = datetime(strcat(info.StartDate, {' '}, info.StartTime), ...
                              'InputFormat', 'dd.MM.yyyy HH.mm.ss');
    end
    start_time = start_time(1);
    ts = start_time + seconds((0:total_epochs-1)' * epoch_len);

    % Output
    output = struct();
    output.EEG1 = EEG1_epochs;
    output.EEG2 = EEG2_epochs;
    output.fs   = fs;
    output.ts   = ts;
end


%%
%Artifact rejection function P2P and Butterworth nadpass filter
function cleaned_data = preprocessing_artifact_p2p_matrix(preproc)
    % ======= CONFIGURATION =======
    fs = preproc.fs;
    EEG1 = preproc.EEG1;
    EEG2 = preproc.EEG2;
    ts = preproc.ts;
    percentile_threshold = 95; %change this to whatever you need

    %Filter
    [b, a] = butter(2, [0.5 100] / (fs / 2), 'bandpass');
    n_epochs = size(EEG1, 1);
    epoch_samples = size(EEG1, 2);

    % Initialize output matrices with NaN
    cleaned_EEG1 = NaN(n_epochs, epoch_samples);
    cleaned_EEG2 = NaN(n_epochs, epoch_samples);
    filtered_EEG1 = NaN(n_epochs, epoch_samples);
    filtered_EEG2 = NaN(n_epochs, epoch_samples);

    % Arrays to store metrics for thresholding
    all_max_ch1 = NaN(n_epochs, 1);
    all_p2p_ch1 = NaN(n_epochs, 1);
    all_max_ch2 = NaN(n_epochs, 1);
    all_p2p_ch2 = NaN(n_epochs, 1);
    rejection_reasons = strings(n_epochs, 1);

    % ======= LOOP THROUGH EPOCHS =======
    for i = 1:n_epochs
        epoch1 = double(EEG1(i, :));
        epoch2 = double(EEG2(i, :));

     %Apply filter
        epoch_filt1 = filtfilt(b, a, epoch1);
        epoch_filt2 = filtfilt(b, a, epoch2);
        
     %Store filtered data
        filtered_EEG1(i, :) = epoch_filt1;
        filtered_EEG2(i, :) = epoch_filt2;

     % Calculate metrics for channel 1
        all_max_ch1(i) = max(abs(epoch_filt1));
        all_p2p_ch1(i) = max(epoch_filt1) - min(epoch_filt1);
        
     % Calculate metrics for channel 2
        all_max_ch2(i) = max(abs(epoch_filt2));
        all_p2p_ch2(i) = max(epoch_filt2) - min(epoch_filt2);
    end

    % ======= DEFINE THRESHOLDS - SEPARATE FOR EACH CHANNEL =======
    max_thresh_ch1 = 399; % µV
    p2p_thresh_ch1 = prctile(all_p2p_ch1, percentile_threshold);
    
    max_thresh_ch2 = 399; % µV
    p2p_thresh_ch2 = prctile(all_p2p_ch2, percentile_threshold);

    % ======= SECOND PASS: APPLY THRESHOLDS =======
    keep_idx_ch1 = true(n_epochs, 1);
    keep_idx_ch2 = true(n_epochs, 1);
    
    % Start with all filtered data
    cleaned_EEG1 = filtered_EEG1;
    cleaned_EEG2 = filtered_EEG2;

    for i = 1:n_epochs
        reason_ch1 = [];
        reason_ch2 = [];

        %Check thresholds for channel 1
        if all_max_ch1(i) > max_thresh_ch1
            reason_ch1 = [reason_ch1, "ch1_high_max"];
        end
        if all_p2p_ch1(i) > p2p_thresh_ch1
            reason_ch1 = [reason_ch1, "ch1_high_p2p"];
        end
        
        % Check thresholds for channel 2
        if all_max_ch2(i) > max_thresh_ch2
            reason_ch2 = [reason_ch2, "ch2_high_max"];
        end
        if all_p2p_ch2(i) > p2p_thresh_ch2
            reason_ch2 = [reason_ch2, "ch2_high_p2p"];
        end
        
        %Apply rejection logic
        if ~isempty(reason_ch1)
            cleaned_EEG1(i, :) = NaN;
            keep_idx_ch1(i) = false;
        end

        if ~isempty(reason_ch2)
            cleaned_EEG2(i, :) = NaN;
            keep_idx_ch2(i) = false;
        end

        % Combine rejection reasons
        all_reasons = [reason_ch1, reason_ch2];
        if ~isempty(all_reasons)
            rejection_reasons(i) = strjoin(all_reasons, " & ");
        end
    end

 % Overall keep index (epoch is kept if BOTH channels are clean)
    keep_idx = keep_idx_ch1 & keep_idx_ch2;

 % ======= CREATE OUTPUT STRUCTURE =======
    cleaned_data = struct();
    cleaned_data.EEG1 = cleaned_EEG1;  % Filtered data with NaN for rejected epochs
    cleaned_data.EEG2 = cleaned_EEG2;  % Filtered data with NaN for rejected epochs
    cleaned_data.EEG1_original = EEG1; % Original unfiltered data (all epochs)
    cleaned_data.EEG2_original = EEG2; % Original unfiltered data (all epochs)
    cleaned_data.fs = fs;
    cleaned_data.ts = ts;
    cleaned_data.keep_idx = keep_idx;
    cleaned_data.keep_idx_ch1 = keep_idx_ch1;  % Individual channel rejection info
    cleaned_data.keep_idx_ch2 = keep_idx_ch2;  % Individual channel rejection info
    cleaned_data.rejected_idx = find(~keep_idx);
    cleaned_data.rejection_reasons = rejection_reasons;
    cleaned_data.n_clean_epochs = sum(keep_idx);
    cleaned_data.n_rejected_epochs = sum(~keep_idx);
    cleaned_data.rejection_rate = sum(~keep_idx) / n_epochs * 100;
    
    % Summary statistics - separate thresholds for each channel
    cleaned_data.thresholds = struct();
    cleaned_data.thresholds.ch1.max_thresh = max_thresh_ch1;
    cleaned_data.thresholds.ch1.p2p_thresh = p2p_thresh_ch1;
    cleaned_data.thresholds.ch2.max_thresh = max_thresh_ch2;
    cleaned_data.thresholds.ch2.p2p_thresh = p2p_thresh_ch2;
    
    % Display summary
    fprintf('Artifact Rejection Summary:\n');
    fprintf('Total epochs: %d\n', n_epochs);
    fprintf('Clean epochs (both channels): %d\n', cleaned_data.n_clean_epochs);
    fprintf('Rejected epochs: %d (%.1f%%)\n', cleaned_data.n_rejected_epochs, cleaned_data.rejection_rate);
    fprintf('Channel 1 clean epochs: %d (%.1f%%)\n', sum(keep_idx_ch1), sum(keep_idx_ch1)/n_epochs*100);
    fprintf('Channel 2 clean epochs: %d (%.1f%%)\n', sum(keep_idx_ch2), sum(keep_idx_ch2)/n_epochs*100);
    fprintf('Ch1 Thresholds - Max: %.2f, P2P: %.2f\n', max_thresh_ch1, p2p_thresh_ch1);
    fprintf('Ch2 Thresholds - Max: %.2f, P2P: %.2f\n', max_thresh_ch2, p2p_thresh_ch2);
end

%%
% === Function: PSD Calculation ===
function psd_struct = compute_full_psd(cleaned_data)
    fs = cleaned_data.fs;
    n_epochs = size(cleaned_data.EEG1, 1);
    desired_freqs = (0:1:40)';
    n_freqs = length(desired_freqs);
    
    % Initialize matrices with NaN
    EEG1_psd_matrix = NaN(n_epochs, n_freqs);
    EEG2_psd_matrix = NaN(n_epochs, n_freqs);

    % Keep track of which epochs were processed
    processed_epochs = false(n_epochs, 1);

    for i = 1:n_epochs
         % Check if this epoch contains NaN (was rejected)
        if any(isnan(cleaned_data.EEG1(i, :))) || any(isnan(cleaned_data.EEG2(i, :)))
            % This epoch was rejected - leave as NaN (already initialized)
            processed_epochs(i) = false;
        else
        x1 = cleaned_data.EEG1(i, :)';
        x2 = cleaned_data.EEG2(i, :)';

        [pxx1, f] = pwelch(x1, [], [], [], fs);
        [pxx2, ~] = pwelch(x2, [], [], [], fs);

        % Interpolate to desired 1-Hz resolution from 0 to 40 Hz
        pxx1_interp = interp1(f, pxx1, desired_freqs, 'linear', 'extrap');
        pxx2_interp = interp1(f, pxx2, desired_freqs, 'linear', 'extrap');
            
        % Store in matrix
        EEG1_psd_matrix(i, :) = pxx1_interp;
        EEG2_psd_matrix(i, :) = pxx2_interp;
        processed_epochs(i) = true;
        end
    end
% Create output structure
    psd_struct = struct();
    psd_struct.freq = desired_freqs;
    psd_struct.EEG1_psd = EEG1_psd_matrix;  % Matrix format: epochs x frequencies
    psd_struct.EEG2_psd = EEG2_psd_matrix;  % Matrix format: epochs x frequencies
    psd_struct.processed_epochs = processed_epochs;
    psd_struct.n_clean_epochs = sum(processed_epochs);
    psd_struct.n_rejected_epochs = sum(~processed_epochs);
    
    % Display summary
    fprintf('PSD Computation Summary:\n');
    fprintf('Total epochs: %d\n', n_epochs);
    fprintf('Processed epochs: %d\n', psd_struct.n_clean_epochs);
    fprintf('Skipped (NaN) epochs: %d\n', psd_struct.n_rejected_epochs);
end
%%
%Fundtion to add mouse ID, treatment and create power spectra file
function aggregated_data = process_all_psd_files(data_dirs)
    % Define treatment groups
    TBI_animals = {'Mouse#3','Mouse#4','Mouse#7','Mouse#6', 'Mouse#12','Mouse#13','Mouse#14','Mouse#15','Mouse#16','Mouse#17'};
    Sham_animals = {'Mouse#1','Mouse#2','Mouse#5','Mouse#10','Mouse#11','Mouse#18','Mouse#19','Mouse#20'};
    
    % Allow single directory string to be passed as well
    if ischar(data_dirs) || isstring(data_dirs)
        data_dirs = {char(data_dirs)};
    end

    aggregated_data = table();

    % Loop through each directory
    for dir_idx = 1:length(data_dirs)
        current_dir = data_dirs{dir_idx};
        psd_files = dir(fullfile(current_dir, '*_psd.tsv'));
        fprintf('Found %d PSD files in: %s\n\n', length(psd_files), current_dir);

        % Loop through each PSD file in this directory
        for file_idx = 1:length(psd_files)
            try
                psd_filename = psd_files(file_idx).name;
                fprintf('Processing file %d/%d: %s\n', file_idx, length(psd_files), psd_filename);

                % Pass the specific directory for this file
                file_data = process_single_file(current_dir, psd_filename, TBI_animals, Sham_animals);

                if ~isempty(file_data)
                    aggregated_data = [aggregated_data; file_data];
                end

            catch ME
                fprintf('  ERROR processing %s: %s\n\n', psd_filename, ME.message);
                continue;
            end
        end
    end
end
%%
function file_aggregated_data = process_single_file(data_dir, psd_filename, TBI_animals, Sham_animals)
    psd_filepath = fullfile(data_dir, psd_filename);
    
    % Generate corresponding scores filename
    scores_filename = strrep(psd_filename, '_psd.tsv', '_scores.tsv');
    scores_filepath = fullfile(data_dir, scores_filename);
    
    % Check if corresponding scores file exists
    if ~exist(scores_filepath, 'file')
        fprintf('  WARNING: Scores file not found: %s\n', scores_filename);
        fprintf('  Skipping this file...\n\n');
        file_aggregated_data = table();
        return;
    end
    
    % Load and process PSD data
    psd_table = load_psd_data(psd_filepath);
    
    % Load and process scores data
    scores_tt = load_scores_data(scores_filepath);
    
    % Match timestamps and assign classes
    psd_table = match_timestamps(psd_table, scores_tt);
    
    % Extract mouse ID and assign treatment
    mouse_id = extract_mouse_id(psd_filename);
    treatment = assign_treatment(mouse_id, TBI_animals, Sham_animals);
    
    % Add mouse ID and treatment columns
    psd_table.mouse_id = repmat({mouse_id}, height(psd_table), 1);
    psd_table.treatment = repmat(treatment, height(psd_table), 1);
    
    % Aggregate power spectral data
    file_aggregated_data = aggregate_power_spectra(psd_table, mouse_id, treatment);
    
    fprintf('  Mouse ID: %s, Treatment: %s\n', mouse_id, treatment);
    fprintf('  Successfully processed: %d PSD records\n\n', height(psd_table));
end

%%
function psd_table = load_psd_data(psd_filepath)
    psd_table = readtable(psd_filepath, 'FileType', 'text', 'Delimiter', '\t');
    
    % Check datetime column in PSD file (assuming it's 'StartTime')
    if ~ismember('StartTime', psd_table.Properties.VariableNames)
        error('PSD file missing "StartTime" datetime column');
    end
    
    % Convert PSD timestamp to datetime if not already
    if ~isdatetime(psd_table.StartTime)
        psd_table.StartTime = datetime(psd_table.StartTime, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
    else
        % Remove timezone if it exists
        psd_table.StartTime = datetime(psd_table.StartTime, 'TimeZone', '');
    end
    
    % Sort by time
    psd_table = sortrows(psd_table, 'StartTime');
end

%%
function scores_tt = load_scores_data(scores_filepath)
    opts = detectImportOptions(scores_filepath, 'FileType', 'text', 'Delimiter', '\t');

    % Explicitly set date format to avoid MM/dd vs dd/MM ambiguity ----
    date_var = opts.VariableNames{1};  % column 1 is the date column
    opts = setvaropts(opts, date_var, 'InputFormat', 'MM/dd/yyyy');
    % ----
    
    score_table = readtable(scores_filepath, opts);
    
    % Combine date (col 1) and time (col 2) into datetime
    score_dates = string(score_table{:,1});
    score_times = string(score_table{:,2});
    score_dt_str = strcat(score_dates, {' '}, score_times);
    
    % Parse datetime
    scores_dt = datetime(score_dt_str, 'InputFormat', 'MM/dd/yyyy HH:mm:ss');
    
    % Round scores datetime to nearest 10 seconds (to fix offset)
    scores_seconds = second(scores_dt);
    rounded_seconds = 10 * floor(scores_seconds / 10);
    
    % Create new datetime with rounded seconds
    scores_dt_rounded = datetime(year(scores_dt), month(scores_dt), day(scores_dt), ...
        hour(scores_dt), minute(scores_dt), rounded_seconds);
    
    % Extract class column (usually column 5)
    if width(score_table) < 5
        error('Scores file missing class column (5th column)');
    end
    class_codes = score_table{:,5};
    
    % Map numeric class codes
    class_map = containers.Map({'255','1','2','3'}, {'Unscored','Wake','NREM','REM'});
    
    % Convert class codes to strings properly
    if isnumeric(class_codes)
        class_codes_str = cellstr(string(class_codes));
    end
    
    % Round and convert to proper format
    class_codes_str = cellfun(@(x) num2str(round(str2double(x))), class_codes_str, 'UniformOutput', false);
    
    % Map class codes to labels
    class_labels = cell(size(class_codes_str));
    for i = 1:numel(class_codes_str)
        key = class_codes_str{i};
        if isKey(class_map, key)
            class_labels{i} = class_map(key);
        else
            class_labels{i} = 'Unknown';
        end
    end
    
    % Create timetable for scores
    scores_tt = timetable(scores_dt_rounded, class_labels, 'VariableNames', {'class'});
    scores_tt = sortrows(scores_tt);
end

%%
function psd_table = match_timestamps(psd_table, scores_tt)
    % Initialize class column for PSD
    psd_table.class = repmat({'Unscored'}, height(psd_table), 1);
    
    % Match timestamps using a tolerance-based approach
    tol_seconds = 5; % tolerance in seconds
    
    for i = 1:height(psd_table)
        psd_time = psd_table.StartTime(i);
        
        % Find scores within tolerance
        time_diff = abs(seconds(scores_tt.scores_dt_rounded - psd_time));
        within_tolerance = time_diff <= tol_seconds;
        
        if any(within_tolerance)
            % Find the closest match within tolerance
            [~, closest_idx] = min(time_diff(within_tolerance));
            valid_indices = find(within_tolerance);
            matched_idx = valid_indices(closest_idx);
            
            % Assign the class label
            psd_table.class(i) = scores_tt.class(matched_idx);
        end
    end
    psd_table = assign_lights_phase(psd_table);
end

%%
function mouse_id = extract_mouse_id(psd_filename)
    % Extract mouse ID from filename
    mouse_id_match = regexp(psd_filename, '^(Mouse[^_]+)', 'tokens');
    if ~isempty(mouse_id_match)
        mouse_id = mouse_id_match{1}{1};
    else
        mouse_id = 'Unknown';
        fprintf('  WARNING: Could not extract mouse ID from filename: %s\n', psd_filename);
    end
end

%%
function treatment = assign_treatment(mouse_id, TBI_animals, Sham_animals)
    if ismember(mouse_id, TBI_animals)
        treatment = "TBI";
    elseif ismember(mouse_id, Sham_animals)
        treatment = "Sham";
    else
        treatment = "Unknown";
        fprintf('  WARNING: Mouse ID "%s" not found in treatment groups\n', mouse_id);
    end
end

%%
function aggregated_data = aggregate_power_spectra(psd_table, mouse_id, treatment)
    % Find frequency columns (0-30 Hz)
    [freq_cols, freq_values] = find_frequency_columns(psd_table);
    
    if isempty(freq_cols)
        fprintf('  WARNING: No frequency columns found (0-30 Hz)\n');
        aggregated_data = table();
        return;
    end
    
    fprintf('  Found %d frequency columns (0-30 Hz)\n', length(freq_cols));
    
    % Separate EEG channels
    [eeg1_cols, eeg1_freqs, eeg2_cols, eeg2_freqs] = separate_eeg_channels(freq_cols, freq_values);
    
    % Create aggregated data for each sleep state, light phase, and EEG channel
    sleep_states = {'Wake', 'NREM', 'REM'};
    lights_phases = {'Light', 'Dark'};
    aggregated_data = table();
    
    for state_idx = 1:length(sleep_states)
        current_state = sleep_states{state_idx};

        for phase_idx = 1:length(lights_phases)
            current_phase = lights_phases{phase_idx};

        %Mask for this sleep state AND lights phase combination
        state_mask = strcmp(psd_table.class, current_state) & ...
                     strcmp(psd_table.lights_phase, current_phase);
        
        if ~any(state_mask)
            continue; % Skip if no data for this state
        end
        
        % Process EEG1 data
        if ~isempty(eeg1_cols)
            agg_record = create_aggregated_record(psd_table, state_mask, eeg1_cols, eeg1_freqs, ...
                mouse_id, treatment, current_state, current_phase, 'EEG1');
            aggregated_data = [aggregated_data; agg_record];
        end
        
        % Process EEG2 data
        if ~isempty(eeg2_cols)
            agg_record = create_aggregated_record(psd_table, state_mask, eeg2_cols, eeg2_freqs, ...
                mouse_id, treatment, current_state, current_phase, 'EEG2');
            aggregated_data = [aggregated_data; agg_record];
            end
        end
    end
end

%%
function [freq_cols, freq_values] = find_frequency_columns(psd_table)
    freq_cols = {};
    freq_values = [];
    
    % Find columns that represent frequencies (0-30 Hz)
    for col_idx = 1:width(psd_table)
        col_name = psd_table.Properties.VariableNames{col_idx};
        
        % Initialize freq_val for this iteration
        freq_val = NaN;
        
        % Try different naming patterns for frequency columns
        freq_match = [];
        
        % Pattern 1: Hz_X_Y format (e.g., Hz_1_5 for 1.5 Hz)
        freq_match = regexp(col_name, '^EEG[12]_f_(\d+)Hz$', 'tokens');
        if ~isempty(freq_match)
            freq_val = str2double(freq_match{1}{1});
        else
            % Pattern 2: Direct frequency as column name (e.g., '1.5', '2.0')
            freq_match = regexp(col_name, '^(\d+\.?\d*)$', 'tokens');
            if ~isempty(freq_match)
                freq_val = str2double(freq_match{1}{1});
            else
                % Pattern 3: Freq_X format (e.g., Freq_1_5)
                freq_match = regexp(col_name, '^Freq_(\d+)_(\d+)$', 'tokens');
                if ~isempty(freq_match)
                    freq_val = str2double(freq_match{1}{1}) + str2double(freq_match{1}{2})/10;
                end
            end
        end
        
        % Only include frequencies between 0 and 30 Hz
        if ~isnan(freq_val) && freq_val >= 0 && freq_val <= 30
            freq_cols{end+1} = col_name;
            freq_values(end+1) = freq_val;
        end
    end
end

%%
function [eeg1_cols, eeg1_freqs, eeg2_cols, eeg2_freqs] = separate_eeg_channels(freq_cols, freq_values)
    eeg1_cols = {};
    eeg2_cols = {};
    eeg1_freqs = [];
    eeg2_freqs = [];
    
    for i = 1:length(freq_cols)
        col_name = freq_cols{i};
        freq_val = freq_values(i);
        
        % Check if column name indicates EEG1 or EEG2
        if contains(col_name, 'EEG1', 'IgnoreCase', true) || contains(col_name, 'Ch1', 'IgnoreCase', true)
            eeg1_cols{end+1} = col_name;
            eeg1_freqs(end+1) = freq_val;
        elseif contains(col_name, 'EEG2', 'IgnoreCase', true) || contains(col_name, 'Ch2', 'IgnoreCase', true)
            eeg2_cols{end+1} = col_name;
            eeg2_freqs(end+1) = freq_val;
        else
            % If no clear channel indicator, assume it's a general frequency column
            eeg1_cols{end+1} = col_name;
            eeg1_freqs(end+1) = freq_val;
        end
    end
end

%%
function agg_record = create_aggregated_record(psd_table, state_mask, eeg_cols, eeg_freqs, ...
    mouse_id, treatment, current_state, current_phase, eeg_channel)
    
    eeg_data = psd_table(state_mask, eeg_cols);
    if height(eeg_data) > 0
        % Calculate mean power for each frequency
        mean_powers = mean(table2array(eeg_data), 1, 'omitnan');
        
        % Create aggregated record
        agg_record = table();
        agg_record.Mouse_ID = {mouse_id};
        agg_record.Treatment = treatment;
        agg_record.Sleep_State = {current_state};
        agg_record.Lights_Phase = {current_phase};
        agg_record.EEG_Channel = {eeg_channel};
        agg_record.N_Epochs = sum(state_mask);
        
        % Add frequency columns
        for freq_idx = 1:length(eeg_freqs)
            freq_col_name = sprintf('Hz_%.1f', eeg_freqs(freq_idx));
            freq_col_name = strrep(freq_col_name, '.', '_'); % Replace . with _ for valid column names
            agg_record.(freq_col_name) = mean_powers(freq_idx);
        end
    else
        agg_record = table();
    end
end

%%
function save_to_excel(aggregated_data, excel_filename)
    if isempty(aggregated_data)
        fprintf('WARNING: No aggregated data was created. Check your frequency column naming.\n');
        return;
    end
    
    % Create normalized version of the data
    aggregated_data_normalized = create_normalized_data(aggregated_data);
    
    % Save both versions to Excel with different sheets
    writetable(aggregated_data, excel_filename, 'Sheet', 'Raw_Power');
    writetable(aggregated_data_normalized, excel_filename, 'Sheet', 'Normalized_Power');
    
    fprintf('Power spectral data saved to: %s\n', excel_filename);
    fprintf('  - Sheet 1: Raw_Power (absolute power values)\n');
    fprintf('  - Sheet 2: Normalized_Power (normalized by total power per mouse/state/channel)\n');
    fprintf('Total aggregated records: %d\n', height(aggregated_data));
    
    % Display breakdown
    display_data_breakdown(aggregated_data);
end

%%
function aggregated_data_normalized = create_normalized_data(aggregated_data)
    % Create normalized version of the data
    aggregated_data_normalized = aggregated_data;
    
    % Get frequency columns (Hz_X_X format)
    freq_columns = aggregated_data.Properties.VariableNames;
    hz_columns = freq_columns(contains(freq_columns, '_f_') & contains(freq_columns, 'Hz'));
    
    if ~isempty(hz_columns)
        fprintf('Creating normalized power spectral data...\n');
        
        % Normalize each row individually (each row = one mouse/state/channel combination)
        for row_idx = 1:height(aggregated_data)
            % Get power values for this row (all frequencies)
            power_values = table2array(aggregated_data(row_idx, hz_columns));
            
            % Calculate total power across all frequencies for this mouse/state/channel
            total_power = sum(power_values);
            
            % Normalize each frequency by the total power
            if total_power > 0
                normalized_values = power_values ./ total_power;
                
                % Update the normalized table
                aggregated_data_normalized(row_idx, hz_columns) = array2table(normalized_values);
            else
                fprintf('  WARNING: Zero total power for row %d (%s-%s-%s)\n', ...
                    row_idx, aggregated_data.Mouse_ID{row_idx}, ...
                    aggregated_data.Sleep_State{row_idx}, aggregated_data.EEG_Channel{row_idx});
            end
        end
        
        % Verify normalization (each row should sum to 1.0)
        normalized_sums = sum(table2array(aggregated_data_normalized(:, hz_columns)), 2);
        fprintf('  Normalization check: mean row sum = %.6f (should be 1.0)\n', mean(normalized_sums));
    end
end

%%
function display_data_breakdown(aggregated_data)
    % Display breakdown by treatment and state
    fprintf('\nData breakdown:\n');
    unique_treatments = unique(aggregated_data.Treatment);
    unique_states = unique(aggregated_data.Sleep_State);
    unique_phases = unique(aggregated_data.Lights_Phase);
    unique_channels = unique(aggregated_data.EEG_Channel);
    
    for t = 1:length(unique_treatments)
        for s = 1:length(unique_states)
            for p = 1:length(unique_phases)
                for c = 1:length(unique_channels)
                    mask = strcmp(aggregated_data.Treatment, unique_treatments{t}) & ...
                           strcmp(aggregated_data.Sleep_State, unique_states{s}) & ...
                           strcmp(aggregated_data.Lights_Phase, unique_phases{p}) & ...
                           strcmp(aggregated_data.EEG_Channel, unique_channels{c});
                    count = sum(mask);
                    if count > 0
                        fprintf('  %s - %s - %s: %d mice\n', unique_treatments{t}, unique_states{s}, unique_phases{p}, unique_channels{c}, count);
                    end
                end
            end
        end
    end
    
    % Show frequency range covered
    freq_columns = aggregated_data.Properties.VariableNames;
    hz_columns = freq_columns(startsWith(freq_columns, 'Hz_'));
    if ~isempty(hz_columns)
        fprintf('\nFrequency range: %s to %s\n', hz_columns{1}, hz_columns{end});
        fprintf('Total frequencies analyzed: %d\n', length(hz_columns));
    end
end

%%
%Function to create the appropriate files for PLS
% Script to create separate Excel files for EEG1 and EEG2 (Raw and Normalized)

function create_separate_eeg_files(aggregated_data)
    if isempty(aggregated_data)
        fprintf('WARNING: No aggregated data available for creating separate EEG files.\n');
        return;
    end
    
    fprintf('\nCreating separate EEG Excel files...\n');
    
    % Create normalized version of the data
    aggregated_data_normalized = create_normalized_data_for_separate_files(aggregated_data);
    
    % Create files for each EEG channel and data type
    create_eeg_file(aggregated_data, 'EEG1', 'raw', 'EEG1_Raw_Power.xlsx');
    create_eeg_file(aggregated_data_normalized, 'EEG1', 'normalized', 'EEG1_Normalized_Power.xlsx');
    create_eeg_file(aggregated_data, 'EEG2', 'raw', 'EEG2_Raw_Power.xlsx');
    create_eeg_file(aggregated_data_normalized, 'EEG2', 'normalized', 'EEG2_Normalized_Power.xlsx');
    
    fprintf('All EEG files created successfully!\n\n');
end

function create_eeg_file(data, eeg_channel, data_type, filename)
    % Filter data for specific EEG channel
    eeg_data = data(strcmp(data.EEG_Channel, eeg_channel), :);
    
    if isempty(eeg_data)
        fprintf('WARNING: No data found for %s. Skipping %s\n', eeg_channel, filename);
        return;
    end
    
    % Get frequency columns (Hz_X_X format)
    freq_columns = eeg_data.Properties.VariableNames;
    hz_columns = freq_columns(startsWith(freq_columns, 'Hz_'));
    
    if isempty(hz_columns)
        fprintf('WARNING: No frequency columns found for %s. Skipping %s\n', eeg_channel, filename);
        return;
    end
    
    % Sort frequency columns numerically
    hz_columns = sort_frequency_columns(hz_columns);
    
    % Define grouping orders
    sleep_states = {'NREM', 'REM', 'Wake'};
    lights_phases = {'Light', 'Dark'};
    treatment_order = {'TBI', 'Sham'};
    
    % Pre-initialize final_table with correct columns to avoid vertcat mismatch
    empty_freq_row = array2table(zeros(0, length(hz_columns)), 'VariableNames', hz_columns);
    final_table = empty_freq_row;
    mouse_ids = {};
    treatments = {};
    sleep_state_labels = {};
    phase_labels = {};
    
    for state_idx = 1:length(sleep_states)
        current_state = sleep_states{state_idx};

        for phase_idx = 1:length(lights_phases)
            current_phase = lights_phases{phase_idx};
            
            % Filter for this sleep state and lights phase
            group_data = eeg_data(strcmp(eeg_data.Sleep_State,  current_state) & ...
                                  strcmp(eeg_data.Lights_Phase, current_phase), :);
        
        if isempty(group_data)
            fprintf('  No data for %s state in %s\n', current_state, current_phase, eeg_channel);
            continue;
        end
        
        % Process each treatment group in order (TBI first, then Sham)
        for treat_idx = 1:length(treatment_order)
            current_treatment = treatment_order{treat_idx};
            treatment_data = group_data(strcmp(group_data.Treatment, current_treatment), :);
            
            if isempty(treatment_data)
                continue; % Skip if no data for this treatment in this state
            end
                               
            % Collect metadata for each row
            for mouse_idx = 1:height(treatment_data)
                mouse_ids{end+1} = treatment_data.Mouse_ID{mouse_idx};
                treatments{end+1} = char(treatment_data.Treatment(mouse_idx));
                sleep_state_labels{end+1} = current_state;
                phase_labels{end+1} = current_phase;
            end
            
            % Extract and append only the sorted hz_columns to avoid mismatch
                freq_data = treatment_data(:, hz_columns);
                final_table = [final_table; freq_data];
            end
        end
    end
    
    if isempty(final_table)
        fprintf('WARNING: No data assembled for %s. Skipping %s\n', eeg_channel, filename);
        return;
    end

    % Add identification columns at the beginning
    final_table = [table(mouse_ids', treatments', sleep_state_labels', phase_labels', ...
        'VariableNames', {'Mouse_ID', 'Treatment', 'Sleep_State', 'Lights_Phase'}), final_table];
    
    % Write to Excel
    writetable(final_table, filename);
    
    fprintf('Created %s (%s %s): %d animals x %d frequencies\n', ...
        filename, eeg_channel, data_type, height(final_table), length(hz_columns));
    
    % Display breakdown by sleep state
    display_file_breakdown(final_table, eeg_channel, data_type);
end

%%
function sorted_columns = sort_frequency_columns(hz_columns)
    % Extract frequency values from column names for proper sorting
    freq_values = zeros(size(hz_columns));
    
    for i = 1:length(hz_columns)
        col_name = hz_columns{i};
        % Extract frequency value from Hz_X_X format
        freq_match = regexp(col_name, 'Hz_(\d+)_(\d+)', 'tokens');
        if ~isempty(freq_match)
            freq_values(i) = str2double(freq_match{1}{1}) + str2double(freq_match{1}{2})/10;
        else
            % Try simpler Hz_X format
            freq_match = regexp(col_name, 'Hz_(\d+)', 'tokens');
            if ~isempty(freq_match)
                freq_values(i) = str2double(freq_match{1}{1});
            end
        end
    end
    
    % Sort by frequency value
    [~, sort_idx] = sort(freq_values);
    sorted_columns = hz_columns(sort_idx);
end

%%
function aggregated_data_normalized = create_normalized_data_for_separate_files(aggregated_data)
    % Create normalized version of the data for separate files
    aggregated_data_normalized = aggregated_data;
    
    % Get frequency columns (Hz_X_X format)
    freq_columns = aggregated_data.Properties.VariableNames;
    hz_columns = freq_columns(startsWith(freq_columns, 'Hz_'));
    
    if ~isempty(hz_columns)
        fprintf('Creating normalized power spectral data for separate files...\n');
        
        % Normalize each row individually
        for row_idx = 1:height(aggregated_data)
            % Get power values for this row
            power_values = table2array(aggregated_data(row_idx, hz_columns));
            
            % Calculate total power across all frequencies
            total_power = sum(power_values);
            
            % Normalize each frequency by the total power
            if total_power > 0
                normalized_values = power_values ./ total_power;
                aggregated_data_normalized(row_idx, hz_columns) = array2table(normalized_values);
            end
        end
    end
end

%%
function psd_table = assign_lights_phase(psd_table)
    % Lights are anchored to PST year-round (lab does not adjust for DST):
    %   Lights OFF (dark phase starts):  09:00 PST
    %   Lights ON  (light phase starts): 21:00 PST
    %
    % Recordings follow wall-clock Pacific time (PST or PDT depending on date).
    % During PDT, recorded timestamps are 1 hour ahead of PST, so the
    % equivalent lights-off/on in PDT wall-clock time becomes 10:00 / 22:00.
    %
    % Strategy: convert recorded wall-clock time back to PST by subtracting
    % the DST offset when active, then compare against 09:00 / 21:00 PST.
    %
    % DST (US Pacific):
    %   PDT starts: 2:00 AM on the 2nd Sunday of March   (e.g. March 13 2022)
    %   PDT ends:   2:00 AM on the 1st Sunday of November (e.g. Nov 6 2022)

    lights_off_pst = 9;   % 09:00 PST
    lights_on_pst  = 21;  % 21:00 PST

    n = height(psd_table);
    phase = repmat({'Unknown'}, n, 1);

    for i = 1:n
        t_wall = psd_table.StartTime(i);  % recorded wall-clock time (PST or PDT)
        yr = year(t_wall);

        % Compute DST boundaries for this year
        dst_start = nth_weekday_of_month(yr, 3, 1, 2) + hours(2); % 2nd Sun March 02:00
        dst_end   = nth_weekday_of_month(yr, 11, 1, 1) + hours(2); % 1st Sun November 02:00

        % Determine DST offset: PDT = UTC-7, PST = UTC-8
        % Wall-clock during PDT is 1 hour ahead of PST
        if t_wall >= dst_start && t_wall < dst_end
            dst_correction = 1;  % subtract 1 hour to convert PDT wall-clock -> PST equivalent
        else
            dst_correction = 0;  % already in PST, no correction needed
        end

        % Convert to PST-equivalent time
        t_pst = t_wall - hours(dst_correction);

        % Fractional hour for exact boundary comparison
        pst_hour = hour(t_pst) + minute(t_pst)/60 + second(t_pst)/3600;

        % Assign phase based on PST hour
        if pst_hour >= lights_off_pst && pst_hour < lights_on_pst
            phase{i} = 'Dark';
        else
            phase{i} = 'Light';
        end
    end

    psd_table.lights_phase = phase;
end

% === Helper: Find nth weekday of a given month/year ===
function dt = nth_weekday_of_month(yr, mo, weekday_num, n)
    % weekday_num: 1=Sun, 2=Mon, ..., 7=Sat (MATLAB convention)
    % n: which occurrence (1=first, 2=second, etc.)

    first_day = datetime(yr, mo, 1);
    first_dow = weekday(first_day);

    days_ahead = mod(weekday_num - first_dow, 7);
    first_occurrence = first_day + days(days_ahead);

    dt = first_occurrence + days(7 * (n - 1));
end

%%
function display_file_breakdown(final_table, eeg_channel, data_type)
    % Display breakdown by treatment and sleep state
    unique_treatments = unique(final_table.Treatment);
    unique_states = unique(final_table.Sleep_State);
    unique_phase = unique(final_table.Lights_Phase);
    
    fprintf('  Breakdown for %s %s:\n', eeg_channel, data_type);
    for t = 1:length(unique_treatments)
        for s = 1:length(unique_states)
            for p = 1:length(unique_phase)
                mask = strcmp(final_table.Treatment, unique_treatments{t}) & ...
                       strcmp(final_table.Sleep_State, unique_states{s}) & ...
                       strcmp(final_table.Lights_Phase, unique_phase{p});
                count = sum(mask);
                if count > 0
                    fprintf('    %s - %s: %d animals\n', unique_treatments{t}, unique_states{s}, unique_phase{p}, count);
                end
            end
        end
    end
end
