function preprocessing_with_mse_daylight(edf_path, tsv_path, output_base)
    % ======= CONFIGURATION =======
    fs = 400;
    [b, a] = butter(2, [0.5 50] / (fs / 2), 'bandpass');
    percentile_threshold = 90;
    unique_stages = ["Wake", "NREM", "REM"];

    % MSE parameters
    mse_scales = 1:40;  % Scale factors for MSE
    mse_r = 0.5;        % Tolerance for pattern matching (similarity criterion)
    mse_m = 2;          % Pattern length

    % ======= IDENTIFIERS =======
    [edf_folder, edf_file, ~] = fileparts(edf_path);
    condition = string(split(edf_folder, filesep));
    condition = condition(end);
    fprintf("Processing: %s - %s - EEG: %s\n", edf_file, condition);

    % ======= LOAD DATA =======
    data = edfread(edf_path);
    sleep_states = readtable(tsv_path, 'FileType', 'text');

    % Assuming you have two EEG channels - modify as needed
    eeg1_cells = data.SignalLabel1_EEGEEG1A_B;
    eeg2_cells = data.SignalLabel2_EEGEEG2A_B;

    % ======= PREPARE SLEEP STAGE INFO =======
    epoch_numbers = sleep_states.Var4;
    time = sleep_states.Var2;
    date = sleep_states.Var1;
    sleep_stage_codes = sleep_states.Var5;
     
    % ======= FILTER OUT NaT/NaN ROWS FROM THE .TSV FILE =======
    valid_rows = ~(strcmp(string(date), 'NaT') | strcmp(string(time), 'NaN') | isnan(epoch_numbers) | isnan(sleep_stage_codes));
    epoch_numbers = epoch_numbers(valid_rows);
    time = time(valid_rows);
    date = date(valid_rows);
    sleep_stage_codes = sleep_stage_codes(valid_rows);


    % ======= FILTER OUT INVALID SLEEP STAGE CODES =======
    valid_idx = sleep_stage_codes ~= 255;
    epoch_numbers = epoch_numbers(valid_idx);
    sleep_stage_codes = sleep_stage_codes(valid_idx);
    time = time(valid_idx);
    date = date(valid_idx);

    sleep_labels = strings(length(sleep_stage_codes), 1);
    sleep_labels(sleep_stage_codes == 1) = "Wake";
    sleep_labels(sleep_stage_codes == 2) = "NREM";
    sleep_labels(sleep_stage_codes == 3) = "REM";

    epoch_indices = epoch_numbers / 10 + 1;
    sleep_info = table (epoch_indices, epoch_numbers, sleep_stage_codes, sleep_labels, time, date,...
        'VariableNames', {'Index', 'Epoch', 'StageCode', 'StageLabel', 'time','date'});

    % Get total number of epochs from EEG data
    total_epochs = length(eeg1_cells);

    % ======= LIGHTING PERIOD CLASSIFICATION WITH DST =======
    fprintf('Classifying lighting periods with DST consideration...\n');

    % Get the first date (recording start date)
    if ~isempty(sleep_info.date)
        start_date_str = string(sleep_info.date(1));
        fprintf('Recording start date: %s\n', start_date_str);
        
        % Parse the start date to determine if it falls within DST period
        % DST period: March 10, 2024 to November 3, 2024
        try
            % Convert string date to datetime for comparison
            start_date_dt = datetime(start_date_str, 'InputFormat', 'dd-MMM-yyyy');
            
            % Define DST boundaries for 2022
            dst_start = datetime('13-Mar-2022', 'InputFormat', 'dd-MMM-yyyy');
            dst_end = datetime('06-Nov-2022', 'InputFormat', 'dd-MMM-yyyy');
            
            % Check if recording falls within DST period
            is_dst_period = (start_date_dt >= dst_start) && (start_date_dt <= dst_end);
            
            if is_dst_period
                lights_off_start = 10;  % 10 AM
                lights_off_end = 22;    % 10 PM
                lights_on_end = 10;     % 10 AM (next day)
                fprintf('DST period detected: Using 10AM-10PM schedule\n');
            else
                lights_off_start = 9;   % 9 AM
                lights_off_end = 21;    % 9 PM
                lights_on_end = 9;      % 9 AM (next day)
                fprintf('Standard time period: Using 9AM-9PM schedule\n');
            end
            
        catch date_error
            % If date parsing fails, default to standard time
            fprintf('Warning: Could not parse date %s, defaulting to standard time schedule\n', start_date_str);
            lights_off_start = 9;   % 9 AM
            lights_off_end = 21;    % 9 PM
            lights_on_end = 9;      % 9 AM (next day)
            is_dst_period = false;
        end
    else
        error('No date information found in sleep data');
    end

    % Initialize lighting period array
    lighting_periods = strings(total_epochs, 1);
    lighting_periods(:) = "Unknown";  % Default value

    % Classify each epoch's lighting period
    for epoch_idx = 1:total_epochs
        % Find sleep stage info for this epoch
        sleep_row = find(sleep_info.Index == epoch_idx, 1);
        if ~isempty(sleep_row)
            current_date = string(sleep_info.date(sleep_row));
            current_time = string(sleep_info.time(sleep_row));

            % Parse hour from time string (24-hour format)
            time_parts = split(current_time, ':');
            if length(time_parts) >= 1
                hour_of_day = str2double(time_parts(1));

                % Determine lighting period using DST-aware logic
                is_first_day = strcmp(current_date, start_date_str);

                if is_first_day
                    % First day: lights_off_start to lights_off_end = Lights Off
                    % lights_off_end to midnight = Lights On
                    if hour_of_day >= lights_off_start && hour_of_day < lights_off_end
                        lighting_periods(epoch_idx) = "Lights_Off";
                    else
                        lighting_periods(epoch_idx) = "Lights_On";
                    end
                else
                    % Second day: midnight to lights_on_end = Lights On
                    if hour_of_day >= 0 && hour_of_day < lights_on_end
                        lighting_periods(epoch_idx) = "Lights_On";
                    else
                        lighting_periods(epoch_idx) = "Lights_Off";  % Shouldn't happen in normal recordings
                    end
                end

                % Debug first few classifications
                if epoch_idx <= 10
                    dst_status = "";
                    if is_dst_period
                        dst_status = " (DST)";
                    else
                        dst_status = " (Standard)";
                    end
                    fprintf('  Epoch %d: %s %s (Hour %d)%s -> %s\n', ...
                            epoch_idx, current_date, current_time, hour_of_day, dst_status, lighting_periods(epoch_idx));
                end
            else
                fprintf('Warning: Could not parse time for epoch %d: %s\n', epoch_idx, current_time);
                lighting_periods(epoch_idx) = "Unknown";
            end
        else
            % No sleep info for this epoch
            lighting_periods(epoch_idx) = "Unknown";
        end
    end

    % Print lighting period summary
    lights_off_count = sum(lighting_periods == "Lights_Off");
    lights_on_count = sum(lighting_periods == "Lights_On");
    unknown_count = sum(lighting_periods == "Unknown");

    fprintf('\nLighting Period Classification Summary:\n');
    if is_dst_period
        fprintf('  Schedule used: 10AM-10PM (DST period)\n');
    else
        fprintf('  Schedule used: 9AM-9PM (Standard time)\n');
    end
    fprintf('  Lights Off epochs: %d\n', lights_off_count);
    fprintf('  Lights On epochs: %d\n', lights_on_count);
    fprintf('  Unknown epochs: %d\n', unknown_count);
    fprintf('  Total epochs: %d\n', total_epochs);

    % ======= STAGE-SPECIFIC THRESHOLDS =======
    stage_stats = struct();
    for j = 1:length(unique_stages)
        stage = unique_stages(j);
        idx = strcmp(sleep_info.StageLabel, stage);
        indices = sleep_info.Index(idx);

        % Only process indices that exist in the data
        indices = indices(indices <= total_epochs);

        if ~isempty(indices)
            epoch_means = zeros(length(indices), 1);
            p2p_values = zeros(length(indices), 1);
            % FIXED: Removed duplicate for loop
            for k = 1:length(indices)
                epoch = double(eeg1_cells{indices(k)});
                epoch_filt = filtfilt(b, a, epoch);
                epoch_means(k) = mean(abs(epoch_filt));
                p2p_values(k) = max(epoch_filt) - min(epoch_filt);
            end
            stage_stats.(stage).threshold = prctile(epoch_means, percentile_threshold);
            stage_stats.(stage).p2p_threshold = prctile(p2p_values, percentile_threshold);
        else
            % Default thresholds if no epochs for this stage
            stage_stats.(stage).threshold = inf;
            stage_stats.(stage).p2p_threshold = inf;
        end
    end

    % ======= INITIALIZE OUTPUT ARRAYS =======
    EEG1 = cell(total_epochs, 1);
    EEG2 = cell(total_epochs, 1);
    Sleep_Stage = strings(total_epochs, 1);
    Lighting_Period = strings(total_epochs, 1);  % NEW: Add lighting period array
    MSE_EEG1 = cell(total_epochs, 1);
    MSE_EEG2 = cell(total_epochs, 1);
    Epoch_Number = (1:total_epochs)';

    % Initialize with null/empty values
    for i = 1:total_epochs
        EEG1{i} = [];
        EEG2{i} = [];
        Sleep_Stage(i) = "Unknown";
        Lighting_Period(i) = lighting_periods(i);  % Set lighting period
        MSE_EEG1{i} = [];
        MSE_EEG2{i} = [];
    end

    % ======= PROCESS ALL EPOCHS IN ORDER =======
    rejected_count = 0;
    processed_count = 0;

    fprintf('\nProcessing %d epochs...\n', total_epochs);

    for epoch_idx = 1:total_epochs
        % Find sleep stage for this epoch
        sleep_row = find(sleep_info.Index == epoch_idx, 1);

        if ~isempty(sleep_row)
            current_stage = sleep_info.StageLabel(sleep_row);
            Sleep_Stage(epoch_idx) = current_stage;

            % Get thresholds for this stage
            if isfield(stage_stats, char(current_stage))
                threshold = stage_stats.(char(current_stage)).threshold;
                p2p_threshold = stage_stats.(char(current_stage)).p2p_threshold;
            else
                threshold = inf;
                p2p_threshold = inf;
            end
        else
            current_stage = "Unknown";
            Sleep_Stage(epoch_idx) = current_stage;
            threshold = inf;
            p2p_threshold = inf;
        end

        % Process EEG data
        try
            % EEG Channel 1
            epoch1 = double(eeg1_cells{epoch_idx});
            epoch1_filt = filtfilt(b, a, epoch1);

            % EEG Channel 2
            epoch2 = double(eeg2_cells{epoch_idx});
            epoch2_filt = filtfilt(b, a, epoch2);

            % Artifact detection
            epoch1_mean = mean(abs(epoch1_filt));
            epoch1_p2p = max(epoch1_filt) - min(epoch1_filt);
            epoch1_max = max(abs(epoch1_filt));

            % Check if epoch passes quality criteria
            if epoch1_mean <= threshold && epoch1_max <= 400 && epoch1_p2p <= p2p_threshold
                % Keep epoch and calculate MSE
                EEG1{epoch_idx} = epoch1_filt;
                EEG2{epoch_idx} = epoch2_filt;

                % Calculate MSE for both channels
                try
                    [mse1, ~] = get_multiple_mse_curves_matlab(epoch1_filt(:), mse_m, mse_r, mse_scales);
                    [mse2, ~] = get_multiple_mse_curves_matlab(epoch2_filt(:), mse_m, mse_r, mse_scales);
                    MSE_EEG1{epoch_idx} = mse1';  % Transpose to row vector
                    MSE_EEG2{epoch_idx} = mse2';  % Transpose to row vector
                catch mse_error
                    fprintf('MSE calculation failed for epoch %d: %s\n', epoch_idx, mse_error.message);
                    MSE_EEG1{epoch_idx} = NaN(1, length(mse_scales));
                    MSE_EEG2{epoch_idx} = NaN(1, length(mse_scales));
                end

                processed_count = processed_count + 1;
            else
                % Reject epoch - set to empty/null
                EEG1{epoch_idx} = [];
                EEG2{epoch_idx} = [];
                MSE_EEG1{epoch_idx} = [];
                MSE_EEG2{epoch_idx} = [];
                rejected_count = rejected_count + 1;
            end

        catch processing_error
            fprintf('Error processing epoch %d: %s\n', epoch_idx, processing_error.message);
            % Set to empty on error
            EEG1{epoch_idx} = [];
            EEG2{epoch_idx} = [];
            MSE_EEG1{epoch_idx} = [];
            MSE_EEG2{epoch_idx} = [];
            rejected_count = rejected_count + 1;
        end

        % Progress indicator with lighting info
        if mod(epoch_idx, 100) == 0
            fprintf('Processed %d/%d epochs (Current: %s - %s)\n', epoch_idx, total_epochs, ...
                    Sleep_Stage(epoch_idx), Lighting_Period(epoch_idx));
        end
    end

    % ======= ROBUST SAVE FOR LARGE FILES =======
    if ~exist(output_base, 'dir')
        mkdir(output_base);
    end

    [~, edf_name, ~] = fileparts(edf_path);
    output_filename = sprintf('%s_%s_dataset.mat', edf_name, condition);
    output_path = fullfile(output_base, output_filename);

    % Create the dataset structure
    dataset = struct();
    dataset.EEG1 = EEG1;
    dataset.EEG2 = EEG2;
    dataset.Sleep_Stage = Sleep_Stage;
    dataset.Lighting_Period = Lighting_Period;
    dataset.MSE_EEG1 = MSE_EEG1;
    dataset.MSE_EEG2 = MSE_EEG2;
    dataset.Epoch_Number = Epoch_Number;
    dataset.time = time;
    dataset.date = date;
    dataset.sampling_frequency = fs;
    dataset.mse_scales = mse_scales;
    dataset.processing_info = struct('processed_epochs', processed_count, ...
                                   'rejected_epochs', rejected_count, ...
                                   'total_epochs', total_epochs, ...
                                   'lights_off_epochs', lights_off_count, ...
                                   'lights_on_epochs', lights_on_count, ...
                                   'unknown_lighting_epochs', unknown_count, ...
                                   'dst_period', is_dst_period, ...
                                   'lights_off_start_hour', lights_off_start, ...
                                   'lights_off_end_hour', lights_off_end);

    % ROBUST SAVE WITH MULTIPLE ATTEMPTS AND VERIFICATION
    fprintf('\nSaving dataset...\n');
    save_successful = false;
    max_attempts = 3;
    
    % Estimate file size
    s = whos('dataset');
    estimated_size_mb = s.bytes / (1024^2);
    fprintf('Estimated dataset size: %.1f MB\n', estimated_size_mb);
    
    for attempt = 1:max_attempts
        fprintf('Save attempt %d/%d...\n', attempt, max_attempts);
        
        try
            % For large files (>2GB), use -v7.3, otherwise use -v7
            if estimated_size_mb > 1000  % > 1GB, use -v7.3
                fprintf('Using -v7.3 format for large file...\n');
                save(output_path, 'dataset', '-v7.3');
            else
                fprintf('Using -v7 format...\n');
                save(output_path, 'dataset', '-v7');
            end
            
            % Verify the save worked by trying to load it
            fprintf('Verifying saved file...\n');
            test_load = load(output_path, 'dataset');
            
            % Check if the loaded data structure is valid
            if isstruct(test_load.dataset) && ...
               isfield(test_load.dataset, 'EEG1') && ...
               isfield(test_load.dataset, 'MSE_EEG1') && ...
               length(test_load.dataset.EEG1) == total_epochs
                
                fprintf('File saved and verified successfully!\n');
                save_successful = true;
                break;
            else
                error('Loaded dataset structure is invalid');
            end
            
        catch save_error
            fprintf('Save attempt %d failed: %s\n', attempt, save_error.message);
            
            % Clean up failed file
            if exist(output_path, 'file')
                delete(output_path);
            end
            
            if attempt < max_attempts
                fprintf('Retrying with different method...\n');
                pause(2); % Brief pause before retry
            end
        end
    end
    
    if ~save_successful
        error('All save attempts failed. Check disk space and permissions on: %s', output_path);
    end

    fprintf('\n======= PROCESSING COMPLETE =======\n');
    fprintf('Total epochs: %d\n', total_epochs);
    fprintf('Processed epochs: %d\n', processed_count);
    fprintf('Rejected epochs: %d\n', rejected_count);
    fprintf('Lights Off epochs: %d\n', lights_off_count);
    fprintf('Lights On epochs: %d\n', lights_on_count);
    fprintf('Unknown lighting epochs: %d\n', unknown_count);
    if is_dst_period
        fprintf('DST schedule applied: %d:00-%d:00 Lights Off\n', lights_off_start, lights_off_end);
    else
        fprintf('Standard schedule applied: %d:00-%d:00 Lights Off\n', lights_off_start, lights_off_end);
    end
    fprintf('Dataset saved to: %s\n', output_path);
end

% ======= MSE CALCULATION FUNCTION =======
function [mse, scales]= get_multiple_mse_curves_matlab(data,m,r,scales,override)
% $$ Written by Natasha Kovacevic - based on M. Costa's C code - simplified - skipped all sorts of command line options  - here we simultaneously get mse curves across trials or voxels
% $$ Usage:   [mse, scales]= get_multiple_mse_curves_matlab(data,m,r,scales)
% $$ Inputs:
% $$   data = time sereis data in (time variable) format, e.g (time trial) or (time voxel)
% $$   m = pattern length (default 2)
% $$   r = similarity criterion (defualt 0.5)
% $$   scales = vector of scales (default [1:(floor(size(data,1)/50))])
% $$   override = 0 or 1 (default is 0) - in some rare cases we want to oveeride this min50-rule
% $$ Outputs:
% $$   mse = entropy in (scale variable) format
% $$   scales = return vector of scales actually used by the program

    num_tpts = size(data,1);
    num_vars = size(data,2);
    if ~exist('m','var'), m=2; end
    if ~exist('r','var'), r=0.5; end
    if ~exist('scales','var'), scales = [1:(floor(num_tpts/50))] ; end

    % make sure that scales are positive integers
    scales = sort(unique(scales)); % sort scales in increasing order
    if sum(scales == round(scales)) ~= numel(scales)
      error('scales vector must contain positive integers');
    end

    if exist('override','var')
      if ~override
        % make sure that scales are  < (num_tpts/50)
        scales = sort(intersect(scales,[1:(floor(size(data,1)/50))]));
      end
    end

    % mormalize data and r
    sd = std(data,[],1);
    r_new = r * sd;

    mse = zeros(numel(scales),num_vars); %(scale var)
    for s=1:numel(scales)
      sc = scales(s);

      % coarse grind time series at this scale
      num_cg_tpts = floor(num_tpts/sc);
      y = zeros(num_cg_tpts, num_vars);
      for t = 1:num_cg_tpts
        y(t,:) = mean(data((t-1)*sc + [1:sc],:),1);
      end

      % calculate sample entropy of coarse ground time series y
      nlin_sc = num_cg_tpts - m;
      cont = zeros(m+1,num_vars);
      for var = 1:num_vars
        for i = 1:nlin_sc
          for l = (i+1):nlin_sc % self-matches are not counted
            k = 0;
            while ((k < m) & (abs(y(i+k,var) - y(l+k,var)) <= r_new(var)))
              k = k + 1;
              cont(k,var) = cont(k,var) + 1;
            end
            if ((k == m) & (abs(y(i+m,var) - y(l+m,var)) <= r_new(var)))
              cont(m+1,var) = cont(m+1,var) + 1;
            end
          end
        end
      end

      % calculate mse at this scale
      for var = 1:num_vars
        if (cont(m+1,var) == 0 | cont(m,var) == 0)
          mse(s,var) = -log(1/((nlin_sc)*(nlin_sc -1)));
        else
          mse(s,var) = -log(cont(m+1,var)/cont(m,var));
        end
      end

    end % for s=1:numel(scales)
end
