function preprocessing_with_mse(edf_path, tsv_path, output_base)
    % ======= CONFIGURATION =======
    fs = 400;
    [b, a] = butter(2, [0.5 50] / (fs / 2), 'bandpass');
    percentile_threshold = 90;
    unique_stages = ["Wake", "NREM", "REM"];

    % MSE parameters
    mse_scales = 1:50;  % Scale factors for MSE
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
    MSE_EEG1 = cell(total_epochs, 1);
    MSE_EEG2 = cell(total_epochs, 1);
    Epoch_Number = (1:total_epochs)';

    % Initialize with null/empty values
    for i = 1:total_epochs
        EEG1{i} = [];
        EEG2{i} = [];
        Sleep_Stage(i) = "Unknown";
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

        % Progress indicator
        if mod(epoch_idx, 100) == 0
            fprintf('Processed %d/%d epochs (Current: %s)\n', epoch_idx, total_epochs, Sleep_Stage(epoch_idx));
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
    dataset.MSE_EEG1 = MSE_EEG1;
    dataset.MSE_EEG2 = MSE_EEG2;
    dataset.Epoch_Number = Epoch_Number;
    dataset.time = time;
    dataset.date = date;
    dataset.sampling_frequency = fs;
    dataset.mse_scales = mse_scales;
    dataset.processing_info = struct('processed_epochs', processed_count, ...
                                   'rejected_epochs', rejected_count, ...
                                   'total_epochs', total_epochs);

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
