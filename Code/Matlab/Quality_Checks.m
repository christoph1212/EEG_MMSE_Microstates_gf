%% Quality Assessment

dir_Root = "E:/Complexity/";
dir_Log = "E:/Complexity/Data/Log/";

%% Channel Difference between Dreszer et al. and us

UJ = load("C:\Users\Christoph Frühlinger\Downloads\UJ.mat");
original = UJ.UJ(1);
clear UJ

eeglab nogui
EEG = pop_loadset('filename','first_run_eyes_open_sub-AA06WI11_task-Resting_run-1_eeg.set', ...
    'filepath','C:\\Users\\Christoph Frühlinger\\Nextcloud\\PhD\\Forschung\\Complexity\\Data\\Preprocessed_updated\\');

chan_diff1_idx = ~ismember(upper({original.chanlocs.labels}), upper({EEG.chanlocs.labels}));
chan_diff1 = {original.chanlocs(chan_diff1_idx).labels};
disp('We miss channels:')
fprintf('%s\n', chan_diff1{:})

chan_diff2_idx = ~ismember(upper({EEG.chanlocs.labels}), upper({original.chanlocs.labels}));
chan_diff2 = {EEG.chanlocs(chan_diff2_idx).labels};
disp('We have extra channels:')
fprintf('%s\n', chan_diff2{:})

%% MMSE Difference between old and updated Preproc Pipeline

% old = readmatrix('..\..\Data\MMSEData\MMSE_C_26_seconds_first_run_eyes_open_sub-AA06WI11_task-Resting_run-1_eeg.set.csv');
% old = old(1,2:end);
% 
% new = readmatrix('..\..\Data\MMSEData_updated\MMSE_C_26_seconds_first_run_eyes_open_sub-AA06WI11_task-Resting_run-1_eeg.set.csv');
% new = new(1,2:end);
% 
% figure()
% plot(old)
% hold on
% plot(new)
% hold off
% legend(["old", "new"]);

%% Check Error-Files

% List of Error Files
Preproc_Log = "E:\Complexity\Data\Log\Preproc";
error_files = dir(fullfile(Preproc_Log, 'Error_PreProc_*.txt'));

% empty string array for messages
messages = strings(0);

% Loop over files
for error_file = 1:length(error_files)
    
    % read lines
    lines = readlines(fullfile(Preproc_Log, error_files(error_file).name));

    % Loop over lines and collect relevant lines
    for i = 1:length(lines)
        % remove leading or tailing white-spaces
        line = strtrim(lines(i));
        % Actual error message starts two lines after this message
        if startsWith(line, 'The returned Error Message is:')            
            messages(end+1) = strtrim(lines(i + 2));
        end
    end
end

% Count number of unique errors
[unique_errors, ~, idx] = unique(messages);
counts = accumarray(idx, 1);

% Sort errors
[sorted_counts, sort_idx] = sort(counts, 'descend');
sorted_errors = unique_errors(sort_idx);

% Do you also want to save the errors as a file? (0|1)
save_errors = 1;

% Print our errors in command window
fprintf('Error Messages:\n');
for i = 1:length(sorted_errors)
    fprintf('%3d × %s\n', sorted_counts(i), sorted_errors(i));
end

fprintf('Files with not enough Data for Preprocessing:\n')
fprintf('%3d x in the first run\n', length(dir(fullfile(Preproc_Log, 'Error_NotEnoughData_for_Preproc_first*.txt'))))
fprintf('%3d x in the second run\n', length(dir(fullfile(Preproc_Log, 'Error_NotEnoughData_for_Preproc_second*.txt'))))
fprintf('%3d x in the third run\n', length(dir(fullfile(Preproc_Log, 'Error_NotEnoughData_for_Preproc_third*.txt'))))

% Create table to save as csv File in Log Folder
if save_errors == 1
    % Add errors and counts to table
    Error = strrep(sorted_errors(:), ',', ';');
    err_table = table(Error(:), num2cell(sorted_counts(:)), 'VariableNames', {'Error', 'Count'});
    % Add 'NotEnoughData' Errors
    err_table(end+1,:) = {"Error_NotEnoughData_for_Preproc_first_run", {length(dir(fullfile(Preproc_Log, 'Error_NotEnoughData_for_Preproc_first*.txt')))}};
    err_table(end+1,:) = {"Error_NotEnoughData_for_Preproc_second_run", {length(dir(fullfile(Preproc_Log, 'Error_NotEnoughData_for_Preproc_second*.txt')))}};
    err_table(end+1,:) = {"Error_NotEnoughData_for_Preproc_third_run", {length(dir(fullfile(Preproc_Log, 'Error_NotEnoughData_for_Preproc_third*.txt')))}};
    % Sort errors according to counts
    err_table = sortrows(err_table, 'Count', 'descend');
    % Save File
    ErrorFilename = strcat(Preproc_Log, '\Errors.csv');
    writetable(err_table, ErrorFilename);
end


%% Descriptives after Snipping

% Set path and load files
log_path = fullfile(dir_Log, 'Preproc');
log_files = dir(fullfile(log_path, 'Log_*.csv'));

% Initialize empty table
all_logs = table();

for i = 1:length(log_files)
    file = fullfile(log_path, log_files(i).name);
    T = readtable(file);
    all_logs = [all_logs; T];
end

% Adapt Table
all_logs.Properties.VariableNames{'FileName'} = 'ID';
splits = split(string(all_logs.ID), '_');
all_logs.ID = splits(:, 1);

all_logs.Condition = string(all_logs.Condition);

% Save Table as csv File
LogFilename = strcat(Preproc_Log, '\All_Logs.csv');
writetable(all_logs, LogFilename);

% Check complete Participants
[unique_files, ~, idx] = unique(all_logs.ID);
counts = accumarray(idx, 1);
complete_IDs = unique_files(counts == 6);
n_complete_IDs = length(complete_IDs);
fprintf("%s Participants have complete data sets\n", num2str(n_complete_IDs));

% Check Files < 40 seconds Epochs
short_epochs = sum(all_logs.SnippletLength < 40) / height(all_logs) * 100;
fprintf("%.2f %% of the files have epochs shorter than 40s\n", short_epochs);

% Check complete Participants with 40 seconds Epoch length
epoch_40s = all_logs(all_logs.SnippletLength == 40, :);
[unique_files, ~, idx] = unique(epoch_40s.ID);
counts = accumarray(idx, 1);
complete_IDs_40s = unique_files(counts == 6);
n_complete_IDs_40s = length(complete_IDs_40s);
fprintf("%s Participants have complete data sets with 40s epochs\n", num2str(n_complete_IDs_40s));
epoch_40s_complete = epoch_40s(ismember(epoch_40s.ID, complete_IDs_40s),:);

% Other descriptives
summary(all_logs)
summary(epoch_40s)
summary(epoch_40s_complete)

figure()
subplot(2,3,1)
boxplot(all_logs.BadICs)
title_text = sprintf("Bad ICs\nMean = %.2f\nSD = %.2f", ...
    mean(all_logs.BadICs), std(all_logs.BadICs));
title(title_text)

subplot(2,3,2)
boxplot(all_logs.InterpolatedChannels)
title_text = sprintf("Interpolated Channels\nMean = %.2f\nSD = %.2f", ...
    mean(all_logs.InterpolatedChannels), std(all_logs.InterpolatedChannels));
title(title_text)

subplot(2,3,3)
boxplot(all_logs.SnippletLength)
title_text = sprintf("Epoch Length\nMean = %.2f\nSD = %.2f", ...
    mean(all_logs.SnippletLength), std(all_logs.SnippletLength));
title(title_text)

subplot(2,3,4)
histogram(all_logs.BadICs)

subplot(2,3,5)
histogram(all_logs.InterpolatedChannels)

subplot(2,3,6)
histogram(all_logs.SnippletLength)

%% Check Raw File Folders

% Get Directories
dir_Raw  = strcat(dir_Root, "Data/RawData/task-Resting");

subdir = dir(dir_Raw);
subdir = subdir([subdir.isdir] & ~ismember({subdir.name}, {'.', '..'}));

for i = 1:length(subdir)
    subID = subdir(i).name;
    eegDir = fullfile(dir_Raw, subID, 'eeg');

    if isfolder(eegDir)
        files = dir(fullfile(eegDir, '*'));
        files = files(~[files.isdir]);

        numFiles = numel(files);
        
        % 3 Sessions with 4 Files each
        if numFiles < 12
            fprintf('%s: %d Files\n', subID, numFiles);
        end
    else
        fprintf('%s: eeg folder not found\n', subID);
    end
end

%% Check Channel Order
% Get epoched data
data_path = "E:\Complexity\Data\Snipplet";
files = dir(fullfile(data_path, '*.set'));

template_labels = [];
mismatched_files = {};

for i = 1:length(files)
    filename = fullfile(files(i).folder, files(i).name);
    evalc("EEG = pop_loadset('filename', files(i).name, 'filepath', files(i).folder);");

    % Get channel labels
    current_labels = {EEG.chanlocs.labels};

    % Save first as template
    if isempty(template_labels)
        template_labels = current_labels;
        continue
    end

    % Compare
    if isequal(template_labels, current_labels)
        % fprintf('Datensatz %s: OK\n', files(i).name);
    else
        fprintf('File %s: MISMATCH\n', files(i).name);
        mismatched_files{end+1} = files(i).name;
    end
end

if isempty(mismatched_files)
    disp('All Files are sorted equally.');
else
    disp('These Files have mismatched Channels:');
    disp(mismatched_files');
end

%% MMSE Vectors

% Set path and load files
vector_path = strcat(dir_Root, "Data/MMSEData/");
vector_files = dir(fullfile(vector_path, 'MMSE_*.csv'));

% Initialize empty table
all_vectors = table();

% Loop through vector files
for i = 1:length(vector_files)
    file = fullfile(vector_path, vector_files(i).name);
    T = readtable(file);
    splits = strsplit(T.name{1}, '_');
    T.ID = string(splits{9});
    T.Condition = string([splits{5}, '_', splits{6}, '_', splits{7}, '_', splits{8}]);
    T.Set = string(splits{2});
    T.Length = string(splits{3});
    all_vectors = [all_vectors; T(:, [14:17, 2:13])];
end

% make certain columns categorical
all_vectors.ID = categorical(all_vectors.ID);
all_vectors.Condition = categorical(all_vectors.Condition);
all_vectors.Set = categorical(all_vectors.Set);
all_vectors.Length = double(all_vectors.Length);

% Save Table as csv File
VectorFilename = strcat(dir_Log, 'MMSE', '/All_Vectors.csv');
writetable(all_vectors, VectorFilename);

% exclude rows with NaN or Inf
mask = all(isfinite(all_vectors{:, 5:end}), 2);
all_vectors_clean = all_vectors(mask, :);

fprintf("%d vector files removed due to NaN or Inf\n", height(all_vectors) - height(all_vectors_clean));

summary(all_vectors_clean)

% Plotting
conditions = categories(all_vectors.Condition);
sets = categories(all_vectors.Set);
plotidx = 1;

figure()
x = 1:12;
mmse_cols = {'mmse_1', 'mmse_2', 'mmse_3', 'mmse_4', 'mmse_5', 'mmse_6', ...
             'mmse_7', 'mmse_8', 'mmse_9', 'mmse_10', 'mmse_11', 'mmse_12'};

for i_row = 1:numel(conditions)

    for i_col = 1:numel(sets)
        
        % Get subset of data
        subset = all_vectors_clean(all_vectors_clean.Condition == conditions{i_row} ...
            & all_vectors_clean.Set == sets{i_col}, :);

        subplot(3, 2, plotidx);
        hold on

        if ~isempty(subset)
            % Calculate mean of MMSE Vector
            mean_mmse = mean(subset{:, mmse_cols}, 1, 'omitnan');
            plot(x, mean_mmse, '.-', 'DisplayName',char(sets(i_col)))
            title_string = strsplit(conditions{i_row}, '_');
            title_string = [title_string{1}, ' ', title_string{2}, ' ', title_string{3}, ' ', title_string{4}];
            title(sprintf("%s (\\itN\\rm = %d)", title_string, height(subset)));
        end

    end

    lgd = legend('show', 'Location', 'eastoutside');
    title(lgd,'Channel Set');
    xlabel('Scaling Vector');
    ylabel('MMSE');
    xticks(x);
    grid on
    hold off
    plotidx = plotidx + 1;

end
