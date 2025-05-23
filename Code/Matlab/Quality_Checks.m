%% Quality Assessment

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

old = readmatrix('..\..\Data\MMSEData\MMSE_C_26_seconds_first_run_eyes_open_sub-AA06WI11_task-Resting_run-1_eeg.set.csv');
old = old(1,2:end);

new = readmatrix('..\..\Data\MMSEData_updated\MMSE_C_26_seconds_first_run_eyes_open_sub-AA06WI11_task-Resting_run-1_eeg.set.csv');
new = new(1,2:end);

figure()
plot(old)
hold on
plot(new)
hold off
legend(["old", "new"]);

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
    err_table = table(Error(:), num2cell(sorted_counts(:)), ...
        'VariableNames', {'Error', 'Count'});
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

%% Check Raw File Folders

dir_Root = "E:/Complexity/";                                % path to project
dir_Raw  = strcat(dir_Root, "Data/RawData/task-Resting");   % path to raw data

subdir = dir(dir_Raw);
subdir = subdir([subdir.isdir] & ~ismember({subdir.name}, {'.', '..'}));

for i = 1:length(subdir)
    subID = subdir(i).name;
    eegDir = fullfile(dir_Raw, subID, 'eeg');

    if isfolder(eegDir)
        files = dir(fullfile(eegDir, '*'));
        % Ignoriere '.' und '..'
        files = files(~[files.isdir]);

        numFiles = numel(files);
        if numFiles < 12
            fprintf('%s: %d Files\n', subID, numFiles);
        end
    else
        fprintf('%s: eeg folder not found\n', subID);
    end
end