%% Quality Assessment

dir_Root    = "E:/Complexity/";
dir_Log     = "E:/Complexity/Data/Log/";
Preproc_Log = fullfile(dir_Log, 'Preproc');

%% Check Error-Files
% Do you also want to save the errors as a file? (0|1)
save_errors = 1;

% List of Error Files
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
log_files = dir(fullfile(Preproc_Log, 'Log_*.csv'));

% Initialize empty table
all_logs = table();

for i = 1:length(log_files)
    file = fullfile(Preproc_Log, log_files(i).name);
    T = readtable(file);
    all_logs = [all_logs; T];
end

% Adapt Table
all_logs.Properties.VariableNames{'FileName'} = 'ID';
splits = split(string(all_logs.ID), '_');
all_logs.ID = splits(:, 1);
all_logs.Condition = string(all_logs.Condition);

% make certain columns categorical
all_logs.ID = categorical(all_logs.ID);
all_logs.Condition = categorical(all_logs.Condition);

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

% Plot Distributions
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

% Initialize empty lists
template_labels = [];
mismatched_files = {};

for i = 1:length(files)
    % Load file
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

% Print out Results
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

% Check Rank Correlation between conditions

% define conditions and sets
cond1 = "first_run_eyes_open";

otherconds = ["first_run_eyes_closed", "second_run_eyes_open", "third_run_eyes_open"];
sets = ["C", "F", "FL", "FR", "ML", "MR", "P", "PL", "PR"];

% loop over sets and condition to compare
for i_cond = 1:numel(otherconds)
    for i_set = 1:numel(sets)
        cond2 = otherconds(i_cond);
        targetSet = sets(i_set);
        
        % get MMSE data
        mmseVars = startsWith(all_vectors_clean.Properties.VariableNames, 'mmse_');
        mmseCols = all_vectors_clean.Properties.VariableNames(mmseVars);
        
        % Filter for condition and set
        T1 = all_vectors_clean(all_vectors_clean.Condition == cond1 & all_vectors_clean.Set == targetSet, :);
        T2 = all_vectors_clean(all_vectors_clean.Condition == cond2 & all_vectors_clean.Set == targetSet, :);
        
        % get common IDs and filter
        commonIDs = intersect(T1.ID, T2.ID);
        
        T1 = T1(ismember(T1.ID, commonIDs), :);
        T2 = T2(ismember(T2.ID, commonIDs), :);
        
        % sort rows according to ID
        T1 = sortrows(T1, 'ID');
        T2 = sortrows(T2, 'ID');
        
        % get condition-specific MMSE data
        X1 = T1{:, mmseVars};  % Condition 1
        X2 = T2{:, mmseVars};  % Condition 2
        
        % Überprüfen
        assert(isequal(T1.ID, T2.ID), 'IDs stimmen nicht überein!');
        
        % calculate Spearman correlation and save coefficient and p-value
        numMMSE = size(X1, 2);
        rhoValues = zeros(1, numMMSE);
        pValues  = zeros(1, numMMSE);
        
        for i = 1:numMMSE
            [r, p] = corr(X1(:, i), X2(:, i), 'Type', 'Spearman');
            rhoValues(i) = r;
            pValues(i)  = p;
        end
        
        % save in table and print
        resultTable = table(repelem(targetSet, numMMSE)', repelem(cond1, numMMSE)', repelem(cond2, numMMSE)', mmseCols', rhoValues', pValues', ...
            'VariableNames', {'Set', 'Condition1', 'Condition2', 'MMSE_Feature', 'SpearmanRho', 'pValue'});
        
        disp(resultTable);
        
        % Scatterplots
        mmseNames = all_vectors_clean.Properties.VariableNames(mmseVars);
        figure('Name',sprintf('%s, %s', cond2, targetSet));
        for i = 1:numMMSE
            x = T1{:, mmseNames{i}};
            y = T2{:, mmseNames{i}};
            rho = corr(x, y, 'Type', 'Spearman');
            subplot(3, 4, i);  % 3x4 Layout für 12 Plots
            scatter(x, y, 30, 'filled');
            lsline;
            title(sprintf('%s\nρ = %.2f', mmseNames{i}, rho), 'Interpreter', 'none');
            xlabel(cond1, 'Interpreter', 'none');
            ylabel(cond2, 'Interpreter', 'none');
        end
        sgtitle(sprintf('Retest-Reliability for each MMSE Vector (Set %s)', targetSet));
    end
end

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

%% MMSE Features

% Set path and load files
feature_path = strcat(dir_Root, "Data/MMSEFeatures/");
feature_files = dir(fullfile(feature_path, 'features_*.csv'));

% Initialize empty table
all_features = table();

% Loop through vector files
for i = 1:length(feature_files)
    fprintf('%d/%d\n', i, length(feature_files))
    file = fullfile(feature_path, feature_files(i).name);
    T = readtable(file, "Delimiter", ';');
    splits = strsplit(T.name{1}, '_');
    T.ID = string(splits{9});
    T.Condition = string([splits{5}, '_', splits{6}, '_', splits{7}, '_', splits{8}]);
    T.Set = string(splits{2});
    T.Length = string(splits{3});
    all_features = [all_features; T(:, [5:8, 2:4])];
end

% make certain columns categorical
all_features.ID = categorical(all_features.ID);
all_features.Condition = categorical(all_features.Condition);
all_features.Set = categorical(all_features.Set);
all_features.Length = double(all_features.Length);

% Save Table as csv File
FeatureFilename = strcat(dir_Log, 'MMSE', '/All_Features.csv');
writetable(all_features, FeatureFilename);

% exclude rows with NaN or Inf
mask = all(isfinite(all_features{:, 5:end}), 2);
all_features_clean = all_features(mask, :);

fprintf("%d vector files removed due to NaN or Inf\n", height(all_features) - height(all_features_clean));

summary(all_features_clean)

% Plot each Feature in separate Figure
conditions = categories(all_features_clean.Condition);
sets = categories(all_features_clean.Set);
Features = ["AUC", "Max-Slope", "Avg-Entropy"];

for i_feat = 1:numel(Features)
    
    figure(i_feat)

    % Initialize empty lists
    data = [];
    group_labels = {};
    plotidx = 1;

    for i_row = 1:numel(conditions)
        
        condition = conditions{i_row};
        subplot(3, 2, plotidx);
        hold on
        
        % Initialize empty lists per condition
        data_condition = [];
        group_labels_condition = {};
        
        for i_col = 1:numel(sets)
           
            % Get subset of data
            subset = all_features_clean(all_features_clean.Condition == condition & ...
                all_features_clean.Set == sets{i_col}, :);
            
            if ~isempty(subset)
                col_data = table2array(subset(:,i_feat+4));
                n = height(subset);
                data_condition = [data_condition; col_data];
                
                % get Labels
                group_label = repmat({char(sets{i_col})}, n, 1);
                group_labels_condition = [group_labels_condition; group_label];
            end
        end
        
        % Plot Boxplot
        if ~isempty(data_condition)
            boxplot(data_condition, group_labels_condition)
        end
                
        title_str = strrep(char(condition), '_', ' ');
        title(sprintf("%s", title_str))
        xlabel('Channel Set')
        ylabel(Features(i_feat))
        grid on
        
        hold off
        plotidx = plotidx + 1;
    end
end
