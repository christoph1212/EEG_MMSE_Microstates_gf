function MMSE_feature_diff(dir_Root, dir_Log)
%% Calculation of MMSE features difference scores
% Outputs:
%   Creates a .csv file of the MMSE features and difference scores for the
%   following channels sets: F-P, FL-PL, FR-PR, FL-FR, PL-PR, and ML-MR
%
% Inputs:
%   dir_Root:   String pointing to the project's parent folder
%   dir_Log:    String pointing to Log-Files directory
%
% This script was created by: Christoph Fruehlinger (June 2025)
%
% This script was adapted by: Christoph Fruehlinger (June 2025)

%% Load or create file

fprintf("*******************************\nStarting Difference Calculation\n*******************************\n")

FeatureFilename = strcat(dir_Log, 'MMSE', '/All_Features.csv');

if isfile(FeatureFilename)

    fprintf("File already exists.\nLoading...\n")

    all_features = readtable(FeatureFilename);

else

    fprintf("Creating file. Please wait...\n")

    % Set path and load files
    feature_path = strcat(dir_Root, "Data/MMSEFeatures/");
    feature_files = dir(fullfile(feature_path, 'features_*.csv'));
    
    % Initialize empty table
    all_features = table();
    
    % Loop through vector files
    for i = 1:length(feature_files)
        fprintf('Loading file %d/%d\n', i, length(feature_files))
        file = fullfile(feature_path, feature_files(i).name);
        T = readtable(file, "Delimiter", ';');
        splits = strsplit(T.name{1}, '_');
        T.ID = string(splits{9});
        T.Condition = string([splits{5}, '_', splits{6}, '_', splits{7}, '_', splits{8}]);
        T.Set = string(splits{2});
        T.Length = string(splits{3});
        all_features = [all_features; T(:, [5:8, 2:4])];
    end

end

% make certain columns categorical
all_features.ID = categorical(all_features.ID);
all_features.Condition = categorical(all_features.Condition);
all_features.Set = categorical(all_features.Set);
all_features.Length = double(all_features.Length);

% exclude rows with NaN or Inf
mask = all(isfinite(all_features{:, 5:end}), 2);
all_features_clean = all_features(mask, :);

% Calculate Difference Scores
conditions = unique(all_features_clean.Condition);
sets = {["F", "P"], ["FL", "PL"], ["FR", "PR"], ["FL", "FR"], ["PL", "PR"], ["ML", "MR"]};

fprintf("Calculating Difference Scores")

for i_cond = 1:numel(conditions)

    for i_sets = 1:numel(sets)

        % Get Condition and Sets 
        cond = conditions(i_cond);
        set1 = sets{1, i_sets}(1);
        set2 = sets{1, i_sets}(2);

        % Filter
        filtered_1 = all_features_clean(all_features_clean.Condition == cond & all_features_clean.Set == set1,:);
        filtered_2 = all_features_clean(all_features_clean.Condition == cond & all_features_clean.Set == set2,:);

        % Check if ID order is equal
        if ~isequal(filtered_1.ID, filtered_2.ID)

            commonIDs = intersect(filtered_1.ID, filtered_2.ID);
        
            filtered_1 = filtered_1(ismember(filtered_1.ID, commonIDs), :);
            filtered_2 = filtered_2(ismember(filtered_2.ID, commonIDs), :);
        
            % sort rows according to ID
            filtered_1 = sortrows(filtered_1, 'ID');
            filtered_2 = sortrows(filtered_2, 'ID');

        end

        % Calculate Difference and create new table
        Diffname = strcat(set1, "-", set2);
        Diff = table(repelem(Diffname, height(filtered_1))', 'VariableNames', "Set");
	    Diff_table = [filtered_1(:,1:2), Diff, filtered_1(:,4), filtered_1(:,5:7) - filtered_2(:,5:7)];
        Diff_table.Set = categorical(Diff_table.Set);

        % Append Difference Table to feature table
        all_features_clean = [all_features_clean; Diff_table];
    end
end

% Save File
SaveFileName = strcat(dir_Root, 'Data/MMSEFeatures/All_Features_Diff.csv');
writetable(all_features_clean, SaveFileName);

fprintf("\nSaved file to: %s\n\n", SaveFileName)

fprintf("*******************************\nFinished Difference Calculation\n*******************************\n")