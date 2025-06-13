function combine_data(dir_Root, dir_Log)
%% Combine MMSE Features and Vectors and add demographics and gf
% Outputs:
%   Creates a .csv file of the MMSE features, MMSE vectors, demographic and
%   gf data for futher analysis
%
% Inputs:
%   dir_Root:   String pointing to the project's parent folder
%   dir_Log:    String pointing to Log-Files directory
%
% This script was created by: Christoph Fruehlinger (June 2025)
%
% This script was adapted by: Christoph Fruehlinger (June 2025)

%% Load and create relevant files

fprintf("******************\nStarting Combining\n******************\n")

% MMSE Features
FeatureFilename = strcat(dir_Root, 'Data/MMSEFeatures/All_Features_Diff.csv');
features = readtable(FeatureFilename);

features.ID = categorical(features.ID);
features.Condition = categorical(features.Condition);
features.Set = categorical(features.Set);
features.Length = double(features.Length);

% MMSE Vectors
VectorFilename = strcat(dir_Log, 'MMSE', '/All_Vectors.csv');

if isfile(VectorFilename)

    all_vectors = readtable(VectorFilename);

else
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

    % Save Table as csv File
    writetable(all_vectors, VectorFilename);
end

% make certain columns categorical
all_vectors.ID = categorical(all_vectors.ID);
all_vectors.Condition = categorical(all_vectors.Condition);
all_vectors.Set = categorical(all_vectors.Set);
all_vectors.Length = double(all_vectors.Length);

% exclude rows with NaN or Inf
mask = all(isfinite(all_vectors{:, 5:end}), 2);
all_vectors_clean = all_vectors(mask, :);

mmse_data = outerjoin(features, all_vectors_clean, ...
    'Keys', {'ID', 'Condition', 'Set', 'Length'}, ...
    'MergeKeys', true);

%% Get fluid scores
fluid_path = strcat(dir_Root, 'Data/FluidData');

Participants = dir(fullfile(fluid_path, 'task-IST'));
Participants = {Participants(~ismember({Participants.name}, {'.', '..'})).name};

% Initiate Table
gf_data = table();
gf_data.ID = Participants';

% Loop through Fluid Intelligence Data
fluid_cor = readtable(fullfile(fluid_path, "IST_fluid_A.xlsx"));
fluid_cor = table2array(fluid_cor(:,2));

for i_sub = 1:length(Participants)

    file = dir(fullfile(fluid_path, 'task-IST', Participants{i_sub}, 'beh', '*Fluid_beh.csv'));
    
    if isempty(file)
        gf_data.gf_score(i_sub) = NaN;
        continue
    end
    
    opts = detectImportOptions(fullfile(file.folder, file.name));
    opts.SelectedVariableNames = "ratings";
    gf = readtable(fullfile(file.folder, file.name), opts);
    gf = table2array(gf(2:end-1, :));

    if isempty(gf)
        gf_data.gf_score(i_sub) = NaN;
        continue
    end

    gf_data.gf_score(i_sub) = sum(strcmp(fluid_cor, gf));

end

% Standardize gf scores
gf_data.gf_score = zscore(gf_data.gf_score);

%% Add sex to table
demographics_path = strcat(dir_Root, 'Data/SocioDemographics.txt');
demographics_data = readtable(demographics_path);
demographics_data = demographics_data(:,1:3);
demographics_data.Gender = string(demographics_data.Gender);
demographics_data.Gender(strcmp(demographics_data.Gender, "1")) = "female";
demographics_data.Gender(strcmp(demographics_data.Gender, "2")) = "male";
demographics_data.Gender = categorical(demographics_data.Gender);

%% Combine all tables
dem_gf_data = innerjoin(gf_data, demographics_data, "Keys", "ID");
dem_gf_data.ID = categorical(dem_gf_data.ID);
data = outerjoin(dem_gf_data,mmse_data, "Keys", "ID", "MergeKeys", true);

%% Save File
SaveFileName = strcat(dir_Root, 'Data/MMSEData/MMSE_data_full.csv');
writetable(data, SaveFileName);

fprintf("\nSaved file to: %s\n\n", SaveFileName)

fprintf("*******************\nFinisched Combining\n*******************\n")