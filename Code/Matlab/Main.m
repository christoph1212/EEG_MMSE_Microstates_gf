%% Main Analysis Script
%
% This script is the main file to run:
% 1. Preprocessing
% 2. Epoch Extraction
% 3. MMSE Complexity Analysis
%
% Important Note: We adapted the ADJUST.m file to run without a GUI and
% creation of an output file. Therefore, please make sure that you cloned
% the GitHub repository [LINK]. Adaptions are marked in the function file.
%
% Make sure to set the folders according to your workspace. The scripts
% were used in the study "Spatiotemporal complexity patterns of 
% resting-state  bioelectrical activity explain fluid intelligence: Sex
% matters" by Dreszer et al. (2020) and are openly available on GitHub
% [https://github.com/IS-UMK/complexity/tree/master]. They were adapted to
% the CoScience data set by ...
%
% Last modified: April 2025 by Christoph Fruehlinger

%% Housekeeping
clear
clc
close

%% Preprocessing

% Setup
dir_Root    = "E:/Complexity/";                   % path to project
dir_Raw     = strcat(dir_Root, "Data/RawData/");  % path to raw data
dir_Log     = strcat(dir_Root, "Data/Log/");      % path where log data should be stored
Overwrite   = 1;                                  % 0 = no; 1 = yes

% Start EEGLAB
dir_eeglab  = "E:/Complexity/Code/Matlab/eeglab2025.0.0";
addpath(dir_eeglab);
eeglab nogui

% Check for necessary Plugins
if ~ismember('trimOutlier', {PLUGINLIST.plugin})
    error("Error: Plugin 'trimOutlier' is not installed. Please install or " + ...
        "ideally clone 'eeglab2025.0.0' folder from GitHub repository.")
end

if ~ismember('Adjust', {PLUGINLIST.plugin})
    error("Error: Plugin 'Adjust' is not installed. Please install or " + ...
        "ideally clone 'eeglab2025.0.0' folder from GitHub repository.")
end

% Actual Preprocessing
Preproc(dir_Raw, dir_Root, dir_Log, Overwrite)

%% Epoch Extraction

Snipplet(dir_Root, dir_Log, Overwrite)

%% MMSE Complexity Analysis

% MMSE Vectors
MMSE_silent(dir_Root, dir_Log, Overwrite)

% MMSE Feature Extraction
dir_MMSE = strcat(dir_Root, "Data/MMSEData/");
dir_feat = strcat(dir_Root, "Data/MMSEFeatures/");

if ~isfolder(dir_MMSE)
    error("Folder 'Data/MMSEData/' not found. Please run MMSE_silent().")
end

if ~isfolder(dir_feat)
    mkdir(dir_feat)
end

% Get relevant files
files = dir(fullfile(dir_MMSE, '**/*.csv'));
files = {files.name};

% Start parallel processing
delete(gcp('nocreate'));
parpool("Processes");

parfor i_file = 1:length(files)

    file = files{i_file};
    input_file = strcat(dir_MMSE, file);
    output_file = strcat(dir_feat, "features_", file);

    mmse_to_feats(input_file, output_file);

end