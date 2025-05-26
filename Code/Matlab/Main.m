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

%% Setup

dir_Root    = "E:/Complexity/";                   % path to project
dir_Raw     = strcat(dir_Root, "Data/RawData/");  % path to raw data
dir_Log     = strcat(dir_Root, "Data/Log/");      % path where log data should be stored
Overwrite   = 1;                                  % 0 = no; 1 = yes

% Start EEGLAB
dir_eeglab  = "E:/Complexity/Code/Matlab/eeglab2025.0.0";
addpath(dir_eeglab);
eeglab nogui
clc

% Check for necessary Plugins
if ~ismember('trimOutlier', {PLUGINLIST.plugin})
    error("Error: Plugin 'trimOutlier' is not installed. Please install or " + ...
        "ideally clone 'eeglab2025.0.0' folder from GitHub repository.")
end

if ~ismember('Adjust', {PLUGINLIST.plugin})
    error("Error: Plugin 'Adjust' is not installed. Please install or " + ...
        "ideally clone 'eeglab2025.0.0' folder from GitHub repository.")
end

fprintf("Let's get started\n")
fprintf(['Alright, all Plugin are available.\nYour folder settings:\n\n' ...
    'Root-Folder: %s\nRaw-Folder: %s\nLog-Folder: %s\nOverwrite: %d\n'], ...
    dir_Root, dir_Raw, dir_Log, Overwrite)

%% Preprocessing

Preproc(dir_Raw, dir_Root, dir_Log, Overwrite)

%% Epoch Extraction

Snipplet(dir_Root, dir_Log, Overwrite)

%% MMSE Complexity Analysis

% MMSE Vectors
MMSE_silent(dir_Root, dir_Log, Overwrite)

% MMSE Feature Extraction
MMSE_features(dir_Root, dir_Log, Overwrite)
