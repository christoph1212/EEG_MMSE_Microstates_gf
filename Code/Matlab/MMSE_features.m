function MMSE_features(dir_Root, dir_Log, Overwrite)
%% Calculation of MMSE features
% Outputs:
%   Creates a .csv file of the MMSE features. One per subject, condition 
%   and channel set
%
% Inputs:
%   dir_Root:   String pointing to the project's parent folder
%   dir_Log:    String pointing to Log-Files directory
%   Overwrite:  Numeric (0|1). Should calculation be recalculated? 
%               Default: 0
%
% This script was created by: Christoph Fruehlinger (May 2025)
%
% This script was adapted by: Christoph Fruehlinger (May 2025)

%% get from function input
if nargin < 3
    Overwrite = 0;
end

fprintf("*********************************\nStarting MMSE Feature Calculation\n*********************************\n")

%% Directory where file should be saved
dir_MMSE = strcat(dir_Root, "Data/MMSEData/");
dir_feat = strcat(dir_Root, "Data/MMSEFeatures/");
dir_Log_Feat = strcat(dir_Log, 'MMSE/');

if ~isfolder(dir_feat)
    mkdir(dir_feat)
end

%% Prepare List of Files to be Processed
Vector_files = dir(fullfile(dir_MMSE, '**/*.csv'));

Feature_files = dir(dir_feat);

%% Increase calculation speed by running multiple subjects in parallel
delete(gcp('nocreate'));
parpool("Processes");

parfor i_File = 1:length(Vector_files)

    try
        
        if (sum(contains({Feature_files.name}, Vector_files(i_File).name)) == 9) && (Overwrite == 0)
            continue
        end

        file = Vector_files(i_File).name;
        input_file = strcat(dir_MMSE, file);
        output_file = strcat(dir_feat, "features_", file);
    
        mmse_to_feats(input_file, output_file);

    catch e

        % If error ocurrs, create ErrorMessage
        ErrorMessage = string(e.message);
        for ierrors = 1:length(e.stack)
            ErrorMessage = strcat(ErrorMessage, "//", num2str(e.stack(ierrors).name), ", Line: ",  num2str(e.stack(ierrors).line));
        end
        % make error log
        fprintf('File: %s;\nError with Execution: %s.\n', Vector_files(i_File).name, ErrorMessage);
        ErrorFile = strcat(dir_Log_Feat, 'Error_Features_', Vector_files(i_File).name , '.txt' );
        ErrorFile = strrep(ErrorFile, '.set.csv', '');
        fid1 = fopen( ErrorFile, 'wt' );
        fprintf(fid1, 'Error-Subject: %s \nFeature Error: \n%s \n', Vector_files(i_File).name,  ErrorMessage);
        fclose(fid1);

    end

end

fprintf("*********************************\nFinished MMSE Feature Calculation\n*********************************\n")
