function MMSE_silent(dir_Root, dir_Log, Overwrite)
%% Calculation of MMSE vectors
% Outputs:
%   Creates a .csv file of the MMSE vectors to be processed further, one 
%   per subject, condition and channel set
%
% Inputs:
%   dir_Root:   String pointing to the project's parent folder
%   dir_Log:    String pointing to Log-Files directory
%   Overwrite:  Numeric (0|1). Should calculation be recalculated? 
%               Default: 0
%
% This script was created by: Name (Date)
%
% This script was adapted by: Christoph Fruehlinger (May 2025)

%% get from function input
if nargin < 3
    Overwrite = 0;
end

fprintf("********************************\nStarting MMSE Vector Calculation\n********************************\n")

%% Directory where file should be saved
dir_MMSE = strcat(dir_Root, 'Data/MMSEData/');
dir_Snipplet = strcat(dir_Root, 'Data/Snipplet/');
dir_Log_MMSE =  strcat(dir_Log, 'MMSE/');

if ~isfolder(dir_MMSE)
    mkdir(dir_MMSE)
end
if ~isfolder(dir_Log)
    mkdir(dir_Log)
end
if ~isfolder(dir_Log_MMSE)
    mkdir(dir_Log_MMSE)
end

%% Prepare List of Files to be Processed
Snipplet_List = dir(strcat(dir_Snipplet, '*.set'));  %get list of files in *.set format

fprintf('\nMMSE calculation for %d Files. \n', length(Snipplet_List));

%% Define Channel Sets for integrated MMSE vector calculation
elecAnF = {'f7','f8','f3','f4'};
elecAnFL = {'fp1','f7','f3','fc3'};
elecAnFR = {'fp2','f8','f4','fc4'};
elecAnC = {'fz','cz','pz','oz'};
elecAnP = {'p3','p4','p7','p8'};
elecAnPL = {'p7','p3','o1','po3'};
elecAnPR = {'p8','p4','o2','po4'};
elecAnML = {'t7','c3','cp5','cp1'};
elecAnMR = {'t8','c4','cp6','cp2'};

% combine different channel sets into one struct used in the loop
ChannelSets = struct('sets',{elecAnF,elecAnFL,elecAnFR,elecAnC,elecAnP,elecAnPL,elecAnPR,elecAnML,elecAnMR}, ...
    'setnames', {'F','FL','FR','C', 'P', 'PL', 'PR', 'ML', 'MR'});


%% Increase calculation speed by running multiple subjects in parallel
delete(gcp('nocreate')); % make sure that previous pooling is closed
parpool("Processes");

%% Looped preprocessing
MMSE_Files = dir(dir_MMSE);

parfor i_File = 1:length(Snipplet_List)
    try
        
        if (sum(contains({MMSE_Files.name}, Snipplet_List(i_File).name)) == 9) && (Overwrite == 0)
            continue
        end
        
        EEG = pop_loadset(char(strcat(dir_Snipplet, Snipplet_List(i_File).name)));               
        
        %% Calculate MMSE vectors
        fprintf('\nFile: %s; Calculating MMSE vectors.\n', Snipplet_List(i_File).name);

        % Loop over Channel Sets
        for elec = 1:length(ChannelSets)
            
            ElecName =  strcat("MMSE_", ChannelSets(elec).setnames, "_", Snipplet_List(i_File).name);
            if ~isfile (strcat(dir_MMSE, ElecName))
                elecAn = ChannelSets(elec).sets;
                channelset = ChannelSets(elec).setnames;
                try
                    MMSE_set_values = mmse_load_validate_save(EEG,elecAn,channelset,char(Snipplet_List(i_File).name),char(dir_MMSE));

                catch e
                    % If error ocurrs, create ErrorMessage
                    ErrorMessage = string(e.message);
                    for ierrors = 1:length(e.stack)
                        ErrorMessage = strcat(ErrorMessage, "//", e.stack(ierrors).name, ", Line: ",  num2str(e.stack(ierrors).line));
                    end
                    % make error log
                    fprintf('File: %s; \nError with MMSE: %s.\n', Snipplet_List(i_File).name, ErrorMessage);
                    ErrorFile = strcat(dir_Log_MMSE, 'Error_MMSE_', Snipplet_List(i_File).name, '_', ChannelSets(elec).setnames, '.txt' );
                    ErrorFile = strrep(ErrorFile, '.set', '');
                    fid1 = fopen( ErrorFile, 'wt' );
                    fprintf(fid1, 'Error-Subject: %s \nMMSE Error\n%s \n', Snipplet_List(i_File).name,  ErrorMessage);
                    fclose(fid1);
                    
                end
            end
            fprintf('Finished Subject %s.\n',Snipplet_List(i_File).name);
        end

    catch e
        % If error ocurrs, create ErrorMessage
        ErrorMessage = string(e.message);
        for ierrors = 1:length(e.stack)
            ErrorMessage = strcat(ErrorMessage, "//", num2str(e.stack(ierrors).name), ", Line: ",  num2str(e.stack(ierrors).line));
        end
        % make error log
        fprintf('File: %s;\nError with Execution: %s.\n',Snipplet_List(i_File).name, ErrorMessage);
        ErrorFile = strcat(dir_Log_MMSE, 'Error_MMSE_', Snipplet_List(i_File).name , '.txt' );
        ErrorFile = strrep(ErrorFile, '.set', '');
        fid1 = fopen( ErrorFile, 'wt' );
        fprintf(fid1, 'Error-Subject: %s \nExecution Error: \n%s \n', Snipplet_List(i_File).name,  ErrorMessage);
        fclose(fid1);
    end
    
end

fprintf("********************************\nFinished MMSE Vector Calculation\n********************************\n")
