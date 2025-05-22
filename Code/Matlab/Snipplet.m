function Snipplet(dir_Root, dir_Log, Overwrite)
%% Extraction of epochs from preprocessed EEG data
% Outputs:
%   Creates a .set file for each preprocessed file of either:
%   - 40s continuous epochs or 
%   - longest continuous epoch of the preprocessed data
%
% EEG-Channel order will be harmonized.
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

fprintf("*****************\nStarting Snipping\n*****************\n")
%% Create relevant directories
dir_Preproc = strcat(dir_Root, 'Data/Preprocessed/');
dir_Snipplet = strcat(dir_Root, 'Data/Snipplet/');
dir_Log = strcat(dir_Log, 'Preproc/');

if ~isfolder(dir_Snipplet)
    mkdir(dir_Snipplet)
end

%% Prepare List of Files to be Processed
PreProcFiles = dir(fullfile(dir_Preproc, '**/*.set'));  %get list of files in *.set format in any subfolder

%% Prepare Template to sort Channels
EEG_template = struct([]);
template_file = [PreProcFiles(1).folder '\' PreProcFiles(1).name];
evalc("EEG_template = pop_loadset(template_file);");

template_chans = {EEG_template.chanlocs.labels};

%% First step to increase calculation speed- run multiple subjects in parallel
delete(gcp('nocreate')); % make sure that previous pooling is closed
parpool("Processes");

Files_Snipplet = dir(dir_Snipplet);

%% Looped epoching
parfor i_Sub = 1:length(PreProcFiles)

    % Check if Subject has been preprocessed 
    if (sum(contains({Files_Snipplet.name}, PreProcFiles(i_Sub).name)) == 1) && (Overwrite == 0)
        continue
    end

    FileName = PreProcFiles(i_Sub).name;
    run_silent(PreProcFiles, dir_Snipplet, dir_Log, i_Sub, template_chans);

end

fprintf("*****************\nFinished Snipping\n*****************\n")

end

function run_silent(PreProcFiles, dir_Snipplet, dir_Log, i_Sub, template_chans)

try

    fprintf('Snipping: %s. \n', PreProcFiles(i_Sub).name);
    
    % Initate variables, and load preprocessed EEG file
    EEG = struct([]);
    InputFile = [PreProcFiles(i_Sub).folder,'/',PreProcFiles(i_Sub).name];
    evalc("EEG = pop_loadset(InputFile);");

    % Sort Channels according to Template
    current_labels = {EEG.chanlocs.labels};

    [~, sort_idx] = ismember(template_chans, current_labels);

    EEG.data = EEG.data(sort_idx, :, :);
    EEG.chanlocs = EEG.chanlocs(sort_idx);

    % LogFile to append Snipplet length
    startIdx = strfind(PreProcFiles(i_Sub).name, 'sub-');
    endIdx = strfind(PreProcFiles(i_Sub).name, '_eeg.set');
    LogFileName = ['Log_' PreProcFiles(i_Sub).name(startIdx:endIdx + length('eeg')) '.csv'];
    LogFile = strcat(dir_Log, LogFileName);
       
    %% extract 40s of consecutive timepoints
    % Check length of longest clean epoch
    Boundaries = [EEG.event.latency]';
    Start_Boundary = [1; ceil(Boundaries)];
    End_Boundary = [floor(Boundaries); EEG.pnts];
    LengthClean = End_Boundary - Start_Boundary;    
    % which interval is larger than 40s?
    PossibleWindows = find(LengthClean>40*EEG.srate);

    if isempty(PossibleWindows)
        % extract longest clean epoch
        [MaxV, PickedWindow] = max(LengthClean);
        PickedStart = Start_Boundary(PickedWindow);
        evalc("EEG = pop_select(EEG,'point',[PickedStart PickedStart+MaxV] );");
        FileName =[num2str(floor(MaxV/EEG.srate)), '_seconds_', char(PreProcFiles(i_Sub).name)];
        evalc("EEG = pop_saveset(EEG, 'filename', FileName, 'filepath', char(dir_Snipplet), 'savemode', ['onefile']);");

        % Save Snipplet length to Log File
        Log_table = readtable(LogFile);
        Log_table.SnippletLength = num2str(floor(MaxV/EEG.srate));
        writetable(Log_table, LogFile);

    else
        PickedWindow = PossibleWindows(randi(length(PossibleWindows)));
        % in that window, start randomly at a point that allows 40s
        PossibleStarts = Start_Boundary(PickedWindow) : (End_Boundary(PickedWindow) - (40*EEG.srate));
        PickedStart = PossibleStarts(randi(length(PossibleStarts)));

        % Take that Start Trigger and add 40 s
        evalc("EEG = pop_select(EEG,'point',[PickedStart PickedStart+(40*EEG.srate)] );");
        FileName = ['40_seconds_', char(PreProcFiles(i_Sub).name)];
        evalc("EEG = pop_saveset(EEG, 'filename', FileName, 'filepath', char(dir_Snipplet), 'savemode', ['onefile']);");
        
        % Save Snipplet length to Log File
        Log_table = readtable(LogFile);
        Log_table.SnippletLength = 40;
        writetable(Log_table, LogFile);
    end
    
catch e
    % If error ocurrs, create ErrorMessage(concatenated for all nested errors).
    ErrorMessage = string(e.message);
    for ierrors = 1:length(e.stack)
        ErrorMessage = strcat(ErrorMessage, "//", num2str(e.stack(ierrors).name), ", Line: ",  num2str(e.stack(ierrors).line));
    end
    % make error log
    fprintf('Error Snipping File: %s;\n%s.\n', PreProcFiles(i_Sub).name, ErrorMessage);
    ErrorFile = strcat(dir_Log, 'Error_Snipplet_', PreProcFiles(i_Sub).name, '.txt' );
    fid1 = fopen( ErrorFile, 'wt' );
    fprintf(fid1, 'Error-Subject: %s \nSnipping Error: \n%s \n', PreProcFiles(i_Sub).name,  ErrorMessage);
    fclose(fid1);
    
end

end
