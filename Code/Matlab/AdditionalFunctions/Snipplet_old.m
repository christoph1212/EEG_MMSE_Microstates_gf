function Snipplet_old(IndexSubjects, SubsetSize, dir_Raw, dir_Root, dir_Log, Parpoolsize, Overwrite)
%% Extraction of epochs from preprocessed EEG data
% Outputs:
%	Creates a .set file of either:
%      - 40s continuous epochs or 
%      - longest continuous epoch of the preprocessed data
%  one per subject and condition and run (six files per subject)

% Inputs:
%	IndexSubjects: integer >0. subsets are created based on List of Subjects, of size Subsetsize
%	SubsetSize: integer >0. size of subset to be calculated on one node. Default: 32
% 	dir_Raw: String pointing to the folder where Raw Data is Default: "/work/bay2875/RawData/task-Resting/"
%	dir_Root: String pointing to the folder where Output Data will be saved to. Default: "/work/bay2875/Resting_ComplexityIntelligence/"
%	Overwrite: Numeric (0|1). Should calculation be recalculated? Default: 0

%% get from function input
if nargin<7
    Overwrite = "0";
end
if nargin<6
    Parpoolsize = "8";
end
if nargin<5
    dir_Log = "/work/bay2875/Resting_Complexity/Logs/PreProc/";
end
if nargin<4
    dir_Root = "/work/bay2875/Resting_Complexity/";
end
if nargin<3
    dir_Raw = "/work/bay2875/RawData/task-Resting/";
end
if nargin<2
    SubsetSize = "32";
end
if nargin<1
    IndexSubjects = "1";
end

if isstring(IndexSubjects)
    IndexSubjects = str2num(IndexSubjects);
end
if isstring(SubsetSize)
    SubsetSize = str2num(SubsetSize);
end
if isstring(Overwrite)
    Overwrite = str2num(Overwrite);
end
if isstring(Parpoolsize)
    Parpoolsize = str2num(Parpoolsize);
end

Dummy = "";


%% Directory where file should be saved
dir_Preproc = strcat(dir_Root, 'Data/Preprocessed/');
dir_Snipplet = strcat(dir_Root, 'Data/Snipplet/');
dir_MMSE = strcat(dir_Root, 'Data/MMSEData/');
dir_LogData = strcat(dir_Log, 'Incomplete/');

if ~isfolder(dir_Preproc)
    mkdir(dir_Preproc)
end
if ~isfolder(dir_MMSE)
    mkdir(dir_MMSE)
end
if ~isfolder(dir_Log)
    mkdir(dir_Log)
end
if ~isfolder(dir_LogData)
    mkdir(dir_LogData)
end
if ~isfolder(dir_Snipplet)
    mkdir(dir_Snipplet)
end

%% Prepare List of Files to be Processed
PreProcFiles = dir(fullfile(dir_Preproc, '**/*.set'));  %get list of files in *.set format in any subfolder

%% To parallelize across nodes, subsets of Subjects are created
IndexSubjects = ((IndexSubjects-1)*SubsetSize+1): IndexSubjects*SubsetSize;
if max(IndexSubjects) > length(PreProcFiles)
    IndexSubjects = IndexSubjects(ismember(IndexSubjects, 1:length(PreProcFiles)));
end

if length(IndexSubjects) == 0
    fprintf('"Index does not contain any subjects.\n');
    return
end
PreProcFiles = PreProcFiles(IndexSubjects);





%% for writing in logfile
fprintf('\nAttempting to calculate %d Subjects. \n', length(IndexSubjects));
fprintf('InputFolder is %s. \nOutputFolder is %s. \nLogFolder is %s. \n\n', dir_Raw , dir_Preproc, dir_Log);


%% Define Channel Sets for integrated MMSE vector calculation
%% Define channel sets
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
ChannelSets = struct('sets',{elecAnF,elecAnFL,elecAnFR,elecAnC,elecAnP,elecAnPL,elecAnPR,elecAnML,elecAnMR}, 'setnames', {'F','FL','FR','C', 'P', 'PL', 'PR', 'ML', 'MR'});
AllChannels = {ChannelSets.sets}; AllChannels = sort(unique(horzcat(AllChannels{:})));

ElectrodeLocs = ["'C3'" , "'-90'" , "'0.267'" , "'4.55e-17'" , "'0.743'" , "'0.669'" , "'90'" , "'42'" ; ...
    "'C4'" , "'90'" , "'0.267'" , "'4.55e-17'" , "'-0.743'" , "'0.669'" , "'-90'" , "'42'" ; ...
    "'CP1'" , "'-135'" , "'0.181'" , "'-0.381'" , "'0.381'" , "'0.842'" , "'135'" , "'57.4'" ; ...
    "'CP2'" , "'135'" , "'0.181'" , "'-0.381'" , "'-0.381'" , "'0.842'" , "'-135'" , "'57.4'" ; ...
    "'CP5'" , "'-111'" , "'0.408'" , "'-0.344'" , "'0.895'" , "'0.284'" , "'111'" , "'16.5'" ; ...
    "'CP6'" , "'111'" , "'0.408'" , "'-0.344'" , "'-0.895'" , "'0.284'" , "'-111'" , "'16.5'" ; ...
    "'CZ'" , "'0'" , "'0'" , "'6.12e-17'" , "'0'" , "'1'" , "'0'" , "'90'" ; ...
    "'F3'" , "'-39.9'" , "'0.344'" , "'0.677'" , "'0.566'" , "'0.469'" , "'39.9'" , "'28'" ; ...
    "'F4'" , "'39.9'" , "'0.344'" , "'0.677'" , "'-0.566'" , "'0.469'" , "'-39.9'" , "'28'" ; ...
    "'F7'" , "'-53.9'" , "'0.528'" , "'0.587'" , "'0.805'" , "'-0.088'" , "'53.9'" , "'-5.05'" ; ...
    "'F8'" , "'53.9'" , "'0.528'" , "'0.587'" , "'-0.805'" , "'-0.088'" , "'-53.9'" , "'-5.05'" ; ...
    "'FC3'" , "'-62.4'" , "'0.288'" , "'0.365'" , "'0.697'" , "'0.617'" , "'62.4'" , "'38.1'" ; ...
    "'FC4'" , "'62.4'" , "'0.288'" , "'0.365'" , "'-0.697'" , "'0.617'" , "'-62.4'" , "'38.1'" ; ...
    "'FP1'" , "'-17.9'" , "'0.515'" , "'0.951'" , "'0.307'" , "'-0.0471'" , "'17.9'" , "'-2.7'" ; ...
    "'FP2'" , "'17.9'" , "'0.515'" , "'0.951'" , "'-0.307'" , "'-0.0471'" , "'-17.9'" , "'-2.7'" ; ...
    "'FZ'" , "'0'" , "'0.253'" , "'0.714'" , "'0'" , "'0.7'" , "'0'" , "'44.4'" ; ...
    "'O1'" , "'-162'" , "'0.515'" , "'-0.95'" , "'0.309'" , "'-0.0471'" , "'162'" , "'-2.7'" ; ...
    "'O2'" , "'162'" , "'0.515'" , "'-0.95'" , "'-0.309'" , "'-0.0471'" , "'-162'" , "'-2.7'" ; ...
    "'OZ'" , "'180'" , "'0.507'" , "'-1'" , "'-1.22e-16'" , "'-0.0209'" , "'-180'" , "'-1.2'" ; ...
    "'P3'" , "'-140'" , "'0.344'" , "'-0.676'" , "'0.568'" , "'0.469'" , "'140'" , "'28'" ; ...
    "'P4'" , "'140'" , "'0.344'" , "'-0.676'" , "'-0.568'" , "'0.469'" , "'-140'" , "'28'" ; ...
    "'P7'" , "'-126'" , "'0.528'" , "'-0.586'" , "'0.806'" , "'-0.088'" , "'126'" , "'-5.05'" ; ...
    "'P8'" , "'126'" , "'0.528'" , "'-0.586'" , "'-0.806'" , "'-0.088'" , "'-126'" , "'-5.05'" ; ...
    "'PO3'" , "'-158'" , "'0.421'" , "'-0.899'" , "'0.363'" , "'0.245'" , "'158'" , "'14.2'" ; ...
    "'PO4'" , "'158'" , "'0.421'" , "'-0.899'" , "'-0.363'" , "'0.245'" , "'-158'" , "'14.2'" ; ...
    "'PZ'" , "'180'" , "'0.253'" , "'-0.714'" , "'-8.75e-17'" , "'0.7'" , "'-180'" , "'44.4'" ; ...
    "'T7'" , "'-90'" , "'0.533'" , "'6.09e-17'" , "'0.995'" , "'-0.104'" , "'90'" , "'-5.97'" ; ...
    "'T8'" , "'90'" , "'0.533'" , "'6.09e-17'" , "'-0.995'" , "'-0.104'" , "'-90'" , "'-5.97'" ];


%% Infos on Triggers and Conditions
SplitStruct = struct('Trigger', {11 12 21 22 31 32}, 'Condition', {'first_run_eyes_open' 'first_run_eyes_closed' 'second_run_eyes_open' 'second_run_eyes_closed' 'third_run_eyes_open' 'third_run_eyes_closed'});


%% First step to increase calculation speed- run multiple subjects in parallel
delete(gcp('nocreate')); % make sure that previous pooling is closed
parpool(Parpoolsize);


FileName = "";
InputFile = "";
Cond_FileName = "";

Files_Snipplet = dir(dir_Snipplet);
%% Looped preprocessing
parfor i_Sub = 1:length(IndexSubjects)
    % Check if Subject has been preprocessed 
    if (sum(contains({Files_Snipplet.name},PreProcFiles(i_Sub).name)) == 1) && (Overwrite == 0)
        continue
    end
    FileName = PreProcFiles(i_Sub).name;
    run_silent(Overwrite, IndexSubjects, SubsetSize, PreProcFiles, dir_Preproc, dir_Snipplet, dir_MMSE, dir_Log,dir_LogData, ChannelSets, SplitStruct, FileName, InputFile, Cond_FileName, i_Sub, Dummy, ElectrodeLocs, AllChannels);
end

fprintf('\nEnd Snipping. Starting MMSE calculation... \n\n');

end

function run_silent(Overwrite, IndexSubjects, SubsetSize, PreProcFiles, dir_Preproc, dir_Snipplet, dir_MMSE, dir_Log, dir_LogData, ChannelSets, SplitStruct, FileName, InputFile, Cond_FileName, i_Sub, Dummy, ElectrodeLocs, AllChannels)

try %if error occurs still continue with next file, but make log
    % Write to logfile
    fprintf('Currently Snipping File %s. \n', PreProcFiles(i_Sub).name);
    
    % Initate variables, and load EEG lab File
    EEG = struct([]);
    InputFile = [PreProcFiles(i_Sub).folder,'/',PreProcFiles(i_Sub).name];
    evalc("EEG = pop_loadset(InputFile);");
    
    
    %% extract 40s of consecutive timepoints
    Boundaries = [EEG.event.latency]';
    Start_Boundary = [1; ceil(Boundaries)];
    End_Boundary = [floor(Boundaries); EEG.pnts];
    LengthClean = End_Boundary - Start_Boundary;
    % which interval is larger than 40s?
    PossibleWindows = find(LengthClean>40*EEG.srate);
    if isempty(PossibleWindows)
        [MaxV, PickedWindow] = max(LengthClean);
         PickedStart = Start_Boundary(PickedWindow);
         EEG = pop_select(EEG,'point',[PickedStart PickedStart+MaxV] );
         FileName =[num2str(floor(MaxV/EEG.srate)), '_seconds_', char(PreProcFiles(i_Sub).name)];
         evalc("EEG = pop_saveset(EEG, 'filename', FileName, 'filepath', char(dir_Snipplet), 'savemode', ['onefile']);");

    else
        PickedWindow = PossibleWindows(randi(length(PossibleWindows)));
        % in that window, start randomly at a point that allows 40s
        PossibleStarts = Start_Boundary(PickedWindow) : (End_Boundary(PickedWindow) - (40*EEG.srate));
        PickedStart = PossibleStarts(randi(length(PossibleStarts)));
        % Take that Start Trigger and add 40 s
        evalc("EEG = pop_select(EEG,'point',[PickedStart PickedStart+(40*EEG.srate)] );");
        FileName = ['40_seconds_', char(PreProcFiles(i_Sub).name)];
        evalc("EEG = pop_saveset(EEG, 'filename', FileName, 'filepath', char(dir_Snipplet), 'savemode', ['onefile']);");
    end
    
catch e
    % If error ocurrs, create ErrorMessage(concatenated for all nested errors).
    ErrorMessage = string(e.message);
    for ierrors = 1:length(e.stack)
        ErrorMessage = strcat(ErrorMessage, "//", num2str(e.stack(ierrors).name), ", Line: ",  num2str(e.stack(ierrors).line));
    end
    % make error log
    fprintf('Problem with executing File %s. \n',FileName);
    fprintf('The Error Message is: \n %s \n',ErrorMessage);
    ErrorFile = strcat(dir_Log, 'Error_Snipplet_', FileName, '.txt' );
    fid1 = fopen( ErrorFile, 'wt' );
    fprintf(fid1, 'Error-Subject= %s \n The returned Error Message is: \n  \n %s \n', InputFile,  ErrorMessage);
    fclose(fid1);
    
end
end
