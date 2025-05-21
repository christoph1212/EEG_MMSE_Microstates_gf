function MMSE_silent(IndexSubset, SubsetSize, dir_Raw, dir_Root, dir_Log, Parpoolsize, Overwrite)
%% Preprocessing EEG files and inital calculation of MMSE vectors
% for a subset of Subjects, run preprocessing and MMSE vector calculation in parallel
% Outputs:
%	Creates a .set file of the preprocessed data, one per subject and condition and run (six files per subject)
%	Creates a .csv file of the MMSE vectors to be processed further, one per subject and condition and run

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
    dir_Log = "/work/bay2875/Logs/Resting_ComplexityIntelligence/PreProc/";
end
if nargin<4
    dir_Root = "/work/bay2875/Resting_ComplyIntelligence/";
end
if nargin<3
    dir_Raw = "/work/bay2875/RawData/task-Resting/";
end
if nargin<2
    SubsetSize = "32";
end
if nargin<1
    IndexSubset = "1";
end

IndexSubset = str2num(IndexSubset);
SubsetSize = str2num(SubsetSize);
Overwrite = str2num(Overwrite);
Parpoolsize = str2num(Parpoolsize);


%% Directory where file should be saved
dir_PreProc = strcat(dir_Root, 'PreprocessedData/');
dir_MMSE = strcat(dir_Root, 'MMSEData/');
dir_Snipplet = strcat(dir_Root, 'Snipplet/');
dir_Log_MMSE =  strcat(dir_Log, 'MMSE/');
mkdir(dir_PreProc)
mkdir(dir_MMSE)
mkdir(dir_Log)
mkdir(dir_Snipplet)
mkdir(dir_Log_MMSE)

%% Prepare List of Files to be Processed
Snipplet_List = dir(strcat(dir_Snipplet, '*.set'));  %get list of files in *.set format

%% To parallelize across nodes, subsets of Subjects are created
IndexSubset = ((IndexSubset-1)*SubsetSize+1): IndexSubset*SubsetSize;
if max(IndexSubset) > length(Snipplet_List)
    IndexSubset = IndexSubset(ismember(IndexSubset, 1:length(Snipplet_List)));
end

if length(IndexSubset) == 0
    fprintf('"Index does not contain any subjects.\n');
    return
end
Snipplet_List = Snipplet_List(IndexSubset);


%% for writing in logfile
fprintf('Attempting to calculate %d Subjects. \n', length(IndexSubset));
fprintf('InputFolder is %s. \n OutputFolder is %s. \n LogFolder is %s. \n', dir_Raw , dir_PreProc, dir_Log);


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



%% Infos on Triggers and Conditions
SplitStruct = struct('Trigger', {11 12 21 22 31 32}, 'Condition', {'first_run_eyes_open' 'first_run_eyes_closed' 'second_run_eyes_open' 'second_run_eyes_closed' 'third_run_eyes_open' 'third_run_eyes_closed'});


%% First step to increase calculation speed- run multiple subjects in parallel
%change PARFOR
%delete(gcp('nocreate')); % make sure that previous pooling is closed
%parpool(Parpoolsize);

FileName = "";
InputFile = "";
Cond_FileName = "";



%% Looped preprocessing
%change PARFOR
MMSE_Files = dir(dir_MMSE);




for i_File = 1:length(IndexSubset)
    try
	EEG = struct([]);
        if (sum(contains({MMSE_Files.name},Snipplet_List(i_File).name)) == 9) && (Overwrite == 0)
            continue
        end
        
        InputFile = char(strcat(dir_Snipplet,Snipplet_List(i_File).name));
        evalc("EEG = pop_loadset(InputFile);");
        
        
        
        %% Delete finally! Still problems with interpolation
        
        AllChannels = {ChannelSets.sets}; AllChannels = sort(unique(horzcat(AllChannels{:})));
        ElectrodeLocs = {'C3' , -90 , 0.267 , 4.55e-17 , 0.743 , 0.669 , 90 , 42 ; ...
            'C4' , 90 , 0.267 , 4.55e-17 , -0.743 , 0.669 , -90 , 42 ; ...
            'CP1' , -135 , 0.181 , -0.381 , 0.381 , 0.842 , 135 , 57.4 ; ...
            'CP2' , 135 , 0.181 , -0.381 , -0.381 , 0.842 , -135 , 57.4 ; ...
            'CP5' , -111 , 0.408 , -0.344 , 0.895 , 0.284 , 111 , 16.5 ; ...
            'CP6' , 111 , 0.408 , -0.344 , -0.895 , 0.284 , -111 , 16.5 ; ...
            'CZ' , 0 , 0 , 6.12e-17 , 0 , 1 , 0 , 90 ; ...
            'F3' , -39.9 , 0.344 , 0.677 , 0.566 , 0.469 , 39.9 , 28 ; ...
            'F4' , 39.9 , 0.344 , 0.677 , -0.566 , 0.469 , -39.9 , 28 ; ...
            'F7' , -53.9 , 0.528 , 0.587 , 0.805 , -0.088 , 53.9 , -5.05 ; ...
            'F8' , 53.9 , 0.528 , 0.587 , -0.805 , -0.088 , -53.9 , -5.05 ; ...
            'FC3' , -62.4 , 0.288 , 0.365 , 0.697 , 0.617 , 62.4 , 38.1 ; ...
            'FC4' , 62.4 , 0.288 , 0.365 , -0.697 , 0.617 , -62.4 , 38.1 ; ...
            'FP1' , -17.9 , 0.515 , 0.951 , 0.307 , -0.0471 , 17.9 , -2.7 ; ...
            'FP2' , 17.9 , 0.515 , 0.951 , -0.307 , -0.0471 , -17.9 , -2.7 ; ...
            'FZ' , 0 , 0.253 , 0.714 , 0 , 0.7 , 0 , 44.4 ; ...
            'O1' , -162 , 0.515 , -0.95 , 0.309 , -0.0471 , 162 , -2.7 ; ...
            'O2' , 162 , 0.515 , -0.95 , -0.309 , -0.0471 , -162 , -2.7 ; ...
            'OZ' , 180 , 0.507 , -1 , -1.22e-16 , -0.0209 , -180 , -1.2 ; ...
            'P3' , -140 , 0.344 , -0.676 , 0.568 , 0.469 , 140 , 28 ; ...
            'P4' , 140 , 0.344 , -0.676 , -0.568 , 0.469 , -140 , 28 ; ...
            'P7' , -126 , 0.528 , -0.586 , 0.806 , -0.088 , 126 , -5.05 ; ...
            'P8' , 126 , 0.528 , -0.586 , -0.806 , -0.088 , -126 , -5.05 ; ...
            'PO3' , -158 , 0.421 , -0.899 , 0.363 , 0.245 , 158 , 14.2 ; ...
            'PO4' , 158 , 0.421 , -0.899 , -0.363 , 0.245 , -158 , 14.2 ; ...
            'PZ' , 180 , 0.253 , -0.714 , -8.75e-17 , 0.7 , -180 , 44.4 ; ...
            'T7' , -90 , 0.533 , 6.09e-17 , 0.995 , -0.104 , 90 , -5.97 ; ...
            'T8' , 90 , 0.533 , 6.09e-17 , -0.995 , -0.104 , -90 , -5.97 };
        
        
        
        EEG.chaninfo = rmfield(EEG.chaninfo, 'removedchans');
        EEG.chaninfo = rmfield(EEG.chaninfo, 'nodatchans');
        evalc("EEG = eeg_checkset( EEG );");

        toInterpolate = [];
        if ~all(ismember(upper(AllChannels), upper({EEG.chanlocs.labels})))
            % if not, interpolate
            missing = find(~(ismember(upper(AllChannels), upper({EEG.chanlocs.labels}))));
            
            nrChan = EEG.nbchan;
            ToInterpolate = [];
            for i = missing
                nrChan = nrChan+1,
                ToInterpolate = [ToInterpolate, nrChan];
                EEG.data(nrChan,:) = 0;
                EEG.nbchan = size(EEG.data,1);
                EEG.chanlocs(nrChan).labels = ElectrodeLocs{i, 1};
                EEG.chanlocs(nrChan).theta = ElectrodeLocs{i, 2};
                EEG.chanlocs(nrChan).radius = ElectrodeLocs{i, 3};
                EEG.chanlocs(nrChan).X = ElectrodeLocs{i, 4};
                EEG.chanlocs(nrChan).Y = ElectrodeLocs{i, 5};
                EEG.chanlocs(nrChan).Z = ElectrodeLocs{i, 6};
                EEG.chanlocs(nrChan).sph_theta = ElectrodeLocs{i, 7};
                EEG.chanlocs(nrChan).sph_phi = ElectrodeLocs{i, 8};
                EEG.chanlocs(nrChan).sph_radius = 1;
                EEG.chanlocs(nrChan).type = 'EEG';
                EEG.chanlocs(nrChan).ref = 'average';
                
            end
        evalc("EEG = eeg_checkset( EEG );");
        evalc("EEG = pop_interp(EEG, ToInterpolate, 'spherical');");
            
            
        end
        
        evalc("EEG = eeg_checkset( EEG );");
        
        
        %% Calculate MMSE vectors
        % Write to logfile
        fprintf('File %s: ; Calculating MMSE vectors.\n', Snipplet_List(i_File).name);
        for elec = 1:length(ChannelSets)
            
            ElecName =  strcat("MMSE_", ChannelSets(elec).setnames, "_", Snipplet_List(i_File).name);
            if ~isfile (strcat(dir_MMSE, ElecName))
                elecAn = ChannelSets(elec).sets;
                channelset = ChannelSets(elec).setnames;
                try
                    MMSE_set_values = mmse_load_validate_save(EEG,elecAn,channelset,char(Snipplet_List(i_File).name),char(dir_MMSE));
                catch e
                    
                    
                    ErrorMessage = string(e.message);
                    for ierrors = 1:length(e.stack)
                        ErrorMessage = strcat(ErrorMessage, "//", e.stack(ierrors).name, ", Line: ",  num2str(e.stack(ierrors).line));
                    end
                    % make error log
                    ErrorFile =strsplit(InputFile, "/");
                    ErrorFile = ErrorFile{end};
                    fprintf('File %s: ;Error with MMSE: %s.\n', Snipplet_List(i_File).name, e.message);
                    fprintf('The Error Message is: \n %s \n',ErrorMessage);
                    ErrorFile = strcat(dir_Log_MMSE, 'ErrorRunning_MMSE', ErrorFile, '_', Cond_FileName , '_', num2str(elec), '.txt' );
                    ErrorFile = strrep(ErrorFile, '.set', '');
                    fid1 = fopen( ErrorFile, 'wt' );
                    fprintf(fid1, 'Error-Subject= %s \n The returned Error Message is: \n  \n %s \n', InputFile,  ErrorMessage);
                    fclose(fid1);
                    
                end
            end
            fprintf('Finished Subject %s.\n',Snipplet_List(i_File).name);
        end
    catch e
        % If error ocurrs, create ErrorMessage(concatenated for all nested errors).
        ErrorMessage = string(e.message);
        for ierrors = 1:length(e.stack)
            ErrorMessage = strcat(ErrorMessage, "//", num2str(e.stack(ierrors).name), ", Line: ",  num2str(e.stack(ierrors).line));
        end
        % make error log
        ErrorFile =strsplit(InputFile, "/");
        ErrorFile = ErrorFile{end};
        fprintf('Problem with executing File %s. \n',ErrorFile);
        fprintf('The Error Message is: \n %s \n',ErrorMessage);
        ErrorFile = strcat(dir_Log, ErrorFile, '_', Snipplet_List(i_File).name , '.txt' );
        ErrorFile = strrep(ErrorFile, '.set', '');
        fid1 = fopen( ErrorFile, 'wt' );
        fprintf(fid1, 'Error-Subject= %s \n The returned Error Message is: \n  \n %s \n', InputFile,  ErrorMessage);
        fclose(fid1);
    end
    
end
