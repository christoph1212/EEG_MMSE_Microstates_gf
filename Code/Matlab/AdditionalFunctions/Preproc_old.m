function Preproc_old(IndexSubjects, SubsetSize, dir_Raw, dir_Root, dir_Log, Parpoolsize, Overwrite)
%% Preprocessing EEG files
% for a subset of Subjects, run preprocessing 
% Output:
%	Creates a .set file of the preprocessed data, one per subject and condition and run (six files per subject)

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
Dummy = "";  % ??
%% Prepare List of Files to be Processed
Raw_Files = dir(fullfile(dir_Raw, '**/*.set'));  % get list of files in *.set format in any subfolder

%% To parallelize across nodes, subsets of Subjects are created
IndexSubjects = ((IndexSubjects-1)*SubsetSize+1) : IndexSubjects*SubsetSize;
if max(IndexSubjects) > length(Raw_Files)
    IndexSubjects = IndexSubjects(ismember(IndexSubjects, 1:length(Raw_Files)));
end

if length(IndexSubjects) == 0
    fprintf('"Index does not contain any subjects.\n');
    return
end
Raw_Files = Raw_Files(IndexSubjects);

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


%% Infos on Triggers and Conditions
SplitStruct = struct('Trigger', {11, 12 21 22 31 32}, 'Condition', {'first_run_eyes_open' 'first_run_eyes_closed' 'second_run_eyes_open' 'second_run_eyes_closed' 'third_run_eyes_open' 'third_run_eyes_closed'});


%% First step to increase calculation speed- run multiple subjects in parallel
delete(gcp('nocreate')); % make sure that previous pooling is closed
parpool(Parpoolsize);


FileName = "";
InputFile = "";
Cond_FileName = "";


Files_PreProc= dir(dir_Preproc);
%% Looped preprocessing
parfor i_Sub = 1:length(IndexSubjects)
    % Check if Subject has been preprocessed 
    if (sum(contains({Files_PreProc.name},Raw_Files(i_Sub).name)) == 2) && (Overwrite == 0)
        continue
    end
    run_silent(Overwrite, IndexSubjects, SubsetSize, Raw_Files, dir_Preproc, dir_Snipplet, dir_MMSE, dir_Log,dir_LogData, ChannelSets, SplitStruct, FileName, InputFile, Cond_FileName, i_Sub, Dummy, ElectrodeLocs, AllChannels);
end

fprintf('End Preprocessing. Starting snipping... \n\n');

end

function run_silent(Overwrite, IndexSubjects, SubsetSize, Raw_Files, dir_Preproc, dir_Snipplet, dir_MMSE, dir_Log, dir_LogData, ChannelSets, SplitStruct, FileName, InputFile, Cond_FileName, i_Sub, Dummy, ElectrodeLocs, AllChannels)

try %if error occurs still continue with next file, but make log
    % Write to logfile
    %fprintf('Currently Standardizing Subject %s: Select Channels, Reference, Sampling Rate. \n', Raw_Files(i_Sub).name);
    
    % Initate variables, and load EEG lab File
    EEG = struct([]);
    temp_out = [];
    InputFile = [Raw_Files(i_Sub).folder,'/',Raw_Files(i_Sub).name];
    Cond_FileName = "General";
    evalc("EEG = pop_loadset(InputFile);");
    
    %% Standardize: remove unwanted channels (some Labs recorded with some extra
    % electrodes that should not be included)
    Common_Channels =  {'FP1', 'FP2', 'AF7', 'AF8', 'AF3', 'AF4', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'FT7', 'FT8', 'FC5', 'FC6', 'FC3', 'FC4', ...
        'FC1', 'FC2', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'T7', 'T8', 'TP7', 'TP8', 'CP5', 'CP6', 'CP3', 'CP4', 'CP1', 'CP2', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', ...
        'P7', 'P8', 'PO7', 'PO8', 'PO3', 'PO4', 'O1', 'O2', 'OZ', 'POZ', 'PZ', 'CPZ', 'CZ', 'FCZ', 'FZ', ... % 'AFZ', 'FPZ' are used as grounds in some labs
        'VOGabove', 'VOGbelow', 'HOGl', 'HOGr','MASTl', 'MASTr'};
    
    Common_Channels = Common_Channels(ismember(Common_Channels, {EEG.chanlocs.labels})); % some labs miss e.g. VOGabove
    evalc("EEG = pop_select( EEG, 'channel',Common_Channels);");
    
    %% Standardize: rereference to FCZ (labs recorded with different online Refs,
    % is standadized across labs to FCZ here)
    if strcmp('CMS/DRL', EEG.Info_Lab.Reference)
        evalc("EEG = pop_reref(EEG, 'FCZ');");
    elseif strcmp('CZ', EEG.Info_Lab.Reference)
        evalc("EEG = pop_reref(EEG, 'FCZ', 'refloc', struct('labels',{'CZ'},'type',{'EEG'}, 'ref', [], 'urchan', [], 'theta',{0},'radius',{0},'X',{6.12e-17},'Y',{0},'Z',{1},'sph_theta',{0},'sph_phi',{90},'sph_radius',{1}));");
    elseif strcmp('FCZ', EEG.Info_Lab.Reference)
        % do nothing
    else
        msg = 'An error occured while rereferencing to FCz. EEG.Info_Lab.Reference does not contain CMS/DRL, Cz or FCz.';
        error(msg)
    end
    
    %% Bookkeeping electrodes
    EEG.cc.d0001__dnSamp.chanlocs_003_removed_unwanted = EEG.chanlocs;
    temp_comStr = ['% removed unwanted channels'];
    EEG = eegh(temp_comStr,EEG);
    %append new reference channel 'FCz' to chanlocs (as no data channel)
    % CF changed
    %RefInfo = {EEG.nbchan 'FCZ' 1 0 67.2 -0 0.127 0.388 0 0.922 'EEG' '' '' 0 };
    %{number label sph_radius sph_theta sph_phi theta radius X Y Z type ref urchan datachan}
    RefInfo = {EEG.nbchan 'labels' 'FCZ' 'theta' 0 'radius' 0.127 'X' 0.388 'Y' 0 'Z' 0.922 'sph_theta' 0 'sph_phi' 67.2 'sph_radius' 1 'type' 'EEG' 'datachan' 0 };
    %{number label theta radius X Y Z sph_theta sph_phi sph_radius type datachan}
    evalc("EEG = pop_chanedit(EEG, 'append', RefInfo);");
    % Chanlocs without Ref:
    EEG.cc.d0001__dnSamp.chanlocs_005_no_ref = EEG.chanlocs;
    % Chanloc of Reference
    EEG.cc.d0001__dnSamp.chanlocs_006_ref = EEG.chaninfo.nodatchans;
    
    %% Preparation: Downsample to 250/256 Hz/s
    EEG.cc.d0001__dnSamp.dnSampFreq = EEG.srate/2;
    evalc("[EEG,temp_out.comStr] = pop_resample(EEG,EEG.cc.d0001__dnSamp.dnSampFreq);");
    EEG = eegh(temp_out.comStr,EEG);
    
    
    %% Split into different Conditions
    % Prepare Matrix identfiying all points (to keep it continous, not
    % epoched)
    EEG_Complete = EEG;
    % check if all Conditions need to be run or not
    
    if contains(Raw_Files(i_Sub).name, 'run-1')
        Rel_SplitStruct.Trigger = [SplitStruct(1:2).Trigger];
        Rel_SplitStruct.Condition = {SplitStruct(1:2).Condition};

    elseif contains(Raw_Files(i_Sub).name, 'run-2')
        Rel_SplitStruct.Trigger = [SplitStruct(3:4).Trigger];
        Rel_SplitStruct.Condition = {SplitStruct(3:4).Condition};

    else contains(Raw_Files(i_Sub).name, 'run-3')
        Rel_SplitStruct.Trigger = [SplitStruct(5:6).Trigger];
        Rel_SplitStruct.Condition = {SplitStruct(5:6).Condition};
    end
    
    for i_cond = 1:length(Rel_SplitStruct.Trigger)
        try
            Cond_FileName = Rel_SplitStruct.Condition{i_cond};
            FileName = [Rel_SplitStruct.Condition{i_cond},'_',Raw_Files(i_Sub).name];
            % test that file has not been calculated yet
            if isfile(fullfile(dir_Preproc, FileName)) &&  Overwrite == 0
                fprintf('Previously finished Subject %s.\n',Raw_Files(i_Sub).name);
                continue
            end
            
            % Write to logfile
            fprintf('\nSubject: %s; Condition: %s; Filtering ICA AC.\n', Raw_Files(i_Sub).name, Cond_FileName);
            
            % Check if file is not empty or not big enough
            if EEG.pnts < EEG.srate
                fprintf('Subject %s: ; Condition: %s; Dataset empty.\n', Raw_Files(i_Sub).name, Cond_FileName);
                
                MissingDataFile = strcat(dir_LogData, "Error_NotEnoughData_for_Preproc_", FileName, '.txt');
                fid1 = fopen( MissingDataFile, 'wt' );
                fprintf(fid1, 'Missing Data-File= %s \n Data only includes %i Points. \n', FileName,  EEG.pnts);
                fclose(fid1);
                continue
                
                
            end
            
            evalc("EEG = pop_epoch( EEG_Complete, {num2str( Rel_SplitStruct.Trigger(i_cond))}, [0  60], 'epochinfo', 'yes');");
            evalc("EEG = eeg_epoch2continuous(EEG);");
            
            %% Demeaning: Remove mean from each channel
            for chanNum = 1:size(EEG.data,1)
                EEG.data(chanNum,:) = single(double(EEG.data(chanNum,:))-mean(double(EEG.data(chanNum, :))));
            end
            
            %% Filtering: Highpassfilter of 1 Hz
            EEG.cc.d0101__hiPass.hiPassFreq = 1;
            evalc("[EEG,temp_out.comStr,temp_out.hiPassCoefs] = pop_eegfiltnew(EEG,[],EEG.cc.d0101__hiPass.hiPassFreq,[],logical(1),[],0);");
            EEG.cc.d0101__hiPass.hiPassCoefs = temp_out.hiPassCoefs;
            EEG = eegh(temp_out.comStr,EEG);
            
            %% Cleaning: Remove bad Channels
            EEG.cc.d0201__remBad.freqlims =  [  0.00, 5.00;  5.00, 40.00   ]; % Hz
            EEG.cc.d0201__remBad.stdthresh = [ -5.00, 5.00; -2.50,  2.50   ]; % SD
            evalc("[EEG,temp_out.allrmchan,temp_out.specdata,temp_out.specfreqs,temp_out.comStr] = pop_rejchanspec(EEG,'freqlims',EEG.cc.d0201__remBad.freqlims,'stdthresh',EEG.cc.d0201__remBad.stdthresh);");
            %remember chanlocs wiInthout these removed bad channels (=only good ones)
            EEG.cc.d0201__remBad.chanlocs_006_BadChansRem = EEG.chanlocs;
            %remember chanlocs of these  bad channels (=only bad ones)
            EEG.cc.d0201__remBad.chanlocs_006_BadChans = EEG.cc.d0001__dnSamp.chanlocs_005_no_ref(temp_out.allrmchan);
            
            
            %% Cleaning: Bad Segments using trimOutlier I/III
            EEG.cc.d0201__remBad.channelSdLowerBound = -Inf; % SD
            EEG.cc.d0201__remBad.channelSdUpperBound =  Inf; % SD (consider using 35)
            EEG.cc.d0201__remBad.amplitudeThreshold  =  444; % uV (use 444 or 500 or 555)
            EEG.cc.d0201__remBad.pointSpreadWidth    = 2000; % samples (use 2000 or 4000)
            evalc("EEG = trimOutlier(EEG,EEG.cc.d0201__remBad.channelSdLowerBound,EEG.cc.d0201__remBad.channelSdUpperBound, EEG.cc.d0201__remBad.amplitudeThreshold,EEG.cc.d0201__remBad.pointSpreadWidth);");
            if EEG.pnts == 0
               error('All Datapoints deemed outlier with first trimOutlier')
            end
            
            %% Filtering: LowPassFilter of 40 Hz
            EEG.cc.d0301__loPass.loPassFreq = 40;
            evalc("[EEG,temp_out.comStr,temp_out.loPassCoefs] = pop_eegfiltnew(EEG,[],EEG.cc.d0301__loPass.loPassFreq,[],logical(0),[],0);");
            EEG.cc.d0301__loPass.loPassCoefs = temp_out.loPassCoefs;
            EEG = eegh(temp_out.comStr,EEG);
            
            
            %% Referencing: Average (to only good channels and FCZ)
            evalc("[EEG,temp_comStr] = pop_reref(EEG,[],'refloc',struct('labels',{'FCZ'},'type',{'EEG'},'ref', [], 'urchan', [], 'theta',{0},'radius',{0.127},'X',{0.388},'Y',{0},'Z',{0.922},'sph_theta',{0},'sph_phi',{67.2},'sph_radius',{1}, 'sph_theta_besa', [], 'sph_phi_besa', []));");
            EEG = eegh(temp_comStr,EEG);
            
            
            %% Cleaning: Bad Segments using trimOutlier II/III
            EEG.cc.d0501__remBad.channelSdLowerBound = -Inf; % SD
            EEG.cc.d0501__remBad.channelSdUpperBound =  Inf; % SD
            EEG.cc.d0501__remBad.amplitudeThreshold  =  222; % uV
            EEG.cc.d0501__remBad.pointSpreadWidth    = 2000; % samples (use 2000 or 4000)
            evalc("EEG = trimOutlier(EEG,EEG.cc.d0501__remBad.channelSdLowerBound,EEG.cc.d0501__remBad.channelSdUpperBound,EEG.cc.d0501__remBad.amplitudeThreshold,EEG.cc.d0501__remBad.pointSpreadWidth);");
            if EEG.pnts == 0
               error('All Datapoints deemed outlier with second trimOutlier')
            end
            
            %% Cleaning: Correction of occular artefacts using ICA & Adjust
            % find rank of ICA through function or the number of channels
            EEG.cc.d0601__pcaICA.rankData = rank(double(EEG.data'));
            EEG.cc.d0601__pcaICA.rankChan = length(EEG.cc.d0201__remBad.chanlocs_006_BadChansRem)-1; % Dreszer et al forgot to add -1 (to correct for AV)
            EEG.cc.d0601__pcaICA.rankUsed = min([EEG.cc.d0601__pcaICA.rankData,EEG.cc.d0601__pcaICA.rankChan]);
            % perform ICA
            evalc("EEG = pop_runica(EEG,'icatype','runica','pca',EEG.cc.d0601__pcaICA.rankUsed,'extended',1,'interupt','off');");
            % Identify bad components with ADJUST
            evalc("[badIC , ~] = ADJUST(EEG)");
            % Remove bad Components
            EEG.cc.clean_ICA_Mask = ones(size(EEG.icaact,1), 1);
            evalc("EEG = pop_subcomp(EEG, badIC, 0);");
            EEG.cc.clean_ICA_Mask(badIC) = 0;
            
            %% Cleaning: Interpolate bad Channels
            %% Saving with same Channel number across all
            evalc("EEG = pop_select( EEG, 'nochannel', {'VOGabove', 'HOGl'});");
            % check that all relevant channels are in the set
            EEG.chaninfo = rmfield(EEG.chaninfo, 'removedchans');
            EEG.chaninfo = rmfield(EEG.chaninfo, 'nodatchans');
            EEG = eeg_checkset( EEG );
            ToInterpolate = []; % necessary here?
            if ~all( contains(upper(AllChannels), upper({EEG.chanlocs.labels})))
                % if not, interpolate
                missing = find(~(contains(upper(AllChannels), upper({EEG.chanlocs.labels}))));
                nrChan = EEG.nbchan;
                ToInterpolate = [];
                for i = missing
                    nrChan = nrChan+1;
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


            
            %% Cleaning: Bad Segments using trimOutlier III/III
            EEG.cc.d0701__adjust.channelSdLowerBound = -Inf; % SD
            EEG.cc.d0701__adjust.channelSdUpperBound =  Inf; % SD
            EEG.cc.d0701__adjust.amplitudeThreshold  =  111; % uV (111)
            EEG.cc.d0701__adjust.pointSpreadWidth    = 2000; % samples (use 4000 or 2000)
            evalc("EEG = trimOutlier(EEG, EEG.cc.d0701__adjust.channelSdLowerBound, EEG.cc.d0701__adjust.channelSdUpperBound, EEG.cc.d0701__adjust.amplitudeThreshold, EEG.cc.d0701__adjust.pointSpreadWidth);");
            
            if EEG.pnts == 0
               error('All Datapoints deemed outlier with third trimOutlier')
            end
            
            evalc("EEG = pop_saveset(EEG, 'filename', FileName, 'filepath', char(dir_Preproc), 'savemode', ['onefile']);");
            
        catch e
            ErrorMessage = string(e.message);
            for ierrors = 1:length(e.stack)
                ErrorMessage = strcat(ErrorMessage, "//", e.stack(ierrors).name, ", Line: ",  num2str(e.stack(ierrors).line));
            end
            
            fprintf('"Subject= %s Problem with one Condition: %s.\n', FileName, ErrorMessage);
            
            ErrorFile = strcat(dir_Log, 'Error_PreProc', '_', strrep(FileName, '.set', '.txt' ));
            fid1 = fopen( ErrorFile, 'wt' );
            fprintf(fid1, 'Error-Subject= %s \n The returned Error Message is: \n  \n %s \n', FileName,  ErrorMessage);
            fclose(fid1);
            
        end
        fprintf('Finished Subject %s.\n',Raw_Files(i_Sub).name);
    end
    
    

    
    
catch e
    % If error ocurrs, create ErrorMessage(concatenated for all nested errors).
    ErrorMessage = string(e.message);
    for ierrors = 1:length(e.stack)
        ErrorMessage = strcat(ErrorMessage, "//", e.stack(ierrors).name, ", Line: ",  num2str(e.stack(ierrors).line));
    end
    % make error log
    ErrorFile =strsplit(InputFile, "/");
    ErrorFile = ErrorFile{end};
    fprintf('Problem with executing File %s. \n',ErrorFile);
    fprintf('The Error Message is: \n %s \n',ErrorMessage);
    ErrorFile = strcat(dir_Log, 'Error_PreProc', ErrorFile, '_', Cond_FileName , '.txt' );
    ErrorFile = strrep(ErrorFile, '.set', '');
    fid1 = fopen( ErrorFile, 'wt' );
    fprintf(fid1, 'Error-Subject= %s \n The returned Error Message is: \n  \n %s \n', InputFile,  ErrorMessage);
    fclose(fid1);
    
end
end