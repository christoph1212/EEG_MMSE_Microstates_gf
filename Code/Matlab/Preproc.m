function Preproc(dir_Raw, dir_Root, dir_Log, Overwrite)
%% Preprocessing EEG files
% Run preprocessing 
% Output:
%   (1) a .set file of the preprocessed data, one per subject and condition
%       and run (six files per subject)
%   (2) a .csv Log-file of the preprocessing
%   (3) a potential .txt Error-file
%
% Inputs:
%   dir_Raw:    String pointing to the folder where Raw Data is
%   dir_Root:   String pointing to the project's parent folder
%   dir_Log:    String pointing to where Log-Files should be saved
%   Overwrite:  Numeric (0|1). Should calculation be recalculated? 
%               Default: 0
%
% This script was created by: Name (Date), and is based on the
% preprocessing pipeline by Dreszer et al. (2020):
% [https://github.com/IS-UMK/complexity/tree/master/Preprocessing]
%
% This script was adapted by: Christoph Fruehlinger (May 2025)

%% get from function input
if nargin < 4
    Overwrite = 0;
end

fprintf("**********************\nStarting Preprocessing\n**********************\n")

%% Prepare List of Files to be Processed
Raw_Files = dir(fullfile(dir_Raw, '**/*.set'));  % get list of files in *.set format in any subfolder

%% Directory where file should be saved
dir_Preproc = strcat(dir_Root, 'Data/Preprocessed/');
dir_Log = strcat(dir_Log, 'Preproc/');

if ~isfolder(dir_Preproc)
    mkdir(dir_Preproc)
end

if ~isfolder(dir_Log)
    mkdir(dir_Log)
end

fprintf('\nPreprocessing %d Files. \n', length(Raw_Files));
fprintf('InputFolder is %s. \nOutputFolder is %s. \nLogFolder is %s. \n\n', dir_Raw , dir_Preproc, dir_Log);

%% Infos on Triggers and Conditions
SplitStruct = struct('Trigger', {11, 12 21 22 31 32}, ...
    'Condition', {'first_run_eyes_open' 'first_run_eyes_closed' 'second_run_eyes_open' 'second_run_eyes_closed' 'third_run_eyes_open' 'third_run_eyes_closed'});

%% Increase calculation speed by running multiple subjects in parallel
delete(gcp('nocreate')); % make sure that previous pooling is closed
parpool("Processes");

FileName = "";
InputFile = "";
Cond_FileName = "";

Files_PreProc = dir(dir_Preproc);

%% Looped preprocessing
parfor i_Sub = 1:length(Raw_Files)

    % Check if Subject has been preprocessed 
    if (sum(contains({Files_PreProc.name},Raw_Files(i_Sub).name)) == 2) && (Overwrite == 0)
        continue
    end
    % ignore subs with unknown reference
    if contains(Raw_Files(i_Sub).name, "sub-SS08EL29") || contains(Raw_Files(i_Sub).name, "sub-EL07EL18")
        continue
    end

    run_silent(Overwrite, Raw_Files, dir_Preproc, dir_Log, SplitStruct, FileName, InputFile, Cond_FileName, i_Sub);

end

fprintf("**********************\nFinished Preprocessing\n**********************\n")

end

function run_silent(Overwrite, Raw_Files, dir_Preproc, dir_Log, SplitStruct, FileName, InputFile, Cond_FileName, i_Sub)

try

    fprintf('Standardizing File: %s; Select Channels, Reference, Sampling Rate. \n', Raw_Files(i_Sub).name);
    
    % Initate variables, and load EEG lab File
    EEG = struct([]);
    temp_out = [];
    InputFile = [Raw_Files(i_Sub).folder,'/',Raw_Files(i_Sub).name];
    Cond_FileName = "General";
    evalc("EEG = pop_loadset(InputFile);");
    
    %% Standardize: remove unwanted channels (some Labs recorded with some extra electrodes that should not be included)
    % AFZ and FPZ are used as grounds in some labs
    Common_Channels =  {'FP1', 'FP2', 'AF7', 'AF8', 'AF3', 'AF4', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', ...
        'F7', 'F8', 'FT7', 'FT8', 'FC5', 'FC6', 'FC3', 'FC4', 'FC1', 'FC2', 'C1', 'C2', 'C3', 'C4', ...
        'C5', 'C6', 'T7', 'T8', 'TP7', 'TP8', 'CP5', 'CP6', 'CP3', 'CP4', 'CP1', 'CP2', 'P1', 'P2', ...
        'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'PO7', 'PO8', 'PO3', 'PO4', 'O1', 'O2', 'OZ', 'POZ', 'PZ', ...
        'CPZ', 'CZ', 'FCZ', 'FZ', 'VOGabove', 'VOGbelow', 'HOGl', 'HOGr','MASTl', 'MASTr'};
    
    Common_Channels = Common_Channels(ismember(Common_Channels, {EEG.chanlocs.labels}));
    evalc("EEG = pop_select( EEG, 'channel', Common_Channels);");
    
    %% Standardize: rereference all Files to FCZ (differs between Labs)
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
    temp_comStr = '% removed unwanted channels';
    EEG = eegh(temp_comStr,EEG);

    % append new reference channel 'FCz' to chanlocs (as no data channel)
    RefInfo = {EEG.nbchan 'labels' 'FCZ' 'theta' 0 'radius' 0.127 'X' 0.388 'Y' 0 'Z' 0.922 ...
        'sph_theta' 0 'sph_phi' 67.2 'sph_radius' 1 'type' 'EEG' 'datachan' 0 };
    evalc("EEG = pop_chanedit(EEG, 'append', RefInfo);");

    % Chanlocs without Ref - save chanlocs for interpolation:
    evalc("EEG_interp = pop_select( EEG, 'chantype','EEG');");
    evalc("EEG_interp = pop_select( EEG_interp, 'nochannel', {'MASTl', 'MASTr'});");

    % some datasets miss OZ
    if length(EEG_interp.chanlocs) ~= 58 
        if length(EEG_interp.chanlocs) == 57 && ~ismember('OZ', upper({EEG_interp.chanlocs.labels}))
            % Add OZ to chanlocs interpolation template
            Oz.labels     = 'OZ';
            Oz.theta      = 180;
            Oz.radius     = 0.5000;
            Oz.X          = -1;
            Oz.Y          = -1.2246e-16;
            Oz.Z          = 0;
            Oz.sph_theta  = -180;
            Oz.sph_phi    = 0;
            Oz.sph_radius = 1;
            Oz.type       = 'EEG';
            Oz.ref        = '';
            Oz.urchan     = [];

            EEG_interp.chanlocs(end+1) = Oz;
            EEG_interp.nbchan = EEG_interp.nbchan + 1;
        else
            msg = ['Less than 57 channels. Number of channels in dataset is ' num2str(length(EEG_interp.chanlocs)) '.'];
            error(msg);
        end
    end

    EEG.cc.d0001__dnSamp.chanlocs_005_no_ref = EEG_interp.chanlocs;
    % Chanloc of Reference
    EEG.cc.d0001__dnSamp.chanlocs_006_ref = EEG.chaninfo.nodatchans;
    
    %% Downsample to 250/256 Hz (Sampling Rate differs between Labs)
    EEG.cc.d0001__dnSamp.dnSampFreq = EEG.srate/2;
    evalc("[EEG,temp_out.comStr] = pop_resample(EEG,EEG.cc.d0001__dnSamp.dnSampFreq);");
    EEG = eegh(temp_out.comStr,EEG);
        
    %% Split into different Conditions (Eyes Open / Eyes Closed)
    % Prepare Matrix identifying all points (to keep it continous, not
    % epoched)
    EEG_Complete = EEG;

    % Select triggers and conditions
    if contains(Raw_Files(i_Sub).name, 'run-1')
        Rel_SplitStruct.Trigger = [SplitStruct(1:2).Trigger];
        Rel_SplitStruct.Condition = {SplitStruct(1:2).Condition};

    elseif contains(Raw_Files(i_Sub).name, 'run-2')
        Rel_SplitStruct.Trigger = [SplitStruct(3:4).Trigger];
        Rel_SplitStruct.Condition = {SplitStruct(3:4).Condition};

    elseif contains(Raw_Files(i_Sub).name, 'run-3')
        Rel_SplitStruct.Trigger = [SplitStruct(5:6).Trigger];
        Rel_SplitStruct.Condition = {SplitStruct(5:6).Condition};
    end
    
    for i_cond = 1:length(Rel_SplitStruct.Trigger)
        try
            Cond_FileName = Rel_SplitStruct.Condition{i_cond};
            FileName = [Rel_SplitStruct.Condition{i_cond},'_',Raw_Files(i_Sub).name];

            % test that file has been calculated yet
            if isfile(fullfile(dir_Preproc, FileName)) &&  Overwrite == 0
                fprintf('Previously finished Subject %s.\n',Raw_Files(i_Sub).name);
                continue
            end

            fprintf('Processing File: %s; Condition: %s\n', Raw_Files(i_Sub).name, Cond_FileName);
            
            % Check if file is not empty
            if EEG.pnts < EEG.srate
                fprintf('Dataset empty (%s: Condition: %s)\n', Raw_Files(i_Sub).name, Cond_FileName);
                
                MissingDataFile = strcat(dir_Log, "Error_NotEnoughData_for_Preproc_", FileName, '.txt');
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
            EEG.cc.d0201__remBad.chanlocs_before_rem = EEG.chanlocs;
            EEG.cc.d0201__remBad.freqlims =  [  0.00, 5.00;  5.00, 40.00   ]; % Hz
            EEG.cc.d0201__remBad.stdthresh = [ -5.00, 5.00; -2.50,  2.50   ]; % SD
            evalc("[EEG,temp_out.allrmchan,temp_out.specdata,temp_out.specfreqs,temp_out.comStr] = pop_rejchanspec(EEG,'freqlims',EEG.cc.d0201__remBad.freqlims,'stdthresh',EEG.cc.d0201__remBad.stdthresh);");
            
            %remember chanlocs wiInthout these removed bad channels (=only good ones)
            EEG.cc.d0201__remBad.chanlocs_006_BadChansRem = EEG.chanlocs;
            % remember chanlocs of these bad channels (=only bad ones)
            % EEG.cc.d0201__remBad.chanlocs_006_BadChans = EEG.cc.d0001__dnSamp.chanlocs_005_no_ref(temp_out.allrmchan);
            EEG.cc.d0201__remBad.removed_chans = {EEG.cc.d0201__remBad.chanlocs_before_rem(temp_out.allrmchan).labels};            
            
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
            % (adapted ADJUST.m file - make sure to use the file from the Github Repository)
            evalc("[badIC , ~] = ADJUST(EEG)");

            % Remove bad Components
            EEG.cc.clean_ICA_Mask = ones(size(EEG.icaact,1), 1);
            evalc("EEG = pop_subcomp(EEG, badIC, 0);");
            EEG.cc.clean_ICA_Mask(badIC) = 0;
            
            %% Cleaning: Interpolate bad Channels
            evalc("EEG = pop_select( EEG, 'chantype','EEG');");
            evalc("EEG = pop_select( EEG, 'nochannel', {'MASTl', 'MASTr'});");
            evalc("EEG = pop_interp(EEG, EEG.cc.d0001__dnSamp.chanlocs_005_no_ref, 'spherical');");           
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
            
            %% Check Channel Number
            if length(EEG.chanlocs) ~= 59
                msg = ['Incorrect Channel Number: ' num2str(length(EEG.chanlocs))];
                error(msg)
            end

            %% Save Data set
            evalc("EEG = pop_saveset(EEG, 'filename', FileName, 'filepath', char(dir_Preproc), 'savemode', ['onefile']);");

            % Save Log Data
            if isempty(EEG.cc.d0201__remBad.removed_chans)
                EEG.cc.d0201__remBad.removed_chans = 0;
            end

            Log_table = table({Raw_Files(i_Sub).name(1:end-4)}, {Cond_FileName}, length(EEG.chanlocs), length(badIC), length(EEG.cc.d0201__remBad.removed_chans), ...
                'VariableNames', {'File Name', 'Condition', 'Number Channels', 'BadICs', 'Interpolated Channels'});

            ID = strsplit(Raw_Files(i_Sub).name, '_');
            ID = ID{1};
            log_filename = strcat(dir_Log, 'Log_', ID, '_', Cond_FileName, '.csv');
            writetable(Log_table, log_filename);
            
        catch e
            % If error ocurrs, create ErrorMessage
            ErrorMessage = string(e.message);
            for ierrors = 1:length(e.stack)
                ErrorMessage = strcat(ErrorMessage, "//", e.stack(ierrors).name, ", Line: ",  num2str(e.stack(ierrors).line));
            end
            
            fprintf('"Error in File: %s; %s.\n', FileName, ErrorMessage);
            
            ErrorFile = strcat(dir_Log, 'Error_PreProc', '_', strrep(FileName, '.set', '.txt' ));
            fid1 = fopen( ErrorFile, 'wt' );
            fprintf(fid1, 'Error-File: %s \nThe returned Error Message is: \n\n%s \n', FileName,  ErrorMessage);
            fclose(fid1);
            
        end
        
    end
    fprintf('\nFinished File: %s.\n\n',Raw_Files(i_Sub).name);
    
catch e
    % If error ocurrs, create ErrorMessage
    ErrorMessage = string(e.message);
    for ierrors = 1:length(e.stack)
        ErrorMessage = strcat(ErrorMessage, "//", e.stack(ierrors).name, ", Line: ",  num2str(e.stack(ierrors).line));
    end
    % make error log
    ErrorFile = strsplit(InputFile, "/");
    ErrorFile = ErrorFile{end};
    fprintf('Problem executing File: %s\n',ErrorFile);
    fprintf('The Error Message is: \n%s \n',ErrorMessage);
    ErrorFile = strcat(dir_Log, 'Error_PreProc', ErrorFile, '_', Cond_FileName , '.txt' );
    ErrorFile = strrep(ErrorFile, '.set', '');
    fid1 = fopen( ErrorFile, 'wt' );
    fprintf(fid1, 'Error-Subject: %s \nThe returned Error Message is: \n\n%s \n', InputFile,  ErrorMessage);
    fclose(fid1);
    
end

end