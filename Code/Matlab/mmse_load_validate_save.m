function x = mmse_load_validate_save(file,elecAn,channelset,FileName,outputfolder)     
    
     
    numberSamples = size(file.data,2);
    disp(['CC: number of samples: ',num2str(numberSamples)]); 
    
    
    disp(['CC: ',datestr(now),': starting']);        
        
    % checks minimal length of data and whether
    % all required electrodes are present for the current subject
    labelsTemp = lower({file.chanlocs.labels});
    labelsIndex = find(ismember(labelsTemp,elecAn));
    if length(labelsIndex) ~= length(elecAn)
      msg = 'Electrodes handed over with elecAn are not equal to electrodes in dataset';
      error(msg);
    end        
    disp(['CC: electrodes indices: ',num2str(labelsIndex)]);                
        
    % actual MMSE computation
    EEGData = file.data(labelsIndex,1:numberSamples)';
    file.MMSE.MMSEValues = cc_mmse_prepare_computation(EEGData);

    disp(['CC: ',datestr(now),': finished']);    
        
    x=struct;
    x.name=['MMSE','_',channelset,'_',FileName];
    x.mmse=file.MMSE.MMSEValues(:,1)'; 
    filename_temp=[x.name,'.csv'];
    writetable(struct2table(x),[outputfolder,'/',filename_temp],'Delimiter','\t'); %exports table to the current folder
    
end