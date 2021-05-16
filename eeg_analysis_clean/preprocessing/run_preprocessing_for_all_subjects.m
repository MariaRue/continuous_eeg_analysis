%subjectList = [16, 18:21, 24, 26, 27, 28, 29, 31, 32, 33, 34, 35, 39, 40, 41, 42, 43, 47, 50, 51, 52, 54, 55, 57, 58];
subjectList = [62:64,66,68,70];

reference = 'LMRM';
options = continuous_RDK_set_options('iMac');
eeglab;

forceEEGLABConvert = 0;
forceSPMConvert = 0;
forceDownsample = 0;
forceRereference =0;
forceFilter = 0;
forceArtefactRejection = 0;
forceNanArtefacts = 0;

csdFlag = 1; 

% transform raw data to set files that can be read by spm
for subject = 1:length(subjectList)
    subID = subjectList(subject);
    [details,paths] =  conrdk_subjects(subID,options,reference,0);
    
    
    if forceEEGLABConvert
        for session = 1:length(details.preproc.sessionIDs)
            sessionID = details.sessionIDs(session);
            
            create_set_file(paths.(reference).preprocessing.rawData.EEG(session).sessionList,paths.(reference).preprocessing.EEGlab.EEG(sessionID).sessionList);
            
        end
    end
    
    for session = 1%:length(details.sessionIDs)
        sessionID = details.preproc.sessionIDs(session);
        
        cd(paths.preprocessed.SPM)
        
        if forceSPMConvert
            
            D{session} = convert_eeglab_SPM(paths.preprocessed.setFile.path,paths.preprocessing.SPMFiles.converted(sessionID).sessionList,details.setFiles(sessionID).names,options.scriptdir);
            
        else
            
            D{session} = spm_eeg_load(details.preproc.convertedFiles(session).names);
            
            
        end
        
        
        
        if forceDownsample
            
            D{session} =  downsample_EEG(D{session},options);
            
        else
            
            D{session} = spm_eeg_load(details.preproc.downsampledFiles(session).names);
            
        end
        
        
        if forceRereference
            
            D{session} = rereference_EEG(D{session},reference);
        else
            D{session} = spm_eeg_load(details.preproc.referencedFiles(session).names);
            
        end
        
        if forceFilter
        
        D{session} = filter_EEG(D{session},options);
        
        else 
          D{session} = spm_eeg_load(details.preproc.filteredFiles(session).names);   
            
        end
        
        if forceArtefactRejection
        D{session} = reject_eyeblinks_artefacts(D{session},subID,session,options,details.preproc.eyeblink.threshold(1,sessionID));
        else 
            
            
        end 
        
        
        if forceNanArtefacts
        
        make_copy_for_fieldtrip(D{session},0);
        
        else 
            
        end 
        
        
       if csdFlag
           
        D{session} = spm_eeg_load(details.preproc.eyeCorrection(session).names);      
        D{session} = csd_transform_EEG(D{session});
        
        make_copy_for_fieldtrip(D{session},1);

        end 
    end
    cd (options.scriptdir)
end


%