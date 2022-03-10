%subjectList = [ 16, 18:21, 24, 26, 27, 28, 29, 31, 32, 33, 34, 35, 39, 40, 41, 42, 43, 47, 50, 51, 52, 54, 55, 57, 58];
subjectList = [ 16, 18:21, 24, 26, 27, 28, 29, 31, 33, 34, 35, 42, 43, 47, 50, 51, 52, 54, 55, 57, 58];
%subjectList = [ 19,21,24,26,27,31,33,35,39,40,41,47,50,55];
subjectList = [62:64,66,68,70]; % vertical motion 





reference = 'LMRM';
options = continuous_RDK_set_options('iMac');


forceEEGLABConvert = 0;
forceSPMConvert = 0;
forceDownsample = 0;
forceRereference = 0;
forceFilter = 0;
forceEyeblink = 0;
forceArtefactRejection = 0;
forceNanArtefacts = 0;

csdFlag = 1;

% transform raw data to set files that can be read by spm
for subject = 1:length(subjectList)
    subID = subjectList(subject);
    [details,paths] =  conrdk_subjects(subID,options,reference,0);
    
    
    
    if forceEEGLABConvert
        eeglab;
        for session = 1:length(details.preproc.sessionIDs)
            sessionID = details.preproc.sessionIDs(session);
            
            create_set_file(paths.(reference).preprocessing.rawData.EEG(session).sessionList,paths.(reference).preprocessing.EEGlab.EEG(sessionID).sessionList);
            
        end
        
    end
    %
    for session = 1:length(details.preproc.sessionIDs)
        sessionID = details.preproc.sessionIDs(session);
        
        options.preproc.windowForEyeblinkdetection = details.preproc.windowEyblinkdetection(:,session);
        
        
        
        cd(paths.preprocessed.SPM)
        %
        if forceSPMConvert
            
            D{sessionID} = convert_eeglab_SPM(paths.preprocessed.setFile.path,paths.preprocessing.SPMFiles.converted(sessionID).sessionList,details.setFiles(sessionID).names,options.scriptdir);
            
        else
            
            D{sessionID} = spm_eeg_load(details.preproc.convertedFiles(sessionID).names);
            
            
        end
        
        
        
        if forceDownsample
            
            D{sessionID} =  downsample_EEG(D{sessionID},options);
            
        else
            
            D{sessionID} = spm_eeg_load(details.preproc.downsampledFiles(sessionID).names);
            
        end
        
        
        if forceRereference
            
            D{sessionID} = rereference_EEG(D{sessionID},reference);
        else
            D{sessionID} = spm_eeg_load(details.preproc.referencedFiles(sessionID).names);
            
        end
        
        if forceFilter
            
            D{sessionID} = filter_EEG(D{sessionID},options);
            
        else
            D{sessionID} = spm_eeg_load(details.preproc.filteredFiles(sessionID).names);
            
            if isfield (D{sessionID},'sconfounds')
                D{sessionID} = rmfield(D{sessionID},'sconfounds');
            end
        end
        
        
        if forceEyeblink
            
            [ D{sessionID} ] = tnueeg_sensor_locations( D{sessionID} );
            
            unNew = repmat({'uV'}, 1, 64);
            D{sessionID}= units(D{sessionID}, 1:64 , unNew);
            %-- eye blink detection ----------------------------------------------%
            % subject-specific EB detection thresholds
            options.preproc.eyeblinkthreshold = details.preproc.eyeblink.threshold(1,sessionID);
            
            [Dm, trialStats.numEyeblinks] = tnueeg_eyeblink_detection_spm(D{sessionID}, options);
            savefig(details.preproc.ebdetectfig(sessionID).names);
            fprintf('\nEye blink detection done.\n\n');
            
            %-- eye blink epoching and confounds ----------------------------------------------------------%
            Db = tnueeg_epoch_blinks(Dm, options);
            fprintf('\nEpoching to eye blinks done.\n\n');
            
            % preclean the EB epochs to get better confound estimates
            % Dbc = dprst_preclean_eyeblink_epochs(Db, details, options);
            
            Dbc = tnueeg_get_spatial_confounds(Db, options);
            savefig(details.preproc.ebspatialfig(sessionID).names);
            fprintf('\nSpatial confounds done.\n\n');
            
            %-- eye blink correction & headmodel ---------------------------------%
            % if options.preproc.preclean.doBadChannels
            %     badPrecleanChannels = badchannels(Dbc);
            %     if ~isempty(badPrecleanChannels)
            %         fprintf('\nMarking %s channels as bad due to precleaning.', ...
            %             num2str(numel(badPrecleanChannels)));
            %         D = badchannels(D, badPrecleanChannels, ones(1, numel(badPrecleanChannels)));
            %     end
            % end
            
            
            if subID == 19 && (sessionID == 5 || sessionID == 2)
                
                
                Dbc = rmfield(Dbc,'sconfounds');
                newsconfounds.label = chanlabels(D{sessionID-1});
                newsconfounds.coeff = D{sessionID-1}.sconfounds;
                newsconfounds.bad = zeros(61,1);
                %
                Dbc = sconfounds(Dbc, newsconfounds);
                save(Dbc);
            elseif subID == 21 && (sessionID == 1 || sessionID == 2)
                Dbc = rmfield(Dbc,'sconfounds');
                newsconfounds.label = chanlabels(D{3});
                newsconfounds.coeff = D{3}.sconfounds;
                newsconfounds.bad = zeros(61,1);
                %
                Dbc = sconfounds(Dbc, newsconfounds);
                save(Dbc);
                
            elseif subID == 24 && sessionID == 1
                Dbc = rmfield(Dbc,'sconfounds');
                newsconfounds.label = chanlabels(D{2});
                newsconfounds.coeff = D{2}.sconfounds;
                newsconfounds.bad = zeros(61,1);
                %
                Dbc = sconfounds(Dbc, newsconfounds);
                save(Dbc);
                
            elseif subID == 27 && sessionID == 4
                Dbc = rmfield(Dbc,'sconfounds');
                newsconfounds.label = chanlabels(D{3});
                newsconfounds.coeff = D{3}.sconfounds;
                newsconfounds.bad = zeros(61,1);
                %
                Dbc = sconfounds(Dbc, newsconfounds);
                save(Dbc);
                
            elseif subID == 33 && (sessionID == 2 || sessionID == 6)
                Dbc = rmfield(Dbc,'sconfounds');
                newsconfounds.label = chanlabels(D{1});
                newsconfounds.coeff = D{1}.sconfounds;
                newsconfounds.bad = zeros(61,1);
                %
                Dbc = sconfounds(Dbc, newsconfounds);
                save(Dbc);
                
            end
            
            D{sessionID} = tnueeg_add_spatial_confounds(D{sessionID}, Dbc, options);
            fprintf('\nAdding spatial confounds done.\n\n');
            
            
            
            D{sessionID} = tnueeg_eyeblink_correction_simple(D{sessionID}, options);
            fprintf('\nSimple eye blink correction done.\n\n');
        else
            
            D{sessionID} = spm_eeg_load(details.preproc.eyeCorrection(sessionID).names);
        end
        
        
        if forceNanArtefacts
            
            D{sessionID} = reject_eyeblinks_artefacts(D{sessionID},subID,sessionID,options,details.preproc.eyeblink.threshold(1,sessionID));
            
            make_copy_for_fieldtrip(D{sessionID},0);
            
            
            
        end
        %
        %
        if csdFlag
            
            
            D{sessionID} = csd_transform_EEG(D{sessionID});
            
            make_copy_for_fieldtrip(D{sessionID},csdFlag);
            
        end
    end
    cd (options.scriptdir)
end


%