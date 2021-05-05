function crdm_eeg_preprocess( sid, options )
%CRDM_EEG_PREPROCESS Preprocess raw EEG data for one subject of the CRDM
%task (Maria's version).
%   IN:     sid         - subject id, e.g. 64
%           options     - the struct that holds all analysis options
%   OUT:    -


[details, paths] =  crdm_eeg_subjects(sid, options);

%{
forceEEGLABConvert = 0;
forceSPMConvert = 0;
forceDownsample = 0;
forceRereference = 0;
forceFilter = 0;
forceEyeblink = 0;
forceArtefactRejection = 0;
forceNanArtefacts = 0;
%}

% Loop through all sessions for this subject id
for session = 1:length(details.preproc.sessionIDs)
    sessionID = details.preproc.sessionIDs(session);
	
    % step 1: transform raw data to set files that can be read by spm
    if forceEEGLABConvert
        eeglab;
        create_set_file(paths.(options.preproc.reference).preprocessing.rawData.EEG(session).sessionList, ...
            paths.(options.preproc.reference).preprocessing.EEGlab.EEG(sessionID).sessionList);
        fprintf('\nCreated raw set file.\n\n');
    end
    
    % step 2: read data into SPM 
    cd(paths.preprocessed.SPM)
    if forceSPMConvert
        D{sessionID} = convert_eeglab_SPM(paths.preprocessed.setFile.path, ...
            paths.preprocessing.SPMFiles.converted(sessionID).sessionList, ...
            details.setFiles(sessionID).names, options.scriptdir);
        fprintf('\nConverted to SPM format.\n\n');
    else
        D{sessionID} = spm_eeg_load(details.preproc.convertedFiles(sessionID).names);
    end

    % step 3: downsample
    if forceDownsample
        D{sessionID} =  downsample_EEG(D{sessionID},options);
        fprintf('\nDownsampled.\n\n');
    else
        D{sessionID} = spm_eeg_load(details.preproc.downsampledFiles(sessionID).names);
    end

    % step 4: downsample
    if forceRereference
        D{sessionID} = rereference_EEG(D{sessionID},options.preproc.reference);
        fprintf('\nRereferencing done.\n\n');
    else
        D{sessionID} = spm_eeg_load(details.preproc.referencedFiles(sessionID).names);
    end

    % step 5: bandpass filter
    if forceFilter
        D{sessionID} = filter_EEG(D{sessionID},options);
        fprintf('\nFiltering done.\n\n');
    else
        D{sessionID} = spm_eeg_load(details.preproc.filteredFiles(sessionID).names);

        if isfield (D{sessionID}, 'sconfounds')
            D{sessionID} = rmfield(D{sessionID}, 'sconfounds');
        end
    end

    % step 6: eyeblink correction
    if forceEyeblink
        %{
        % this is most likely not needed
        % make sure we have some 3D sensor locations
        [ D{sessionID} ] = tnueeg_sensor_locations( D{sessionID} );
        unNew = repmat({'uV'}, 1, 64);
        D{sessionID}= units(D{sessionID}, 1:64 , unNew);
        %}
        
        %-- eye blink detection -----------------------------------------------%
        % subject-specific EB detection thresholds
        options.preproc.eyeblinkthreshold = details.preproc.eyeblink.threshold(1,sessionID);
        options.preproc.windowForEyeblinkdetection = details.preproc.windowEyblinkdetection(:, session);

        [Dm, trialStats.numEyeblinks] = tnueeg_eyeblink_detection_spm(D{sessionID}, options);
        savefig(details.preproc.ebdetectfig(sessionID).names);
        fprintf('\nEye blink detection done.\n\n');

        %-- eye blink epoching and confounds ----------------------------------%
        Db = tnueeg_epoch_blinks(Dm, options);
        fprintf('\nEpoching to eye blinks done.\n\n');

        % preclean the EB epochs to get better confound estimates
        % Dbc = dprst_preclean_eyeblink_epochs(Db, details, options);

        Dbc = tnueeg_get_spatial_confounds(Db, options);
        savefig(details.preproc.ebspatialfig(sessionID).names);
        fprintf('\nSpatial confounds done.\n\n');

        %-- eye blink correction ----------------------------------------------%
        % if options.preproc.preclean.doBadChannels
        %     badPrecleanChannels = badchannels(Dbc);
        %     if ~isempty(badPrecleanChannels)
        %         fprintf('\nMarking %s channels as bad due to precleaning.', ...
        %             num2str(numel(badPrecleanChannels)));
        %         D = badchannels(D, badPrecleanChannels, ones(1, numel(badPrecleanChannels)));
        %     end
        % end


        if sid == 19 && (sessionID == 5 || sessionID == 2)
            Dbc = rmfield(Dbc,'sconfounds');
            newsconfounds.label = chanlabels(D{sessionID-1});
            newsconfounds.coeff = D{sessionID-1}.sconfounds;
            newsconfounds.bad = zeros(61,1);
            %
            Dbc = sconfounds(Dbc, newsconfounds);
            save(Dbc);
            
        elseif sid == 21 && (sessionID == 1 || sessionID == 2)
            Dbc = rmfield(Dbc,'sconfounds');
            newsconfounds.label = chanlabels(D{3});
            newsconfounds.coeff = D{3}.sconfounds;
            newsconfounds.bad = zeros(61,1);
            %
            Dbc = sconfounds(Dbc, newsconfounds);
            save(Dbc);

        elseif sid == 24 && sessionID == 1
            Dbc = rmfield(Dbc,'sconfounds');
            newsconfounds.label = chanlabels(D{2});
            newsconfounds.coeff = D{2}.sconfounds;
            newsconfounds.bad = zeros(61,1);
            %
            Dbc = sconfounds(Dbc, newsconfounds);
            save(Dbc);

        elseif sid == 27 && sessionID == 4
            Dbc = rmfield(Dbc,'sconfounds');
            newsconfounds.label = chanlabels(D{3});
            newsconfounds.coeff = D{3}.sconfounds;
            newsconfounds.bad = zeros(61,1);
            %
            Dbc = sconfounds(Dbc, newsconfounds);
            save(Dbc);

        elseif sid == 33 && (sessionID == 2 || sessionID == 6)
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

    % step 7: mark additional (remaining) artefacts as NaN
    if forceNanArtefacts
        D{sessionID} = crdm_eeg_reject_remaining_artefacts(D{sessionID}, options);
        make_copy_for_fieldtrip(D{sessionID},0);
        fprintf('\nRejected additional artefacts.\n\n');
    end
    
    % step 8: CSD transformation
    if options.preproc.csdFlag
        D{sessionID} = csd_transform_EEG(D{sessionID});
        make_copy_for_fieldtrip(D{sessionID}, options.preproc.csdFlag);
        fprintf('\nCSD transformation done.\n\n');
    end
end
cd(options.scriptdir)

end