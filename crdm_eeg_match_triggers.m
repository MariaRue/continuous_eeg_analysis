function crdm_eeg_match_triggers( subID, options)
%CRDM_EEG_MATCH_TRIGGERS Tries to find the optimal match of triggers in the
% EEG data and triggers recorded with the behavioural data for each
% stimulus

[details, paths] = crdm_eeg_subjects(subID, options);
    
for sessionCount = 1: length(details.sessionIDs)

    session = details.sessionIDs(sessionCount);
    D = spm_eeg_load(...
        paths.(options.preproc.reference).continuousPreproc(sessionCount).sessionList);

    % load stimulus triggers
    load(paths.behaviour(sessionCount).sessionList, 'B');
    bhvTriggersAllBlocks = B.trigger_vals;

    % load sequence of conditions for a session
    load(paths.stimulus(sessionCount).sessionList, 'S');
    conditionID = cellfun(@str2double, S.block_ID_cells);

    % find start and end of each block in the eeg events
    eeg_events = D.events; % all events
    eeg_events = eeg_events(strmatch('trigger', {eeg_events.type})); % only trigger events

    sob = find([eeg_events.value]== 11); % index of start of block trigger
    eob = find([eeg_events.value]==210); % index of end of block trigger

    if length(sob) < 4
        % if first sob is missing, remove first eob
        blocksIDs = [2 3 4];
        eob(1) = [];
    else
        blocksIDs = [1 2 3 4];
    end
    
    % sanity check 
    if numel(sob) ~= numel(eob)
        error('Unequal number of blocks');
    end

    for iBlock = 1: length(blocksIDs) % loop through blocks
        
        eegTriggers = eeg_events(sob(iBlock) : eob(iBlock));
        bhvTriggers = bhvTriggersAllBlocks{blocksIDs(iBlock)};

        [~, flagBlockMatched] = ...
            crdm_eeg_compare_eeg_behav_triggers(eegTriggers, bhvTriggers);

        if flagBlockMatched

            [blockStartEEGIdx, blockEndEEGIdx] = ...
                match_eeg_data_with_behavioural_triggers(...
                eegTriggers, bhvTriggers, D);

            % select EEG data and corresponding artefacts that matches stimulus
            EEGDat{session}{conditionID(blocksIDs(iBlock))} = ...
                D(:, blockStartEEGIdx:blockEndEEGIdx, 1);
            badSamples{session}{conditionID(blocksIDs(iBlock))} = ...
                D.badsamples(:, blockStartEEGIdx:blockEndEEGIdx, 1);

        else
            % don't save anything because EEG data doesn't match
            % stimulus
            EEGDat {session} {conditionID(blocksIDs(iBlock))} = {};
            badSamples{session} {conditionID(blocksIDs(iBlock))} = {};

            fprintf('Haven''t matched session %0.0f, block %0.0f\n', ...
                session, blocksIDs(iBlock));

        end % block matched?

    end % blocks

end %session

%%%%%%%
save(paths.(options.preproc.reference).matchedEEG.saveName, 'EEGDat', 'badSamples');
    

end
