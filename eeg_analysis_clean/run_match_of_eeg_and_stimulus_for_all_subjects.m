% this script runs across subjects and tries to find the optimal match of
% EEG samples and triggers recorded with the behavioural data for each
% stimulus
options = continuous_RDK_set_options('iMac');

subjectList = [16, 18:21, 24, 26, 27, 28, 29, 31, 33, 34, 35,  42, 43, 47, 50, 51, 52, 54, 55, 57, 58];
%subjectList = [ 42, 43, 47, 50, 51, 52, 54, 55, 57, 58];
%subjectList = [62:64,66,68,70];

csdFlag = 0; % 1 for csd transformed data
reference = 'LMRM';
for subject = 1:length(subjectList)
    
    subID = subjectList(subject);
    disp('subject: ')
    disp(subID)
    [details,paths] =  conrdk_subjects( subID,options,reference,csdFlag);
    
    for sessionCount = 1:length(details.sessionIDs)
        
        session = details.sessionIDs(sessionCount);
        
   
        D = spm_eeg_load(paths.(reference).continuousPreproc(sessionCount).sessionList);
       
        D.trialonset

        % load stimulus triggers
        bhv = load(paths.behaviour(sessionCount).sessionList,'B');
        bhvTriggersAllBlocks = bhv.B.trigger_vals;
       
        % load sequence of conditions for a session
        stim = load(paths.stimulus(sessionCount).sessionList,'S');
        conditionID = cellfun(@str2double,stim.S.block_ID_cells);
       
        
        % find start and end of each block in the eeg events
        eeg_events = D.events; %events in EEG data
        eeg_events = eeg_events(strmatch('trigger', {eeg_events.type})); %triggers in EEG data
        
        eob = find([eeg_events.value]==210); %end of block trigger
        sob = find([eeg_events.value]== 11); % start of block
        
        % if first sob is missing
        if length(sob) < 4
            
            nBlocks = [2 3 4];
            eob(1) = [];
         
        else

            nBlocks = [1 2 3 4];
            
        end
        
        for block = 1:length(nBlocks) % loop through blocks
            
            [eegTriggerVals,flagBlockMatched]=compare_eeg_and_behavioural_triggers(eeg_events(sob(block):eob(block)),bhvTriggersAllBlocks{nBlocks(block)});
            
            if flagBlockMatched
                
                [blockStartEEGIdx,blockEndEEGIdx] = match_eeg_data_with_behavioural_triggers(eeg_events(sob(block):eob(block)),bhvTriggersAllBlocks{nBlocks(block)},D);
                
                % select EEG data and corresponding artefacts that matches stimulus
                EEGDat{session}{conditionID(nBlocks(block))} = D(:,blockStartEEGIdx:blockEndEEGIdx,1);
                badSamples{session}{conditionID(nBlocks(block))} = D.badsamples(:,blockStartEEGIdx:blockEndEEGIdx,1);
             
            else
                % don't save anything because EEG data doesn't match
                % stimulus
                EEGDat {session} {conditionID(nBlocks(block))} = {};
                badSamples{session} {conditionID(nBlocks(block))} = {};
                
                fprintf('Haven''t matched session %0.0f, block %0.0f\n',session,nBlocks(block));
                
            end % block matched?
            
        end % blocks
        
    end %session
    
    %%%%%%%
    save(paths.(reference).matchedEEG.saveName,'EEGDat','badSamples');
    
end % subject

