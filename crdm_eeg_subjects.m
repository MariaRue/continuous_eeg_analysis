function [ details, paths ] =  crdm_eeg_subjects(subID, options )
%CRDM_EEG_SUBJECTS Function that sets all subject-specific options and paths
%   IN:     subID       - subject id, e.g. 64
%           options     - the struct that holds all analysis options
%   OUT:    details     - a struct that holds all options/details for sid
%           paths       - a struct that holds all files/paths for sid


% which sessions do we have for this sid
switch subID
    case 55
        details.sessionIDs = 1:5;
        details.preproc.sessionIDs = 1:5;
    case 40
        details.sessionIDs = [2,3,4,5];
        details.preproc.sessionIDs = 1:6;
    case 47
        details.sessionIDs = [1,2,3,5];
        details.preproc.sessionIDs = 1:6;
    case 42
        details.sessionIDs = [1,2,3,4,6]; % for some reason there is a problem with the triggers in session 5 where the end trigger of a block (210) is missing
        details.preproc.sessionIDs = 1:6;
    case 21
        details.sessionIDs = [1,2,3,4,5,6]; 
        details.preproc.sessionIDs = [3,1,2,4,5,6];
    case 24 
        details.sessionIDs = [1,2,3,4,5,6]; 
        details.preproc.sessionIDs = [2,1,3:6];
    case 70 
        details.sessionIDs = [1,2,3,5,6]; 
        details.preproc.sessionIDs = 1:6;
    otherwise
        details.sessionIDs = 1:6;
        details.preproc.sessionIDs = 1:6;
end

% path to preprocessed data to match EEG with stim stream for continuous GLM
paths.referenceType.LMRM = fullfile(options.path.EEG.analysis, ...
    'preprocessedData', 'EEG', 'LMRM', sprintf('sub%03.0f', subID));
paths.referenceType.averageReference = fullfile(options.path.EEG.subjects, ...
    sprintf('sub%03.0f', subID), 'eeg', 'average_reference');

paths.preprocessed.SPM = fullfile(options.path.EEG.raw, 'preprocessedData', ...
    'EEG', options.preproc.reference, sprintf('sub%03.0f', subID));
if ~exist( paths.preprocessed.SPM,'dir' )
    mkdir( paths.preprocessed.SPM )
end
paths.rawData.SPM = fullfile(options.path.EEG.raw, 'rawData', 'experiment', ...
    sprintf('sub%03.0f', subID), 'eeg');
paths.preprocessed.setFile.path = fullfile(options.path.EEG.raw, ...
    'preprocessedData', 'EEG', 'setFiles');


if options.preproc.csdFlag
    paths.(options.preproc.reference).matchedEEG.saveName = fullfile(options.path.EEG.analysis,'convGLM','matchedEegData','LMRM',[sprintf('sub%03.0f',subID),'_csdEEGdat.mat']);
    paths.(options.preproc.reference).subjectLevelGLM.jumps_absolute.saveName = fullfile(options.path.EEG.analysis,'convGLM','betasGLMData','LMRM',sprintf('sub%03.0f_betas_jumps_absolute_csd.mat',subID));
    paths.(options.preproc.reference).subjectLevelGLM.response.saveName = fullfile(options.path.EEG.analysis,'convGLM','betasGLMData','LMRM',sprintf('sub%03.0f_betas_response_csd.mat',subID));
    paths.(options.preproc.reference).subjectLevelGLM.coherence_responses.saveName = fullfile(options.path.EEG.analysis,'convGLM','betasGLMData','LMRM',sprintf('sub%03.0f_betas_coherence_responses_csd.mat',subID));
    paths.(options.preproc.reference).subjectLevelGLM.vertical_jumps_absolute.saveName = fullfile(options.path.EEG.analysis,'convGLM','betasGLMData','LMRM',sprintf('sub%03.0f_betas_vertical_jumps_absolute_csd.mat',subID));
    
    paths.(options.preproc.reference).singleTrial.appendedData.trialStart = fullfile(options.path.EEG.analysis,'conventionalEEGAnalysis',options.preproc.reference,['csd_trial_start_locked_EEG_dat_',sprintf('sub%03.0f',subID),'.mat']);
    paths.(options.preproc.reference).singleTrial.appendedData.buttonPress = fullfile(options.path.EEG.analysis,'conventionalEEGAnalysis',options.preproc.reference,['csd_response_locked_EEG_dat_',sprintf('sub%03.0f',subID),'.mat']);
else
    paths.(options.preproc.reference).matchedEEG.saveName = fullfile(options.path.EEG.analysis,'convGLM','matchedEegData','LMRM',[sprintf('sub%03.0f',subID),'_EEGdat.mat']);
    paths.(options.preproc.reference).subjectLevelGLM.jumps_absolute.saveName = fullfile(options.path.EEG.analysis,'convGLM','betasGLMData','LMRM',sprintf('sub%03.0f_betas_jumps_absolute.mat',subID));
    paths.(options.preproc.reference).subjectLevelGLM.response.saveName = fullfile(options.path.EEG.analysis,'convGLM','betasGLMData','LMRM',sprintf('sub%03.0f_betas_response.mat',subID));
    
    paths.(options.preproc.reference).subjectLevelGLM.coherence_responses.saveName = fullfile(options.path.EEG.analysis,'convGLM','betasGLMData','LMRM',sprintf('sub%03.0f_betas_coherence_responses.mat',subID));
    paths.(options.preproc.reference).subjectLevelGLM.vertical_jumps_absolute.saveName = fullfile(options.path.EEG.analysis,'convGLM','betasGLMData','LMRM',sprintf('sub%03.0f_betas_vertical_jumps_absolute.mat',subID));
    
    %paths.(options.preproc.reference).singleTrial.appendedData.trialStart = fullfile(options.path.EEG.analysis,'test_April_2020',['trial_start_locked_EEG_dat_',sprintf('sub%03.0f',id),'.mat']);
    %     paths.(options.preproc.reference).singleTrial.appendedData.buttonPress = fullfile(options.path.EEG.analysis,['response_locked_EEG_dat',sprintf('sub%03.0f',id),'.mat']);
    %
    paths.(options.preproc.reference).singleTrial.appendedData.trialStart = fullfile(options.path.EEG.analysis,'conventionalEEGAnalysis',options.preproc.reference,['trial_start_locked_EEG_dat_',sprintf('sub%03.0f',subID),'.mat']);
    %paths.(options.preproc.reference).singleTrial.appendedData.buttonPress = fullfile(options.path.EEG.analysis,'conventionalEEGAnalysis',options.preproc.reference,['response_locked_EEG_dat',sprintf('sub%03.0f',id),'.mat']);
    paths.(options.preproc.reference).singleTrial.appendedData.buttonPress = fullfile('/Volumes/crdkData','conventionalEEGAnalysis',options.preproc.reference,['response_locked_EEG_dat',sprintf('sub%03.0f',subID),'.mat']);    
end



for sessID = 1:length(details.sessionIDs)
    if options.preproc.csdFlag
        paths.(options.preproc.reference).continuousPreproc(sessID).sessionList = fullfile(paths.referenceType.(options.preproc.reference),sprintf('csdcfMdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',subID,details.sessionIDs(sessID)));
        paths.(options.preproc.reference).singleTrial.preproc(sessID).preprocSessionList = fullfile(paths.preprocessed.SPM,sprintf('csdnanart_csdcfMdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',subID,details.sessionIDs(sessID)));
    else
        paths.(options.preproc.reference).continuousPreproc(sessID).sessionList = fullfile(paths.referenceType.(options.preproc.reference),sprintf('cfMdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',subID,details.sessionIDs(sessID)));
        %
        paths.(options.preproc.reference).singleTrial.preproc(sessID).preprocSessionList = fullfile(paths.preprocessed.SPM,sprintf('nanart_cfMdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',subID,details.sessionIDs(sessID)));
        % paths.(options.preproc.reference).singleTrial.preproc(session).preprocSessionList = fullfile(paths.referenceType.(options.preproc.reference),sprintf('fMdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',id,details.sessionIDs(session)));
    end

    % paths.(options.preproc.reference).continuousPreproc(session).sessionList = fullfile(paths.referenceType.averageReference,sprintf('fMdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',id,details.sessionIDs(session)));
    % paths.(options.preproc.reference).csd.continuousPreproc(session).sessionList = fullfile(paths.referenceType.averageReference,sprintf('cfMdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',id,details.sessionIDs(session)));
    paths.behaviour(sessID).sessionList = fullfile(options.path.EEG.subjects,sprintf('sub%03.0f',subID),'behaviour',sprintf('sub%03.0f_sess%03.0f_behav.mat',subID,details.sessionIDs(sessID)));
    paths.stimulus(sessID).sessionList = fullfile(options.path.EEG.subjects,sprintf('sub%03.0f',subID),'stim',sprintf('sub%03.0f_sess%03.0f_stim.mat',subID,details.sessionIDs(sessID)));
end

for sessID = 1: length(details.preproc.sessionIDs)
    details.setFiles(sessID).names = sprintf('sub%03.0f_sess%03.0f_eeg.set',subID,details.preproc.sessionIDs(sessID));
    details.preproc.convertedFiles(sessID).names = sprintf('spmeeg_sub%03.0f_sess%03.0f_eeg.mat',subID,details.preproc.sessionIDs(sessID));
    details.preproc.downsampledFiles(sessID).names = sprintf('dspmeeg_sub%03.0f_sess%03.0f_eeg.mat',subID,details.preproc.sessionIDs(sessID));
    details.preproc.referencedFiles(sessID).names = sprintf('Mdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',subID,details.preproc.sessionIDs(sessID));
    details.preproc.filteredFiles(sessID).names = sprintf('fMdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',subID,details.preproc.sessionIDs(sessID));
    details.preproc.eyeCorrection(sessID).names = sprintf('cfMdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',subID,details.preproc.sessionIDs(sessID));

    %%%%%%
    details.preproc.ebdetectfig(sessID).names = fullfile(paths.preprocessed.SPM,sprintf('sub%03.0f_sess%03.0f_eyeblinkDetection.fig',subID,details.preproc.sessionIDs(sessID)));
    details.preproc.ebspatialfig(sessID).names = fullfile(paths.preprocessed.SPM,sprintf('sub%03.0f_sess%03.0f_eyeblinkConfounds.fig',subID,details.preproc.sessionIDs(sessID)));
    
    paths.(options.preproc.reference).preprocessing.rawData.EEG(sessID).sessionList = fullfile(options.path.EEG.raw,'rawData','experiment',sprintf('sub%03.0f',subID),'eeg',sprintf('sub%03.0f_sess%03.0f_eeg.cdt',subID,details.preproc.sessionIDs(sessID)));
    paths.(options.preproc.reference).preprocessing.EEGlab.EEG(sessID).sessionList = fullfile('/Volumes/crdkData/preprocessedData/EEG/setFiles',sprintf('sub%03.0f_sess%03.0f_eeg.set',subID,details.preproc.sessionIDs(sessID)));
    
    %%%% spm files
    paths.preprocessing.SPMFiles.converted(sessID).sessionList = fullfile(paths.preprocessed.SPM,sprintf('spmeeg_sub%03.0f_sess%03.0f_eeg.mat',subID,details.preproc.sessionIDs(sessID)));
    paths.preprocessing.SPMFiles.downsampled(sessID).sessionList = fullfile(paths.preprocessed.SPM,sprintf('dspmeeg_sub%03.0f_sess%03.0f_eeg.mat',subID,details.preproc.sessionIDs(sessID)));
    paths.preprocessing.SPMFiles.rereferenced(sessID).sessionList = fullfile(paths.preprocessed.SPM,sprintf('Mdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',subID,details.preproc.sessionIDs(sessID)));
    paths.preprocessing.SPMFiles.filtered(sessID).sessionList = fullfile(paths.preprocessed.SPM,sprintf('fMdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',subID,details.preproc.sessionIDs(sessID)));
    paths.preprocessing.SPMFiles.eyeCorrection(sessID).sessionList = fullfile(paths.preprocessed.SPM,sprintf('cfMdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',subID,details.preproc.sessionIDs(sessID)));
    paths.preprocessing.SPMFiles.csd(sessID).sessionList = fullfile(paths.preprocessed.SPM,sprintf('csdfMdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',subID,details.preproc.sessionIDs(sessID)));
end


%%--preprocessing eyeblink thresholds-------------------------------------%

details.preproc.eyeblink.threshold = nan(2,6);
switch subID
    
    case 16
        details.preproc.eyeblink.threshold(1,1) = 4;
        details.preproc.eyeblink.threshold(2,1) = 2;

        % session 2 huge artefact at start - won't allow me to set an appropriate
        % eyeblink threshold - this needs to be removed before the
        % eyeblink detection is run again
        
        % for now
        details.preproc.eyeblink.threshold(1,2) = 0.9;
        details.preproc.eyeblink.threshold(2,[2,5:6]) = 2;
        
        
        details.preproc.eyeblink.threshold(1,3) = 3;
        details.preproc.eyeblink.threshold(2,3:4) = 1.5;
        
        details.preproc.eyeblink.threshold(1,4) = 2;
        details.preproc.eyeblink.threshold(1,5:6) = 1.5; % not sure
        % data towards end quite noise - not sure what to make of this
        
    case 18
        % sess 2 good example of more data thrown out during artefact
        % rejection - is that correct?
        details.preproc.eyeblink.threshold(1,1:6) = 1.5;
        details.preproc.eyeblink.threshold(2,[1:3,5:6]) = 1.5;
        details.preproc.eyeblink.threshold(2,4) = 2;

    case 19
        % session 2/3 eyeblinks weird? are these other artefacts?
        details.preproc.eyeblink.threshold(1,1) = 15;
        details.preproc.eyeblink.threshold(2,1:6) = 1.5;
        
        details.preproc.eyeblink.threshold(1,2) = 5;
        details.preproc.eyeblink.threshold(1,3:6) = 15;
        
    case 20
        details.preproc.eyeblink.threshold(1,1:6) = 2;
        details.preproc.eyeblink.threshold(2,1:2) = 3;
        
        details.preproc.eyeblink.threshold(2,3:4) = 1.5;
        details.preproc.eyeblink.threshold(2,5:6) = 5;
        
    case 21
        details.preproc.eyeblink.threshold(1,2:6) = 4;
        details.preproc.eyeblink.threshold(1,1) = 4;
        
        details.preproc.eyeblink.threshold(2,1:3) = 3;
        details.preproc.eyeblink.threshold(2,4:6) = 2;
        
    case 24
        details.preproc.eyeblink.threshold(1,1:4) = 5;
        
        % with around 10 eyeblinks per minute if reduced we
        % catch each little bump. Right now we might miss one
        % every now and then        
        details.preproc.eyeblink.threshold(1,5:6) = 5;
        details.preproc.eyeblink.threshold(2,1:6) = 2;
        
    case 26
        details.preproc.eyeblink.threshold(1,1:6) = 9;
        details.preproc.eyeblink.threshold(2,1:5) = 2;
        details.preproc.eyeblink.threshold(1,5:6) = 19;
        % sess 2 super weird for threshold for eyeblinks -
        % missing big peaks - why?
        % same for 3 and 4, 5 threshold for normal artefacts for 5
        % ok?
        details.preproc.eyeblink.threshold(2,6) = 3.5;
        
    case 27
        details.preproc.eyeblink.threshold(1,1:6) = 4;
        
    case 28
        details.preproc.eyeblink.threshold(1,1:6) = 3;
        details.preproc.eyeblink.threshold(2,1:6) = 3.5;
        
    case 32
        details.preproc.eyeblink.threshold(1,1) = 3;
        details.preproc.eyeblink.threshold(2,1:6) = 3.5;
        details.preproc.eyeblink.threshold(1,2:6) = 4;
        
    case 33
        details.preproc.eyeblink.threshold(1,3) = 2;
        details.preproc.eyeblink.threshold(1,[1:2,4:6]) = 3.5;
        
    case 34
        details.preproc.eyeblink.threshold(1,1:6) = 2;
        details.preproc.eyeblink.threshold(2,1:6) = 3.5;

    case 35
        details.preproc.eyeblink.threshold(1,1:6) = 4;
        details.preproc.eyeblink.threshold(1,3) = 6;
        details.preproc.eyeblink.threshold(2,1:6) = 3.5;
        % good example where I think we don't catch artefacts in sessions 4
        % or 5
        
    case 39
        details.preproc.eyeblink.threshold(1,1:6) = 2.5;
        
    case 40
        % eyeblink detection - is that correct - maybe the signal is just
        % distorted because of big peaks? - after inspecting channel, I
        % believe that there are much less eyeblinks - but not sure
        details.preproc.eyeblink.threshold(1,1:6) = 4;
        details.preproc.eyeblink.threshold(2,1:6) = 3.5;
        
    case 41
        details.preproc.eyeblink.threshold(1,1:6) = 2;
        details.preproc.eyeblink.threshold(2,1:6) = 3.5;
        
    case 47
        details.preproc.eyeblink.threshold(1,1:6) = 6;
        details.preproc.eyeblink.threshold(2,1:6) = 3.5;
        % allowed of spikes specifically in last sessions - maybe check it
        % out?

    case 50
        details.preproc.eyeblink.threshold(1,1:6) = 6.5;
        
    case 51
        details.preproc.eyeblink.threshold(1,1:6) = 2;
        details.preproc.eyeblink.threshold(2,1:6) = 3.5;
    case 52
        
        details.preproc.eyeblink.threshold(1,1:6) = 2;
        details.preproc.eyeblink.threshold(2,1:6) = 3.5;
        
    case 55
        details.preproc.eyeblink.threshold(1,[3,4,5]) = 5;
        details.preproc.eyeblink.threshold(1,[1,2,6]) = 2;
        details.preproc.eyeblink.threshold(2,1:6) = 3.5;
        
        % high num of eyeblinks per minute... like  over 25 per minute
    otherwise
        details.preproc.eyeblink.threshold(1,1:6) = 3.5;
        details.preproc.eyeblink.threshold(2,1:6) = 3.5;
end

switch subID
    
    case 18
        for i = 1:6
            if i ~= 4 || i ~= 1
                details.preproc.windowEyblinkdetection(1:2,i) = [10 900];
            else
                details.preproc.windowEyblinkdetection(1:2,i) = [1400 2000];
            end
        end
        
        
    case 19
        for i = 1:6
            if i ~= 5 || i ~= 1 || i ~= 2
                details.preproc.windowEyblinkdetection(1:2,i) = [10 900];
            else
                details.preproc.windowEyblinkdetection(1:2,i) = [2500 3000];
            end
        end
        
    case 21
        for i = 1:6
            details.preproc.windowEyblinkdetection(1:2,i) = [1400 2000];
        end
        
    case 26
        for i = 1:6
            if i <= 4
                details.preproc.windowEyblinkdetection(1:2,i) = [10 900];
            else
                details.preproc.windowEyblinkdetection(1:2,i) = [1400 2000];
            end
        end
        
        
    case 24
        for i = 1:6
            if i ~= 5
                details.preproc.windowEyblinkdetection(1:2,i) = [10 900];
            else
                details.preproc.windowEyblinkdetection(1:2,i) = [1400 2000];
            end
        end
        
    case 31
        for i = 1:6
            if i ~= 1
                details.preproc.windowEyblinkdetection(1:2,i) = [10 900];
            else
                details.preproc.windowEyblinkdetection(1:2,i) = [1400 2000];
            end
        end
        
    case 35
        for i = 1:6
            if i ~= 5
                details.preproc.windowEyblinkdetection(1:2,i) = [10 900];
            else
                details.preproc.windowEyblinkdetection(1:2,i) = [1400 2000];
            end
        end
       
    case 39
        for i = 1:6
            if  i == 6
                details.preproc.windowEyblinkdetection(1:2,i) = [1400 2000];
                
            elseif i == 3
                
                details.preproc.windowEyblinkdetection(1:2,i) = [2500 3000];
            else
                
                details.preproc.windowEyblinkdetection(1:2,i) = [10 900];
            end
        end
        
    case 40
        for i = 1:6
            if  i ~= 2
                details.preproc.windowEyblinkdetection(1:2,i) = [1400 2000];
                
            elseif i == 1
                details.preproc.windowEyblinkdetection(1:2,i) = [2500 3000];
                
            else
                
                details.preproc.windowEyblinkdetection(1:2,i) = [10 900];
            end
        end
         
    case 47
        for i = 1:6
            if  i == 2  || i == 3
                details.preproc.windowEyblinkdetection(1:2,i) = [5 800];
                
            elseif i == 1
                
                details.preproc.windowEyblinkdetection(1:2,i) = [2500 3000];
            else
                
                details.preproc.windowEyblinkdetection(1:2,i) = [1400 2000];
            end
        end
        
    case 50
        for i = 1:6
            if  i == 5
                details.preproc.windowEyblinkdetection(1:2,i) = [2500 3000];
            else
                details.preproc.windowEyblinkdetection(1:2,i) = [5 800];
            end
        end
        
        
	case 55
        for i = 1:6
            if  i == 1
                details.preproc.windowEyblinkdetection(1:2,i) = [1400 2000];
            else
                details.preproc.windowEyblinkdetection(1:2,i) = [5 800];
            end
        end
        
    otherwise
        for i = 1:length(details.preproc.sessionIDs)
            details.preproc.windowEyblinkdetection(1:2,i) = [10 900];
        end
end

% bad channels before EB confound estimation
switch subID
    case 19
        for sessID = 1:length(details.preproc.sessionIDs)
            if sessID == 2 
        details.preclean(sessID).badchannels = [30];
            else
             details.preclean(sessID).badchannels = [];    
            end
        end 
                
    otherwise
        details.preclean.badchannels = [];
end



end