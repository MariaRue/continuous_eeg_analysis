function [details,paths] =  conrdk_subjects(id,options,reference,csdFlag)
%DPRST_SUBJECTS Function that sets all filenames and paths
%   IN:     EITHER (for quering both general and subject-specific paths:
%           id                  - the subject number as a string, e.g. '0001'
%           options (optional)  - the struct that holds all analysis options
%           OR (for only quering general paths & files):
%           options - the struct that holds all analysis options
%   OUT:    details     - a struct that holds all filenames and paths for
%                       the subject with subject number id
%           paths       - a struct that holds all paths and config files
%                       that are not subject-specific



% path to preprocessed data to match EEG with stim stream for continuous
% GLM


paths.referenceType.LMRM = fullfile(options.path.EEG.analysis,'preprocessedData','EEG','LMRM',sprintf('sub%03.0f',id));
paths.referenceType.averageReference = fullfile(options.path.EEG.subjects,sprintf('sub%03.0f',id),'eeg','average_reference');


if id == 55
    
    details.sessionIDs = 1:5;
    details.preproc.sessionIDs = 1:5;
    
elseif id == 40
    details.sessionIDs = [2,3,4,5];
    details.preproc.sessionIDs = 1:6;
elseif id == 47
    details.sessionIDs = [1,2,3,5];
    details.preproc.sessionIDs = 1:6;
elseif id == 42
    details.sessionIDs = [1,2,3,4,6]; % for some reason there is a problem with the triggers in session 5 where the end trigger of a block (210) is missing
    details.preproc.sessionIDs = 1:6;
elseif id == 21
    details.sessionIDs = [1,2,3,4,5,6]; 
    details.preproc.sessionIDs = [3,1,2,4,5,6];
elseif id == 24 
     details.sessionIDs = [1,2,3,4,5,6]; 
    details.preproc.sessionIDs = [2,1,3:6];
    elseif id == 70 
    details.sessionIDs = [1,2,3,5,6]; 
    details.preproc.sessionIDs = 1:6;
else
    details.sessionIDs = 1:6;
    details.preproc.sessionIDs = 1:6;
end




paths.preprocessed.SPM = fullfile(options.path.EEG.analysis,'preprocessedData','EEG',reference,sprintf('sub%03.0f',id));
if ~exist( paths.preprocessed.SPM,'dir' )
    mkdir( paths.preprocessed.SPM )
end
paths.rawData.SPM = fullfile(options.path.EEG.raw,'experiment',sprintf('sub%03.0f',id),'eeg');
paths.preprocessed.setFile.path = fullfile(options.path.EEG.analysis,'preprocessedData','EEG','setFiles');


if csdFlag
    paths.(reference).matchedEEG.saveName = fullfile(options.path.EEG.analysis,'convGLM','matchedEegData','LMRM',[sprintf('sub%03.0f',id),'_csdEEGdat.mat']);
    paths.(reference).subjectLevelGLM.jumps_absolute.saveName = fullfile(options.path.EEG.analysis,'convGLM','betasGLMData','LMRM',sprintf('sub%03.0f_betas_jumps_absolute_csd.mat',id));
    paths.(reference).subjectLevelGLM.response.saveName = fullfile(options.path.EEG.analysis,'convGLM','betasGLMData','LMRM',sprintf('sub%03.0f_betas_response_csd.mat',id));
    paths.(reference).subjectLevelGLM.coherence_responses.saveName = fullfile(options.path.EEG.analysis,'convGLM','betasGLMData','LMRM',sprintf('sub%03.0f_betas_coherence_responses_csd.mat',id));
    paths.(reference).subjectLevelGLM.vertical_jumps_absolute.saveName = fullfile(options.path.EEG.analysis,'convGLM','betasGLMData','LMRM',sprintf('sub%03.0f_betas_vertical_jumps_absolute_csd.mat',id));
    
    
    paths.(reference).singleTrial.appendedData.trialStart = fullfile(options.path.EEG.analysis,'conventionalEEGAnalysis',reference,['csd_trial_start_locked_EEG_dat_',sprintf('sub%03.0f',id),'.mat']);
    paths.(reference).singleTrial.appendedData.buttonPress = fullfile(options.path.EEG.analysis,'conventionalEEGAnalysis',reference,['csd_response_locked_EEG_dat_',sprintf('sub%03.0f',id),'.mat']);
    
    
else
    
    
    paths.(reference).matchedEEG.saveName = fullfile(options.path.EEG.analysis,'convGLM','matchedEegData','LMRM',[sprintf('sub%03.0f',id),'_EEGdat.mat']);
    paths.(reference).matchedEEG_tf.saveName = fullfile(options.path.EEG.analysis,'convGLM','matchedEegData','LMRM',[sprintf('sub%03.0f',id),'_EEGdat_tf.mat']);
    paths.(reference).subjectLevelGLM.jumps_absolute.saveName = fullfile(options.path.EEG.analysis,'convGLM','betasGLMData','LMRM',sprintf('sub%03.0f_betas_jumps_absolute.mat',id));
    paths.(reference).subjectLevelGLM.response.saveName = fullfile(options.path.EEG.analysis,'convGLM','betasGLMData','LMRM',sprintf('sub%03.0f_betas_response.mat',id));
    
    paths.(reference).subjectLevelGLM.coherence_responses.saveName = fullfile(options.path.EEG.analysis,'convGLM','betasGLMData','LMRM',sprintf('sub%03.0f_betas_coherence_responses.mat',id));
    paths.(reference).subjectLevelGLM.vertical_jumps_absolute.saveName = fullfile(options.path.EEG.analysis,'convGLM','betasGLMData','LMRM',sprintf('sub%03.0f_betas_vertical_jumps_absolute.mat',id));
    paths.(reference).subjectLevelGLM.all_regressors.saveName = fullfile(options.path.EEG.analysis,'convGLM','betasGLMData','LMRM',sprintf('sub%03.0f_betas_all_regressors.mat',id));
    paths.(reference).subjectLevelGLM.all_regressors_tf.saveName = fullfile(options.path.EEG.analysis,'convGLM','betasGLMData','LMRM',sprintf('sub%03.0f_betas_all_regressors_tfdecomp.mat',id));
    paths.(reference).subjectLevelGLM.all_regressors_with_reg.saveName = fullfile(options.path.EEG.analysis,'convGLM','betasGLMData','LMRM',sprintf('sub%03.0f_betas_all_regressors_regularised.mat',id));
     
    
    paths.(reference).singleTrial.appendedData.trialStart = fullfile(options.path.EEG.analysis,'conventionalEEGAnalysis',reference,['trial_start_locked_EEG_dat_',sprintf('sub%03.0f',id),'.mat']);
    paths.(reference).singleTrial.appendedData.buttonPress = fullfile('/Volumes/crdkData','conventionalEEGAnalysis',reference,['response_locked_EEG_dat',sprintf('sub%03.0f',id),'.mat']);
    
end



for session = 1:length(details.sessionIDs)
    if csdFlag
        paths.(reference).continuousPreproc(session).sessionList = fullfile(paths.referenceType.(reference),sprintf('csdcfMdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',id,details.sessionIDs(session)));
        paths.(reference).singleTrial.preproc(session).preprocSessionList = fullfile(paths.preprocessed.SPM,sprintf('csdnanart_csdcfMdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',id,details.sessionIDs(session)));
    
    else



        paths.(reference).continuousPreproc(session).sessionList = fullfile(paths.referenceType.(reference),sprintf('cfMdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',id,details.sessionIDs(session)));
        paths.(reference).singleTrial.preproc(session).preprocSessionList = fullfile(paths.preprocessed.SPM,sprintf('nanart_cfMdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',id,details.sessionIDs(session)));

    end

    paths.behaviour(session).sessionList = fullfile(options.path.EEG.subjects,sprintf('sub%03.0f',id),'behaviour',sprintf('sub%03.0f_sess%03.0f_behav.mat',id,details.sessionIDs(session)));
    paths.stimulus(session).sessionList = fullfile(options.path.EEG.subjects,sprintf('sub%03.0f',id),'stim',sprintf('sub%03.0f_sess%03.0f_stim.mat',id,details.sessionIDs(session)));



end

for session = 1:length(details.preproc.sessionIDs)
    details.setFiles(session).names = sprintf('sub%03.0f_sess%03.0f_eeg.set',id,details.preproc.sessionIDs(session));
    details.preproc.convertedFiles(session).names = sprintf('spmeeg_sub%03.0f_sess%03.0f_eeg.mat',id,details.preproc.sessionIDs(session));
    details.preproc.downsampledFiles(session).names = sprintf('dspmeeg_sub%03.0f_sess%03.0f_eeg.mat',id,details.preproc.sessionIDs(session));
    details.preproc.referencedFiles(session).names = sprintf('Mdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',id,details.preproc.sessionIDs(session));
    details.preproc.filteredFiles(session).names = sprintf('fMdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',id,details.preproc.sessionIDs(session));
    details.preproc.eyeCorrection(session).names = sprintf('cfMdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',id,details.preproc.sessionIDs(session));
    details.preproc.nanCorrection(session).names = sprintf('nanart_cfMdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',id,details.preproc.sessionIDs(session));
    details.preproc.tfDecomp(session).names = sprintf('tf_nanart_cfMdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',id,details.preproc.sessionIDs(session));

    
    %%%%%%
    details.preproc.ebdetectfig(session).names = fullfile(paths.preprocessed.SPM,sprintf('sub%03.0f_sess%03.0f_eyeblinkDetection.fig',id,details.preproc.sessionIDs(session)));
    details.preproc.ebspatialfig(session).names = fullfile(paths.preprocessed.SPM,sprintf('sub%03.0f_sess%03.0f_eyeblinkConfounds.fig',id,details.preproc.sessionIDs(session)));
    
    
    paths.(reference).preprocessing.rawData.EEG(session).sessionList = fullfile(options.path.EEG.raw,'experiment',sprintf('sub%03.0f',id),'eeg',sprintf('sub%03.0f_sess%03.0f_eeg.cdt',id,details.preproc.sessionIDs(session)));
    paths.(reference).preprocessing.EEGlab.EEG(session).sessionList = fullfile('/Volumes/crdkData/preprocessedData/EEG/setFiles',sprintf('sub%03.0f_sess%03.0f_eeg.set',id,details.preproc.sessionIDs(session)));
    
    
    
    
    %%%% spm files
    paths.preprocessing.SPMFiles.converted(session).sessionList = fullfile(paths.preprocessed.SPM,sprintf('spmeeg_sub%03.0f_sess%03.0f_eeg.mat',id,details.preproc.sessionIDs(session)));
    paths.preprocessing.SPMFiles.downsampled(session).sessionList = fullfile(paths.preprocessed.SPM,sprintf('dspmeeg_sub%03.0f_sess%03.0f_eeg.mat',id,details.preproc.sessionIDs(session)));
    paths.preprocessing.SPMFiles.rereferenced(session).sessionList = fullfile(paths.preprocessed.SPM,sprintf('Mdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',id,details.preproc.sessionIDs(session)));
    paths.preprocessing.SPMFiles.filtered(session).sessionList = fullfile(paths.preprocessed.SPM,sprintf('fMdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',id,details.preproc.sessionIDs(session)));
    paths.preprocessing.SPMFiles.eyeCorrection(session).sessionList = fullfile(paths.preprocessed.SPM,sprintf('cfMdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',id,details.preproc.sessionIDs(session)));
    paths.preprocessing.SPMFiles.csd(session).sessionList = fullfile(paths.preprocessed.SPM,sprintf('csdfMdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',id,details.preproc.sessionIDs(session)));
    
    
end







%%--preprocessing eyeblink thresholds-------------------------------------%

details.preproc.eyeblink.threshold = nan(2,6);
switch id
    
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

switch id
    
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
switch id
    case 19
        for session = 1:length(details.preproc.sessionIDs)
            if session == 2 
        details.preclean(session).badchannels = [30];
            else
             details.preclean(session).badchannels = [];    
            end
            
        end 
                
    otherwise
        details.preclean.badchannels = [];
end



end