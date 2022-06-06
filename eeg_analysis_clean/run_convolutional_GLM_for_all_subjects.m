glmFlag = 'all_regressors';
doReg = false; %apply regularisation?
lambda =    1.7600e+09; %4 * value used in Gonclaves et al 2014 (as our regularisation matrix is same as theirs, divided by 4)
vertical = 0;
options = continuous_RDK_set_options('LTHiMac');
% load all_responses,stim_streams, mean_stim_streams
load_name = fullfile(options.path.preproc.behaviour,'behav_data_all_subjs_all3.mat'); % load behavioural data

if vertical 
 load(load_name,'all_responses','mean_stim_streams','stim_streams','vertical_stim_streams');  
else 
  load(load_name,'all_responses','mean_stim_streams','stim_streams');  
end 

% subject list
subjectList = [16, 18:21, 24, 26, 27, 28, 29, 31,  33, 34, 35, 42, 43, 47, 50, 51, 52, 54, 55, 57, 58]; %32 taken out
subjectList = subjectList(1); %let's just try for the first subject...
%subjectList = [62:64,66,68,70];

csdFlag = 0; % 1 for csd transformed data
reference = 'LMRM';

% loop through subjects
for subject = 1:length(subjectList)
    subID = subjectList(subject);
  
    disp('subject: ')
    disp(subID)
    [details,paths] =  conrdk_subjects( subID,options,reference,csdFlag);



    subjectMatchedEEGdata = load(paths.(reference).matchedEEG.saveName);
    
    
    
    for sessionCount = 1:length(details.sessionIDs) % loop through sessions
                session = details.sessionIDs(sessionCount);
        % load D structure from preprocessed data to get chanlabels that need
        % to be saved at a later point - I think this can be deleted in a
        % future version in which we save the chan labels with the EEGdat file
        % in the previous step where we match EEG with stimulus stream
        

        preprocSubjectSessionData = spm_eeg_load(paths.(reference).continuousPreproc(sessionCount).sessionList);
        chanlabels = preprocSubjectSessionData.chanlabels; 
        
        
        VEOG_indx = selectchannels(preprocSubjectSessionData,'VEOG'); % we need this to remove bad eye blinks in run_subject_level_convolutional_glm - in theory that is the same for each subject at the moment but not sure what to make with this or how to do it differently
        
        
        % if condition 4 could not be matched we only have 3 conditions,
        % therefore we need to check the number of conditions in a session 
        nConditions = length(subjectMatchedEEGdata.EEGDat{session});
        
        % loop through conditions
        for condition = 1:nConditions % DO WE HAVE CONDITIONS ALREADY SORTED? YES!
            
            if ~isempty(subjectMatchedEEGdata.EEGDat{session}{condition})
                
                if ~vertical 
                    vertical_stim_streams = []; 
                end 
                
            % create all regressor types
            all_regressors = select_and_prepare_regressor_data_for_subject_level_glm(subID,session,condition,all_responses,stim_streams, mean_stim_streams,vertical,vertical_stim_streams);
            
            % create designmatrix
            laggedDesignMatrix = create_subject_level_design_matrix_for_convolutional_glm(all_regressors,options,glmFlag);
            
        
            %identify breakpoints (important if regularisation is switched on)
            breakPoints = [];
            for reg = 2:length(laggedDesignMatrix.regressor_indices)
                breakPoints = [breakPoints laggedDesignMatrix.regressor_indices(reg).dm_row_idx(1)];
            end
            
            % run glm for each channel of a condition
            betas_per_condition = run_subject_level_convolutional_glm(laggedDesignMatrix.regressors_matrix, subjectMatchedEEGdata.EEGDat{session}{condition}, subjectMatchedEEGdata.badSamples{session}{condition}, VEOG_indx, chanlabels, doReg, lambda, breakPoints);
            
            % loop through regressors to save betas accordingly 
            for regressor = 1:length(laggedDesignMatrix.regressor_indices)
                betas_subject{regressor}(:,:,sessionCount,condition) = betas_per_condition(:,laggedDesignMatrix.regressor_indices(regressor).dm_row_idx); %betas (indexed by regressors): channels * lags * sessions * blocks
            end
            
            end
            
        end
    end
    
    % save betas
    
    if doReg %regularisation was switched on...
        [pp,nn,ee] = fileparts(paths.(reference).subjectLevelGLM.(glmFlag).saveName);
        nn = [nn '_regularised'];
        save(fullfile(pp,[nn ee]),'betas_subject' ,'chanlabels','lambda');
    else
        save(paths.(reference).subjectLevelGLM.(glmFlag).saveName,'betas_subject' ,'chanlabels');
    end
end