options = continuous_RDK_set_options('iMac');
% load all_responses,stim_streams, mean_stim_streams
load_name = fullfile(options.path.preproc.behaviour,'behav_data_all_subjs_all3.mat'); % load behavioural data

load(load_name,'all_responses','mean_stim_streams','stim_streams');  

subjectList = [16, 18:21, 24, 26, 27, 28, 29, 31,  33, 34, 35, 42, 43, 47, 50, 51, 52, 54, 55, 57, 58];

reference = 'LMRM';
csdFlag = 0;
vertical = 0;
for subject = 1:length(subjectList)
    subID = subjectList(subject);
  
    
    disp('subject: ')
    disp(subID)
    [details,paths] =  conrdk_subjects( subID,options,reference,csdFlag);
    
    corrMatrixSession =[];
      for sessionCount = 1:length(details.sessionIDs) % loop through sessions
                session = details.sessionIDs(sessionCount);
           
                    % if condition 4 could not be matched we only have 3 conditions,
        % therefore we need to check the number of conditions in a session 
        nConditions = 4;
        
        % loop through conditions
        for condition = 1:nConditions % DO WE HAVE CONDITIONS ALREADY SORTED? YES!
           
            vertical_stim_streams = []; 
            
            all_regressors = select_and_prepare_regressor_data_for_subject_level_glm(subID,session,condition,all_responses,stim_streams, mean_stim_streams,vertical,vertical_stim_streams);
            
            
            regMatrix = [all_regressors.coherence_jump, all_regressors.coherence_jump_level,all_regressors.prediction_error,all_regressors.absoluted_stimulus];
            corrMatrixSession(:,:,condition,sessionCount) = corr(regMatrix);
           
        end
      end
     
      corrMatrixSubject(:,:,:,subject) = mean(corrMatrixSession,4);
end 

corrMatrixCondition = mean(corrMatrixSubject,4);
corrMatrixAll = round(mean(corrMatrixCondition,3).^2,2);

%seAll = std(mean(corrMatrixSubject,3)); 