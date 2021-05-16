timelockCondition = 'trialStart';
reference = 'LMRM'; 
csdFlag = 1;
subjectList = [16, 18:21, 24, 26, 27, 28, 29, 31, 33, 34, 35,  42, 43, 47, 50, 51, 52, 54, 55, 57, 58];
%subjectList = 58;
subjectList = [62:64,66,68,70];
options = continuous_RDK_set_options('iMac');
load_name = fullfile(options.path.preproc.behaviour,'behav_data_all_subjs_all3.mat'); % load behavioural data
load(load_name)

for subjects = 1:length(subjectList)
    subID = subjectList(subjects); 
    
    [details,paths] =  conrdk_subjects(subID,options,reference,csdFlag); 
    
    data = []; 
    
    for session = 1:length(details.sessionIDs)
        
        
            sessionID = details.sessionIDs(session); 
            dataName =  paths.(reference).singleTrial.preproc(session).preprocSessionList;
            
            responses = all_responses(all_responses(:,10)==sessionID & all_responses(:,12)==subID,:);
            stream_sj = unique(all_responses(all_responses(:,12)==subID,11));

            meanStimSession = mean_stim_streams_org{stream_sj,sessionID};

    
            data{session} = create_single_trials_session(responses,meanStimSession,options,timelockCondition,dataName);


    end
    
    
    cfg = [];
    
    dataAppend =  ft_appenddata(cfg,data{:});
      
    save(paths.(reference).singleTrial.appendedData.(timelockCondition),'dataAppend')
    
    
end