timelockCondition = 'buttonPress';
reference = 'LMRM';
csdFlag = 0;
subjectList = [16, 18:21, 24, 26, 27, 28, 29, 31, 33, 34, 35,  42, 43, 47, 50, 51, 52, 54, 55, 57, 58];
%subjectList = 58;
%subjectList = [62:64,66,68,70];
% for ft analysis a few subjects (19, 21, 24, 27, 33) had to be removed because the
% preprocessing pipeline did not work - some eyeblink confound issue? 
subjectList = [16, 18, 20, 26, 28, 29, 31, 34, 35,  42, 43, 47, 50, 51, 52, 54, 55, 57, 58];
options = continuous_RDK_set_options('iMac');
load_name = fullfile(options.path.preproc.behaviour,'behav_data_all_subjs_all3.mat'); % load behavioural data
load(load_name)
whichAnalysis = 'tf_analysis';

for subjects = 1:length(subjectList)
    subID = subjectList(subjects);
    
    [details,paths] =  conrdk_subjects(subID,options,reference,csdFlag);
    
    data = [];
    
    for session = 1:length(details.sessionIDs)
        
        
        sessionID = details.sessionIDs(session);
        switch whichAnalysis
            case 'tf_analysis'
                dataName = paths.(reference).singleTrial.preproc(session).tf_analysis.preprocSessionList;
            case 'timeseries_analysis'
                dataName =  paths.(reference).singleTrial.preproc(session).preprocSessionList;
        end
        
        responses = all_responses(all_responses(:,10)==sessionID & all_responses(:,12)==subID,:);
        stream_sj = unique(all_responses(all_responses(:,12)==subID,11));
        
        meanStimSession = mean_stim_streams_org{stream_sj,sessionID};
        
        data{session} = create_single_trials_session(responses,meanStimSession,options,timelockCondition,dataName);
        
        
    end
    
    
    
    cfg = [];
    
    dataAppend =  ft_appenddata(cfg,data{:});
    
    switch whichAnalysis
        case 'tf_analysis'
            
            % pre-process trials - demean and detrend for fourier analysis and spectral leakage problems
            [dataDetrendDemean] = detrend_demean_trials_for_tf_analysis(dataAppend);
            
            % redefine trial - shorter lenght 
            dataRedefined = redefine_trial_to_shorter_length_for_tf_analysis(dataDetrendDemean, options);

     
            tfdataAppend = compute_tf_on_single_trial_data(dataRedefined,options);
            save(paths.(reference).singleTrial.tf_analysis.appendedData.(timelockCondition),'tfdataAppend')
            
        case 'timeseries_analysis'
            save(paths.(reference).singleTrial.appendedData.(timelockCondition),'dataAppend')
    end
    
    
    
end