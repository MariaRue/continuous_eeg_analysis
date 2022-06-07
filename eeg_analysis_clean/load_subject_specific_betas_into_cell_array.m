function  [betas_all_subjects, chanlabels] = load_subject_specific_betas_into_cell_array (subjectList,options, reference, glmFlag, csdFlag, nS)

for subject = 1:nS
    subID = subjectList(subject);
    
    
    disp('subject: ')
    disp(subID)
    [details,paths] =  conrdk_subjects( subID,options,reference,csdFlag);

    % load betas for one subject
    betaSubject = load( paths.(reference).subjectLevelGLM.(glmFlag).saveName);
    
    % get chanlabes = needed later for preparing data for permtest 
    
    chanlabels = betaSubject.chanlabels; 
    
    for regressor = 1:length(options.subjectLevelGLM.(glmFlag).regressors)
        
        for sessionCount = 1:length(details.sessionIDs) % loop through sessions
            session = details.sessionIDs(sessionCount);
            
            for condition = 1:4
               
                betas_all_subjects{regressor}(subject,:,:,session,condition) = betaSubject.betas_subject{regressor}(:,:,session,condition);
               
            end
        end
    end
end

end 