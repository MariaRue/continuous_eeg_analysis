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
    %{
    out = validate_betas(betaSubject);
    if out < 0
        keyboard;
    else
        fprintf('validated with %0.0f missing blocks\n',out);
    end
    %}
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

function out = validate_betas(betaSubject)
out = 0;
for i = 1:length(betaSubject.betas_subject)
    if all(size(betaSubject.betas_subject{i})~=[ 64   251     6     4]) %first, check dimensions of betas
        out = -i;
    elseif any(isnan(betaSubject.betas_subject{i}(:))) %then check we don't have any nans
        out = -999;
    else % otherwise, count the number of missing blocks
        tmp = permute(betaSubject.betas_subject{i},[4 3 2 1]);
        if i == 1 
            out = sum(sum(all(tmp(:,:,:)==0,3)));
        end
    end
end


end

        
        
    