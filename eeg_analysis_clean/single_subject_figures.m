clear all;
addpath(genpath(pwd))

glmFlag = 'jumps_absolute';

options = continuous_RDK_set_options('iMac');

% subject list
subjectList = [16 18:21, 24, 26, 27, 28, 29, 31,  33, 34, 35, 42, 43, 47, 50, 51, 52, 54, 55, 57, 58]; %32 taken out
%subjectList = [62:64,66,68,70]; % vertical motion only 

csdFlag = 0; % 1 for csd transformed data
reference = 'LMRM';

for subject = 1:length(subjectList)
    subID = subjectList(subject);
    
    
    disp('subject: ')
    disp(subID)
    [details,paths] =  conrdk_subjects( subID,options,reference,csdFlag);
    
    % load betas for one subject
    betaSubject = load( paths.(reference).subjectLevelGLM.(glmFlag).saveName);
    
    for regressor = 1:length(options.subjectLevelGLM.(glmFlag).regressors)
        
        for sessionCount = 1:length(details.sessionIDs) % loop through sessions
            session = details.sessionIDs(sessionCount);
            
            for condition = 1:4
                
                betas_all_subjects{regressor}(subject,:,:,session,condition) = betaSubject.betas_subject{regressor}(:,:,session,condition);
                
                
            end
        end
    end
end
%%
for r = 1:length(options.subjectLevelGLM.(glmFlag).regressors) % loop over regressors
    betas_all_subjects_sessavg{r} = squeeze(mean(betas_all_subjects{r},4));
end
%% 


channel_ID = find(strcmp( betaSubject.chanlabels,'CPZ'));


RegressorsToPlot = betas_all_subjects_sessavg{3}; % coherence jump = 1, PE = 3 



for sj = 1:24
    
topo =  conv2(ones(10,1)/10,1,squeeze(RegressorsToPlot(sj,channel_ID,:,1)),'same');
subplot(6,4,sj)
plot(options.subjectLevelGLM.(glmFlag).regressors(3).timeBins, topo)
tidyfig;

end 