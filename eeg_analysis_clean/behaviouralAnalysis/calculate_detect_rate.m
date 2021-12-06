function [detectRateSubject, detectRate, SEDetectRate] = calculate_detect_rate(options,all_responses)


% WHAT NAME TO USE FOR TRIAL???

for subject = 1:options.totalNumberofSubjects
     details = conrdk_subjects_behaviour(subject);
     
    for condition = 1:4
      
        
  
        
        for session = 1:length(details.sessionIDs)

            % extract trials for subject, session and condition
            allTrialsIdx = all_responses(:,9) == condition & all_responses(:,10) == details.sessionIDs(session) & all_responses(:,11) == subject;
            allTrialResponses = all_responses(allTrialsIdx,[4,3,6,7]);
            
            
         
            for coherence = 1:length(options.trialCoherence) % count number of correct responses per coherence/condition and calculate rate of correct responses over all trials
         
                numberOfTrials = sum(abs(allTrialResponses(:,1)) == options.trialCoherence(coherence));
                
                hitTrials =  sum(abs(allTrialResponses(:,1)) == options.trialCoherence(coherence) & allTrialResponses(:,4) ==1);
                 
                detectRateSession(coherence,details.sessionIDs(session)) = hitTrials/numberOfTrials;
   
            end
            
            
            
            
        end
        
        % average detect Rate across sessions for each subject
        detectRateSubject(1,condition,subject) = nanmean(detectRateSession(1,:));
        detectRateSubject(2,condition,subject) = nanmean(detectRateSession(2,:));
        detectRateSubject(3,condition,subject) = nanmean(detectRateSession(3,:));
        
    end
    
    
    
end

% group level detect rate 
for condition = 1:4
    for coherence = 1:length(options.trialCoherence)
        
        
       detectRate(coherence,condition) = nanmean(squeeze(detectRateSubject(coherence,condition,:)));
        SEDetectRate(coherence,condition) = std(squeeze(detectRateSubject(coherence,condition,:)))/sqrt(options.totalNumberofSubjects);
        
    end
end
