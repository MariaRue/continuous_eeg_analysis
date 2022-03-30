function detectRate = calculate_detect_rate(responses,SubjectList,nS,coherence, mean_stim_streams)

for subject = 1:nS
    
    subjectID = SubjectList(subject);
    
    for condition = 1:4
        sessions = 6;
        
        if subjectID == 26
            sessions = 5; 
        end 
        
        
        for sessionID = 1:sessions

            all_trials_idx = responses(:,9) == condition & responses(:,10) == sessionID & responses(:,11) == subjectID;
            all_trials = responses(all_trials_idx,[4,3,6,7]);
            
            
         
            for coh = 1:3
         
                num_of_trials = sum(abs(all_trials(:,1)) == coherence(coh));
                
                num_of_trials_1 = sum(abs(mean_stim_streams{subjectID,sessionID}(2:end,condition)) == coherence(coh) & abs(mean_stim_streams{subjectID,sessionID}(1:end-1,condition)) ==0 );

%                
                

                hit_trials =  sum(abs(all_trials(:,1)) == coherence(coh) & all_trials(:,4) ==1);
                 
                
                detect_rate_se(coh,sessionID) = hit_trials/num_of_trials;
                
                ntr(subjectID,condition,sessionID,coh) = num_of_trials; %from the response matrix
                ntr1(subjectID,condition,sessionID,coh) = num_of_trials_1; %from the stimulus stream
                
                if num_of_trials~=num_of_trials_1
                    keyboard;
                end
                
            end
            
            
            
            
        end
        
        subjectLevel(1,condition,subject) = nanmean(detect_rate_se(1,:));
        subjectLevel(2,condition,subject) = nanmean(detect_rate_se(2,:));
        subjectLevel(3,condition,subject) = nanmean(detect_rate_se(3,:));
        
    end
    
    
    
end

for condition = 1:4
    for coh = 1:3
        
        
        groupLevel(coh,condition) = nanmean(squeeze(subjectLevel(coh,condition,:)));
        groupLevelSe(coh,condition) = std(squeeze(subjectLevel(coh,condition,:)))/sqrt(nS);
        
    end
end

detectRate.subjectLevel = subjectLevel; 
detectRate.groupLevel = groupLevel; 
detectRate.groupLevelSe = groupLevelSe; 




end 