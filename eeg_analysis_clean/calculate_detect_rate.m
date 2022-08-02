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

%% now perform 3-way ANOVA for factors of FREQUENCY, LENGTH, COHERENCE

slDataforANOVA = permute(subjectLevel,[3 2 1]);

varNames = {'Y1' 'Y2' 'Y3' 'Y4' 'Y5' 'Y6' 'Y7' 'Y8' 'Y9' 'Y10' 'Y11' 'Y12'};
t = array2table(slDataforANOVA(:,:),'VariableNames',varNames);

factorNames = {'Freq','Length','Coherence'};
within = table({'F';'F';'R';'R';'F';'F';'R';'R';'F';'F';'R';'R'},...
               {'S';'L';'S';'L';'S';'L';'S';'L';'S';'L';'S';'L'},...
               {'30';'30';'30';'30';'40';'40';'40';'40';'50';'50';'50';'50'},'VariableNames',factorNames); %F = frequent, R = Rare, S = short, L = long

% fit the repeated measures model
rm = fitrm(t,'Y1-Y12~1','WithinDesign',within);
[ranovatbl] = ranova(rm, 'WithinModel','Freq*Length*Coherence');

%% compress into output structure

detectRate.subjectLevel = subjectLevel; 
detectRate.groupLevel = groupLevel; 
detectRate.groupLevelSe = groupLevelSe; 
detectRate.ranovatbl = ranovatbl;



end 