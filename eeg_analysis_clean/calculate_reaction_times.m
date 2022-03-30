function ReactionTimes = calculate_reaction_times(responses,SubjectList,nS,coherence) 

for subject = 1:nS
    subjectID = SubjectList(subject);
    for con = 1:4
        for coh = 1:3
            
    idx_rts = abs(responses(:,4)) == coherence(coh) & responses(:,9) == con & responses(:,7) == 1 & responses(:,11) == subjectID & responses(:,2) <= 3.5; % idx to all trials with a  a given coherence, condition , correct response,and a specific subject and only reactiones <= 3.5s 
    
   
    ReactionTimes.subjectLevel(coh,con,subject)  = nanmedian(responses(idx_rts,2));
        end 
    end 

end 

ReactionTimes.groupLevel = nanmedian(ReactionTimes.subjectLevel,3);



end 