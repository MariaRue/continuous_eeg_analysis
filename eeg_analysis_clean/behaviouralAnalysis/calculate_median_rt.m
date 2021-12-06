function [medianRtGroupLevel, medianRtSubject] = calculate_median_rt(options, all_responses)


for subject = 1:options.totalNumberofSubjects
    for condition = 1:4
        for coherence = 1:length(options.trialCoherence)
            
    % idx to all trials with a  a given coherence, condition , correct response,and a specific subject 
    % - only get rts < 3.5s (length of short trial condition)         
    idxRts = abs(all_responses(:,4)) == options.trialCoherence(coherence) & all_responses(:,9) == condition & all_responses(:,7) == 1 & all_responses(:,11) == subject & all_responses(:,2) <= 3.5; 
    
   % calculate median RT per subject 
    medianRtSubject(coherence,condition,subject)  = nanmedian(all_responses(idxRts,2));
    
        end 
    end 

end 

medianRtGroupLevel = nanmedian(medianRtSubject,3);



end 