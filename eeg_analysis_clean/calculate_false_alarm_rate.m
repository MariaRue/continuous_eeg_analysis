function falseAlarmRate =  calculate_false_alarm_rate(responses,SubjectList,nS, mean_stim_streams)


for subject = 1:nS
    sessions = 6; 
    
    subjectID =  SubjectList(subject);
    if subjectID == 26 
        sessions = 5; 
    end 
    for sessionID = 1:sessions
        for condition = 1:4
            
            
            time_ITI_sec =   sum(mean_stim_streams{subjectID,sessionID}(:,condition) == 0)/100;
            
            
            FA_num = sum(responses(:,7) == 2 & responses(:,9) == condition & responses(:,10) == sessionID & responses(:,11) == subjectID);
            
            FA_rate(sessionID, subjectID, condition) = FA_num/time_ITI_sec;
            
            FA_rate(sessionID, subjectID, condition) = FA_rate(sessionID, subjectID, condition)*60; %convert to minutes
        end
    end
    
    for condition = 1:4
    % mean per subject
subjectLevel(condition,subject) = mean(squeeze(FA_rate(:,subjectID,condition))); 
    end
end

%%%% This was for excluding outliers - but don't think we do that anymore 
% mean_FA_rate_total = mean(subjectLevel(:)); 
% sd_FA_rate_total = 2 * std(subjectLevel(:));
% 
% [~,sj_id] = find(subjectLevel>= mean_FA_rate_total + sd_FA_rate_total | subjectLevel<= mean_FA_rate_total - sd_FA_rate_total );
% subjectLevel(:,sj_id) = [];

% mean across subjects
for condition = 1:4
    
groupLevel(condition) = mean(subjectLevel(condition,:)); 
%se_FA_rate_all(condition) = std(subjectLevel(condition,:))/sqrt(nS);

end 

falseAlarmRate.subjectLevel = subjectLevel; 
falseAlarmRate.groupLevel = groupLevel; 



% y = [[FA_rate_sj(1,:)';FA_rate_sj(2,:)'],[FA_rate_sj(3,:)';FA_rate_sj(4,:)']];
% 
% [p,tbl,stats] = anova2(y,length(FA_rate_sj(1,:)));





end 