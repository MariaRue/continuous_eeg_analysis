function [FalseAlarmRatePerSubject, FalseAlarmRateGroupLevel] = calculate_false_alarm_rate(options,mean_stim_streams,all_responses)

for subject = 1:options.totalNumberofSubjects
    
    details = conrdk_subjects_behaviour(subject);
    
    
    FalseAlarmRate = [];
    
    for session = 1:length(details.sessionIDs)
        for condition = 1:4
            
            % calculate seconds spent in baseline per session and condition
            timeSpentInBaseline =   sum(mean_stim_streams{subject,details.sessionIDs(session)}(:,condition) == 0)/options.SampleFrequency;
            
            % count false alarms
            NumberOfFalseAlarms = sum(all_responses(:,7) == 2 & all_responses(:,9) == condition & all_responses(:,10) == details.sessionIDs(session) & all_responses(:,11) == subject);
            
            % calculate false alarm rate per second
            FalseAlarmRate(details.sessionIDs(session), condition) =  NumberOfFalseAlarms/timeSpentInBaseline;
            
            
        end
        
        % average across sessions
        FalseAlarmRatePerSubject(subject,:) = squeeze(mean(FalseAlarmRate,1));
    end
    
    
    
    
end

% calculate overall mean and standard deviation to detect and remove
% outliers
TotalMeanFalseAlarmRate = mean(FalseAlarmRatePerSubject(:));
SDFalseAlarmRate = 2 * std(FalseAlarmRatePerSubject(:));

% remove outliers??? not done for detectrion rate - should I do this? SEs
% look reasonable for detect rate
[sj_id,~] = find(FalseAlarmRatePerSubject >= TotalMeanFalseAlarmRate + SDFalseAlarmRate | FalseAlarmRatePerSubject <= TotalMeanFalseAlarmRate - SDFalseAlarmRate );
FalseAlarmRatePerSubject(sj_id,:) = [];

% mean Group Level per condition
FalseAlarmRateGroupLevel = mean(FalseAlarmRatePerSubject,1);

end 