addpath(genpath('/Users/maria/Documents/Matlab/cbrewer'))
addpath('/Users/maria/Documents/MATLAB/raacampbell-shadedErrorBar-9b19a7b')
EEGpreproc = '/Volumes/crdkData/preprocessedData/behaviour'; % path to behav data all subjs

condition = {'Tr frequent TR short', 'Tr frequent Tr long','Tr rare Tr short','Tr rare Tr long'};
% plot each subject separately for conditions
cl = cbrewer('div','RdBu', 12);
cl = cl([2 4 12 10],:);

%% load data
load_name = fullfile(EEGpreproc,'behav_data_all_subjs_all3');
load(load_name)

lags = 500; % frames back in time that are leading up to FA

totalNumberOfSubjects = max(all_responses(:,11)); % number of subjects

coherence = [0.3 0.4 0.5];

%%
for subject = 1:totalNumberOfSubjects
    totalNumberOfSessions = 6;
    if subject == 26
        totalNumberOfSessions = 5;
    end
    
    LengthOfIntertrialInterval = cell(4,1);
    for session = 1:totalNumberOfSessions
        
        for condition = 1:4 % loop conditions
            
            
            StartOfIntertrialIntervals{session,condition} = find( (mean_stim_streams{subject,session}(1:end-1,condition) ~= 0 & mean_stim_streams{subject,session}(2:end,condition) == 0));
            EndOfIntertrialIntervals{session,condition} = find( (mean_stim_streams{subject,session}(1:end-1,condition) == 0 & mean_stim_streams{subject,session}(2:end,condition) ~= 0));
            
            StartOfIntertrialIntervals{session,condition} = [1; StartOfIntertrialIntervals{session,condition}];
            EndOfIntertrialIntervals{session,condition} = [EndOfIntertrialIntervals{session,condition}; length(mean_stim_streams{subject,session}(:,condition))];
            
            LengthOfIntertrialInterval{condition} = [LengthOfIntertrialInterval{condition}; (EndOfIntertrialIntervals{session,condition} - StartOfIntertrialIntervals{session,condition})];
            
        end
        
    end
    
    for condition = 1:4
        MeanLengthOfInterval(condition) = mean(LengthOfIntertrialInterval{condition}).*0.5;
    end
    
   
    
    for session = 1:totalNumberOfSessions
        for condition = 1:4
            FalseAlarmList = all_responses((all_responses(:,7) == 2 & all_responses(:,9) == condition & all_responses(:,10) == session & all_responses(:,11) == subject),6);
            
            
            for falseAlarm = 1:length(FalseAlarmList)
                
                % find interval
                
                IntervalStart = StartOfIntertrialIntervals{session,condition}(FalseAlarmList(falseAlarm) < EndOfIntertrialIntervals{session,condition} & FalseAlarmList(falseAlarm) > StartOfIntertrialIntervals{session,condition});
                
                
                % early or late?
                
                if FalseAlarmList(falseAlarm) > (IntervalStart + MeanLengthOfInterval(condition))
                    
                    
                    FalseAlarmType{subject,session,condition}(falseAlarm) = 1;
                    
                    
                else
                    
                    FalseAlarmType{subject,session,condition}(falseAlarm) = 0;
                    
                end
                
        
              
                
                
            end
            
                    %calculate time in early and late interval periods 
                
                EarlyLength(subject,session,condition) = length(StartOfIntertrialIntervals{session,condition}) .* MeanLengthOfInterval(condition);
                
                LengthOfIntervals{session,condition} = EndOfIntertrialIntervals{session,condition} - StartOfIntertrialIntervals{session,condition};
                
                LateLength(subject,session,condition) = sum(LengthOfIntervals{session,condition} - MeanLengthOfInterval(condition));
            
            
            
            sumFalseAlarmsLate(subject,session,condition) = sum(FalseAlarmType{subject,session,condition})/(LateLength(subject,session,condition)/100);
            sumFalseAlarmsEarly(subject,session,condition) = sum(FalseAlarmType{subject,session,condition}==0)/(EarlyLength(subject,session,condition)/100);
            
        end
    end
    
    
end



% mean per subject

subjectMeanFalseAlarmsLate = squeeze(mean(sumFalseAlarmsLate,2));

subjectMeanFalseAlarmsEarly = squeeze(mean(sumFalseAlarmsEarly,2));


GroupMeanFalseAlarmsLate = squeeze(mean(subjectMeanFalseAlarmsLate,1));
GroupMeanFalseAlarmsEarly = squeeze(mean(subjectMeanFalseAlarmsEarly,1));


ToPlot = [GroupMeanFalseAlarmsEarly(1), GroupMeanFalseAlarmsLate(1); GroupMeanFalseAlarmsEarly(2), GroupMeanFalseAlarmsLate(2); GroupMeanFalseAlarmsEarly(3), GroupMeanFalseAlarmsLate(3);GroupMeanFalseAlarmsEarly(4), GroupMeanFalseAlarmsLate(4)];

X = categorical({'freq short','freq long','rare short','rare long'});
X = reordercats(X,{'freq short','freq long','rare short','rare long'});
bar(X,ToPlot)

ylabel('false alarms/s')
legend('early','late')


% length in early and late trials per subject 



