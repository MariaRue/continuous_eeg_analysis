addpath(genpath('/Users/maria/Documents/Matlab/cbrewer'))
addpath('/Users/maria/Documents/MATLAB/raacampbell-shadedErrorBar-9b19a7b')
EEGpreproc = '/Volumes/crdkData/preprocessedData/behaviour';  % path to behav data all subjs


condition = {'Tr frequent TR short', 'Tr frequent Tr short','Tr rare Tr short', 'Tr rare Tr long'};

% plot each subject separately for conditions
cl = cbrewer('div','RdBu', 12);
cl = cl([2 4 12 10],:);
%% load data
load_name = fullfile(EEGpreproc,'behav_data_all_subjs_all3');
load(load_name)

lags = 500; % frames back in time that are leading up to FA

which_responses = 'false alarms';  % calculating integration kernel for either FA button presses or button presses during trials options: 'false alarms' or 'trials', '3sec_rts',

with_coherence = 'without coherence level'; % if we want to know the coherence levels for the trials version
totalNumberOfSubjects = max(all_responses(:,11)); % number of subjects


coherence = [0.3 0.4 0.5];

counter = 0;

%% loop through subjects and find button presses
mean_coherences = [];
for subject = 1 :totalNumberOfSubjects
    
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
        MeanLengthOfInterval(subject, condition) = mean(LengthOfIntertrialInterval{condition}).*0.5;
    end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    for condition = 1 :4
        combinedLate = [];
        combinedEarly = [];
        
        totalNumberOfSessions = 6;
    if subject == 26
        totalNumberOfSessions = 5;
    end
    
        for session = 1:totalNumberOfSessions
            % select only stimstreams from all sessions that belong to specific
            % block
            
            % select all stim streams that belong to one subject
            stim_streams_sj = [];
            stim_streams_sj = stim_streams{subject,session}(:,condition);
            
            
            % select trigger streams that belong to one subject
            %             trigger_streams_sj = [];
            %             trigger_streams_sj = trigger_streams{sj,se}(:,bl);
            
            mean_streams = [];
            mean_streams = mean_stim_streams{subject,session}(:,condition);
            
            
            responses = all_responses((all_responses(:,9)== condition & all_responses(:,10) == session & all_responses(:,11) == subject),:);
            
            
         
            
            
            % find false alarms (called triggers)
            
            
            triggers_right = [];
            triggers_left = [];
            
            
     
            
            % this is with triggers from the response matrix
            %
            triggers_right = responses((responses(:,7) == 2 & responses(:,3) == 1),6);
            triggers_left = responses((responses(:,7) == 2 & responses(:,3) == 0),6);
            
            
             for falseAlarm = 1:length(triggers_right)
                
                % find interval
                
                IntervalStart = StartOfIntertrialIntervals{session,condition}(triggers_right(falseAlarm) < EndOfIntertrialIntervals{session,condition} & triggers_right(falseAlarm) > StartOfIntertrialIntervals{session,condition});
                
                % early or late?
                
                if triggers_right(falseAlarm) > (IntervalStart + MeanLengthOfInterval(condition))
                    
                    
                    triggers_right(falseAlarm,2) = 1;
                    
                    
                else
                    
                    triggers_right(falseAlarm,2) = 0;
                    
                end
  
                
             end
            
             
             
             
             
             for falseAlarm = 1:length(triggers_left)
                
                % find interval
                
                IntervalStart = StartOfIntertrialIntervals{session,condition}(triggers_left(falseAlarm) < EndOfIntertrialIntervals{session,condition} & triggers_left(falseAlarm) > StartOfIntertrialIntervals{session,condition});
                
                % early or late?
                
                if triggers_left(falseAlarm) > (IntervalStart + MeanLengthOfInterval(condition))
                    
                    triggers_left(falseAlarm,2) = 1;
                    
                    
                else
                    
                     triggers_left(falseAlarm,2) = 0;
                    
                end
  
                
             end
             
            
             
            % only choose triggers that are bigger than the lags we go
            % back
            
            
            triggers_right(triggers_right(:,1)<=lags,:) = [];
            triggers_left(triggers_left(:,1)<=lags,:) = [];
             
             % loop through triggers for right and left button presses and
                    % select coherences from stim_streams_bl
                    matrix_right_early = [];
                    matrix_right_late = [];
                    countLate = 1;
                    countEarly = 1;
                    for i = 1:length(triggers_right(:,1))
                        
                        if triggers_right(i,2) == 1
                        
                        
                        matrix_right_Late(:,i) = stim_streams_sj(triggers_right(i,1) - lags : triggers_right(i,1));
                        countLate = countLate +1;
                        
                        else 
                            
                         matrix_right_Early(:,i) = stim_streams_sj(triggers_right(i,1) - lags : triggers_right(i,1));
                            countEarly  =  countEarly +1;
                        end 
                    end
                    
                    matrix_left_early = [];
                    matrix_left_late = [];
                    countLate = 1;
                    countEarly = 1;
                  
                    for i = 1:length(triggers_left(:,1))
                        if triggers_left(i,2) == 1
                        matrix_left_Late(:,countLate) = stim_streams_sj(triggers_left(i,1) - lags : triggers_left(i,1)) .* -1 ; % to be able to combine with right button presses * -1
                        
                        countLate = countLate+1;
                        else 
                        matrix_left_Early(:,countEarly) = stim_streams_sj(triggers_left(i,1) - lags : triggers_left(i,1)) .* -1 ; % to be able to combine with right button presses * -1
                        countEarly = countEarly+1;
                        end 
                        
                    end
                    
 
                    
                    
                    combinedLate = [combinedLate,matrix_right_Late, matrix_left_Late];
                    combinedEarly = [combinedEarly,matrix_right_Early, matrix_left_Early];
             

            
        end
      
        
                mean_coherences_Late(:,condition,subject) = nanmean(combinedLate,2);
                sem_coherence_Late(:,condition,subject) = nanstd(combinedLate')/sqrt(size(mean_coherences_Late,1));
                %num_false_alarms(subject,condition) = size(combined,2);
                
                
                mean_coherences_Early(:,condition,subject) = nanmean(combinedEarly,2);
                sem_coherence_Early(:,condition,subject) = nanstd(combinedEarly')/sqrt(size(mean_coherences_Early,1));
                
                
                
        
    end
end

        sem_across_subjects_Late = squeeze(std(permute(mean_coherences_Late,[3,1,2]))/sqrt(size(mean_coherences_Late,3)));
        mean_across_subjects_Late = nanmean(mean_coherences_Late,3);
        
        sem_across_subjects_Early = squeeze(std(permute(mean_coherences_Early,[3,1,2]))/sqrt(size(mean_coherences_Early,3)));
        mean_across_subjects_Early = nanmean(mean_coherences_Early,3);
        
        
        
        
        % plot
        
        
        figure
        subplot(2,2,1)
        hold on 
        plot(1:lags+1, mean_across_subjects_Early(:,1),'k','LineWidth',3)
        plot(1:lags+1, mean_across_subjects_Late(:,1),'b','LineWidth',3)
        hold off 
        tidyfig;
        legend('early fas', 'late fas')
        title ('frequent and short') 
        
                subplot(2,2,2)
        hold on 
        plot(1:lags+1, mean_across_subjects_Early(:,2),'k','LineWidth',3)
        plot(1:lags+1, mean_across_subjects_Late(:,2),'b','LineWidth',3)
        hold off 
        title ('frequent and long') 
          tidyfig;
                        subplot(2,2,3)
        hold on 
        plot(1:lags+1, mean_across_subjects_Early(:,3),'k','LineWidth',3)
        plot(1:lags+1, mean_across_subjects_Late(:,3),'b','LineWidth',3)
        hold off 
        title ('rare and short') 
          tidyfig;
        
                                subplot(2,2,4)
        hold on 
        plot(1:lags+1, mean_across_subjects_Early(:,4),'k','LineWidth',3)
        plot(1:lags+1, mean_across_subjects_Late(:,4),'b','LineWidth',3)
        hold off 
        title ('rare and long') 
          tidyfig;