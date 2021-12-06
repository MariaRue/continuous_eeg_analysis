clear all


options = get_options_behaviour_rdk();
load(fullfile(options.path.behaviour,'behav_data_all_subjs_all3'));


counter = 0;

%% loop through subjects and calculate mean interval length
mean_coherences = [];
for subject = 1 : options.totalNumberofSubjects
    details = conrdk_subjects_behaviour(subject);
    
    
    [MeanLengthOfInterval,StartOfIntertrialIntervals,EndOfIntertrialIntervals] = calculate_mean_length_of_interval({mean_stim_streams{subject,:}},details);
    
    
    for condition = 1:4
        
        kernelsLate = []; 
        kernelsEarly = []; 
        for session = details.sessionIDs
            
            
            
            stimulusStream = stim_streams{subject,session}(:,condition);
            responses = all_responses((all_responses(:,9)== condition & all_responses(:,10) == session & all_responses(:,11) == subject),:);
            
            triggersRight = responses((responses(:,7) == 2 & responses(:,3) == 1),6);
            triggersRight(triggersRight <= options.lags) = [];
            triggersLeft = responses((responses(:,7) == 2 & responses(:,3) == 0),6);
            triggersLeft(triggersLeft <= options.lags) = [];
            
            
            if any(triggersRight) 
            [falseAlarmsSortedRight{session}] = divide_fas_in_late_and_early(triggersRight,MeanLengthOfInterval(condition),StartOfIntertrialIntervals{session,condition},EndOfIntertrialIntervals{session,condition});
            else 
               falseAlarmsSortedRight{session} = nan(1,1);
               
            end 
            
            if any(triggersLeft)
            [falseAlarmsSortedLeft{session}] = divide_fas_in_late_and_early(triggersLeft,MeanLengthOfInterval(condition),StartOfIntertrialIntervals{session,condition},EndOfIntertrialIntervals{session,condition});
            else 
                falseAlarmsSortedLeft{session} = nan(1,1);
            end 
            
            
            if  (isempty(triggersLeft) && isempty(triggersRight))
            kernelsLate(:,session) = nan(options.lags+1,1); 
            kernelsEarly(:,session) = nan(options.lags+1,1); 
            else
            [kernelsLate(:,session), kernelsEarly(:,session)] = extract_kernels_from_stimStreams(stimulusStream,options,falseAlarmsSortedLeft{session},falseAlarmsSortedRight{session});
            end 
            
        end
        
        kernelsLateCondition(:,condition,subject) = nanmean(kernelsLate,2); 
        kernelsEarlyCondition(:,condition,subject) = nanmean(kernelsEarly,2); 
        
        %falseAlarmsSortedRightCondition{condition} = transform_cell_2_mat(falseAlarmsSortedRight);
        %falseAlarmsSortedLeftCondition{condition} = transform_cell_2_mat(falseAlarmsSortedLeft);
        
       
    
    end
    
  
   
    
    
    
    
    
end

kernelsLateOverSubjects = nanmean(kernelsLateCondition,3); 
kernelsEarlyOverSubjects = nanmean(kernelsEarlyCondition,3); 


% plot 
for condition = 1:4
    subplot(2,2,condition)
    hold on 
    plot(1:options.lags+1, kernelsLateOverSubjects(:,condition),'k-','LineWidth',3)
    plot(1:options.lags+1, kernelsEarlyOverSubjects(:,condition),'b-','LineWidth',3)
    title(options.conditionLabels{condition})
    hold off 
    
    if condition == 1
        
        legend('late', 'early')
        xlabel('ms 0 = button press') 
        
        
    end 
    
    tidyfig; 
    
end 

% plot 

%frequent vs rare 

% long vs short 
kernelsLateOverSubjectsLong = mean(kernelsLateOverSubjects(:,[2,4]),2);
kernelsLateOverSubjectsShort = mean(kernelsLateOverSubjects(:,[1,3]),2);
kernelsLateOverSubjectsFreq = mean(kernelsLateOverSubjects(:,[1,2]),2);
kernelsLateOverSubjectsRare = mean(kernelsLateOverSubjects(:,[3,4]),2);

kernelsEarlyOverSubjectsLong = mean(kernelsEarlyOverSubjects(:,[2,4]),2);
kernelsEarlyOverSubjectsShort = mean(kernelsEarlyOverSubjects(:,[1,3]),2);
kernelsEarlyOverSubjectsFreq = mean(kernelsEarlyOverSubjects(:,[1,2]),2);
kernelsEarlyOverSubjectsRare = mean(kernelsEarlyOverSubjects(:,[3,4]),2);

figure
subplot(2,2,1)
hold on 
plot(1:options.lags+1,kernelsLateOverSubjectsLong,'k','LineWidth',3)
plot(1:options.lags+1,kernelsLateOverSubjectsShort,'b','LineWidth',3)
legend('late long', 'late short')
 xlabel('ms 0 = button press') 
hold off 
tidyfig;

subplot(2,2,2)
hold on 
plot(1:options.lags+1,kernelsLateOverSubjectsFreq,'k','LineWidth',3)
plot(1:options.lags+1,kernelsLateOverSubjectsRare,'b','LineWidth',3)
legend('late freq', 'late rare')
hold off 
tidyfig;

subplot(2,2,3)
hold on 
plot(1:options.lags+1,kernelsEarlyOverSubjectsLong,'k','LineWidth',3)
plot(1:options.lags+1,kernelsEarlyOverSubjectsShort,'b','LineWidth',3)
legend('early long', 'early short')
hold off 
tidyfig;

subplot(2,2,4)
hold on 
plot(1:options.lags+1,kernelsEarlyOverSubjectsFreq,'k','LineWidth',3)
plot(1:options.lags+1,kernelsEarlyOverSubjectsRare,'b','LineWidth',3)
legend('early freq', 'early rare')
hold off 
tidyfig;