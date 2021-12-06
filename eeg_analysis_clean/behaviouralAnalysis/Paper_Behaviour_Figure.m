% This script produces the plots of the behavioural data shown in figure
% ???
%% Trial has to be replaced with what ever name we use for this
%%
options = get_options_behaviour_rdk(); % sets options such as path to data, colour scheme etc
load(fullfile(options.path.behaviour,'behav_data_all_subjs_all3')); % load the behavioural data

%% calculate the false alarm rate per seconds baseline/interval???

[FalseAlarmRatePerSubject, FalseAlarmRateGroupLevel] = calculate_false_alarm_rate(options,mean_stim_streams,all_responses);


% plot false alarm rate
figure
hold on

for condition = 1:4
    
    bar(condition, FalseAlarmRateGroupLevel(condition), 'FaceColor',options.colours4Conditions(condition,:),'EdgeColor',options.colours4Conditions(condition,:))
    plot([ones(size(FalseAlarmRatePerSubject,1),1).*condition],squeeze(FalseAlarmRatePerSubject(:,condition)),'k.','MarkerSize',7,'LineWidth',3);
    
end

xticks([1 2 3 4])
xticklabels(options.conditionLabels)
ylabel('False Alarm Rate/s')

tidyfig;

% 2x2 rm Anova 
statsAnovaFalseAlarms = rm_Anova2_for_false_alarm_rate(FalseAlarmRatePerSubject);

%% detection rate

[detectRateSubject, detectRate, SEDetectRate] = calculate_detect_rate(options,all_responses);

% plot detection rate 
figure 
hold on

errorbar(options.trialCoherence,detectRate(:,1), SEDetectRate(:,1),'x-','Color',options.colours4Conditions(1,:),'LineWidth',3)
errorbar(options.trialCoherence,detectRate(:,2), SEDetectRate(:,2),'x-','Color',options.colours4Conditions(2,:),'LineWidth',3)
errorbar(options.trialCoherence,detectRate(:,3), SEDetectRate(:,3),'x-','Color',options.colours4Conditions(3,:),'LineWidth',3)
errorbar(options.trialCoherence,detectRate(:,4), SEDetectRate(:,4),'x-','Color',options.colours4Conditions(4,:),'LineWidth',3)
xlim([0.2 0.6])

xticks([0.3 0.4 0.5])
ylabel('hit rate')
xlabel('coherence')
legend(options.conditionLabels)
tidyfig;

% random effects analysis detect rate 
[statsGroupLevelSummaryDetectionRate] = regression_analysis_for_detection_rate(options,detectRateSubject);

%% plot reaction times 

% median rt for each coherence per subject and condition
[medianRtGroupLevel, medianRtSubject] = calculate_median_rt(options, all_responses);


figure 
hold on 

x = 1; % x axis value to separate median rts for different cohrences and conditions
for coherence = 1:length(options.trialCoherence)
    
    for condition = 1:4 
    
      b(condition) = bar(x,medianRtGroupLevel(coherence,condition),'FaceColor',options.colours4Conditions(condition,:),'EdgeColor',options.colours4Conditions(condition,:)); 
      plot([ones(options.totalNumberofSubjects,1).*x],squeeze(medianRtSubject(coherence,condition,:)),'k.','MarkerSize',7,'LineWidth',3);
        
        
    x = x+1;     
    end 
    
    x = x + 1.5; 
end 
hold off

 xticks([2.5, 8, 13.5])
 xticklabels([0.3 0.4 0.5])
 xlabel('coherence')
 ylabel('median rt')
 legend(b(1:4),options.conditionLabels)
 tidyfig;
 
 

% rt analysis regression analysis for each subject and test betas in t-test 

 [statsGroupLevelSummaryRTs] = regression_analysis_for_median_rts(options, all_responses);


%% plot integration kernels 





