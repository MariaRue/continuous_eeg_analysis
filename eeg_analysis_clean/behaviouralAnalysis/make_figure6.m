function make_figure6(plotVariables, options)

subjectListEEG = [16 18:21, 24, 26, 27, 28, 29, 31,  33, 34, 35, 42, 43, 47, 50, 51, 52, 54, 55, 57, 58]; %32 taken out


EEGpreproc = options.path.preproc.behaviour;  % path to behav data all subjs
load_name = fullfile(EEGpreproc,'behav_data_all_subjs_all3');
load(load_name)
completeSubjectListBehaviour = unique(all_responses(:,12));
%get behavioural subjects who are in subjectList from EEG
[~,~,SubjectListBehaviourEEG] = intersect(subjectListEEG,completeSubjectListBehaviour');
nS = length(SubjectListBehaviourEEG); 
coherence = [0.3 0.4 0.5]; 

% detection rate
detectRate = calculate_detect_rate(all_responses,SubjectListBehaviourEEG,nS,coherence, mean_stim_streams);

% false alarm rate 
falseAlarmRate = calculate_false_alarm_rate(all_responses,SubjectListBehaviourEEG,nS, mean_stim_streams); 

% reaction time 
ReactionTimes = calculate_reaction_times(all_responses,SubjectListBehaviourEEG,nS,coherence);

% integrationn kernels 
lags = 500;
[GroupIntegrationKernels, SubjectIntegrationKernels, SignificantTimePoints, ExpParameters] = calculate_integration_kernels(all_responses,SubjectListBehaviourEEG,nS, mean_stim_streams, stim_streams, trigger_streams,lags);



figure (6)
%% detect rate 
subplot(2,2,1)
hold on

errorbar(coherence,detectRate.groupLevel(:,1), detectRate.groupLevelSe(:,1),'x-','Color',plotVariables.figure6.colours(1,:),'LineWidth',plotVariables.figure6.detectRate.LineWidth)
errorbar(coherence,detectRate.groupLevel(:,2), detectRate.groupLevelSe(:,2),'x-','Color',plotVariables.figure6.colours(2,:),'LineWidth',plotVariables.figure6.detectRate.LineWidth)
errorbar(coherence,detectRate.groupLevel(:,3), detectRate.groupLevelSe(:,3),'x-','Color',plotVariables.figure6.colours(3,:),'LineWidth',plotVariables.figure6.detectRate.LineWidth)
errorbar(coherence,detectRate.groupLevel(:,4), detectRate.groupLevelSe(:,4),'x-','Color',plotVariables.figure6.colours(4,:),'LineWidth',plotVariables.figure6.detectRate.LineWidth)
xlim([0.2 0.6])
xlabel('trial??? coherence')  
ylabel('detect rate') 
legend(plotVariables.originalConditions)
tidyfig;
hold off

%% false Alarm rate 
subplot(2,2,2) 
hold on 
for condition = 1:4

 bar(condition, falseAlarmRate.groupLevel(condition), 'FaceColor',plotVariables.figure6.colours(condition,:),'EdgeColor',plotVariables.figure6.colours(condition,:))
 ln = plot([ones(nS,1).*condition],squeeze(falseAlarmRate.subjectLevel(condition,:)),'.','MarkerSize',plotVariables.figure6.falseAlarmRate.MarkerSize,'LineWidth',plotVariables.figure6.falseAlarmRate.LineWidth, 'Color', plotVariables.figure6.falseAlarmRate.Color);
 
end 
ylabel('false alarm rate per [s]')
xticks([1 2 3 4])
xticklabels(plotVariables.originalConditions)

hold off 


tidyfig; 

%% reaction times 
subplot(2,2,3) 
hold on 
x = 1; % bar spacer for plot
for coherence = 1:3
    
    for condition = 1:4 
    
       b(condition) = bar(x,ReactionTimes.groupLevel(coherence,condition),'FaceColor',plotVariables.figure6.colours(condition,:),'EdgeColor',plotVariables.figure6.colours(condition,:)); 
       ln = plot([ones(nS,1).*x],squeeze(ReactionTimes.subjectLevel(coherence,condition,:)),'.','MarkerSize',plotVariables.figure6.reactionTimes.MarkerSize,'LineWidth',plotVariables.figure6.reactionTimes.LineWidth,'Color',plotVariables.figure6.reactionTimes.Color);
        
        
    x = x+1;     
    end 
    
    x = x + 1.5; 
end 
hold off

 xticks([2.5, 8, 13.5])
 xticklabels([0.3 0.4 0.5])
 xlabel('coherence')
 ylabel('median rt < 3.5s trials only???')
 legend(b(1:4),plotVariables.originalConditions)
tidyfig
%% integration kernels 

subplot(2,2,4)
for condition = 1:4

   hold on
   b(condition) = plot(GroupIntegrationKernels.mean(:,condition), 'LineWidth', plotVariables.figure6.IntegrationKernels.LineWidth,'Color',plotVariables.figure6.colours(condition,:));
    hold on
    h = shadedErrorBar(1:lags+1,GroupIntegrationKernels.mean(:,condition),(GroupIntegrationKernels.sem(:,condition)), 'lineprops', 'k-')
    h.patch.FaceColor = plotVariables.figure6.colours(condition,:);
    h.mainLine.Color = plotVariables.figure6.colours(condition,:);
    
    
    title('mean coherence leading to a button press')
    xlim([0 501])
    
    %tidyfig
    xticks([0:100:500])
    xticklabels([-5 -4 -3 -2 -1 0])
    if condition == 1
        xlabel('time to button press [s]')
        ylabel('mean coherence')
    end
end
plot(1:lags+1,zeros(lags+1,1),'k--','LineWidth',plotVariables.figure6.IntegrationKernels.LineWidth)
hold on 

for id = 1 
    
plot(SignificantTimePoints{id},ones(length(SignificantTimePoints{id}),1).* 0.6,'.','MarkerSize',plotVariables.figure6.IntegrationKernels.MarkerSize,'Color',plotVariables.figure6.IntegrationKernels.Color)

   
end
hold off
ylim([-0.1 0.65])

legend(b(1:4),plotVariables.originalConditions)

tidyfig; 
end 