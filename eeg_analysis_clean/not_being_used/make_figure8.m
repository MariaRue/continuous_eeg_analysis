function make_figure8(plotVariables, options)

subjectListEEG = [16 18:21, 24, 26, 27, 28, 29, 31,  33, 34, 35, 42, 43, 47, 50, 51, 52, 54, 55, 57, 58]; %32 taken out


EEGpreproc = options.path.preproc.behaviour;  % path to behav data all subjs
load_name = fullfile(EEGpreproc,'behav_data_all_subjs_all3');
load(load_name)
completeSubjectListBehaviour = unique(all_responses(:,12));
%get behavioural subjects who are in subjectList from EEG
[~,~,SubjectListBehaviourEEG] = intersect(subjectListEEG,completeSubjectListBehaviour');
nS = length(SubjectListBehaviourEEG); 
coherence = [0.3 0.4 0.5]; 


% false alarm rate 
falseAlarmRate = calculate_false_alarm_rate_collapsed(all_responses,SubjectListBehaviourEEG,nS, mean_stim_streams); 


% integration kernels 
% not yet implemented 




%% figure
figure (8)
subplot(2,1,1) 
hold on 


 bar(1, falseAlarmRate.groupLevelLong, 'FaceColor',plotVariables.figure6.colours(2,:),'EdgeColor',plotVariables.figure6.colours(2,:))
 ln = plot([ones(length(falseAlarmRate.subjectLevelLong),1).*1],squeeze(falseAlarmRate.subjectLevelLong),'.','MarkerSize',plotVariables.figure6.falseAlarmRate.MarkerSize,'LineWidth',plotVariables.figure6.falseAlarmRate.LineWidth, 'Color', plotVariables.figure6.falseAlarmRate.Color);
 

 bar(2, falseAlarmRate.groupLevelShort, 'FaceColor',plotVariables.figure6.colours(1,:),'EdgeColor',plotVariables.figure6.colours(1,:))
 ln = plot([ones(length(falseAlarmRate.subjectLevelShort),1).*2],squeeze(falseAlarmRate.subjectLevelShort),'.','MarkerSize',plotVariables.figure6.falseAlarmRate.MarkerSize,'LineWidth',plotVariables.figure6.falseAlarmRate.LineWidth, 'Color', plotVariables.figure6.falseAlarmRate.Color);
 
 
 
ylabel('false alarm rate per [s]')
xticks([1 2 3 4])
xticklabels({'long'; 'short'})

hold off 


tidyfig; 

end 