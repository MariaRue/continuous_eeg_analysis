function make_figure10b(plotVariables, options)



%behavioural figure for trial periods
subjectListEEG = [16 18:21, 24, 26, 27, 28, 29, 31,  33, 34, 35, 42, 43, 47, 50, 51, 52, 54, 55, 57, 58]; %32 taken out

lags = 500;
EEGpreproc = options.path.preproc.behaviour;  % path to behav data all subjs
load_name = fullfile(EEGpreproc,'behav_data_all_subjs_all3');
load(load_name)
completeSubjectListBehaviour = unique(all_responses(:,12));
%get behavioural subjects who are in subjectList from EEG
[~,~,SubjectListBehaviourEEG] = intersect(subjectListEEG,completeSubjectListBehaviour');
nS = length(SubjectListBehaviourEEG);



[ExpParameters,Correlation] = calculate_collapsed_integration_kernels(all_responses,SubjectListBehaviourEEG,nS, mean_stim_streams, stim_streams, trigger_streams,lags);



figure
subplot(2,2,1)
plot(ExpParameters.sorted.short.tau,ExpParameters.sorted.long.tau,'kd','LineWidth',1,'MarkerSize',8)
%title([condition{con},' ','R= ',num2str(round(R_tau_len(1,2),2)),' ','P= ', num2str(P_tau_len(1,2),2)])
line([0 4],[0 4])
xlabel('short')
ylabel('long')
xlim([0 4]);ylim([0 4]);
tidyfig;

subplot(2,2,2)
plot(ExpParameters.sorted.short.amplitude,ExpParameters.sorted.long.amplitude,'kd','LineWidth',1,'MarkerSize',8)
%title([condition{con},' ','R= ',num2str(round(R_amp_len(1,2),2)),' ','P= ', num2str(P_amp_len(1,2),2)])
line([0 1],[0 1])
xlabel('short')
ylabel('long')
xlim([0 1]);ylim([0 1]);
tidyfig;

subplot(2,2,3)
plot(ExpParameters.sorted.frequent.tau,ExpParameters.sorted.rare.tau,'kd','LineWidth',1,'MarkerSize',8)
%title([condition{con},' ','R= ',num2str(round(R_tau_freq(1,2),2)),' ','P= ', num2str(P_tau_freq(1,2),2)])
line([0 4],[0 4])
xlabel('frequent')
ylabel('rare')
xlim([0 4]);ylim([0 4]);
tidyfig;

subplot(2,2,4)
plot(ExpParameters.sorted.frequent.amplitude,ExpParameters.sorted.rare.amplitude,'kd','LineWidth',1,'MarkerSize',8)
%title([condition{con},' ','R= ',num2str(round(R_amp_freq(1,2),2)),' ','P= ', num2str(P_amp_freq(1,2),2)])
line([0 1],[0 1])
xlabel('frequent')
ylabel('rare')
xlim([0 1]);ylim([0 1]);
tidyfig;

%% STATS TESTING
% frequent vs. rare tau
[h,p,ci,stats] = ttest(ExpParameters.sorted.frequent.tau-ExpParameters.sorted.rare.tau);

% frequent vs. rare A
[h,p,ci,stats] = ttest(ExpParameters.sorted.frequent.amplitude-ExpParameters.sorted.rare.amplitude);

% long vs. short tau
[h,p,ci,stats] = ttest(ExpParameters.sorted.long.tau-ExpParameters.sorted.short.tau);

% long vs. short A
[h,p,ci,stats] = ttest(ExpParameters.sorted.long.amplitude-ExpParameters.sorted.short.amplitude);



end