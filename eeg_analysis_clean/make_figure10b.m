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



[EperimentParameters,Correlation] = calculate_collapsed_integration_kernels(all_responses,SubjectListBehaviourEEG,nS, mean_stim_streams, stim_streams, trigger_streams,lags);


keyboard;

figure

subplot(2,2,1)
plot(tau_short,tau_long,'kd','LineWidth',1,'MarkerSize',8)
title([condition{con},' ','R= ',num2str(round(R_tau_len(1,2),2)),' ','P= ', num2str(P_tau_len(1,2),2)])
xlabel('short')
ylabel('long')
tidyfig;

subplot(2,2,2)
plot(amp_short,amp_long,'kd','LineWidth',1,'MarkerSize',8)
title([condition{con},' ','R= ',num2str(round(R_amp_len(1,2),2)),' ','P= ', num2str(P_amp_len(1,2),2)])
xlabel('short')
ylabel('long')
tidyfig;

subplot(2,2,3)
plot(tau_freq,tau_rare,'kd','LineWidth',1,'MarkerSize',8)
title([condition{con},' ','R= ',num2str(round(R_tau_freq(1,2),2)),' ','P= ', num2str(P_tau_freq(1,2),2)])
xlabel('frequent')
ylabel('rare')
tidyfig;

subplot(2,2,4)
plot(amp_freq,amp_rare,'kd','LineWidth',1,'MarkerSize',8)
title([condition{con},' ','R= ',num2str(round(R_amp_freq(1,2),2)),' ','P= ', num2str(P_amp_freq(1,2),2)])
xlabel('frequent')
ylabel('rare')
tidyfig;


end 