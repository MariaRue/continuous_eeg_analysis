timelockFA = load(fullfile(options.path.EEG.analysis,'fa_perm_stat_trial_freq_rare.mat'));
timelock = load(fullfile(options.path.EEG.analysis,'perm_stat_trial_freq_rare_correct_response.mat'));

glmFA = load('glm_stat_false_alarm.mat'); 
glm = load('glm_stat_correct_resp.mat');



figure 
subplot(1,2,1) 
cl = cbrewer('qual','Set1',3);


subplot(1,2,1)
hold on
plot(timelock.stat.time,timelock.stat.stat)
plot(glm.stat.time,glm.stat.stat) 
legend('timelock','glm')
title('correct resp')
hold off

subplot(1,2,2)
hold on 
plot(timelockFA.stat.time,timelockFA.stat.stat)
plot(glmFA.stat.time,glmFA.stat.stat) 
legend('timelock','glm')
title('fa')
%time points in which either or all of the channes cp1 cp2 cpz are
%significant