timelockCondition = 'buttonPress';
reference = 'LMRM';
csdFlag = 0;
rtFlag = 0;

baseline = 0;
faFlag = 0;
subjectList = [16, 18, 20, 26, 28, 29, 31, 34, 35,  42, 43, 47, 50, 51, 52, 54, 55, 57, 58];


numberSubjects = length(subjectList);
options = continuous_RDK_set_options('iMac');

per_condition = 1;

for subjects = 1:length(subjectList)
    
    subID = subjectList(subjects);
    [details,paths] =  conrdk_subjects(subID,options,reference,csdFlag);
   
    data_load = load(paths.(reference).singleTrial.tf_analysis.appendedData.(timelockCondition));
    data{subjects} = data_load.tfdataAppend;
    
    
        % change labels to easycap
    [easy_cap_labels] = change_electrode_labels(data{subjects}.label);
    data{subjects}.label = easy_cap_labels;
    

    % separate left and rightward trials 
    [dataLeft{subjects}, dataRight{subjects}] = select_trials_for_tf_analysis(data{subjects},per_condition);
    
    % baseline correct 
    
     [dataLeftBaseCorr{subjects}, dataRightBaseCorr{subjects}] = baseline_correct_tf_data(dataLeft{subjects}, dataRight{subjects},per_condition);

%      cfg = [];
% %cfg.baseline     = [-0.5 -0.1];
% %cfg.baselinetype = 'absolute';
% %cfg.xlim         = [0.9 1.3];
% cfg.zlim         = [-1e-27 1e-27];
% % cfg.ylim         = [5 40];
% % cfg.marker       = 'on';
%  cfg.layout       = 'easycapM3';
% % cfg.colorbar     = 'yes';
% figure
% ft_topoplotTFR(cfg,dataLeft{subjects});




% 
end



% param for number of subjects has to be inserted for input
% diff wave for left vs right response 
[diffWaveLeftRight] = calculate_diff_wave_tf_data(dataLeftBaseCorr,dataRightBaseCorr,per_condition);

for_diffWaves = 0;

[grandavgLeft, grandavgRight, grandAvgDiffWave] = calculate_grandaverage_tf_data(dataLeftBaseCorr, dataRightBaseCorr, diffWaveLeftRight, per_condition,for_diffWaves);



%% plotting

% for  DIFFWAVE
cfg = [];

cfg.xlim         = [-0.4 -0.2];
cfg.zlim         ='maxabs';
cfg.ylim         = [10 35];
cfg.marker       = 'on';
cfg.layout       = 'easycapM3';
cfg.colorbar     = 'yes';
cfg.marker = 'labels';
figure
ft_topoplotTFR(cfg,grandAvgDiffWave);
% % % 

cfg = [];
cfg.ylim         = [15 35]; 
cfg.xlim         = [-3 0.2];
cfg.zlim         = [-2.5 2.5];
%cfg.baseline     = [-4 -3];
% cfg.baselinetype = 'db';
cfg.channel      = {'C4' 'CP4'};
cfg.layout       = 'easycapM3';
cfg.parameter = 'powspctrm';

figure
ft_singleplotTFR(cfg, grandAvgDiffWave);
ylabel ('freq Hz')
xlabel('time in s 0 = button press')


% for normal data 
% freq stats analysis - t-stats contrast left-right response 
per_condition = 1;
[output] =  tStatistics_left_minus_right(dataLeftBaseCorr, dataRightBaseCorr, length(subjectList), per_condition);

    cfg = [];

cfg.xlim         = [-0.4 -0.2];
cfg.zlim         = 'maxabs';
 cfg.ylim         = [10  35];
cfg.marker       = 'on';
 cfg.layout       = 'easycapM11';
cfg.colorbar     = 'yes';
cfg.parameter   = 'stat';
cfg.marker = 'labels';
figure
ft_topoplotTFR(cfg,output{1});


cfg = [];


% cfg.ylim         = [15 35]; 
cfg.xlim         = [-3 0.2];
cfg.zlim         = [-3 3];

cfg.channel      = {'C3' 'CP3'};
cfg.avgoverchan = 'yes';
cfg.parameter   = 'stat';
cfg.layout       = 'easycapM11';
figure
ft_singleplotTFR(cfg, output);
ylabel ('freq Hz')
xlabel('time in s 0 = button press')
tidyfig;






% per condition 
conditionLabels = {'shortFreq', 'longFreq', 'shortRare', 'longRare'};

for condition = 1:4

cfg = [];


% cfg.ylim         = [15 35]; 
cfg.xlim         = [-3 0.2];
cfg.zlim         =  'maxabs'; %[-3 3];

cfg.channel      = {'C3' 'CP3'};
cfg.avgoverchan = 'yes';
cfg.parameter   = 'stat';
cfg.layout       = 'easycapM11';
figure
ft_singleplotTFR(cfg, output{condition});

 
 
 title(conditionLabels{condition})
 ylabel('frequency')
 xlabel('time (s) 0 = response')
    tidyfig;
end 




% compare conditions 
[dataPerConditionLeft, dataPerConditionRight] = compare_conditions_for_different_frequencies(grandavgLeft, grandavgRight, 0,1);


cl = cbrewer('div','RdBu', 12);
cl = cl([2 4 12 10],:);

subplot(1,2,1)
hold on 
for condition = 1:4

    plot(dataPerConditionLeft{condition}.time, squeeze(dataPerConditionLeft{condition}.powspctrm),'Color',cl(condition,:),'LineWidth',4)
    
end 
legend({'shortFreg','longFreq', 'shortRare', 'longRare'} )
hold off
title('avg over C4 CP4 and frequencies 10 to 15 for left resp')
ylabel('power')
xlabel('time (s) response at 0s')
tidyfig;

subplot(1,2,2)
hold on 
for condition = 1:4

    plot(dataPerConditionRight{condition}.time, squeeze(dataPerConditionRight{condition}.powspctrm),'Color',cl(condition,:),'LineWidth',4)
    
end 
legend({'shortFreg','longFreq', 'shortRare', 'longRare'} )
hold off
title('avg over C3 CP3 and frequencies 28 to 33 for right resp')
ylabel('power')
xlabel('time (s) response at 0s')
tidyfig;


[dataPerConditionLeft, dataPerConditionRight] = compare_conditions_for_different_frequencies(grandavgLeft, grandavgRight, 0,0);


subplot(2,2,1)
hold on 
for condition = 1:2

    plot(dataPerConditionLeft{condition}.time, squeeze(dataPerConditionLeft{condition}.powspctrm),'Color',cl(condition,:),'LineWidth',4)
    
end 
legend({'Freq', 'Rare'} )
hold off
title('avg over C4 CP4 and 13Hz to 30Hz left resp')
ylabel('power')
xlabel('time (s) response at 0s')
tidyfig;

subplot(2,2,2)
hold on 
for condition = 1:2

    plot(dataPerConditionRight{condition}.time, squeeze(dataPerConditionRight{condition}.powspctrm),'Color',cl(condition,:),'LineWidth',4)
    
end 
legend({'Freq', 'Rare'} )
hold off
title('avg over C3 CP3 and 13Hz to 30Hz right resp')
ylabel('power')
xlabel('time (s) response at 0s')
tidyfig;


subplot(2,2,3)
hold on 
for condition = 3:4

    plot(dataPerConditionLeft{condition}.time, squeeze(dataPerConditionLeft{condition}.powspctrm),'Color',cl(condition,:),'LineWidth',4)
    
end 
legend({'short', 'long'} )
hold off

ylabel('power')
xlabel('time (s) response at 0s')
tidyfig;

subplot(2,2,4)
hold on 
for condition = 3:4

    plot(dataPerConditionRight{condition}.time, squeeze(dataPerConditionRight{condition}.powspctrm),'Color',cl(condition,:),'LineWidth',4)
    
end 
legend({'short', 'long'} )
hold off

ylabel('power')
xlabel('time (s) response at 0s')
tidyfig;

% calculate diff waves for different conditionns
[diffWaveLeftAvg diffWaveRightAvg] = calculate_diffWaves_tf_data_between_conditions(dataLeftBaseCorr, dataRightBaseCorr);

figure

subplot(1,2,1)
hold on 
plot(diffWaveLeftAvg{1}.time, squeeze(diffWaveLeftAvg{1}.powspctrm),'Color',cl(1,:),'LineWidth',4)
plot(diffWaveRightAvg{1}.time, squeeze(diffWaveRightAvg{1}.powspctrm),'Color',cl(2,:),'LineWidth',4)
plot(diffWaveRightAvg{1}.time,zeros(size(diffWaveRightAvg{1}.time)),'k-','LineWidth',4)
hold off 
title('frequent - rare') 
legend({'left response', 'right response'})
ylabel('power avg 13-30Hz') 
xlabel('time (s) 0 = response')
tidyfig;

subplot(1,2,2)
hold on 
plot(diffWaveLeftAvg{1}.time, squeeze(diffWaveLeftAvg{2}.powspctrm),'Color',cl(3,:),'LineWidth',4)
plot(diffWaveRightAvg{1}.time, squeeze(diffWaveRightAvg{2}.powspctrm),'Color',cl(4,:),'LineWidth',4)
plot(diffWaveRightAvg{1}.time,zeros(size(diffWaveRightAvg{1}.time)),'k-','LineWidth',4)
hold off 
title('short - long') 
legend({'left response', 'right response'})
ylabel('power avg 13-30Hz') 
xlabel('time (s) 0 = response')
tidyfig;

%% average across central electrodes and left/right 

[diffWave] = average_over_central_electrodes_tf(dataLeftBaseCorr, dataRightBaseCorr);

cl = cbrewer('div','RdBu', 12);
cl = cl([2 4 12 10],:);

subplot(2,2,1)
hold on 
plot(dataPerCondition{1}.time, squeeze(dataPerCondition{1}.powspctrm),'Color',cl(1,:),'LineWidth',4)
plot(dataPerCondition{2}.time, squeeze(dataPerCondition{2}.powspctrm),'Color',cl(2,:),'LineWidth',4)
legend({'frequent' 'rare'})
title('frequent vs rare') 
ylabel('power avg 13-30Hz') 
xlabel('time (s) 0 = response')

subplot(2,2,2)
hold on 
plot(dataPerCondition{3}.time, squeeze(dataPerCondition{3}.powspctrm),'Color',cl(3,:),'LineWidth',4)
plot(dataPerCondition{4}.time, squeeze(dataPerCondition{4}.powspctrm),'Color',cl(4,:),'LineWidth',4)
legend({'short' 'long'})
title('short vs long') 
ylabel('power avg 13-30Hz') 
xlabel('time (s) 0 = response')


subplot(2,2,3)
hold on 
plot(diffWave{1}.time, squeeze(diffWave{1}.powspctrm),'Color',cl(1,:),'LineWidth',4)
plot(diffWave{1}.time,zeros(size(diffWave{1}.time)),'k-','LineWidth',4)
hold off 
title('frequent - rare') 
ylabel('power avg 13-30Hz') 
xlabel('time (s) 0 = response')
tidyfig;

subplot(2,2,4)
hold on 
plot(diffWave{2}.time, squeeze(diffWave{2}.powspctrm),'Color',cl(3,:),'LineWidth',4)
plot(diffWave{2}.time,zeros(size(diffWave{2}.time)),'k-','LineWidth',4)
hold off 
title('short - long') 
ylabel('power avg 13-30Hz') 
xlabel('time (s) 0 = response')
tidyfig;
