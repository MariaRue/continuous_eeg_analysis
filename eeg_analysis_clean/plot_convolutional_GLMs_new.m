clear all;
addpath(genpath(pwd))

glmFlag = 'jumps_absolute';

options = continuous_RDK_set_options('LTHiMac');

% subject list
subjectList = [16 18:21, 24, 26, 27, 28, 29, 31,  33, 34, 35, 42, 43, 47, 50, 51, 52, 54, 55, 57, 58]; %32 taken out
%subjectList = [62:64,66,68,70]; % vertical motion only

csdFlag = 0; % 1 for csd transformed data
if csdFlag %temporary fix for bug in first subjects' CSD transform
    subjectList(1) = [];
end

reference = 'LMRM';
nS = length(subjectList); %number of subjects

electrodesForPermTest = {'CPz', 'CP1', 'CP2'};
regressorIDx = 3;

jumpEvent = 1; %flag that defines specific variables for jump Event regressors or response locked ones 



%% loop over subjects and load in the betas for each subject into cell array betas_all_subjects

[betas_all_subjects, chanlabels] = load_subject_specific_betas_into_cell_array (subjectList,options, reference, glmFlag, csdFlag, nS);

%% average across the sessions for each subject

betas_all_subjects_sessavg = average_betas_across_sessions (betas_all_subjects, glmFlag, options);

%% prepare data for permutation test 

% electrode specification

[selectedData, allDataAvg] ...
    = prepare_data_for_perm_test(betas_all_subjects_sessavg, options.subjectLevelGLM.(glmFlag).regressors, chanlabels, regressorIDx, electrodesForPermTest, nS);

%% calculate difference waves 
[diffWaveFreqVSRareAvg, ~] = calculate_difference_waveform(selectedData.subjectLevel.Freq, selectedData.subjectLevel.Rare);

[diffWaveShortVSLongAvg, ~] = calculate_difference_waveform(selectedData.subjectLevel.Short, selectedData.subjectLevel.Long);

%this is going to be the input for the permutation test for the interaction
% which is freqShort - rareShort  AND freqLong - rareLong 
[diffWaveFreqShortVSRareShortAvg, diffWaveFreqShortVSRareShortSjs] = calculate_difference_waveform(selectedData.subjectLevel.FreqShort, selectedData.subjectLevel.RareShort);

[diffWaveFreqLongVSRareLongAvg, diffWaveFreqLongVSRareLongSjs] = calculate_difference_waveform(selectedData.subjectLevel.FreqLong, selectedData.subjectLevel.RareLong);

% this is the difference wave for the interaction - which means (freqShort
% - rareShort) - (freqLong - rareLong)
[diffWaveInteractionAvg, ~] = calculate_difference_waveform(diffWaveFreqShortVSRareShortSjs, diffWaveFreqLongVSRareLongSjs);


%% run permutation test 


if jumpEvent
[statFrequency] = permutation_testGLM(selectedData.subjectLevel.Freq, selectedData.subjectLevel.Rare, 1);
[statLength] = permutation_testGLM(selectedData.subjectLevel.Long, selectedData.subjectLevel.Short, 1);
[statInteraction] = permutation_testGLM(diffWaveFreqShortVSRareShortSjs, diffWaveFreqLongVSRareLongSjs, 1);

else % response locked  
    
[statFrequency] = permutation_testGLM(selectedData.subjectLevel.Freq, selectedData.subjectLevel.Rare, 0);
[statLength] = permutation_testGLM(selectedData.subjectLevel.Short, selectedData.subjectLevel.Long, 0);

% we do actually have 1 cluster under alpha two-tailed 0.025 with p = 0.018
% but also do have 2 clusters under alpha 0.05 (one tailed) with p = 0.037
% and p = 0.026  - for now I am only plotting the one cluster because that
% was our original criteria - but we might want to include the othes as
% well? 
[statInteraction] = permutation_testGLM(diffWaveFreqShortVSRareShortSjs, diffWaveFreqLongVSRareLongSjs, 0);
end 


% get correct timepoints for plotting 
TimeBinsRegressor = options.subjectLevelGLM.(glmFlag).regressors(regressorIDx).timeBins/1000; %in seconds becaust permutation test also in seconds

[SignificantTimePointsFrequency] = get_significant_labels_for_plotting_ERPs(statFrequency, TimeBinsRegressor);
[SignificantTimePointsLength] = get_significant_labels_for_plotting_ERPs(statLength, TimeBinsRegressor);
[SignificantTimePointsInteraction] = get_significant_labels_for_plotting_ERPs(statInteraction, TimeBinsRegressor);




%% plot ERP averages 



figure;
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [200 100]);
set(gcf, 'Position',  [500, 500, 800, 1290])

cl = cbrewer('qual','Set1',8);

% ERP frequency 
subplot(2,3,1)
hold on 

h = shadedErrorBar(selectedData.Average.Freq.time,selectedData.Average.Freq.avg,diffWaveFreqVSRareAvg.se, 'lineprops', '-k');
h.patch.FaceColor = cl(1,:);
h.mainLine.Color = cl(1,:);
h.mainLine.LineWidth = 2; 

h = shadedErrorBar(selectedData.Average.Rare.time,selectedData.Average.Rare.avg,diffWaveFreqVSRareAvg.se, 'lineprops', '-k');
h.patch.FaceColor = cl(2,:);
h.mainLine.Color = cl(2,:);
h.mainLine.LineWidth = 2; 


% plot significant time points, if there are any 
if any(SignificantTimePointsFrequency)
    
    plot(selectedData.Average.Freq.time(SignificantTimePointsFrequency), ones(length(SignificantTimePointsFrequency),1), '*k')
    
end 

legend({'frequent', 'rare'})
xlabel('time(s)')
ylabel('effect size')

if jumpEvent
    xlim([-0.1 0.8])
  
else % response locked5
    
end 

tidyfig; 

% ERP Length
subplot(2,3,2)
hold on 

h = shadedErrorBar(selectedData.Average.Short.time,selectedData.Average.Short.avg,diffWaveShortVSLongAvg.se, 'lineprops', '-k');
h.patch.FaceColor = cl(3,:);
h.mainLine.Color = cl(3,:);
h.mainLine.LineWidth = 2; 

h = shadedErrorBar(selectedData.Average.Long.time,selectedData.Average.Long.avg,diffWaveShortVSLongAvg.se, 'lineprops', '-k');
h.patch.FaceColor = cl(5,:);
h.mainLine.Color = cl(5,:);
h.mainLine.LineWidth = 2; 

% plot significant time points, if there are any 
if any(SignificantTimePointsLength)
    
    plot(selectedData.Average.Freq.time(SignificantTimePointsLength), ones(length(SignificantTimePointsLength),1), '*k')
    
end 

legend({'short', 'long'})

if jumpEvent
    xlim([-0.1 0.8])

else % response locked
    
end 

tidyfig;

% ERP Interaction (FreqShort-RareShort) and (FreqLong - RareLong) 
%%%%%%%%%%%
subplot(2,3,3) 
hold on 
h = shadedErrorBar(diffWaveFreqShortVSRareShortAvg.time,diffWaveFreqShortVSRareShortAvg.avg,diffWaveFreqShortVSRareShortAvg.se, 'lineprops', '-k');
h.patch.FaceColor = cl(4,:);
h.mainLine.Color = cl(4,:);
h.mainLine.LineWidth = 2; 

h = shadedErrorBar(diffWaveFreqLongVSRareLongAvg.time,diffWaveFreqLongVSRareLongAvg.avg,diffWaveFreqLongVSRareLongAvg.se, 'lineprops', '-k');
h.patch.FaceColor = cl(8,:);
h.mainLine.Color = cl(8,:);
h.mainLine.LineWidth = 2; 

plot(diffWaveShortVSLongAvg.time,zeros(length(diffWaveShortVSLongAvg.time),1),'-k','LineWidth',2)

% plot significant time points, if there are any 
if any(SignificantTimePointsInteraction)
    
    plot(diffWaveFreqLongVSRareLongAvg.time(SignificantTimePointsInteraction), ones(length(SignificantTimePointsInteraction),1), '*k')
    
end 


legend('freqShort - rareShort', 'freqLong - rareLong')
xlim([-0.1 0.8])
tidyfig; 

% ERP DiffWave Freq vs Rare 
subplot(2,3,4)
hold on 
h = shadedErrorBar(diffWaveFreqVSRareAvg.time,diffWaveFreqVSRareAvg.avg,diffWaveFreqVSRareAvg.se, 'lineprops', '-k');
h.patch.FaceColor = [0.5 0.5 0.5];
h.mainLine.Color = [0.5 0.5 0.5];
h.mainLine.LineWidth = 2; 

plot(diffWaveShortVSLongAvg.time,zeros(length(diffWaveShortVSLongAvg.time),1),'-k','LineWidth',2)

% plot significant time points, if there are any 
if any(SignificantTimePointsFrequency)
    
    plot(selectedData.Average.Freq.time(SignificantTimePointsFrequency), ones(length(SignificantTimePointsFrequency),1), '*k')
    
end 

title('DiffWave Freq - Rare')
if jumpEvent
    xlim([-0.1 0.8])
    
else % response locked
    
end 


tidyfig; 

% ERP DiffWave Long Vs Short 
subplot(2,3,5)
hold on 
title('DiffWave Short - Long')
h = shadedErrorBar(diffWaveShortVSLongAvg.time,diffWaveShortVSLongAvg.avg,diffWaveShortVSLongAvg.se, 'lineprops', '-k');
h.patch.FaceColor = [0.5 0.5 0.5];
h.mainLine.Color = [0.5 0.5 0.5];
h.mainLine.LineWidth = 2; 

plot(diffWaveShortVSLongAvg.time,zeros(length(diffWaveShortVSLongAvg.time),1),'-k','LineWidth',2)

% plot significant time points
if any(SignificantTimePointsLength)
    
    plot(selectedData.Average.Freq.time(SignificantTimePointsLength), ones(length(SignificantTimePointsLength),1), '*k')
    
end 

if jumpEvent
    xlim([-0.1 0.8])
    
else % response locked
    
end 

tidyfig;

% ERP DiffWave Interaction
subplot(2,3,6)
hold on 
title('DiffWave Interaction')
h = shadedErrorBar(diffWaveInteractionAvg.time,diffWaveInteractionAvg.avg,diffWaveInteractionAvg.se, 'lineprops', '-k');
h.patch.FaceColor = [0.5 0.5 0.5];
h.mainLine.Color = [0.5 0.5 0.5];
h.mainLine.LineWidth = 2; 

plot(diffWaveInteractionAvg.time,zeros(length(diffWaveInteractionAvg.time),1),'-k','LineWidth',2)

% plot significant time points
if any(SignificantTimePointsInteraction)
    
    plot(diffWaveInteractionAvg.time(SignificantTimePointsInteraction), ones(length(SignificantTimePointsInteraction),1), '*k')
    
end 

if jumpEvent
    xlim([-0.1 0.8])
    
else % response locked
    
end 

tidyfig;


%% plot topoplots 


% ft_struct is the data structure with fields required by the ft topoplot function
cl = cbrewer('div','RdBu',100);


figure

fig_id = 1; % subplot index
start_times = [0.250  ]; %start times of plotting windows, in ms
end_times   = [0.350 ]; %end times of plotting windows, in ms

% 1:length(time_idx) % loop through regressors
    
    %cax_lim = max(abs(prctile(squash(ft_struct.avg(1:60,:)),[1 99])));
    
    %first subplot is displaying title of regressor
    subplot(1,length(start_times)+1,fig_id);
    text(0.5,0.5,options.subjectLevelGLM.(glmFlag).regressors(regressorIDx).name,'FontSize',20);axis off
    fig_id = fig_id + 1;
    
    
    for t = 1:length(start_times) % loop through time windows
        
        
        
        cfg = [];
        cfg.highlight = 'on';
        cfg.highlightchannel = electrodesForPermTest;
        cfg.highlightsymbol = '^';
        cfg.xlim = [start_times(t) end_times(t)];  % time window for which we create a topo plot
        %cfg.zlim = [-cax_lim cax_lim]; %colorbar scale
        cfg.layout = 'easycapM1.mat';
        cfg.colormap = flip(cl);
        subplot(1,length(start_times)+1,fig_id);
        fig_id = fig_id + 1;
        
        ft_topoplotER(cfg,allDataAvg); colorbar
        
        
    end




