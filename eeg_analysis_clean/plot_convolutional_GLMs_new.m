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



[selectedDataFreq, selectedDataRare, selectedDataShort, selectedDataLong,...
    selectedDataAvgFreq, selectedDataAvgRare, selectedDataAvgShort, selectedDataAvgLong, allDataAvg] ...
    = prepare_data_for_perm_test(betas_all_subjects_sessavg, options.subjectLevelGLM.(glmFlag).regressors, chanlabels, regressorIDx, electrodesForPermTest, nS);
%% run permutation test 


if jumpEvent
[statFrequency] = permutation_testGLM(selectedDataFreq, selectedDataRare, 1);
[statLength] = permutation_testGLM(selectedDataLong, selectedDataShort, 1);

else % response locked  
    
[statFrequency] = permutation_testGLM(selectedDataFreq, selectedDataRare, 0);
[statLength] = permutation_testGLM(selectedDataShort, selectedDataLong, 0);
    
end 


% get correct timepoints for plotting 
TimeBinsRegressor = options.subjectLevelGLM.(glmFlag).regressors(regressorIDx).timeBins/1000; %in seconds becaust permutation test also in seconds

[SignificantTimePointsFrequency] = get_significant_labels_for_plotting_ERPs(statFrequency, TimeBinsRegressor);
[SignificantTimePointsLength] = get_significant_labels_for_plotting_ERPs(statLength, TimeBinsRegressor);

%% calculate difference waves 
[diffWaveFreqVSRare] = calculate_difference_waveform(selectedDataFreq, selectedDataRare);

[diffWaveShortVSLong] = calculate_difference_waveform(selectedDataShort, selectedDataLong);


%% plot ERP averages 



figure;
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [200 100]);
set(gcf, 'Position',  [500, 500, 800, 1290])

cl = cbrewer('qual','Set1',5);

% ERP frequency 
subplot(2,2,1)
hold on 

h = shadedErrorBar(selectedDataAvgFreq.time,selectedDataAvgFreq.avg,diffWaveFreqVSRare.se, 'lineprops', '-k');
h.patch.FaceColor = cl(1,:);
h.mainLine.Color = cl(1,:);
h.mainLine.LineWidth = 2; 

h = shadedErrorBar(selectedDataAvgRare.time,selectedDataAvgRare.avg,diffWaveFreqVSRare.se, 'lineprops', '-k');
h.patch.FaceColor = cl(2,:);
h.mainLine.Color = cl(2,:);
h.mainLine.LineWidth = 2; 


% plot significant time points, if there are any 
if any(SignificantTimePointsFrequency)
    
    plot(selectedDataAvgFreq.time(SignificantTimePointsFrequency), ones(length(SignificantTimePointsFrequency),1), '*k')
    
end 

legend({'frequent', 'rare'})
xlabel('time(s)')
ylabel('effect size')

if jumpEvent
    xlim([-0.1 0.8])
  
else % response locked
    
end 

tidyfig; 

% ERP Length
subplot(2,2,2)
hold on 

h = shadedErrorBar(selectedDataAvgShort.time,selectedDataAvgShort.avg,diffWaveShortVSLong.se, 'lineprops', '-k');
h.patch.FaceColor = cl(3,:);
h.mainLine.Color = cl(3,:);
h.mainLine.LineWidth = 2; 

h = shadedErrorBar(selectedDataAvgLong.time,selectedDataAvgLong.avg,diffWaveShortVSLong.se, 'lineprops', '-k');
h.patch.FaceColor = cl(5,:);
h.mainLine.Color = cl(5,:);
h.mainLine.LineWidth = 2; 

% plot significant time points, if there are any 
if any(SignificantTimePointsLength)
    
    plot(selectedDataAvgFreq.time(SignificantTimePointsLength), ones(length(SignificantTimePointsLength),1), '*k')
    
end 

legend({'short', 'long'})

if jumpEvent
    xlim([-0.1 0.8])

else % response locked
    
end 

tidyfig;

% ERP DiffWave Freq vs Rare 
subplot(2,2,3)
hold on 
h = shadedErrorBar(diffWaveFreqVSRare.time,diffWaveFreqVSRare.avg,diffWaveFreqVSRare.se, 'lineprops', '-k');
h.patch.FaceColor = [0.5 0.5 0.5];
h.mainLine.Color = [0.5 0.5 0.5];
h.mainLine.LineWidth = 2; 

plot(diffWaveShortVSLong.time,zeros(length(diffWaveShortVSLong.time),1),'-k','LineWidth',2)

% plot significant time points, if there are any 
if any(SignificantTimePointsFrequency)
    
    plot(selectedDataAvgFreq.time(SignificantTimePointsFrequency), ones(length(SignificantTimePointsFrequency),1), '*k')
    
end 

title('DiffWave Freq - Rare')
if jumpEvent
    xlim([-0.1 0.8])
    
else % response locked
    
end 


tidyfig; 

% ERP DiffWave Long Vs Short 
subplot(2,2,4)
hold on 
title('DiffWave Short - Long')
h = shadedErrorBar(diffWaveShortVSLong.time,diffWaveShortVSLong.avg,diffWaveShortVSLong.se, 'lineprops', '-k');
h.patch.FaceColor = [0.5 0.5 0.5];
h.mainLine.Color = [0.5 0.5 0.5];
h.mainLine.LineWidth = 2; 

plot(diffWaveShortVSLong.time,zeros(length(diffWaveShortVSLong.time),1),'-k','LineWidth',2)

% plot significant time points
if any(SignificantTimePointsLength)
    
    plot(selectedDataAvgFreq.time(SignificantTimePointsLength), ones(length(SignificantTimePointsLength),1), '*k')
    
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




