function make_figure2(plotVariables,options)

glmFlag = 'vertical_jumps_absolute';



% subject list
subjectList = [62:64,66,68,70]; % vertical motion only

csdFlag = 0; % 1 for csd transformed data
if csdFlag %temporary fix for bug in first subjects' CSD transform
    subjectList(1) = [];
end

reference = 'LMRM';
nS = length(subjectList); %number of subjects

electrodesForPermTest = {'CPz', 'CP1', 'CP2'};
regressorIDx = 6; % correct response

jumpEvent = 1; %flag that defines specific variables for jump Event regressors or response locked ones ag


%% prepare data for ploting
HorzJumpRegressorIDx = 1;
[HorzJumpSelectedData, HorzJumpAllDataAvg, ~, ~, ~] = ...
    prepare_data_for_convolutional_GLM_and_topoplots(glmFlag, options, ...
    subjectList, csdFlag, reference, nS, electrodesForPermTest, HorzJumpRegressorIDx, jumpEvent);

VertJumpRegressorIDx = 8;
[VertJumpSelectedData, VertJumpAllDataAvg, ~, ~, ~] = ...
    prepare_data_for_convolutional_GLM_and_topoplots(glmFlag, options, ...
    subjectList, csdFlag, reference, nS, electrodesForPermTest, VertJumpRegressorIDx, jumpEvent);
 

[diffWave.Avg.diffWaveJump, ~] = calculate_difference_waveform(HorzJumpSelectedData.subjectLevel.All, VertJumpSelectedData.subjectLevel.All);

[stats.Jump] = permutation_testGLM(VertJumpSelectedData.subjectLevel.All, HorzJumpSelectedData.subjectLevel.All, jumpEvent);

TimeBinsRegressor = options.subjectLevelGLM.(glmFlag).regressors(HorzJumpRegressorIDx).timeBins/1000; %in seconds because permutation test also in seconds

[SignificantTimePoints.Jump] = get_significant_labels_for_plotting_ERPs(stats.Jump, TimeBinsRegressor);
%%
HorzCohLevelRegressorIDx = 2;
[HorzCohLevelSelectedData, HorzCohLevelAllDataAvg, ~, ~, ~] = ...
    prepare_data_for_convolutional_GLM_and_topoplots(glmFlag, options, ...
    subjectList, csdFlag, reference, nS, electrodesForPermTest, HorzCohLevelRegressorIDx, jumpEvent);

VertCohLevelRegressorIDx = 9;
[VertCohLevelSelectedData, VertCohLevelAllDataAvg, ~, ~, ~] = ...
    prepare_data_for_convolutional_GLM_and_topoplots(glmFlag, options, ...
    subjectList, csdFlag, reference, nS, electrodesForPermTest, VertCohLevelRegressorIDx, jumpEvent);

[diffWave.Avg.diffWaveCohLevel, ~] = calculate_difference_waveform(HorzCohLevelSelectedData.subjectLevel.All, VertCohLevelSelectedData.subjectLevel.All);

[stats.CohLevel] = permutation_testGLM(VertCohLevelSelectedData.subjectLevel.All, HorzCohLevelSelectedData.subjectLevel.All, jumpEvent);

TimeBinsRegressor = options.subjectLevelGLM.(glmFlag).regressors(HorzCohLevelRegressorIDx).timeBins/1000; %in seconds because permutation test also in seconds

[SignificantTimePoints.CohLevel] = get_significant_labels_for_plotting_ERPs(stats.CohLevel, TimeBinsRegressor);

%%
HorzPERegressorIDx = 3;
[HorzPESelectedData, HorzPEAllDataAvg, ~, ~, ~] = ...
    prepare_data_for_convolutional_GLM_and_topoplots(glmFlag, options, ...
    subjectList, csdFlag, reference, nS, electrodesForPermTest, HorzPERegressorIDx, jumpEvent);

% because we evaluated the horizontal regressor from -1.5 to 1.5, but the
% vertical one only from -1 to 1.5 we need to shorten the horizontal
% regressor by 0.5 so that the time axis for both regressors is the same
% for calculating the correct SE for plotting. 

for subject = 1:length(HorzPESelectedData.subjectLevel.All)
HorzPESelectedData.subjectLevel.All{subject}.time(1:50) = [];
HorzPESelectedData.subjectLevel.All{subject}.avg(1:50) = [];
end 
HorzPESelectedData.Average.All.time(1:50) = [];
HorzPESelectedData.Average.All.avg(1:50) = [];

VertPERegressorIDx = 10;
[VertPESelectedData, VertPEAllDataAvg, ~, ~, ~] = ...
    prepare_data_for_convolutional_GLM_and_topoplots(glmFlag, options, ...
    subjectList, csdFlag, reference, nS, electrodesForPermTest, VertPERegressorIDx, jumpEvent);

[diffWave.Avg.diffWavePE, ~] = calculate_difference_waveform(HorzPESelectedData.subjectLevel.All, VertPESelectedData.subjectLevel.All);

[stats.PE] = permutation_testGLM(VertPESelectedData.subjectLevel.All, HorzPESelectedData.subjectLevel.All, jumpEvent);

TimeBinsRegressor = options.subjectLevelGLM.(glmFlag).regressors(VertPERegressorIDx).timeBins/1000; %in seconds because permutation test also in seconds

[SignificantTimePoints.PE] = get_significant_labels_for_plotting_ERPs(stats.PE, TimeBinsRegressor);

%%
HorzAbsoluteStimRegressorIDx = 4;
[HorzAbsoluteStimSelectedData, HorzAbsoluteStimAllDataAvg, ~, ~, ~] = ...
    prepare_data_for_convolutional_GLM_and_topoplots(glmFlag, options, ...
    subjectList, csdFlag, reference, nS, electrodesForPermTest, HorzAbsoluteStimRegressorIDx, jumpEvent);

VertAbsoluteStimRegressorIDx = 11;
[VertAbsoluteStimSelectedData, VertAbsoluteStimAllDataAvg, ~, ~, ~] = ...
    prepare_data_for_convolutional_GLM_and_topoplots(glmFlag, options, ...
    subjectList, csdFlag, reference, nS, electrodesForPermTest, VertAbsoluteStimRegressorIDx, jumpEvent);

[diffWave.Avg.diffWaveAbsoluteStim, ~] = calculate_difference_waveform(HorzAbsoluteStimSelectedData.subjectLevel.All, VertAbsoluteStimSelectedData.subjectLevel.All);


[stats.AbsoluteStim] = permutation_testGLM(VertAbsoluteStimSelectedData.subjectLevel.All, HorzAbsoluteStimSelectedData.subjectLevel.All, jumpEvent);

TimeBinsRegressor = options.subjectLevelGLM.(glmFlag).regressors(HorzAbsoluteStimRegressorIDx).timeBins/1000; %in seconds because permutation test also in seconds

[SignificantTimePoints.AbsoluteStim] = get_significant_labels_for_plotting_ERPs(stats.AbsoluteStim, TimeBinsRegressor);

%%

figure 
set(gcf, 'PaperUnits', 'inches');
%set(gcf, 'PaperSize', plotVariables.figure3.paperSize);
set(gcf, 'Position',  [200, 1600, 1000, 1000]);

%%
subplot(7,4,2)
text(0,0.3,'horizontal motion','FontSize',20);axis off

subplot(7,4,3)
text(0,0.3,'vertical motion','FontSize',20);axis off

%% Jump Regressor
subplot(7,4,5)
text(0,0.5,options.subjectLevelGLM.(glmFlag).regressors(HorzJumpRegressorIDx).name,'FontSize',20);axis off

subplot(7,4,6)
create_topo_plot(HorzJumpAllDataAvg,0.3, 0.4, plotVariables.figure3.topoPlot.JumpRegressor.zlim, electrodesForPermTest, plotVariables.figure3.topoPlot.colour);
tidyfig; 

subplot(7,4,7)
create_topo_plot(VertJumpAllDataAvg,0.3, 0.4, plotVariables.figure3.topoPlot.JumpRegressor.zlim, electrodesForPermTest, plotVariables.figure3.topoPlot.colour);
tidyfig; 

subplot(7,4,8)
hold on

create_ERP_plot(HorzJumpSelectedData.Average.All.time, HorzJumpSelectedData.Average.All.avg, diffWave.Avg.diffWaveJump.se, plotVariables.figure3.ERP.colour(1,:),plotVariables.figure3.LineWidth)
create_ERP_plot(VertJumpSelectedData.Average.All.time, VertJumpSelectedData.Average.All.avg, diffWave.Avg.diffWaveJump.se, plotVariables.figure3.ERP.colour(2,:),plotVariables.figure3.LineWidth)

% plot significant time points, if there are any 
if any(SignificantTimePoints.Jump)

    plot_significant_timepoints(HorzJumpSelectedData.Average.All.time, SignificantTimePoints.Jump, plotVariables.figure3.JumpRegressor.ylim(2))

end 

legend(plotVariables.figure3.legend);

xlim(plotVariables.figure3.xlim)
ylim(plotVariables.figure3.JumpRegressor.ylim)

xlabel(plotVariables.figure3.xlabel)
ylabel(plotVariables.figure3.ylabel)

tidyfig;

hold off

%% Coh Level
subplot(7,4,9)
text(0,0.5,options.subjectLevelGLM.(glmFlag).regressors(HorzCohLevelRegressorIDx).name,'FontSize',20);axis off

subplot(7,4,10)
create_topo_plot(HorzCohLevelAllDataAvg,0.2, 0.3, plotVariables.figure3.topoPlot.CohLevelRegressor.zlim, electrodesForPermTest, plotVariables.figure3.topoPlot.colour);
tidyfig; 

subplot(7,4,11)
create_topo_plot(VertCohLevelAllDataAvg,0.2, 0.3, plotVariables.figure3.topoPlot.CohLevelRegressor.zlim, electrodesForPermTest, plotVariables.figure3.topoPlot.colour);
tidyfig; 

subplot(7,4,12)
hold on

create_ERP_plot(HorzCohLevelSelectedData.Average.All.time, HorzCohLevelSelectedData.Average.All.avg, diffWave.Avg.diffWaveCohLevel.se, plotVariables.figure3.ERP.colour(1,:),plotVariables.figure3.LineWidth)
create_ERP_plot(VertCohLevelSelectedData.Average.All.time, VertCohLevelSelectedData.Average.All.avg, diffWave.Avg.diffWaveCohLevel.se, plotVariables.figure3.ERP.colour(2,:),plotVariables.figure3.LineWidth)

if any(SignificantTimePoints.CohLevel)

    plot_significant_timepoints(HorzCohLevelSelectedData.Average.All.time, SignificantTimePoints.CohLevel, plotVariables.figure3.CohLevelRegressor.ylim(2))

end 

xlim(plotVariables.figure3.xlim)
ylim(plotVariables.figure3.CohLevelRegressor.ylim)

xlabel(plotVariables.figure3.xlabel)
ylabel(plotVariables.figure3.ylabel)

tidyfig;

hold off

%% PE
subplot(7,4,13)
text(0,0.5,options.subjectLevelGLM.(glmFlag).regressors(HorzPERegressorIDx).name,'FontSize',20);axis off

subplot(7,4,14)
create_topo_plot(HorzPEAllDataAvg,0.3, 0.4, plotVariables.figure3.topoPlot.PERegressor.zlim, electrodesForPermTest, plotVariables.figure3.topoPlot.colour);
tidyfig; 

subplot(7,4,15)
create_topo_plot(VertPEAllDataAvg,0.3, 0.4, plotVariables.figure3.topoPlot.PERegressor.zlim, electrodesForPermTest, plotVariables.figure3.topoPlot.colour);
tidyfig; 

subplot(7,4,16)
hold on

create_ERP_plot(HorzPESelectedData.Average.All.time, HorzPESelectedData.Average.All.avg, diffWave.Avg.diffWavePE.se, plotVariables.figure3.ERP.colour(1,:),plotVariables.figure3.LineWidth)
create_ERP_plot(VertPESelectedData.Average.All.time, VertPESelectedData.Average.All.avg, diffWave.Avg.diffWavePE.se, plotVariables.figure3.ERP.colour(2,:),plotVariables.figure3.LineWidth)

if any(SignificantTimePoints.PE)

    plot_significant_timepoints(HorzPESelectedData.Average.All.time, SignificantTimePoints.PE, plotVariables.figure3.PERegressor.ylim(2))

end 

xlim(plotVariables.figure3.xlim)
ylim(plotVariables.figure3.PERegressor.ylim)

xlabel(plotVariables.figure3.xlabel)
ylabel(plotVariables.figure3.ylabel)

tidyfig;

hold off

%% absolute stimulus 

subplot(7,4,21)
text(0,0.5,options.subjectLevelGLM.(glmFlag).regressors(HorzAbsoluteStimRegressorIDx).name,'FontSize',20);axis off

subplot(7,4,18)
create_topo_plot(HorzAbsoluteStimAllDataAvg,0.1, 0.2, plotVariables.figure3.topoPlot.AbsoluteStimRegressor.zlim, electrodesForPermTest, plotVariables.figure3.topoPlot.colour);
tidyfig; 

subplot(7,4,19)
create_topo_plot(VertAbsoluteStimAllDataAvg,0.1, 0.2, plotVariables.figure3.topoPlot.AbsoluteStimRegressor.zlim, electrodesForPermTest, plotVariables.figure3.topoPlot.colour);
tidyfig; 

subplot(7,4,22)
create_topo_plot(HorzAbsoluteStimAllDataAvg,0.35, 0.45, plotVariables.figure3.topoPlot.AbsoluteStimRegressor.zlim, electrodesForPermTest, plotVariables.figure3.topoPlot.colour);
tidyfig; 

subplot(7,4,23)
create_topo_plot(VertAbsoluteStimAllDataAvg,0.35, 0.45, plotVariables.figure3.topoPlot.AbsoluteStimRegressor.zlim, electrodesForPermTest, plotVariables.figure3.topoPlot.colour);
tidyfig; 

subplot(7,4,24)
hold on

create_ERP_plot(HorzAbsoluteStimSelectedData.Average.All.time, HorzAbsoluteStimSelectedData.Average.All.avg, diffWave.Avg.diffWaveAbsoluteStim.se, plotVariables.figure3.ERP.colour(1,:),plotVariables.figure3.LineWidth)
create_ERP_plot(VertAbsoluteStimSelectedData.Average.All.time, VertAbsoluteStimSelectedData.Average.All.avg, diffWave.Avg.diffWaveAbsoluteStim.se, plotVariables.figure3.ERP.colour(2,:),plotVariables.figure3.LineWidth)

if any(SignificantTimePoints.AbsoluteStim)

    plot_significant_timepoints(HorzAbsoluteStimSelectedData.Average.All.time, SignificantTimePoints.AbsoluteStim, plotVariables.figure3.AbsoluteStimRegressor.ylim(2))

end 



xlim(plotVariables.figure3.xlim)
ylim(plotVariables.figure3.AbsoluteStimRegressor.ylim)

xlabel(plotVariables.figure3.xlabel)
ylabel(plotVariables.figure3.ylabel)

tidyfig;

hold off

subplot(7,4,26)
create_topo_plot(HorzAbsoluteStimAllDataAvg,0.45, 0.55, plotVariables.figure3.topoPlot.AbsoluteStimRegressor.zlim, electrodesForPermTest, plotVariables.figure3.topoPlot.colour);
tidyfig; 

subplot(7,4,27)
create_topo_plot(VertAbsoluteStimAllDataAvg,0.45, 0.55, plotVariables.figure3.topoPlot.AbsoluteStimRegressor.zlim, electrodesForPermTest, plotVariables.figure3.topoPlot.colour);
tidyfig; 


end