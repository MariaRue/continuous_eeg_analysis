function make_figure2(plotVariables, options)

%this used to be known as make_figure5, but has been renamed to align with
%the figure naming in the final submission


glmFlag = 'all_regressors';



% subject list
subjectList = [16 18:21, 24, 26, 27, 28, 29, 31,  33, 34, 35, 42, 43, 47, 50, 51, 52, 54, 55, 57, 58]; %32 taken out
%subjectList = [62:64,66,68,70]; % vertical motion only

csdFlag = 0; % 1 for csd transformed data
if csdFlag %temporary fix for bug in first subjects' CSD transform
    subjectList(1) = [];
end

reference = 'LMRM';
nS = length(subjectList); %number of subjects


jumpEvent = 1; %flag that defines specific variables for jump Event regressors or response locked ones ag

%% prepare data for plotting
JumpRegressorIDx = 1;
electrodesForPermTest = {'CPz', 'CP1', 'CP2'};
[JumpSelectedData, JumpAllDataAvg, JumpDiffWave, JumpStats,JumpSignificantTimePoints] = ...
    prepare_data_for_convolutional_GLM_and_topoplots(glmFlag, options, ...
    subjectList, csdFlag, reference, nS, electrodesForPermTest, JumpRegressorIDx, jumpEvent);

%% prepare data for plotting
CohLevelRegressorIDx = 2;
electrodesForPermTest = {'CPz', 'CP1', 'CP2'};
[CohLevelSelectedData, CohLevelAllDataAvg, CohLevelDiffWave, CohLevelStats, CohLevelSignificantTimePoints] = ...
    prepare_data_for_convolutional_GLM_and_topoplots(glmFlag, options, ...
    subjectList, csdFlag, reference, nS, electrodesForPermTest, CohLevelRegressorIDx, jumpEvent);

%% prepare data for plotting
PERegressorIDx = 3;
electrodesForPermTest = {'CPz', 'CP1', 'CP2'};
[PESelectedData, PEAllDataAvg, PEDiffWave, PEStats, PESignificantTimePoints] = ...
    prepare_data_for_convolutional_GLM_and_topoplots(glmFlag, options, ...
    subjectList, csdFlag, reference, nS, electrodesForPermTest, PERegressorIDx, jumpEvent);

%% prepare data for plotting
AbsolutedStimRegressorIDx = 4;
electrodesForPermTest = {'CPz', 'CP1', 'CP2'};
[CPAbsolutedStimSelectedData, CPAbsolutedStimAllDataAvg, CPAbsolutedStimDiffWave, CPAbsolutedStimStats, CPAbsolutedStimSignificantTimePoints] = ...
    prepare_data_for_convolutional_GLM_and_topoplots(glmFlag, options, ...
    subjectList, csdFlag, reference, nS, electrodesForPermTest, AbsolutedStimRegressorIDx, jumpEvent);

%% prepare data for plotting
electrodesForPermTest = {'Pz', 'P1', 'P2'};
[PAbsolutedStimSelectedData, PAbsolutedStimAllDataAvg, PAbsolutedStimDiffWave, PAbsolutedStimStats, PAbsolutedStimSignificantTimePoints] = ...
    prepare_data_for_convolutional_GLM_and_topoplots(glmFlag, options, ...
    subjectList, csdFlag, reference, nS, electrodesForPermTest, AbsolutedStimRegressorIDx, jumpEvent);

%% prepare data for plotting
electrodesForPermTest = {'POz', 'PO3', 'PO4'};
[POAbsolutedStimSelectedData, POAbsolutedStimAllDataAvg, POAbsolutedStimDiffWave, POAbsolutedStimStats, POAbsolutedStimSignificantTimePoints] = ...
    prepare_data_for_convolutional_GLM_and_topoplots(glmFlag, options, ...
    subjectList, csdFlag, reference, nS, electrodesForPermTest, AbsolutedStimRegressorIDx, jumpEvent);

%% plot
figure
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', plotVariables.figure5.paperSize);
set(gcf, 'Position',  plotVariables.figure5.position);

%% ERP Jump regressor
subplot(4,5,1)
hold on
% plot average jump regressor
create_ERP_plot(JumpSelectedData.Average.All.time, JumpSelectedData.Average.All.avg,JumpSelectedData.Average.All.se, plotVariables.figure5.ERP.colour,plotVariables.figure5.LineWidth)

xlim(plotVariables.figure5.xlim)
ylim(plotVariables.figure5.jumpRegressor.ylim)


xlabel(plotVariables.figure5.xlabel)
ylabel(plotVariables.figure5.ERP.ylabel)

tidyfig;

hold off

%% ERP coherence level

subplot(4,5,6)
hold on
% plot coherence level regressor
create_ERP_plot(CohLevelSelectedData.Average.All.time, CohLevelSelectedData.Average.All.avg, CohLevelSelectedData.Average.All.se, plotVariables.figure5.ERP.colour,plotVariables.figure5.LineWidth)

xlim(plotVariables.figure5.xlim)
ylim(plotVariables.figure5.cohLevelRegressor.ylim)

xlabel(plotVariables.figure5.xlabel)
ylabel(plotVariables.figure5.ERP.ylabel)

tidyfig;

hold off

%% ERP PE

subplot(4,5,11)
hold on
% plot PE regressor
create_ERP_plot(PESelectedData.Average.All.time, PESelectedData.Average.All.avg, PESelectedData.Average.All.se, plotVariables.figure5.ERP.colour,plotVariables.figure5.LineWidth)

xlim(plotVariables.figure5.xlim)
ylim(plotVariables.figure5.PERegressor.ylim)

xlabel(plotVariables.figure5.xlabel)
ylabel(plotVariables.figure5.ERP.ylabel)

tidyfig;

hold off

%% ERP absolute stim

subplot(4,5,16)
hold on
% plot absolute stim regressor
create_ERP_plot(CPAbsolutedStimSelectedData.Average.All.time, CPAbsolutedStimSelectedData.Average.All.avg, CPAbsolutedStimSelectedData.Average.All.se, plotVariables.figure5.absoluteStim.colour(1,:),plotVariables.figure5.LineWidth)
create_ERP_plot(PAbsolutedStimSelectedData.Average.All.time, PAbsolutedStimSelectedData.Average.All.avg, PAbsolutedStimSelectedData.Average.All.se, plotVariables.figure5.absoluteStim.colour(2,:),plotVariables.figure5.LineWidth)
create_ERP_plot(POAbsolutedStimSelectedData.Average.All.time, POAbsolutedStimSelectedData.Average.All.avg, POAbsolutedStimSelectedData.Average.All.se, plotVariables.figure5.absoluteStim.colour(3,:),plotVariables.figure5.LineWidth)

xlim(plotVariables.figure5.xlim)
ylim(plotVariables.figure5.absoluteStimRegressor.ylim)

xlabel(plotVariables.figure5.xlabel)
ylabel(plotVariables.figure5.ERP.ylabel)

tidyfig;

hold off


%%  topo plot
%first subplot is displaying title of regressor
electrodesForPermTest = {'CPz', 'CP1', 'CP2'};

subplot(4,5,2);
text(0.5,0.5,options.subjectLevelGLM.(glmFlag).regressors(JumpRegressorIDx).name,'FontSize',20);axis off

subplot(4,5,3);
create_topo_plot(JumpAllDataAvg,0.25, 0.35, plotVariables.figure5.topoPlot.jumpRegressor.zlim, electrodesForPermTest, plotVariables.figure5.topoPlot.colour);
tidyfig; 
%% 
subplot(4,5,7)
text(0.5,0.5,options.subjectLevelGLM.(glmFlag).regressors(CohLevelRegressorIDx).name,'FontSize',20);axis off

subplot(4,5,8);
create_topo_plot(CohLevelAllDataAvg,0.35, 0.45, plotVariables.figure5.topoPlot.cohLevelRegressor.zlim, electrodesForPermTest, plotVariables.figure5.topoPlot.colour);
tidyfig; 
%%
subplot(4,5,12)
text(0.5,0.5,options.subjectLevelGLM.(glmFlag).regressors(PERegressorIDx).name,'FontSize',20);axis off

subplot(4,5,13);
create_topo_plot(PEAllDataAvg,0.25, 0.35, plotVariables.figure5.topoPlot.PERegressor.zlim, electrodesForPermTest, plotVariables.figure5.topoPlot.colour);
tidyfig; 
%%
subplot(4,5,17)
text(0.5,0.5,options.subjectLevelGLM.(glmFlag).regressors(AbsolutedStimRegressorIDx).name,'FontSize',20);axis off

electrodesForPermTest = {'Pz', 'P1', 'P2'};
subplot(4,5,18);
create_topo_plot(PAbsolutedStimAllDataAvg,0.1, 0.2, plotVariables.figure5.topoPlot.absoluteStimRegressor.zlim, electrodesForPermTest, plotVariables.figure5.topoPlot.colour);
tidyfig; 

electrodesForPermTest = {'POz', 'PO3', 'PO4'};
subplot(4,5,19);
create_topo_plot(POAbsolutedStimAllDataAvg,0.25, 0.35, plotVariables.figure5.topoPlot.absoluteStimRegressor.zlim, electrodesForPermTest, plotVariables.figure5.topoPlot.colour);
tidyfig; 

electrodesForPermTest = {'CPz', 'CP1', 'CP2'};
subplot(4,5,20);
create_topo_plot(CPAbsolutedStimAllDataAvg,0.45, 0.55, plotVariables.figure5.topoPlot.absoluteStimRegressor.zlim, electrodesForPermTest, plotVariables.figure5.topoPlot.colour);
tidyfig; 


end
