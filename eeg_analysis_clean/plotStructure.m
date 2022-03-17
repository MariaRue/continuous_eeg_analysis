plotVariables.conditions = {'frequent', 'rare', 'short', 'long'};
plotVariables.diffWavesLabels = {'frequent - rare', 'short - long'};
%% 
plotVariables.figure1.paperSize =  [200 100];
plotVariables.figure1.position = [500, 500, 800, 1290];
plotVariables.figure1.ERP.colour = cbrewer('qual', 'Set1',8);
plotVariables.figure1.LineWidth = 2;
plotVariables.figure1.xlim = [-0.1 0.8];
plotVariables.figure1.ERP.ylim = [-0.3 1];
plotVariables.figure1.xlabel = {'time (s)'};
plotVariables.figure1.ERP.ylabel = {'effect size'};
plotVariables.figure1.diffWave.ylabel = {};
plotVariables.figure1.diffWave.colour = [0.5 0.5 0.5];
plotVariables.figure1.diffWave.ylim = [-0.5 0.4];

plotVariables.figure1.topoPlot.colour = cbrewer('div','RdBu',100);
plotVariables.figure1.topoPlot.zlim = [-0.7 0.7];
%% 
plotVariables.figure2.paperSize =  [200 100];
plotVariables.figure2.position = [500, 500, 800, 1290];
plotVariables.figure2.ERP.colour = cbrewer('qual', 'Set1',8);
plotVariables.figure2.LineWidth = 2;
plotVariables.figure2.xlim = [-4 0.1];
plotVariables.figure2.ERP.ylim = [-5 4];
plotVariables.figure2.xlabel = {'time (s)'};
plotVariables.figure2.ERP.ylabel = {'effect size'};
plotVariables.figure2.diffWave.ylabel = {};
plotVariables.figure2.diffWave.colour = [0.5 0.5 0.5];
plotVariables.figure2.diffWave.ylim = [-2 3];

plotVariables.figure2.topoPlot.colour = cbrewer('div','RdBu',100);
plotVariables.figure2.topoPlot.zlim = [-3 3];
%%
plotVariables.figure3.paperSize =  [200 200];
plotVariables.figure3.position = [200, 1600, 1000, 1000];
plotVariables.figure3.ERP.colour = [0 0 0; 0.3 0.3 0.3];
plotVariables.figure3.topoPlot.colour = cbrewer('div','RdBu',100);
plotVariables.figure3.LineWidth = 2; 
plotVariables.figure3.xlim = [-0.1 0.8];
plotVariables.figure3.ylabel = {'effect size'}; 
plotVariables.figure3.xlabel = {'time (s)'}; 
plotVariables.figure3.JumpRegressor.ylim = [-0.15 0.35]; 
plotVariables.figure3.topoPlot.JumpRegressor.zlim = [-0.25 0.25];
plotVariables.figure3.CohLevelRegressor.ylim = [-0.6 0.8];
plotVariables.figure3.topoPlot.CohLevelRegressor.zlim = [-0.5 0.5];
plotVariables.figure3.PERegressor.ylim = [-0.4 0.8];
plotVariables.figure3.topoPlot.PERegressor.zlim = [-0.5 0.5];
plotVariables.figure3.AbsoluteStimRegressor.ylim = [-0.25 0.2];
plotVariables.figure3.topoPlot.AbsoluteStimRegressor.zlim = [-0.1 0.1]; 
plotVariables.figure3.legend = {'decision relevant', 'decision irrelevant'}; 
%% 

plotVariables.figure4.paperSize =  [200 100];
plotVariables.figure4.position = [500, 500, 700, 500];

plotVariables.figure4.LineWidth = 2; 
plotVariables.figure4.Colour = [0 0 0]; 
plotVariables.figure4.xlim = [-0.1 0.8]; 
plotVariables.figure4.ERP.ylim = [-0.6 0.6]; 
%%
plotVariables.figure5.paperSize =  [500 500];
plotVariables.figure5.position = [200, 1200, 1700, 1000];

plotVariables.figure5.ERP.colour = [0 0 0];
plotVariables.figure5.absoluteStim.colour = [0 0 0; 0.25 0.25 0.25; 0.5 0.5 0.5];
plotVariables.figure5.LineWidth = 2;
plotVariables.figure5.xlim = [-0.1 0.8];
plotVariables.figure5.jumpRegressor.ylim = [-0.05 0.35];
plotVariables.figure5.cohLevelRegressor.ylim = [-0.5 0.5];
plotVariables.figure5.PERegressor.ylim = [-0.3 0.8];
plotVariables.figure5.absoluteStimRegressor.ylim = [-0.2 0.2];
plotVariables.figure5.xlabel = {'time (s)'};
plotVariables.figure5.ERP.ylabel = {'effect size'};

plotVariables.figure5.topoPlot.colour = cbrewer('div','RdBu',100);
plotVariables.figure5.topoPlot.jumpRegressor.zlim = [-0.25 0.25];
plotVariables.figure5.topoPlot.cohLevelRegressor.zlim = [-0.3 0.3];
plotVariables.figure5.topoPlot.PERegressor.zlim = [-0.7 0.7];
plotVariables.figure5.topoPlot.absoluteStimRegressor.zlim = [-0.1 0.1];