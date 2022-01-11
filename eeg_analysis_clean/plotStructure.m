plotVariables.conditions = {'frequent', 'rare', 'short', 'long'}; 
plotVariables.diffWavesLabels = {'frequent - rare', 'short - long'}; 

plotVariables.figure1.paperSize =  [200 100]; 
plotVariables.figure1.position = [500, 500, 800, 1290]; 
plotVariables.figure1.ERP.colour = cbrewer('qual', 'Set1',8);
plotVariables.figure1.LineWidth = 2; 
plotVariables.figure1.xlim = [-0.1 0.8]; 
plotVariables.figure1.ERP.ylim = [-0.3 1.1]; 
plotVariables.figure1.xlabel = {'time (s)'}; 
plotVariables.figure1.ERP.ylabel = {'effect size'}; 
plotVariables.figure1.diffWave.ylabel = {};
plotVariables.figure1.diffWave.colour = [0.5 0.5 0.5]; 
plotVariables.figure1.diffWave.ylim = [-0.5 0.4]; 

plotVariables.figure1.topoPlot.colour = cbrewer('div','RdBu',100);
plotVariables.figure1.topoPlot.zlim = [-0.7 0.7];

plotVariables.figure2.paperSize =  [200 100]; 
plotVariables.figure2.position = [500, 500, 800, 1290]; 
plotVariables.figure2.ERP.colour = cbrewer('qual', 'Set1',8);
plotVariables.figure2.LineWidth = 2; 
plotVariables.figure2.xlim = [-4 0.1]; 
plotVariables.figure2.ERP.ylim = [-5 8]; 
plotVariables.figure2.xlabel = {'time (s)'}; 
plotVariables.figure2.ERP.ylabel = {'effect size'}; 
plotVariables.figure2.diffWave.ylabel = {};
plotVariables.figure2.diffWave.colour = [0.5 0.5 0.5]; 
plotVariables.figure2.diffWave.ylim = [-5 5]; 

plotVariables.figure2.topoPlot.colour = cbrewer('div','RdBu',100);
plotVariables.figure2.topoPlot.zlim = [-3 3];
