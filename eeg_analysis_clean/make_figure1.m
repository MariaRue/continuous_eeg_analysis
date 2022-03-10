function make_figure1(plotVariables, options)




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

electrodesForPermTest = {'CPz', 'CP1', 'CP2'};
regressorIDx = 3;

jumpEvent = 1; %flag that defines specific variables for jump Event regressors or response locked ones ag

%% prepare data for ploting 

[selectedData, allDataAvg, diffWave, stats, SignificantTimePoints] = ...
         prepare_data_for_convolutional_GLM_and_topoplots(glmFlag, options, ...
         subjectList, csdFlag, reference, nS, electrodesForPermTest, regressorIDx, jumpEvent);
     
     
%% plot 

figure (1) 
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', plotVariables.figure1.paperSize);
set(gcf, 'Position',  plotVariables.figure1.position); 

%% ERP frequency 
subplot(3,2,1)
hold on 
% plot frequent ERP
create_ERP_plot(selectedData.Average.Freq.time, selectedData.Average.Freq.avg,diffWave.Avg.diffWaveFreqVSRare.se,plotVariables.figure1.ERP.colour(1,:),plotVariables.figure1.LineWidth)
% plot rare ERP
create_ERP_plot(selectedData.Average.Rare.time, selectedData.Average.Rare.avg,diffWave.Avg.diffWaveFreqVSRare.se,plotVariables.figure1.ERP.colour(2,:),plotVariables.figure1.LineWidth)

% plot significant time points, if there are any 
if any(SignificantTimePoints.Frequency)

    plot_significant_timepoints(selectedData.Average.Freq.time,SignificantTimePoints.Frequency, plotVariables.figure1.ERP.ylim(2))

end 
   


xlim(plotVariables.figure1.xlim)
ylim(plotVariables.figure1.ERP.ylim)
legend(plotVariables.conditions(1:2))
 
xlabel(plotVariables.figure1.xlabel)
ylabel(plotVariables.figure1.ERP.ylabel)

tidyfig; 

hold off 

%% ERP Length 

subplot(3,2,2)
hold on 
% plot short ERP
create_ERP_plot(selectedData.Average.Short.time, selectedData.Average.Short.avg,diffWave.Avg.diffWaveShortVSLong.se,plotVariables.figure1.ERP.colour(3,:),plotVariables.figure1.LineWidth)

% plot long ERP
create_ERP_plot(selectedData.Average.Long.time, selectedData.Average.Long.avg,diffWave.Avg.diffWaveShortVSLong.se,plotVariables.figure1.ERP.colour(5,:),plotVariables.figure1.LineWidth)

% plot significant time points, if there are any 
if any(SignificantTimePoints.Length)

    plot_significant_timepoints(selectedData.Average.Short.time,SignificantTimePoints.Length, plotVariables.figure1.ERP.ylim(2))

end 
   
xlim(plotVariables.figure1.xlim)
ylim(plotVariables.figure1.ERP.ylim)
legend(plotVariables.conditions(3:4))
    
xlabel(plotVariables.figure1.xlabel)
ylabel(plotVariables.figure1.ERP.ylabel)

tidyfig; 

hold off 

%% diffWave frequent vs rare
subplot(3,2,3)
hold on 

create_ERP_plot(diffWave.Avg.diffWaveFreqVSRare.time, diffWave.Avg.diffWaveFreqVSRare.avg,diffWave.Avg.diffWaveFreqVSRare.se,plotVariables.figure1.diffWave.colour,plotVariables.figure1.LineWidth)

if any(SignificantTimePoints.Frequency)

    plot_significant_timepoints(selectedData.Average.Freq.time,SignificantTimePoints.Frequency, plotVariables.figure1.diffWave.ylim(2))

end 

plot_x_axis_line(diffWave.Avg.diffWaveFreqVSRare.time, plotVariables.figure1.LineWidth)

xlim(plotVariables.figure1.xlim)
ylim(plotVariables.figure1.diffWave.ylim)
xlabel(plotVariables.figure1.xlabel)
ylabel(plotVariables.diffWavesLabels(1))
     
tidyfig; 

hold off 


%% diffWave short vs long
subplot(3,2,4)
hold on 

create_ERP_plot(diffWave.Avg.diffWaveShortVSLong.time, diffWave.Avg.diffWaveShortVSLong.avg,diffWave.Avg.diffWaveShortVSLong.se,plotVariables.figure1.diffWave.colour,plotVariables.figure1.LineWidth)

if any(SignificantTimePoints.Length)

    plot_significant_timepoints(selectedData.Average.Short.time,SignificantTimePoints.Length, plotVariables.figure1.diffWave.ylim(2))

end 

plot_x_axis_line(diffWave.Avg.diffWaveShortVSLong.time, plotVariables.figure1.LineWidth)

xlim(plotVariables.figure1.xlim)
ylim(plotVariables.figure1.diffWave.ylim)

     
xlabel(plotVariables.figure1.xlabel)
ylabel(plotVariables.diffWavesLabels(2))
tidyfig; 

hold off 

%%  topo plot 


    %first subplot is displaying title of regressor
    
    subplot(3,2,5);
    text(0.5,0.5,options.subjectLevelGLM.(glmFlag).regressors(regressorIDx).name,'FontSize',20);axis off
    
    
    
    subplot(3,2,6);
    create_topo_plot(allDataAvg,0.25, 0.35, plotVariables.figure1.topoPlot.zlim, electrodesForPermTest, plotVariables.figure1.topoPlot.colour)
    

end 