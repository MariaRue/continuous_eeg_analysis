function make_figure5(plotVariables,options) 

%this used to be called make_figure2.m, but has been renamed to align with
%the figures for final publication

%%% change regressor indices when using the response or response_coherence
%%% GLM!!!

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

electrodesForPermTest = {'CPz', 'CP1', 'CP2', 'C1', 'Cz', 'C2'};
corrResponseRegressorIDx = 6; % correct response 

jumpEvent = 0; %flag that defines specific variables for jump Event regressors or response locked ones ag

 %% prepare data for correct responses

[corrResponseSelectedData, corrResponseAllDataAvg, corrResponseDiffWave, corrResponseStats, corrResponseSignificantTimePoints] = ...
         prepare_data_for_convolutional_GLM_and_topoplots(glmFlag, options, ...
         subjectList, csdFlag, reference, nS, electrodesForPermTest, corrResponseRegressorIDx, jumpEvent);

%% prepare data for false alarms

falseAlarmRegressorIDx = 7;
[falseAlarmSelectedData, falseAlarmAllDataAvg, falseAlarmDiffWave, falseAlarmStats, falseAlarmSignificantTimePoints] = ...
         prepare_data_for_convolutional_GLM_and_topoplots(glmFlag, options, ...
         subjectList, csdFlag, reference, nS, electrodesForPermTest, falseAlarmRegressorIDx, jumpEvent);

    
%% plot 
figure (2) 
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', plotVariables.figure1.paperSize);
set(gcf, 'Position',  plotVariables.figure1.position); 

%% ERP frequency correct response
subplot(4,3,1)
hold on 
% plot frequent ERP
create_ERP_plot(corrResponseSelectedData.Average.Freq.time, corrResponseSelectedData.Average.Freq.avg,corrResponseDiffWave.Avg.diffWaveFreqVSRare.se,plotVariables.figure2.ERP.colour(1,:),plotVariables.figure2.LineWidth)
% plot rare ERP
create_ERP_plot(corrResponseSelectedData.Average.Rare.time, corrResponseSelectedData.Average.Rare.avg,corrResponseDiffWave.Avg.diffWaveFreqVSRare.se,plotVariables.figure2.ERP.colour(2,:),plotVariables.figure2.LineWidth)

% plot significant time points, if there are any 
if any(corrResponseSignificantTimePoints.Frequency)

    plot_significant_timepoints(corrResponseSelectedData.Average.Freq.time,corrResponseSignificantTimePoints.Frequency, plotVariables.figure2.ERP.ylim(2))

end 
   


xlim(plotVariables.figure2.xlim)
ylim(plotVariables.figure2.ERP.ylim)
legend(plotVariables.conditions(1:2))
 
xlabel(plotVariables.figure2.xlabel)
ylabel(plotVariables.figure2.ERP.ylabel)

tidyfig; 

hold off 
     
 %% ERP length correct response
subplot(4,3,2)
hold on 
% plot frequent ERP
create_ERP_plot(corrResponseSelectedData.Average.Short.time, corrResponseSelectedData.Average.Short.avg,corrResponseDiffWave.Avg.diffWaveShortVSLong.se,plotVariables.figure2.ERP.colour(3,:),plotVariables.figure2.LineWidth)
% plot rare ERP
create_ERP_plot(corrResponseSelectedData.Average.Long.time, corrResponseSelectedData.Average.Long.avg,corrResponseDiffWave.Avg.diffWaveShortVSLong.se,plotVariables.figure2.ERP.colour(5,:),plotVariables.figure2.LineWidth)

% plot significant time points, if there are any 
if any(corrResponseSignificantTimePoints.Length)

    plot_significant_timepoints(corrResponseSelectedData.Average.Short.time,corrResponseSignificantTimePoints.Length, plotVariables.figure2.ERP.ylim(2))

end 
   


xlim(plotVariables.figure2.xlim)
ylim(plotVariables.figure2.ERP.ylim)
legend(plotVariables.conditions(3:4))
 
xlabel(plotVariables.figure2.xlabel)
ylabel(plotVariables.figure2.ERP.ylabel)

tidyfig; 

hold off     
%% diffWave frequent vs rare correct responses 
subplot(4,3,4)
hold on 

create_ERP_plot(corrResponseDiffWave.Avg.diffWaveFreqVSRare.time, corrResponseDiffWave.Avg.diffWaveFreqVSRare.avg,corrResponseDiffWave.Avg.diffWaveFreqVSRare.se,plotVariables.figure2.diffWave.colour,plotVariables.figure2.LineWidth)



plot_x_axis_line(corrResponseDiffWave.Avg.diffWaveFreqVSRare.time, plotVariables.figure2.LineWidth)

xlim(plotVariables.figure2.xlim)
ylim(plotVariables.figure2.diffWave.ylim)
xlabel(plotVariables.figure2.xlabel)
ylabel(plotVariables.diffWavesLabels(1))
     
tidyfig; 

hold off 

%% diffWave short vs long correct responses 
subplot(4,3,5)
hold on 

create_ERP_plot(corrResponseDiffWave.Avg.diffWaveShortVSLong.time, corrResponseDiffWave.Avg.diffWaveShortVSLong.avg,corrResponseDiffWave.Avg.diffWaveShortVSLong.se,plotVariables.figure2.diffWave.colour,plotVariables.figure2.LineWidth)



plot_x_axis_line(corrResponseDiffWave.Avg.diffWaveShortVSLong.time, plotVariables.figure2.LineWidth)

xlim(plotVariables.figure2.xlim)
ylim(plotVariables.figure2.diffWave.ylim)
xlabel(plotVariables.figure2.xlabel)
ylabel(plotVariables.diffWavesLabels(1))
     
tidyfig; 

hold off 
%% ERP frequency false alarm
subplot(4,3,7)
hold on 
% plot frequent ERP
create_ERP_plot(falseAlarmSelectedData.Average.Freq.time, falseAlarmSelectedData.Average.Freq.avg,falseAlarmDiffWave.Avg.diffWaveFreqVSRare.se,plotVariables.figure2.ERP.colour(1,:),plotVariables.figure2.LineWidth)
% plot rare ERP
create_ERP_plot(falseAlarmSelectedData.Average.Rare.time, falseAlarmSelectedData.Average.Rare.avg,falseAlarmDiffWave.Avg.diffWaveFreqVSRare.se,plotVariables.figure2.ERP.colour(2,:),plotVariables.figure2.LineWidth)

% plot significant time points, if there are any 
if any(falseAlarmSignificantTimePoints.Frequency)

    plot_significant_timepoints(falseAlarmSelectedData.Average.Freq.time,falseAlarmSignificantTimePoints.Frequency)

end 
   


xlim(plotVariables.figure2.xlim)
ylim(plotVariables.figure2.ERP.ylim)
legend(plotVariables.conditions(1:2))
 
xlabel(plotVariables.figure2.xlabel)
ylabel(plotVariables.figure2.ERP.ylabel)

tidyfig; 

hold off 
 %% ERP length false alarm
subplot(4,3,8)
hold on 
% plot frequent ERP
create_ERP_plot(falseAlarmSelectedData.Average.Short.time, falseAlarmSelectedData.Average.Short.avg,falseAlarmDiffWave.Avg.diffWaveShortVSLong.se,plotVariables.figure2.ERP.colour(3,:),plotVariables.figure2.LineWidth)
% plot rare ERP
create_ERP_plot(falseAlarmSelectedData.Average.Long.time, falseAlarmSelectedData.Average.Long.avg,falseAlarmDiffWave.Avg.diffWaveShortVSLong.se,plotVariables.figure2.ERP.colour(5,:),plotVariables.figure2.LineWidth)

% plot significant time points, if there are any 
if any(falseAlarmSignificantTimePoints.Length)

    plot_significant_timepoints(falseAlarmSelectedData.Average.Short.time,falseAlarmSignificantTimePoints.Length)

end 
   


xlim(plotVariables.figure2.xlim)
ylim(plotVariables.figure2.ERP.ylim)
legend(plotVariables.conditions(3:4))
 
xlabel(plotVariables.figure2.xlabel)
ylabel(plotVariables.figure2.ERP.ylabel)

tidyfig; 

hold off  
%% diffWave frequent vs rare false alarm 
subplot(4,3,10)
hold on 

create_ERP_plot(falseAlarmDiffWave.Avg.diffWaveFreqVSRare.time, falseAlarmDiffWave.Avg.diffWaveFreqVSRare.avg,falseAlarmDiffWave.Avg.diffWaveFreqVSRare.se,plotVariables.figure2.diffWave.colour,plotVariables.figure2.LineWidth)



plot_x_axis_line(falseAlarmDiffWave.Avg.diffWaveFreqVSRare.time, plotVariables.figure2.LineWidth)

xlim(plotVariables.figure2.xlim)
ylim(plotVariables.figure2.diffWave.ylim)
xlabel(plotVariables.figure2.xlabel)
ylabel(plotVariables.diffWavesLabels(1))
     
tidyfig; 

hold off 
%% diffWave short vs long false alarm
subplot(4,3,11)
hold on 

create_ERP_plot(falseAlarmDiffWave.Avg.diffWaveShortVSLong.time, falseAlarmDiffWave.Avg.diffWaveShortVSLong.avg,falseAlarmDiffWave.Avg.diffWaveShortVSLong.se,plotVariables.figure2.diffWave.colour,plotVariables.figure2.LineWidth)



plot_x_axis_line(falseAlarmDiffWave.Avg.diffWaveShortVSLong.time, plotVariables.figure2.LineWidth)

xlim(plotVariables.figure2.xlim)
ylim(plotVariables.figure2.diffWave.ylim)
xlabel(plotVariables.figure2.xlabel)
ylabel(plotVariables.diffWavesLabels(1))
     
tidyfig; 

hold off 

%% topo plot correct responses 

    %first subplot is displaying title of regressor
    
    subplot(4,3,3);
    text(0.5,0.5,options.subjectLevelGLM.(glmFlag).regressors(corrResponseRegressorIDx).name,'FontSize',20);axis off
    
    
    
    subplot(4,3,6);
    create_topo_plot(corrResponseAllDataAvg,-0.5, 0, plotVariables.figure2.topoPlot.zlim, electrodesForPermTest, plotVariables.figure1.topoPlot.colour)
    

%% topo plot false alarm responses 

    %first subplot is displaying title of regressor
    
    subplot(4,3,9);
    text(0.5,0.5,options.subjectLevelGLM.(glmFlag).regressors(falseAlarmRegressorIDx).name,'FontSize',20);axis off
    
    
    
    subplot(4,3,12);
    create_topo_plot(falseAlarmAllDataAvg,-0.5, 0, plotVariables.figure2.topoPlot.zlim, electrodesForPermTest, plotVariables.figure1.topoPlot.colour)
    

end 
