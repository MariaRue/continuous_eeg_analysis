function make_figure4_lilian_posterversion(plotVariables,options)

glmFlag = 'jumps_absolute';



% subject list
subjectList = [16 18:21, 24, 26, 27, 28, 29, 31,  33, 34, 35, 42, 43, 47, 50, 51, 52, 54, 55, 57, 58]; %32 taken out


csdFlag = 0; % 1 for csd transformed data
if csdFlag %temporary fix for bug in first subjects' CSD transform
    subjectList(1) = [];
end

reference = 'LMRM';
nS = length(subjectList); %number of subjects

AbsoluteSimRegressorIdx = 4;

load('elec_field_for_GLM');

electrodesForPermTest =  {'Cz'};

% create fieldtrip structure
timeBins = options.subjectLevelGLM.(glmFlag).regressors(AbsoluteSimRegressorIdx).timeBins/1000; % get timebins in seconds
Fs = 1/(timeBins(2)-timeBins(1)); %sampling rate

%% loop over subjects and load in the betas for each subject into cell array betas_all_subjects

[betas_all_subjects, chanlabels] = load_subject_specific_betas_into_cell_array (subjectList,options, reference, glmFlag, csdFlag, nS);

new_labels = change_electrode_labels(chanlabels);
%% average across the sessions for each subject

betas_all_subjects_sessavg = average_betas_across_sessions (betas_all_subjects, glmFlag, options);


%% now load in the Taus from behavioural analysis (calculated using integration_kernels_for_within_subject_comparisons)

load('behaviouralAnalysis/stored_integration_kernels.mat','allTau','subj_list_behav');

%reduce allTau down to the subjects who are in subjectList
[~,~,ai] = intersect(subjectList,subj_list_behav');
allTau_reduced = allTau(ai,:);

allTau_reducedCollapsed(:,1) = (allTau_reduced(:,1) + allTau_reduced(:,2)) ./2; %frequent average
allTau_reducedCollapsed(:,2) = (allTau_reduced(:,3) + allTau_reduced(:,4)) ./2; %rare average

%and average and demean this for inclusion in Spearman's correlation:
allTau_reducedCollapsed_dm = mean(allTau_reducedCollapsed,2);
allTau_reducedCollapsed_dm = allTau_reducedCollapsed_dm - mean(allTau_reducedCollapsed_dm);

%% Spearmans correlation

B = betas_all_subjects_sessavg{AbsoluteSimRegressorIdx}; %betas of interest
Bcollapsed(:,:,:,1) = (B(:,:,:,1) + B(:,:,:,2)) ./ 2; %frequent average
Bcollapsed(:,:,:,2) = (B(:,:,:,3) + B(:,:,:,4)) ./ 2; %rare average

%first, smooth the betas with a gaussian (FMHM = 50ms)
smoothedBcollapsed = smooth_betas(Bcollapsed,Fs,50,3);

%now, correlate across subjects with tau
[corr_across_subjects] = calculate_spearmans_correlation(smoothedBcollapsed,allTau_reducedCollapsed_dm,nS);

for condition = 1:2
    
    dataCorrelated{condition} = create_fieldTrip_structure(timeBins, new_labels, corr_across_subjects(:,:,condition), elecs);
    
    cfg = [];
    cfg.channel = electrodesForPermTest;
    cfg.avgoverchan = 'yes';
    cfg.nanmean = 'yes';
    
    
    selectedData{condition} = select_electrode_specific_data_for_permtest(cfg, dataCorrelated{condition});
    
    
end

cfg = [];
selectedDataAvg = ft_timelockgrandaverage(cfg, selectedData{:});

%% extract averaged betas in specific timewindow

chanIdx = find(strcmp(selectedData{1}.elec.label,'CZ'));
timeWindow = findc(timeBins,[0.4 0.6]); %time window in samples, 400-600ms
dataForCorrelationPlots = (Bcollapsed(:,chanIdx,timeWindow(1):timeWindow(2),:));
dataForCorrelationPlots = squeeze(mean(dataForCorrelationPlots,3));

%% plot
figure
subplot(2,1,1);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [200 100]);
set(gcf, 'Position',  [500, 500, 700, 500]);

hold on

l(1) = plot(timeBins, selectedData{1}.avg,'-', 'LineWidth', plotVariables.figure4.LineWidth, 'Color', [1 0 0]);
l(2) = plot(timeBins, selectedData{2}.avg,'-', 'LineWidth', plotVariables.figure4.LineWidth, 'Color', [0 0 1]);


ll = line(plotVariables.figure4.xlim,[0 0]);
ll.Color = [0.5 0.5 0.5];

plotVariables.figure4.ylim = [-0.7 0.4];
h = fill([0.4 0.4 0.6 0.6] ,[-0.7 0.4 0.4 -0.7],[0.75 0.75 0.75]);
h.FaceAlpha=0.1;
h.EdgeColor=[0.75 0.75 0.75];
xlim(plotVariables.figure4.xlim);
ylim(plotVariables.figure4.ylim);
xlabel('Time (s)');
ylabel('Spearman''s rho');
leg = legend(l,{'Trial periods frequent' 'Trial periods rare'},'Box','off','Location','southwest');

title(sprintf('Correlation between individual evidence integration decay parameter (tau) \n and EEG correlate of absoluted sensory evidence at sensor Cz'));
tidyfig;

subplot(2,2,3);
scatteroptions = {'MarkerSize' 14 'Color' 'r'};
lineoptions = {'LineWidth' 2 'Color' 'k'};
corroptions = {'type' 'spearman'};
scatter_plus_fit(log(allTau_reducedCollapsed(:,1)),dataForCorrelationPlots(:,1),scatteroptions,lineoptions,[],corroptions);
box off
xlabel('log(Tau), frequent trials');
ylabel(sprintf('Average EEG beta\n(400-600ms)'));
lsline;title('Trial periods frequent');

tidyfig;

subplot(2,2,4);
scatteroptions = {'MarkerSize' 14 'Color' 'b'};
scatter_plus_fit(log(allTau_reducedCollapsed(:,2)),dataForCorrelationPlots(:,2),scatteroptions,lineoptions,[],corroptions);
xlabel('log(Tau), rare trials');
ylabel('');
lsline;title('Trial periods rare');

tidyfig;

box off




