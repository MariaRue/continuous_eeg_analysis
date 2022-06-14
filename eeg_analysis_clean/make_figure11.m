function make_figure11(plotVariables,options)

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

%load('behaviouralAnalysis/stored_integration_kernels.mat','allTau','subj_list_behav');

%behavioural figure for trial periods
lags = 500;
EEGpreproc = options.path.preproc.behaviour;  % path to behav data all subjs
load_name = fullfile(EEGpreproc,'behav_data_all_subjs_all3');
load(load_name)
completeSubjectListBehaviour = unique(all_responses(:,12));
%get behavioural subjects who are in subjectList from EEG
[~,~,SubjectListBehaviourEEG] = intersect(subjectList,completeSubjectListBehaviour');

[~, ~, ~, ExpParameters] =  calculate_integration_kernels(all_responses,SubjectListBehaviourEEG,nS, mean_stim_streams, stim_streams, trigger_streams,lags);

allTau_reduced = squeeze(ExpParameters.parameters(:,2,:))';

%and average and demean this for inclusion in Spearman's correlation:
allTau_reduced_dm = mean(allTau_reduced,2);
allTau_reduced_dm = allTau_reduced_dm - mean(allTau_reduced_dm);

%% Spearmans correlation

%first, smooth the betas with a gaussian (FMHM = 75ms)
smoothedBetas = smooth_betas(betas_all_subjects_sessavg{AbsoluteSimRegressorIdx},Fs,75,3);

%now, correlate across subjects with tau
[corr_across_subjects] = calculate_spearmans_correlation(smoothedBetas,allTau_reduced_dm,nS);

for condition = 1:4
    
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
timeWindow = findc(timeBins,[0.42 0.75]); %time window in samples, 420-750ms
dataForCorrelationPlots = (betas_all_subjects_sessavg{AbsoluteSimRegressorIdx}(:,chanIdx,timeWindow(1):timeWindow(2),:));
dataForCorrelationPlots = squeeze(mean(dataForCorrelationPlots,3));

%% plot
figure
subplot(2,1,1);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [200 100]);
set(gcf, 'Position',  [446  226 1023  876]);

hold on
for condition = 1:4
    ll(condition) = plot(timeBins, selectedData{condition}.avg,'-', 'LineWidth', plotVariables.figure4.LineWidth, 'Color', plotVariables.figure4.Colour(condition,:))
end
leg = legend(ll,options.conditionLabels); set(leg,'box','off');
hold off

xlim(plotVariables.figure4.xlim)
ylim(plotVariables.figure4.ERP.ylim)
xlabel('Time (s)');
ylabel('Spearman''s rho');

for condition = 1:4
    subplot(2,4,4+condition);
    scatteroptions = {'MarkerSize' 14 'Color' plotVariables.figure4.Colour(condition,:)};
    lineoptions = {'LineWidth' 2 'Color' plotVariables.figure4.Colour(condition,:)};
    corroptions = {'type' 'spearman'};
    scatter_plus_fit(log(allTau_reduced(:,condition)),dataForCorrelationPlots(:,condition),scatteroptions,lineoptions,[],corroptions);
    xlabel('log(Tau)'); 
    if condition==1; ylabel('Average EEG beta (420-750ms)'); else; ylabel(''); end
    title(options.conditionLabels{condition});
    box off;
    tidyfig;
end


%%

end