function make_figure10(plotVariables, options)

%behavioural figure for trial periods 
subjectListEEG = [16 18:21, 24, 26, 27, 28, 29, 31,  33, 34, 35, 42, 43, 47, 50, 51, 52, 54, 55, 57, 58]; %32 taken out

lags = 500;
EEGpreproc = options.path.preproc.behaviour;  % path to behav data all subjs
load_name = fullfile(EEGpreproc,'behav_data_all_subjs_all3');
load(load_name)
completeSubjectListBehaviour = unique(all_responses(:,12));
%get behavioural subjects who are in subjectList from EEG
[~,~,SubjectListBehaviourEEG] = intersect(subjectListEEG,completeSubjectListBehaviour');
nS = length(SubjectListBehaviourEEG); 
coherence = [0.3 0.4 0.5]; 

% single subject kernels with exponential fits 
[GroupIntegrationKernels, SubjectIntegrationKernels, SignificantTimePoints, ExParameters] = calculate_integration_kernels(all_responses,SubjectListBehaviourEEG,nS, mean_stim_streams, stim_streams, trigger_streams,lags);

selectedSubjects = [1 2 3];

[Model] = calculate_integration_kernel_based_on_model(ExParameters, selectedSubjects);

keyboard; 

% correlation of exponential fits 

plotVariables.conditions;


end 