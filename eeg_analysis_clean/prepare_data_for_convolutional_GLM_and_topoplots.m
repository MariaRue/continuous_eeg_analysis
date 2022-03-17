
function [selectedData, allDataAvg, diffWave, stats, SignificantTimePoints] = ...
         prepare_data_for_convolutional_GLM_and_topoplots(glmFlag, options, ...
         subjectList, csdFlag, reference, nS, electrodesForPermTest, regressorIDx, jumpEvent)

%% loop over subjects and load in the betas for each subject into cell array betas_all_subjects

[betas_all_subjects, chanlabels] = load_subject_specific_betas_into_cell_array (subjectList,options, reference, glmFlag, csdFlag, nS);

%% average across the sessions for each subject

betas_all_subjects_sessavg = average_betas_across_sessions (betas_all_subjects, glmFlag, options);

%% prepare data for permutation test 

% electrode specification

[selectedData, allDataAvg] ...
    = prepare_data_for_perm_test(betas_all_subjects_sessavg, options.subjectLevelGLM.(glmFlag).regressors, chanlabels, regressorIDx, electrodesForPermTest, nS);

%% calculate difference waves 
[diffWave.Avg.diffWaveFreqVSRare, ~] = calculate_difference_waveform(selectedData.subjectLevel.Freq, selectedData.subjectLevel.Rare);

[diffWave.Avg.diffWaveShortVSLong, ~] = calculate_difference_waveform(selectedData.subjectLevel.Short, selectedData.subjectLevel.Long);

%this is going to be the input for the permutation test for the interaction
% which is freqShort - rareShort  AND freqLong - rareLong 
[diffWave.Avg.diffWaveFreqShortVSRareShort, diffWave.subjectLevel.diffWaveFreqShortVSRareShort] = calculate_difference_waveform(selectedData.subjectLevel.FreqShort, selectedData.subjectLevel.RareShort);

[diffWave.Avg.diffWaveFreqLongVSRareLong, diffWave.subjectLevel.diffWaveFreqLongVSRareLong] = calculate_difference_waveform(selectedData.subjectLevel.FreqLong, selectedData.subjectLevel.RareLong);

% this is the difference wave for the interaction - which means (freqShort
% - rareShort) - (freqLong - rareLong)
[diffWave.Avg.diffWaveInteraction, ~] = calculate_difference_waveform(diffWave.subjectLevel.diffWaveFreqShortVSRareShort, diffWave.subjectLevel.diffWaveFreqLongVSRareLong);



%% run permutation test 



[stats.statFrequency] = permutation_testGLM(selectedData.subjectLevel.Freq, selectedData.subjectLevel.Rare, jumpEvent);
[stats.statLength] = permutation_testGLM(selectedData.subjectLevel.Long, selectedData.subjectLevel.Short, jumpEvent);

% we do actually have 1 cluster under alpha two-tailed 0.025 with p = 0.018
% but also do have 2 clusters under alpha 0.05 (one tailed) with p = 0.037
% and p = 0.026  - for now I am only plotting the one cluster because that
% was our original criteria - but we might want to include the othes as
% well? 


[stats.statInteraction] = permutation_testGLM(diffWave.subjectLevel.diffWaveFreqShortVSRareShort, diffWave.subjectLevel.diffWaveFreqLongVSRareLong, jumpEvent);




% get correct timepoints for plotting 
TimeBinsRegressor = options.subjectLevelGLM.(glmFlag).regressors(regressorIDx).timeBins/1000; %in seconds because permutation test also in seconds

[SignificantTimePoints.Frequency] = get_significant_labels_for_plotting_ERPs(stats.statFrequency, TimeBinsRegressor);
[SignificantTimePoints.Length] = get_significant_labels_for_plotting_ERPs(stats.statLength, TimeBinsRegressor);
[SignificantTimePoints.Interaction] = get_significant_labels_for_plotting_ERPs(stats.statInteraction, TimeBinsRegressor);



end 