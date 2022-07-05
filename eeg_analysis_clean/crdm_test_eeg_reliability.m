function crdm_test_eeg_reliability( options )
%CRDM_TEST_EEG_RELIABILITY Loads the results of a convGLM: a matrix of
%betas per subject, session, and sample for each regressor, and performs
%several reliability tests on them. 

glmFlag = 'jumps_absolute';
subjectList = [16 18:21, 24, 26, 27, 28, 29, 31,  33, 34, 35, 42, 43, 47, 50, 51, 52, 54, 55, 57, 58]; %32 taken out
csdFlag = 0; % 1 for csd transformed data
reference = 'LMRM';
nSubjects = numel(subjectList);

% load betas per regressor, subject, electrode, session, and
% condition/block
[allBetas, chanlabels] = load_subject_specific_betas_into_cell_array(...
    subjectList, options, reference, glmFlag, csdFlag, nSubjects);

% for each regressor, we want to average over the relevant sensors, and
% then compute reliability in different ways:
% - reliability of a single condition across sessions (5min per session)
% - reliability of the session average across sessions (20min per session)
% - reliability of 10min of one condition across sessions (10min per
% session, but only 3 sessions as we have to concatenate 2
% - reliability of 1 contrast over sessions (e.g., frequent vs rare, 10min
% per condition per session)

% starting with the last option
freqIdx = 1:2;
rareIdx = 3:4;

%% jump regressor
jump.idx = 1;
jump.sensors = {'CPZ', 'CP1', 'CP2'};

jump.freq.data = select_data_for_reliability(allBetas{jump.idx}, jump.sensors, freqIdx, chanlabels);
[jump.freq.maxCorr_within, jump.freq.bestLag_within, jump.freq.maxCorr_between, jump.freq.bestLag_between, ...
    jump.freq.R, jump.freq.R_bounds] = crdm_eeg_reliability(jump.freq.data, 150);

jump.rare.data = select_data_for_reliability(allBetas{jump.idx}, jump.sensors, rareIdx, chanlabels);
[jump.rare.maxCorr_within, jump.rare.bestLag_within, jump.rare.maxCorr_between, jump.rare.bestLag_between, ...
    jump.rare.R, jump.rare.R_bounds] = crdm_eeg_reliability(jump.rare.data, 150);

%% coherence level at jump regressor
cohLevel.idx = 2;
cohLevel.sensors = {'CPZ', 'CP1', 'CP2'};
cohLevel.freq.data = select_data_for_reliability(allBetas{cohLevel.idx}, cohLevel.sensors, freqIdx, chanlabels);
[cohLevel.freq.maxCorr_within, cohLevel.freq.bestLag_within, cohLevel.freq.maxCorr_between, cohLevel.freq.bestLag_between, ...
    cohLevel.freq.R, cohLevel.freq.R_bounds] = crdm_eeg_reliability(cohLevel.freq.data, 150);
cohLevel.rare.data = select_data_for_reliability(allBetas{cohLevel.idx}, cohLevel.sensors, rareIdx, chanlabels);
[cohLevel.rare.maxCorr_within, cohLevel.rare.bestLag_within, cohLevel.rare.maxCorr_between, cohLevel.rare.bestLag_between, ...
    cohLevel.rare.R, cohLevel.rare.R_bounds] = crdm_eeg_reliability(cohLevel.rare.data, 150);

%% PE regressor
pe.idx = 3;
pe.sensors = {'CPZ', 'CP1', 'CP2'};
pe.freq.data = select_data_for_reliability(allBetas{pe.idx}, pe.sensors, freqIdx, chanlabels);
[pe.freq.maxCorr_within, pe.freq.bestLag_within, pe.freq.maxCorr_between, pe.freq.bestLag_between, ...
    pe.freq.R, pe.freq.R_bounds] = crdm_eeg_reliability(pe.freq.data, 150);
pe.rare.data = select_data_for_reliability(allBetas{pe.idx}, pe.sensors, rareIdx, chanlabels);
[pe.rare.maxCorr_within, pe.rare.bestLag_within, pe.rare.maxCorr_between, pe.rare.bestLag_between, ...
    pe.rare.R, pe.rare.R_bounds] = crdm_eeg_reliability(pe.rare.data, 150);

%% absoluted stimulus regressor
absCP.idx = 4;
absCP.sensors = {'CPZ', 'CP1', 'CP2'};
absCP.freq.data = select_data_for_reliability(allBetas{absCP.idx}, absCP.sensors, freqIdx, chanlabels);
[absCP.freq.maxCorr_within, absCP.freq.bestLag_within, absCP.freq.maxCorr_between, absCP.freq.bestLag_between, ...
    absCP.freq.R, absCP.freq.R_bounds] = crdm_eeg_reliability(absCP.freq.data, 150);
absCP.rare.data = select_data_for_reliability(allBetas{absCP.idx}, absCP.sensors, rareIdx, chanlabels);
[absCP.rare.maxCorr_within, absCP.rare.bestLag_within, absCP.rare.maxCorr_between, absCP.rare.bestLag_between, ...
    absCP.rare.R, absCP.rare.R_bounds] = crdm_eeg_reliability(absCP.rare.data, 150);

absPA.idx = 4;
absPA.sensors = {'Pz', 'P1', 'P2'};
absPA.freq.data = select_data_for_reliability(allBetas{absPA.idx}, absPA.sensors, freqIdx, chanlabels);
[absPA.freq.maxCorr_within, absPA.freq.bestLag_within, absPA.freq.maxCorr_between, absPA.freq.bestLag_between, ...
    absPA.freq.R, absPA.freq.R_bounds] = crdm_eeg_reliability(absPA.freq.data, 150);
absPA.rare.data = select_data_for_reliability(allBetas{absPA.idx}, absPA.sensors, rareIdx, chanlabels);
[absPA.rare.maxCorr_within, absPA.rare.bestLag_within, absPA.rare.maxCorr_between, absPA.rare.bestLag_between, ...
    absPA.rare.R, absPA.rare.R_bounds] = crdm_eeg_reliability(absPA.rare.data, 150);

absPO.idx = 4;
absPO.sensors = {'POz', 'PO3', 'PO4'};
absPO.freq.data = select_data_for_reliability(allBetas{absPO.idx}, absPO.sensors, freqIdx, chanlabels);
[absPO.freq.maxCorr_within, absPO.freq.bestLag_within, absPO.freq.maxCorr_between, absPO.freq.bestLag_between, ...
    absPO.freq.R, absPO.freq.R_bounds] = crdm_eeg_reliability(absPO.freq.data, 150);
absPO.rare.data = select_data_for_reliability(allBetas{absPO.idx}, absPO.sensors, rareIdx, chanlabels);
[absPO.rare.maxCorr_within, absPO.rare.bestLag_within, absPO.rare.maxCorr_between, absPO.rare.bestLag_between, ...
    absPO.rare.R, absPO.rare.R_bounds] = crdm_eeg_reliability(absPO.rare.data, 150);


% now averaging across all conditions
sessIdx = 1:4;

%% jump regressor
jump.idx = 1;
jump.sensors = {'CPZ', 'CP1', 'CP2'};

jump.sess.data = select_data_for_reliability(allBetas{jump.idx}, jump.sensors, sessIdx, chanlabels);
[jump.sess.maxCorr_within, jump.sess.bestLag_within, jump.sess.maxCorr_between, jump.sess.bestLag_between, ...
    jump.sess.R, jump.sess.R_bounds] = crdm_eeg_reliability(jump.sess.data, 150);

%% coherence level at jump regressor
cohLevel.idx = 2;
cohLevel.sensors = {'CPZ', 'CP1', 'CP2'};
cohLevel.sess.data = select_data_for_reliability(allBetas{cohLevel.idx}, cohLevel.sensors, sessIdx, chanlabels);
[cohLevel.sess.maxCorr_within, cohLevel.sess.bestLag_within, cohLevel.sess.maxCorr_between, cohLevel.sess.bestLag_between, ...
    cohLevel.sess.R, cohLevel.sess.R_bounds] = crdm_eeg_reliability(cohLevel.sess.data, 150);


%% PE regressor
pe.idx = 3;
pe.sensors = {'CPZ', 'CP1', 'CP2'};
pe.sess.data = select_data_for_reliability(allBetas{pe.idx}, pe.sensors, sessIdx, chanlabels);
[pe.sess.maxCorr_within, pe.sess.bestLag_within, pe.sess.maxCorr_between, pe.sess.bestLag_between, ...
    pe.sess.R, pe.sess.R_bounds] = crdm_eeg_reliability(pe.sess.data, 150);


%% absoluted stimulus regressor
absCP.idx = 4;
absCP.sensors = {'CPZ', 'CP1', 'CP2'};
absCP.sess.data = select_data_for_reliability(allBetas{absCP.idx}, absCP.sensors, sessIdx, chanlabels);
[absCP.sess.maxCorr_within, absCP.sess.bestLag_within, absCP.sess.maxCorr_between, absCP.sess.bestLag_between, ...
    absCP.sess.R, absCP.sess.R_bounds] = crdm_eeg_reliability(absCP.sess.data, 150);

absPA.idx = 4;
absPA.sensors = {'Pz', 'P1', 'P2'};
absPA.sess.data = select_data_for_reliability(allBetas{absPA.idx}, absPA.sensors, sessIdx, chanlabels);
[absPA.sess.maxCorr_within, absPA.sess.bestLag_within, absPA.sess.maxCorr_between, absPA.sess.bestLag_between, ...
    absPA.sess.R, absPA.sess.R_bounds] = crdm_eeg_reliability(absPA.sess.data, 150);

absPO.idx = 4;
absPO.sensors = {'POz', 'PO3', 'PO4'};
absPO.sess.data = select_data_for_reliability(allBetas{absPO.idx}, absPO.sensors, sessIdx, chanlabels);
[absPO.sess.maxCorr_within, absPO.sess.bestLag_within, absPO.sess.maxCorr_between, absPO.sess.bestLag_between, ...
    absPO.sess.R, absPO.sess.R_bounds] = crdm_eeg_reliability(absPO.sess.data, 150);

end
