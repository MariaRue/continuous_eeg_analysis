function options = continuous_RDK_set_options(mac)
% set FS for GLM based on preprocessed data!!!!!!!!
Fs = 100; 

options.conditionLabels = {'Trials frequent & short', 'Trials frequent & long', 'Trials rare & short', 'Trials rare & long' };
% set paths for fieldtrip, spm and cbrewer
addpath('/Users/maria/Documents/matlab/continuous_eeg_analysis/spm12');
addpath('/Users/maria/Documents/MATLAB/eeglab14_1_2b/');
addpath('/Users/maria/Documents/MATLAB/fieldtrip'); % fieldtrip tool box to analyse data
 addpath(genpath('/Users/maria/Documents/MATLAB/cbrewer'));
         addpath('/Users/maria/Documents/MATLAB/raacampbell-shadedErrorBar-9b19a7b')

ft_defaults;
%spm('defaults','eeg')
% set paths for stimulus matched EEG data - we still need to define later
% whether we want csd transformed data or not, in addition, these directory
% structures are not very tidy and I need to delete stuff - I have several
% copies of everything and the naming of directories between average and LM
%  RM references is not very consistent 
switch mac
    case 'iMac'
%options.path.EEG.analysis = '/Users/maria/Documents/data/data_preproc';
options.path.EEG.analysis =  '/Volumes/crdkData';
options.path.EEG.subjects = fullfile(options.path.EEG.analysis,'rawData','experiment');
options.path.preproc.behaviour = fullfile('/Volumes/crdkData','preprocessedData','behaviour');
options.path.EEG.raw = '/Volumes/crdkData';
    case 'Pro'
        
        
%options.path.EEG.analysis =  '/Volumes/LaCie/daten_fuer_laptop';
options.path.EEG.analysis =  '/Volumes/crdkData';
options.path.EEG.subjects = '/Volumes/LaCie/daten_fuer_laptop';


end 
options.scriptdir = '/Users/maria/Documents/MATLAB/continuous_eeg_analysis/eeg_analysis_clean';
addpath(genpath(options.scriptdir));
% this has to be completed with folder for subject and if averaged 

options.preproc.fsample = 100; 
options.preproc.bandpass = [0.1 30];
options.preproc.eyeblink.excwin = 500; 
options.preproc.artefact.threshold = 100; 
options.preproc.artefact.excwin = 500; 
options.preproc.artefact.badchanthresh = 1000;

%-- set options for subject level GLM-------------------------------------%
% comment
options.subjectLevelGLM.allRegressors(1).name = 'coherence_jump';
options.subjectLevelGLM.allRegressors(1).nLagsBack = 100;
options.subjectLevelGLM.allRegressors(1).nLagsForward = 150;
options.subjectLevelGLM.allRegressors(1).timeBins = (1000/Fs)*[-options.subjectLevelGLM.allRegressors(1).nLagsBack : options.subjectLevelGLM.allRegressors(1).nLagsForward];

options.subjectLevelGLM.allRegressors(2).name = 'coherence_jump_level';
options.subjectLevelGLM.allRegressors(2).nLagsBack = 100;
options.subjectLevelGLM.allRegressors(2).nLagsForward = 150;
options.subjectLevelGLM.allRegressors(2).timeBins = (1000/Fs)*[-options.subjectLevelGLM.allRegressors(2).nLagsBack : options.subjectLevelGLM.allRegressors(2).nLagsForward];

options.subjectLevelGLM.allRegressors(3).name = 'prediction_error';
options.subjectLevelGLM.allRegressors(3).nLagsBack = 150;
options.subjectLevelGLM.allRegressors(3).nLagsForward = 150;
options.subjectLevelGLM.allRegressors(3).timeBins = (1000/Fs)*[-options.subjectLevelGLM.allRegressors(3).nLagsBack : options.subjectLevelGLM.allRegressors(3).nLagsForward];

options.subjectLevelGLM.allRegressors(4).name = 'absoluted_stimulus';
options.subjectLevelGLM.allRegressors(4).nLagsBack = 150;
options.subjectLevelGLM.allRegressors(4).nLagsForward = 150;
options.subjectLevelGLM.allRegressors(4).timeBins = (1000/Fs)*[-options.subjectLevelGLM.allRegressors(4).nLagsBack : options.subjectLevelGLM.allRegressors(4).nLagsForward];

options.subjectLevelGLM.allRegressors(5).name = 'trial_start';
options.subjectLevelGLM.allRegressors(5).nLagsBack = 50;
options.subjectLevelGLM.allRegressors(5).nLagsForward = 800;
options.subjectLevelGLM.allRegressors(5).timeBins = (1000/Fs)*[-options.subjectLevelGLM.allRegressors(5).nLagsBack : options.subjectLevelGLM.allRegressors(5).nLagsForward];

options.subjectLevelGLM.allRegressors(6).name = 'correct_trial_response';
options.subjectLevelGLM.allRegressors(6).nLagsBack = 500;
options.subjectLevelGLM.allRegressors(6).nLagsForward = 350;
options.subjectLevelGLM.allRegressors(6).timeBins = (1000/Fs)*[-options.subjectLevelGLM.allRegressors(6).nLagsBack : options.subjectLevelGLM.allRegressors(6).nLagsForward];

options.subjectLevelGLM.allRegressors(7).name = 'false_alarm';
options.subjectLevelGLM.allRegressors(7).nLagsBack = 500;
options.subjectLevelGLM.allRegressors(7).nLagsForward = 350;
options.subjectLevelGLM.allRegressors(7).timeBins = (1000/Fs)*[-options.subjectLevelGLM.allRegressors(7).nLagsBack : options.subjectLevelGLM.allRegressors(7).nLagsForward];

options.subjectLevelGLM.allRegressors(8).name = 'difference_waveform_correct_trial_response';
options.subjectLevelGLM.allRegressors(8).nLagsBack = 500;
options.subjectLevelGLM.allRegressors(8).nLagsForward = 350;
options.subjectLevelGLM.allRegressors(8).timeBins = (1000/Fs)*[-options.subjectLevelGLM.allRegressors(8).nLagsBack : options.subjectLevelGLM.allRegressors(8).nLagsForward];

options.subjectLevelGLM.allRegressors(9).name = 'difference_waveform_false_alarm';
options.subjectLevelGLM.allRegressors(9).nLagsBack = 500;
options.subjectLevelGLM.allRegressors(9).nLagsForward = 350;
options.subjectLevelGLM.allRegressors(9).timeBins = (1000/Fs)*[-options.subjectLevelGLM.allRegressors(9).nLagsBack : options.subjectLevelGLM.allRegressors(9).nLagsForward];


options.subjectLevelGLM.allRegressors(10).name = 'coherence_responses';
options.subjectLevelGLM.allRegressors(10).nLagsBack = 500;
options.subjectLevelGLM.allRegressors(10).nLagsForward = 350;
options.subjectLevelGLM.allRegressors(10).timeBins = (1000/Fs)*[-options.subjectLevelGLM.allRegressors(10).nLagsBack : options.subjectLevelGLM.allRegressors(10).nLagsForward];


options.subjectLevelGLM.allRegressors(11).name = 'coherence_jump_vertical';
options.subjectLevelGLM.allRegressors(11).nLagsBack = 100;
options.subjectLevelGLM.allRegressors(11).nLagsForward = 150;
options.subjectLevelGLM.allRegressors(11).timeBins = (1000/Fs)*[-options.subjectLevelGLM.allRegressors(11).nLagsBack : options.subjectLevelGLM.allRegressors(11).nLagsForward];

options.subjectLevelGLM.allRegressors(12).name = 'coherence_jump_level_vertical';
options.subjectLevelGLM.allRegressors(12).nLagsBack = 100;
options.subjectLevelGLM.allRegressors(12).nLagsForward = 150;
options.subjectLevelGLM.allRegressors(12).timeBins = (1000/Fs)*[-options.subjectLevelGLM.allRegressors(12).nLagsBack : options.subjectLevelGLM.allRegressors(12).nLagsForward];

options.subjectLevelGLM.allRegressors(13).name = 'prediction_error_vertical';
options.subjectLevelGLM.allRegressors(13).nLagsBack = 100;
options.subjectLevelGLM.allRegressors(13).nLagsForward = 150;
options.subjectLevelGLM.allRegressors(13).timeBins = (1000/Fs)*[-options.subjectLevelGLM.allRegressors(13).nLagsBack : options.subjectLevelGLM.allRegressors(13).nLagsForward];

options.subjectLevelGLM.allRegressors(14).name = 'absoluted_stimulus_vertical';
options.subjectLevelGLM.allRegressors(14).nLagsBack = 150;
options.subjectLevelGLM.allRegressors(14).nLagsForward = 150;
options.subjectLevelGLM.allRegressors(14).timeBins = (1000/Fs)*[-options.subjectLevelGLM.allRegressors(14).nLagsBack : options.subjectLevelGLM.allRegressors(14).nLagsForward];



%-- set GLM types---------------------------------------------------------%

% these options specify different GLM models for which we select different
% regressors from the list in options.subjectLevelGLM.allRegressors
options.subjectLevelGLM.absolute.regressors = options.subjectLevelGLM.allRegressors(4:7);
options.subjectLevelGLM.absolute.name = 'absolute';

options.subjectLevelGLM.jumps_absolute.regressors = options.subjectLevelGLM.allRegressors(1:9);
options.subjectLevelGLM.jumps_absolute.name = 'jumps_absolute';

options.subjectLevelGLM.response.regressors = options.subjectLevelGLM.allRegressors(5:9);
options.subjectLevelGLM.response.name = 'response';

options.subjectLevelGLM.coherence_responses.regressors = options.subjectLevelGLM.allRegressors(5:10);
options.subjectLevelGLM.coherence_responses.name = 'coherence_responses';

options.subjectLevelGLM.vertical_jumps_absolute.regressors = options.subjectLevelGLM.allRegressors([1:7,11:14]);
options.subjectLevelGLM.coherence_responses.name = 'vertical_jumps_absolute';

%-- set conventional timelocked analysis options-------------------------%
options.singleTrial.buttonPress.preStim = 7; 
options.singleTrial.buttonPress.postStim = 4; 
options.singleTrial.buttonPress.name = 'buttonPress'; 

options.singleTrial.trialStart.preStim = 2; 
options.singleTrial.trialStart.postStim = 8; 
options.singleTrial.trialStart.name = 'trialStart'; 


options.preproc.eyeblinktreatment   = 'ssp'; % 'reject', 'ssp'
options.preproc.eyeblinkchannels    = {'VEOG'};
options.preproc.eyeblinkthreshold   = 5; % for SD thresholding: in standard deviations, for amp in uV
options.preproc.windowForEyeblinkdetection = [10 900]; % first event of interest (and optionally last)

options.preproc.eyeblinkmode        = 'eventbased'; % uses EEG triggers for trial onsets
options.preproc.eyeblinkwindow      = 0.5; % in s around blink events
options.preproc.eyeblinktrialoffset = 0.1; % in s: EBs won't hurt <100ms after tone onset
options.preproc.eyeblinkEOGchannel  = 'VEOG'; % EOG channel (name/idx) to plot
options.preproc.eyebadchanthresh    = 0.4; % prop of bad trials due to EBs

options.preproc.eyeconfoundcomps    = 1;
options.preproc.eyecorrectionchans  = {'Fp1', 'Fz', 'AF8', 'T7', 'Oz'};
options.preproc.preclean.doFilter           = false;
options.preproc.preclean.lowPassFilterFreq  = 10;
options.preproc.preclean.doBadChannels      = true;
options.preproc.preclean.doBadTrials        = false;
options.preproc.preclean.doRejection        = false;
options.preproc.preclean.badtrialthresh     = 500;
options.preproc.preclean.badchanthresh      = 0.5;
options.preproc.preclean.rejectPrefix       = 'cleaned_';
options.preproc.keep = 1;

options.preproc.path.layout = '/Users/maria/Documents/MATLAB/fieldtrip/template/layout/easycapM10.mat';
end
