function options = crdm_eeg_set_analysis_options( compName )
%CRDM_EEG_SET_ANALYSIS_OPTIONS Sets all options for the analysis of EEG
%data from the continuous dot motion task (Maria's task version).
%   IN:     compName        - string indicating on which machine this is
%                           running ('iMac', 'Pro', 'LilianLinux')
%   OUT:    options         - the struct that holds all analysis options

options.user = compName;

%%%---- ENTER YOUR PATHS HERE -----------------------------------------%%%
switch compName
    case 'LilianLinux'
        % This is where we are now (where the code is to be found):
        options.codeDir = fileparts(mfilename('fullpath'));
        % This is the base root for both raw data and analysis:
        options.mainDir = '/media/lil/Elements/CRDM/Maria/EEG';
        % This is the directory where you want the data analysis to happen:
        options.workDir = fullfile(options.maindir, 'prj', 'convGLM');
        % This is the directory where the raw EEG and behavioural data are:
        options.rawDir  = fullfile(options.maindir, 'raw');
        % This is where we are now (where the code is to be found):
        options.codeDir = fileparts(mfilename('fullpath'));
    case 'iMac'
        %options.path.EEG.analysis = '/Users/maria/Documents/data/data_preproc';
        options.path.EEG.analysis =  '/Volumes/crdkData';
        options.path.EEG.subjects = fullfile(options.path.EEG.analysis,'rawData','experiment');
        options.path.preproc.behaviour = fullfile('/Volumes/crdkData','preprocessedData','behaviour');
        options.path.EEG.raw = '/Volumes/crdkData';
        options.scriptdir = '/Users/maria/Documents/MATLAB/continuous_eeg_analysis/eeg_analysis_clean';
        options.preproc.path.layout = '/Users/maria/Documents/MATLAB/fieldtrip/template/layout/easycapM10.mat';
    case 'Pro'
        %options.path.EEG.analysis =  '/Volumes/LaCie/daten_fuer_laptop';
        options.path.EEG.analysis =  '/Volumes/crdkData';
        options.path.EEG.subjects = '/Volumes/LaCie/daten_fuer_laptop';
        options.scriptdir = '/Users/maria/Documents/MATLAB/continuous_eeg_analysis/eeg_analysis_clean';
        options.preproc.path.layout = '/Users/maria/Documents/MATLAB/fieldtrip/template/layout/easycapM10.mat';
end
%%%--------------------------------------------------------------------%%%

% set FS for GLM based on preprocessed data!!!!!!!!
Fs = 100;

options.conditionLabels = {'Trials frequent & short', 'Trials frequent & long', ...
                           'Trials rare & short', 'Trials rare & long' };

%-- preprocessing --------------------------------------------------------%
options.preproc.fsample = 100;
options.preproc.bandpass = [0.1 30];
options.preproc.eyeblink.excwin = 500;
options.preproc.artefact.threshold = 100;
options.preproc.artefact.excwin = 500;
options.preproc.artefact.badchanthresh = 1000;

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



end
