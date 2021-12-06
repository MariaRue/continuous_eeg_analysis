function [EEGdir,EEGdirdata,scriptdir,nSess,nS] = setup_EEG_session(subj_list)

% [scriptdir,nSess,nS] = setup_EEG_session(subj_list)
%
% takes a list of subject identifier numbers and sets up paths/scriptdir
% for other functions

%% specify the datasets to read in by identifying the subid
% in terminal go to the data folder containing sufolders for each
% participant --> ls > txt (make sure that list contains only direct
% subject directories)
%scriptdir = fullfile('/Users/maria/Documents/Matlab/continuous_eeg_analysis/eeg_analysis');
%EEGdir= fullfile('/Users/maria/Documents/data/data.continuous_rdk','data','EEG');
EEGdir = '/Volumes/LaCie/data_preproc';
EEGdirdata = '/Volumes/LaCie/data/EEG/';
user = 'MR';
switch user
    case 'LH'
        [hd,sd] = get_homedir;
        addpath(fullfile(hd,'matlab','hidden_from_matlab','spm12'))
        addpath(fullfile(hd,'matlab','hidden_from_matlab','fieldtrip'))
        scriptdir = fullfile(hd,'projects','continuous_eeg_analysis','eeg_analysis');
        addpath(scriptdir)
    case 'MR'
        scriptdir = fullfile('/Users/maria/Documents/Matlab/continuous_eeg_analysis/eeg_analysis');
        addpath('/Users/maria/Documents/matlab/spm12');
        addpath('/Users/maria/Documents/MATLAB/fieldtrip'); % fieldtrip tool box to analyse data
        addpath('/Users/maria/Documents/MATLAB/continuous_eeg_analysis/eeg_analysis')
        addpath(genpath('/Users/maria/Documents/MATLAB/cbrewer'));
        addpath('/Users/maria/Documents/MATLAB/raacampbell-shadedErrorBar-9b19a7b')
    case 'MR_iMac'
        scriptdir = fullfile('/Users/maria/Documents/Matlab/continuous_eeg_analysis/eeg_analysis');
        addpath('/Users/maria/Documents/matlab/spm12');
        addpath('/Users/maria/Documents/MATLAB/fieldtrip'); % fieldtrip tool box to analyse data
        addpath('/Users/maria/Documents/MATLAB/continuous_eeg_analysis/eeg_analysis')
        addpath(genpath('/Users/maria/Documents/MATLAB/cbrewer'));
        addpath('/Users/maria/Documents/MATLAB/raacampbell-shadedErrorBar-9b19a7b')
        
        EEGdir = '/Users/maria/Documents/data/data_preproc/';
        EEGdirdata = '/Users/maria/Documents/data/data_preproc';
        
end


ft_defaults

nSess = 6; %number of sessions per subject
nS = length(subj_list); %number of subjects
end