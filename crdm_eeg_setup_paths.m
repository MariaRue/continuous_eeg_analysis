function crdm_eeg_setup_paths( options )
%CRDM_EEG_SETUP_PATHS Sets the Matlab path for the CRDM EEG analysis
%   Run this function from the working directory of the analysis (the
%   project folder). Paths will be added automatically based on the 
%   location of this file.
%   IN:     -
%   OUT:    -

% remove all other toolboxes 
restoredefaultpath; 

% add project path with all sub-paths 
pathProject = fileparts(mfilename('fullpath')); 
addpath(genpath(pathProject));

% add toolboxes, paths depending on user
switch options.user
    case 'LilianLinux'
        pathToolboxes = fullfile(pathProject, '..', 'toolboxesForAnalysis');
        % NOTE: NEVER add SPM with subfolders to your path, since it creates 
        % conflicts with Matlab core functions, e.g., uint16 
        addpath(fullfile(pathToolboxes, 'spm12'));
        addpath(fullfile(pathToolboxes, 'eeglab14_1_2b'));
        addpath(fullfile(pathToolboxes, 'fieldtrip'));
        addpath(fullfile(pathToolboxes, 'cbrewer'));
        addpath(fullfile(pathToolboxes, 'shadedErrorBar'));
        
    case 'LTH_iMac'
        pathToolboxes = fullfile(pathProject, '..', 'toolboxesForAnalysis');
        % NOTE: NEVER add SPM with subfolders to your path, since it creates
        % conflicts with Matlab core functions, e.g., uint16
        addpath(fullfile(pathToolboxes, 'spm12'));
        addpath(fullfile(pathToolboxes, 'eeglab14_1_2b'));
        addpath(fullfile(pathToolboxes, 'fieldtrip'));
        addpath(fullfile(pathToolboxes, 'cbrewer'));
        addpath(fullfile(pathToolboxes, 'shadedErrorBar'));
        
    case {'iMac', 'Pro'}
        % set paths for fieldtrip, spm and cbrewer
        addpath('/Users/maria/Documents/matlab/continuous_eeg_analysis/spm12');
        addpath('/Users/maria/Documents/MATLAB/eeglab14_1_2b/');
        addpath('/Users/maria/Documents/MATLAB/fieldtrip'); % fieldtrip tool box to analyse data
        addpath(genpath('/Users/maria/Documents/MATLAB/cbrewer'));
        addpath('/Users/maria/Documents/MATLAB/raacampbell-shadedErrorBar-9b19a7b')
        
end

% initialise toolboxes
ft_defaults;
spm('defaults','eeg')

end

