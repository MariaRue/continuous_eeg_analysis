function [options] = get_options_behaviour_rdk()

addpath(genpath('/Users/maria/Documents/MATLAB/cbrewer')); % sets path for cbrewer (used to define colours for plotting)
addpath('/Users/maria/Documents/MATLAB/raacampbell-shadedErrorBar-9b19a7b') % sets path to code to plot graphs with error bars
options.scriptdir = '/Users/maria/Documents/MATLAB/continuous_eeg_analysis/eeg_analysis_clean'; % path to analysis code
addpath(genpath(options.scriptdir));
options.path.behaviour = fullfile('/Volumes/crdkData','preprocessedData','behaviour'); % path to subject behaviour


options.SampleFrequency = 100; % frequency with which the stimulus is displayed

options.conditionLabels = {'frequent & short', 'frequent & long','rare & short', 'rare & long'}; % labels for plotting

options.trialCoherence = [0.3 0.4 0.5]; % average coherence of signal periods
options.lags = 500; % lags to produce behavioural kernels (in units of samples). 


% load and save color scheme for conditions 
cl = cbrewer('div','RdBu', 12);
options.colours4Conditions = cl([2 4 12 10],:);

options.totalNumberofSubjects = 28; 

end