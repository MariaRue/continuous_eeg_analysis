function [stat] = permutation_testGLM(data1, data2, jumpRegressor)
% this function calculates a permutation test statistic for EEG data based
% on fieldtrip functions
% http://www.fieldtriptoolbox.org/reference/ft_statistics_montecarlo/
% http://www.fieldtriptoolbox.org/tutorial/cluster_permutation_timelock/#within-subjects-experiments




% prepare neighbuours file - finds channel neighbours for spatial
% clustering or interpolation of bad channels - using the distance method,
% neighbours are baded on a minimum neighbourhood distance.



% 2) run neighbours function

% cfg = [];
% cfg.elec = ft_read_sens('standard_1020.elc');
% cfg.template = 'easycapM1_neighb.mat';
% cfg.method = 'distance';
% cfg.channel = data1{1}.label(1:61);
% %
% neighbours = ft_prepare_neighbours(cfg,data1{1});


%in this step we are using ft_timelockstatistics

% these are config fields that are required that are not specific for
% ft_timelockstatistics
cfg = [];
cfg.neighbours = []; %neighbours;

if jumpRegressor
    
    cfg.latency = [0 0.8];
    
else
    cfg.latency = [-2 0]; % 0 to 0.8s for jump event regressors, timelocked to button press -2 to 0s. time interval over which the experimental conditions will be compared
end
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05; % this is just for thresholding
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 0;
cfg.tail = 0;
cfg.clustertail = 0; % why?
cfg.correcttail = 'alpha'; % for a two-sided test this splits alpha to still have a false alarm rate of 0.05
cfg.alpha = 0.05; % % check? because we have correcttail? because it is a 2 side test - and we want to have an alpha of 0.05 for overall test statistic)
cfg.numrandomization = 10000; %10000 for paper, 500 to get reliable estimate
cfg.spmversion = 'spm12';


% build a design matrix


nS = length(data1);

design = zeros(2,nS * 2);
for i = 1:nS
    design(1,i) = i;
    
end

for i = 1:nS
    
    design(1,nS + i) = i;
end

design(2,1:nS) = 1;
design(2,nS+1 : 2*nS) = 2;


cfg.design = design;
cfg.uvar = 1; % subject id - first row in desing matrix
cfg.ivar = 2; % 2nd row is independent variable


[stat] = ft_timelockstatistics(cfg, data1{:}, data2{:});







end