function [dataDetrendDemean] = detrend_demean_trials_for_tf_analysis(dataAppend)


cfg = []; 
cfg.demean = 'yes'; 
cfg.detrend = 'yes'; 

dataDetrendDemean = ft_preprocessing(cfg,dataAppend);


end 