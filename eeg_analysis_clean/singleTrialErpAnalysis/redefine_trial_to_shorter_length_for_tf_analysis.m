function dataRedefined = redefine_trial_to_shorter_length_for_tf_analysis(dataDetrendDemean, options)

cfg = []; 
cfg.toilim = [-4 1-(1/options.preproc.fsample)];

dataRedefined = ft_redefinetrial(cfg, dataDetrendDemean);



end 