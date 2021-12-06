function tfdataAppend = compute_tf_on_single_trial_data(dataAppend,options)


cfg         = [];
cfg.output  = 'pow';
cfg.channel = 'EEG';
cfg.method  = 'mtmconvol';
cfg.taper   = 'hanning';
cfg.foi     = 1:40;
cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
cfg.toi          = -4:0.1:1;                  
cfg.keeptrials = 'yes';

tfdataAppend = ft_freqanalysis(cfg, dataAppend);

end 