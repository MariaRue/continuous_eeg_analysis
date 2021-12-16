function [dat1_vs_dat2_avg] = calculate_difference_waveform(data_cp_all1, data_cp_all2)

for subject = 1:length(data_cp_all1)
    
    cfg = [];
    cfg.operation = 'subtract';
    cfg.parameter = 'avg';
    
    dat1_vs_dat2{subject} = ft_math(cfg,data_cp_all1{subject},data_cp_all2{subject});
    
end

cfg = [];
dat1_vs_dat2_avg = ft_timelockgrandaverage(cfg, dat1_vs_dat2{:} );
% calculate se
[se] = calculate_se(dat1_vs_dat2_avg.var,length(data_cp_all1));
dat1_vs_dat2_avg.se = se;




end