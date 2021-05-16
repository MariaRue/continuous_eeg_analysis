function dataAverage = calculate_timelockGrandAverage_and_se(channels,data,numberSubjects)

cfg = [];
cfg.channels = channels;
dataAverage = ft_timelockgrandaverage(cfg,data{:});
[se] = calculate_se(dataAverage.var,numberSubjects);
dataAverage.se = se;

end
