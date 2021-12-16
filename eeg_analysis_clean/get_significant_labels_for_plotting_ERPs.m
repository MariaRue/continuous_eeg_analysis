function [SignificantTimePoints] = get_significant_labels_for_plotting_ERPs(stat, TimeBinsRegressor)


significantProbabilities = stat.prob < stat.cfg.alpha; 

statTimePoints = stat.time(significantProbabilities); 

SignificantTimePoints = find(ismember(TimeBinsRegressor, statTimePoints)); 

end 