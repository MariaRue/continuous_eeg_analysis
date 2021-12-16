function selectedDataAvg = calculate_average_across_subjects(cfg, selectedData)

selectedDataAvg = ft_timelockgrandaverage(cfg, selectedData{:}); 

[se] = calculate_se_convGLM(selectedData);
selectedDataAvg.se = se; 



end 