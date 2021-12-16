function [selectedDataFreq, selectedDataRare, selectedDataShort, selectedDataLong,...
    selectedDataAvgFreq, selectedDataAvgRare, selectedDataAvgShort, selectedDataAvgLong, allDataAvg] ...
    = prepare_data_for_perm_test(betas, RegressorOptions, chanlabels, regressorIDx, electrodesForPermTest, nS)

% get data sorted by condition

freq_avg_sj = mean(betas{regressorIDx}(:,:,:,1:2),4);

% short averaged
short_avg_sj = mean(betas{regressorIDx}(:,:,:,[1,3]),4);

% rare averaged
rare_avg_sj = mean(betas{regressorIDx}(:,:,:,3:4),4);

% long averaged
long_avg_sj = mean(betas{regressorIDx}(:,:,:,[2,4]),4);

all_avg_sj = mean(betas{regressorIDx}(:,:,:,:),4);


% get elec info for fieldtrip structure
% get correct chanlabels
new_labels = change_electrode_labels(chanlabels);
% get electrode information
load('elec_field_for_GLM');

% create fieldtrip structure
timeBins = RegressorOptions(regressorIDx).timeBins/1000; % get timebins in seconds

for sj = 1:nS
    
    data_ft_freq{sj} = create_fieldTrip_structure(timeBins, new_labels, freq_avg_sj(sj,:,:), elecs);
    
    data_ft_rare{sj} = create_fieldTrip_structure(timeBins, new_labels, rare_avg_sj(sj,:,:), elecs);
    
    data_ft_short{sj} = create_fieldTrip_structure(timeBins, new_labels, short_avg_sj(sj,:,:), elecs);
    
    data_ft_long{sj} = create_fieldTrip_structure(timeBins, new_labels, long_avg_sj(sj,:,:), elecs);
    
    data_ft_all{sj} = create_fieldTrip_structure(timeBins, new_labels, all_avg_sj(sj,:,:), elecs);
    
    
    
    % select electrode specific data for perm-test
    
    cfg = [];
    cfg.channel = electrodesForPermTest;
    cfg.avgoverchan = 'yes';
    cfg.nanmean = 'yes';
    
    selectedDataFreq{sj} = select_electrode_specific_data_for_permtest(cfg, data_ft_freq{sj});
    
    selectedDataRare{sj} = select_electrode_specific_data_for_permtest(cfg, data_ft_rare{sj});
    
    selectedDataShort{sj} = select_electrode_specific_data_for_permtest(cfg, data_ft_short{sj});
    
    selectedDataLong{sj} = select_electrode_specific_data_for_permtest(cfg, data_ft_long{sj});
    
    
    
    
end



% calculate average across subjects per condition and se for plotting

cfg = [];

selectedDataAvgFreq = calculate_average_across_subjects(cfg, selectedDataFreq);

selectedDataAvgRare = calculate_average_across_subjects(cfg, selectedDataRare);

selectedDataAvgShort = calculate_average_across_subjects(cfg, selectedDataShort);

selectedDataAvgLong = calculate_average_across_subjects(cfg, selectedDataLong);

allDataAvg = calculate_average_across_subjects(cfg, data_ft_all);

end