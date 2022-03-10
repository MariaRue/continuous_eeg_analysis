function [selectedData, allDataAvg] ...
    = prepare_data_for_perm_test(betas, RegressorOptions, chanlabels, regressorIDx, electrodesForPermTest, nS)

% data as tested during experiment - this is used for testing for the
% interaction between frequency and length 
try
freqShortSj = squeeze(betas{regressorIDx}(:,:,:,1));
catch 
    keyboard; 
    
end 
freqLongSj = squeeze(betas{regressorIDx}(:,:,:,2));

rareShortSj = squeeze(betas{regressorIDx}(:,:,:,3));

rareLongSj = squeeze(betas{regressorIDx}(:,:,:,4));


% get data sorted by condition

freqAvgSj = mean(betas{regressorIDx}(:,:,:,1:2),4);

% short averaged
shortAvgSj = mean(betas{regressorIDx}(:,:,:,[1,3]),4);

% rare averaged
rareAvgSj = mean(betas{regressorIDx}(:,:,:,3:4),4);

% long averaged
longAvgSj = mean(betas{regressorIDx}(:,:,:,[2,4]),4);



all_avg_sj = mean(betas{regressorIDx}(:,:,:,:),4);


% get elec info for fieldtrip structure
% get correct chanlabels
new_labels = change_electrode_labels(chanlabels);
% get electrode information
load('elec_field_for_GLM');

% create fieldtrip structure
timeBins = RegressorOptions(regressorIDx).timeBins/1000; % get timebins in seconds

for sj = 1:nS
    
    dataFtFreq{sj} = create_fieldTrip_structure(timeBins, new_labels, freqAvgSj(sj,:,:), elecs);
    
    dataFtRare{sj} = create_fieldTrip_structure(timeBins, new_labels, rareAvgSj(sj,:,:), elecs);
    
    dataFtShort{sj} = create_fieldTrip_structure(timeBins, new_labels, shortAvgSj(sj,:,:), elecs);
    
    dataFtLong{sj} = create_fieldTrip_structure(timeBins, new_labels, longAvgSj(sj,:,:), elecs);
    
    dataFtAll{sj} = create_fieldTrip_structure(timeBins, new_labels, all_avg_sj(sj,:,:), elecs);
    
    
    % data as tested in experiment 
    
    dataFtFreqShort{sj} = create_fieldTrip_structure(timeBins, new_labels, freqShortSj(sj,:,:), elecs);
    dataFtFreqLong{sj} = create_fieldTrip_structure(timeBins, new_labels, freqLongSj(sj,:,:), elecs);
    dataFtRareShort{sj} = create_fieldTrip_structure(timeBins, new_labels, rareShortSj(sj,:,:), elecs);
    dataFtRareLong{sj} = create_fieldTrip_structure(timeBins, new_labels, rareLongSj(sj,:,:), elecs);
    
    
    
    % select electrode specific data for perm-test
    
    cfg = [];
    cfg.channel = electrodesForPermTest;
    cfg.avgoverchan = 'yes';
    cfg.nanmean = 'yes';
    
    selectedDataFreq{sj} = select_electrode_specific_data_for_permtest(cfg, dataFtFreq{sj});
    
    selectedDataRare{sj} = select_electrode_specific_data_for_permtest(cfg, dataFtRare{sj});
    
    selectedDataShort{sj} = select_electrode_specific_data_for_permtest(cfg, dataFtShort{sj});
    
    selectedDataLong{sj} = select_electrode_specific_data_for_permtest(cfg, dataFtLong{sj});
    
    selectedDataAll{sj} = select_electrode_specific_data_for_permtest(cfg, dataFtAll{sj});
    
    
    % data as tested in experiment 
    selectedDataFreqShort{sj} = select_electrode_specific_data_for_permtest(cfg, dataFtFreqShort{sj});
    selectedDataFreqLong{sj}  = select_electrode_specific_data_for_permtest(cfg, dataFtFreqLong{sj});
    selectedDataRareShort{sj} = select_electrode_specific_data_for_permtest(cfg, dataFtRareShort{sj});
    selectedDataRareLong{sj}  = select_electrode_specific_data_for_permtest(cfg, dataFtRareLong{sj});
    
    
end


% calculate average across subjects per condition and se for plotting

cfg = [];

selectedDataAvgFreq = calculate_average_across_subjects(cfg, selectedDataFreq);

selectedDataAvgRare = calculate_average_across_subjects(cfg, selectedDataRare);

selectedDataAvgShort = calculate_average_across_subjects(cfg, selectedDataShort);

selectedDataAvgLong = calculate_average_across_subjects(cfg, selectedDataLong);

selectedDataAvgAll = calculate_average_across_subjects(cfg, selectedDataAll);

allDataAvg = calculate_average_across_subjects(cfg, dataFtAll);

  % data as tested in experiment 
  selectedDataAvgFreqShort = calculate_average_across_subjects(cfg, selectedDataFreqShort);
  selectedDataAvgFreqLong = calculate_average_across_subjects(cfg, selectedDataFreqLong);
  selectedDataAvgRareShort = calculate_average_across_subjects(cfg, selectedDataRareShort);
  selectedDataAvgRareLong = calculate_average_across_subjects(cfg, selectedDataRareLong);
  
  
  % output 
  selectedData.subjectLevel.FreqShort = selectedDataFreqShort; 
  selectedData.subjectLevel.FreqLong = selectedDataFreqLong; 
  selectedData.subjectLevel.RareShort = selectedDataRareShort; 
  selectedData.subjectLevel.RareLong = selectedDataRareLong;
  selectedData.subjectLevel.All = selectedDataAll; 
  
  selectedData.subjectLevel.Freq = selectedDataFreq; 
  selectedData.subjectLevel.Rare = selectedDataRare; 
  selectedData.subjectLevel.Short = selectedDataShort; 
  selectedData.subjectLevel.Long = selectedDataLong;
  
  selectedData.Average.FreqShort = selectedDataAvgFreqShort; 
  selectedData.Average.FreqLong = selectedDataAvgFreqLong; 
  selectedData.Average.RareShort = selectedDataAvgRareShort; 
  selectedData.Average.RareLong = selectedDataAvgRareLong;
  selectedData.Average.All = selectedDataAvgAll; 
  
  selectedData.Average.Freq = selectedDataAvgFreq; 
  selectedData.Average.Rare = selectedDataAvgRare; 
  selectedData.Average.Short = selectedDataAvgShort; 
  selectedData.Average.Long = selectedDataAvgLong;
  
  

end