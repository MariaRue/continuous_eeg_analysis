timelockCondition = 'buttonPress';
reference = 'LMRM';
csdFlag = 0;
rtFlag = 0;

baseline = 0;
faFlag = 1;
subjectList = [16, 18:21, 24, 26, 27, 28, 29, 31, 33, 34, 35,  42, 43, 47, 50, 51, 52, 54, 55, 57, 58];
%subjectList = [58];

numberSubjects = length(subjectList);
options = continuous_RDK_set_options('iMac');

for subjects = 1:length(subjectList)
    
    subID = subjectList(subjects);
    [details,paths] =  conrdk_subjects(subID,options,reference,csdFlag);
    
    data_load = load(paths.(reference).singleTrial.appendedData.(timelockCondition));
    data{subjects} = data_load.dataAppend;
    
    
    % change labels to easycap
    [easy_cap_labels] = change_electrode_labels(data{subjects}.label);
    data{subjects}.label = easy_cap_labels;
    
    % average for each coherence level for each subjects
    data_avg =  compute_average_single_subject_level(data{subjects},'coherence',rtFlag,[4 5],baseline);
    coh_30{subjects} = data_avg{1};
    coh_40{subjects} = data_avg{2};
    coh_50{subjects} = data_avg{3};
    
    
    
    %     average for each condition level for each subjects
    data_avg_con =  compute_average_single_subject_level(data{subjects},'condition',rtFlag,[4 5],baseline);
    cond_1{subjects} = data_avg_con{1};
    cond_2{subjects} = data_avg_con{2};
    cond_3{subjects} = data_avg_con{3};
    cond_4{subjects} = data_avg_con{4};
    %
    
    % prepare data for permutation test for freq vs rare conditions
    data_avg_con_frequency =  compute_average_single_subject_level(data{subjects},'perm_test_rare_freq_trial',rtFlag,[-7 -6],baseline);
    cond_tr_freq{subjects} = data_avg_con_frequency{1};
    cond_tr_rare{subjects} = data_avg_con_frequency{2};
    
    
    
    
    % prepare data for permutation test for long vs short conditions
    data_avg_con_length =  compute_average_single_subject_level(data{subjects},'perm_test_short_long_trial',rtFlag,[-2 -1],baseline);
    cond_tr_short{subjects} = data_avg_con_length{1};
    cond_tr_long{subjects} = data_avg_con_length{2};
    
    % conditions rare vs freq
    
    channels = {'CP1', 'CPz', 'CP2'};
    cond_tr_rare_cp{subjects} = select_single_subject_channel_data(channels,cond_tr_rare{subjects});
    [cond_tr_freq_cp{subjects}] = select_single_subject_channel_data(channels,cond_tr_freq{subjects});
    
    [cond_tr_long_cp{subjects}] = select_single_subject_channel_data(channels,cond_tr_long{subjects});
    [cond_tr_short_cp{subjects}] = select_single_subject_channel_data(channels,cond_tr_short{subjects});
    
    cfg = []; 
    avg_all_dat = ft_timelockgrandaverage(cfg,data_avg_con{:});
    cond_tr_cp_all{subjects} = select_single_subject_channel_data(channels, avg_all_dat);
    
%     [coh_03_cp{subjects}] = select_single_subject_channel_data(channels,coh_30{subjects});
%     [coh_04_cp{subjects}] = select_single_subject_channel_data(channels,coh_40{subjects});
%     [coh_05_cp{subjects}] = select_single_subject_channel_data(channels,coh_50{subjects});
% %     
    channels = {'C1', 'Cz', 'C2'};
    cond_tr_rare_c{subjects} = select_single_subject_channel_data(channels,cond_tr_rare{subjects});
    [cond_tr_freq_c{subjects}] = select_single_subject_channel_data(channels,cond_tr_freq{subjects});
    
    [cond_tr_long_c{subjects}] = select_single_subject_channel_data(channels,cond_tr_long{subjects});
    [cond_tr_short_c{subjects}] = select_single_subject_channel_data(channels,cond_tr_short{subjects});
    
%     [coh_03_c{subjects}] = select_single_subject_channel_data(channels,coh_30{subjects});
%     [coh_04_c{subjects}] = select_single_subject_channel_data(channels,coh_40{subjects});
%     [coh_05_c{subjects}] = select_single_subject_channel_data(channels,coh_50{subjects});
%     
    
%     
%     channels = {'FC1', 'FCz', 'FC2'};
%     cond_tr_rare_fc{subjects} = select_single_subject_channel_data(channels,cond_tr_rare{subjects});
%     [cond_tr_freq_fc{subjects}] = select_single_subject_channel_data(channels,cond_tr_freq{subjects});
%     
%     [cond_tr_long_fc{subjects}] = select_single_subject_channel_data(channels,cond_tr_long{subjects});
%     [cond_tr_short_fc{subjects}] = select_single_subject_channel_data(channels,cond_tr_short{subjects});
%     
%     [coh_03_fc{subjects}] = select_single_subject_channel_data(channels,coh_30{subjects});
%     [coh_04_fc{subjects}] = select_single_subject_channel_data(channels,coh_40{subjects});
%     [coh_05_fc{subjects}] = select_single_subject_channel_data(channels,coh_50{subjects});
    
    channels = {'P1', 'Pz', 'P2'};
    cond_tr_rare_p{subjects} = select_single_subject_channel_data(channels,cond_tr_rare{subjects});
    [cond_tr_freq_p{subjects}] = select_single_subject_channel_data(channels,cond_tr_freq{subjects});
    
    [cond_tr_long_p{subjects}] = select_single_subject_channel_data(channels,cond_tr_long{subjects});
    [cond_tr_short_p{subjects}] = select_single_subject_channel_data(channels,cond_tr_short{subjects});
    
    [coh_03_p{subjects}] = select_single_subject_channel_data(channels,coh_30{subjects});
    [coh_04_p{subjects}] = select_single_subject_channel_data(channels,coh_40{subjects});
    [coh_05_p{subjects}] = select_single_subject_channel_data(channels,coh_50{subjects});
    
    
    channels = {'P3', 'P4', 'P5', 'P6'};
    cond_tr_rare_plat{subjects} = select_single_subject_channel_data(channels,cond_tr_rare{subjects});
    [cond_tr_freq_plat{subjects}] = select_single_subject_channel_data(channels,cond_tr_freq{subjects});
    
    [cond_tr_long_plat{subjects}] = select_single_subject_channel_data(channels,cond_tr_long{subjects});
    [cond_tr_short_plat{subjects}] = select_single_subject_channel_data(channels,cond_tr_short{subjects});
    
%     [coh_03_plat{subjects}] = select_single_subject_channel_data(channels,coh_30{subjects});
%     [coh_04_plat{subjects}] = select_single_subject_channel_data(channels,coh_40{subjects});
%     [coh_05_plat{subjects}] = select_single_subject_channel_data(channels,coh_50{subjects});
    
    %%%%%%%% button press false alarms
    
    switch timelockCondition
        case 'buttonPress'
            
            if faFlag
                % average for each condition level for each subjects for fa
                data_avg_con =  compute_average_single_subject_level(data{subjects},'false alarm',rtFlag,[4 5],baseline);
                fa_cond_1{subjects} = data_avg_con{1};
                fa_cond_2{subjects} = data_avg_con{2};
                fa_cond_3{subjects} = data_avg_con{3};
                fa_cond_4{subjects} = data_avg_con{4};
                
                % prepare data for permutation test for freq vs rare conditions false alarms
                data_avg_con =  compute_average_single_subject_level(data{subjects},'perm_test_rare_freq_fa',rtFlag,[-2 -1],baseline);
                fa_cond_tr_freq{subjects} = data_avg_con{1};
                fa_cond_tr_rare{subjects} = data_avg_con{2};
                
                
                
                % prepare data for permutation test for long vs short conditions
                % false alarms
                data_avg_con =  compute_average_single_subject_level(data{subjects},'perm_test_short_long_fa',rtFlag,[-2 -1],baseline);
                fa_cond_tr_short{subjects} = data_avg_con{1};
                fa_cond_tr_long{subjects} = data_avg_con{2};
                
                
                
                channels = {'CP1', 'CPz', 'CP2'};
                [fa_cond_tr_rare_cp{subjects}] = select_single_subject_channel_data(channels,fa_cond_tr_rare{subjects});
                [fa_cond_tr_freq_cp{subjects}] = select_single_subject_channel_data(channels,fa_cond_tr_freq{subjects});
                
                [fa_cond_tr_long_cp{subjects}] = select_single_subject_channel_data(channels,fa_cond_tr_long{subjects});
                [fa_cond_tr_short_cp{subjects}] = select_single_subject_channel_data(channels,fa_cond_tr_short{subjects});
                
                
%                 
%                 channels = {'C1', 'Cz', 'C2'};
                fa_cond_tr_rare_c{subjects} = select_single_subject_channel_data(channels,fa_cond_tr_rare{subjects});
                [fa_cond_tr_freq_c{subjects}] = select_single_subject_channel_data(channels,fa_cond_tr_freq{subjects});
                
                [fa_cond_tr_long_c{subjects}] = select_single_subject_channel_data(channels,fa_cond_tr_long{subjects});
                [fa_cond_tr_short_c{subjects}] = select_single_subject_channel_data(channels,fa_cond_tr_short{subjects});
                
%                 channels = {'FC1', 'FCz', 'FC2'};
%                 fa_cond_tr_rare_fc{subjects} = select_single_subject_channel_data(channels,fa_cond_tr_rare{subjects});
%                 [fa_cond_tr_freq_fc{subjects}] = select_single_subject_channel_data(channels,fa_cond_tr_freq{subjects});
%                 
%                 [fa_cond_tr_long_fc{subjects}] = select_single_subject_channel_data(channels,fa_cond_tr_long{subjects});
%                 [fa_cond_tr_short_fc{subjects}] = select_single_subject_channel_data(channels,fa_cond_tr_short{subjects});
%                 
%                 channels = {'P1', 'Pz', 'P2'};
%                 fa_cond_tr_rare_p{subjects} = select_single_subject_channel_data(channels,fa_cond_tr_rare{subjects});
%                 [fa_cond_tr_freq_p{subjects}] = select_single_subject_channel_data(channels,fa_cond_tr_freq{subjects});
%                 
%                 [fa_cond_tr_long_p{subjects}] = select_single_subject_channel_data(channels,fa_cond_tr_long{subjects});
%                 [fa_cond_tr_short_p{subjects}] = select_single_subject_channel_data(channels,fa_cond_tr_short{subjects});
%                 
%                 channels = {'P3', 'P4', 'P5', 'P6'};
%                 fa_cond_tr_rare_plat{subjects} = select_single_subject_channel_data(channels,fa_cond_tr_rare{subjects});
%                 [fa_cond_tr_freq_plat{subjects}] = select_single_subject_channel_data(channels,fa_cond_tr_freq{subjects});
%                 
%                 [fa_cond_tr_long_plat{subjects}] = select_single_subject_channel_data(channels,fa_cond_tr_long{subjects});
%                 [fa_cond_tr_short_plat{subjects}] = select_single_subject_channel_data(channels,fa_cond_tr_short{subjects});
%                 
%                 
            end
            
    end
    %
end

%lets calculate grand average for each coherence level
cfg = [];
coh_avg_button{1} = ft_timelockgrandaverage(cfg,coh_30{:});
coh_avg_button{1}.cfg = [];
coh_avg_button{2} = ft_timelockgrandaverage(cfg,coh_40{:});
coh_avg_button{2}.cfg = [];
coh_avg_button{3} = ft_timelockgrandaverage(cfg,coh_50{:});
coh_avg_button{3}.cfg = [];
%and for each c1ondition
cond_avg_button{1} = ft_timelockgrandaverage(cfg,cond_1{:});
cond_avg_button{1}.cfg = [];
cond_avg_button{2} = ft_timelockgrandaverage(cfg,cond_2{:});
cond_avg_button{2}.cfg = [];
cond_avg_button{3} = ft_timelockgrandaverage(cfg,cond_3{:});
cond_avg_button{3}.cfg = [];
cond_avg_button{4} = ft_timelockgrandaverage(cfg,cond_4{:});
cond_avg_button{4}.cfg = [];


channels = {'CP1', 'CPz', 'CP2'};
cond_avg_ALL = calculate_timelockGrandAverage_and_se(channels,cond_tr_cp_all,numberSubjects);

keyboard;

switch timelockCondition
    case 'buttonPress'
        if faFlag
            cfg = [];
            fa_cond_avg_button{1} = ft_timelockgrandaverage(cfg,fa_cond_1{:});
            fa_cond_avg_button{1}.cfg = [];
            fa_cond_avg_button{2} = ft_timelockgrandaverage(cfg,fa_cond_2{:});
            fa_cond_avg_button{2}.cfg = [];
            fa_cond_avg_button{3} = ft_timelockgrandaverage(cfg,fa_cond_3{:});
            fa_cond_avg_button{3}.cfg = [];
            fa_cond_avg_button{4} = ft_timelockgrandaverage(cfg,fa_cond_4{:});
            cond_avg_button{4}.cfg = [];
            
            
            
            channels = 'eeg';
            fa_cond_tr_freq_avg = calculate_timelockGrandAverage_and_se(channels,fa_cond_tr_freq,numberSubjects);
            fa_cond_tr_rare_avg = calculate_timelockGrandAverage_and_se(channels,fa_cond_tr_rare,numberSubjects);
            
            fa_cond_tr_short_avg = calculate_timelockGrandAverage_and_se(channels,fa_cond_tr_short,numberSubjects);
            fa_cond_tr_long_avg = calculate_timelockGrandAverage_and_se(channels,fa_cond_tr_long,numberSubjects);
            
            %%-CP--------------------------------------------------------%%%%%
            channels = [];
            [fa_cond_tr_rare_avg_cp] =calculate_timelockGrandAverage_and_se(channels,fa_cond_tr_rare_cp,numberSubjects);
            [fa_cond_tr_freq_avg_cp] =calculate_timelockGrandAverage_and_se(channels,fa_cond_tr_freq_cp,numberSubjects);
            
            [fa_cond_tr_long_avg_cp] =calculate_timelockGrandAverage_and_se(channels,fa_cond_tr_long_cp,numberSubjects);
            [fa_cond_tr_short_avg_cp] =calculate_timelockGrandAverage_and_se(channels,fa_cond_tr_short_cp,numberSubjects);
            
            %%-C--------------------------------------------------------%%%%%
            channels = [];
            [fa_cond_tr_rare_avg_c] =calculate_timelockGrandAverage_and_se(channels,fa_cond_tr_rare_c,numberSubjects);
            [fa_cond_tr_freq_avg_c] =calculate_timelockGrandAverage_and_se(channels,fa_cond_tr_freq_c,numberSubjects);
            
            [fa_cond_tr_long_avg_c] =calculate_timelockGrandAverage_and_se(channels,fa_cond_tr_long_c,numberSubjects);
            [fa_cond_tr_short_avg_c] =calculate_timelockGrandAverage_and_se(channels,fa_cond_tr_short_c,numberSubjects);
            
            %%-FC--------------------------------------------------------%%%%%
%             channels = [];
%             [fa_cond_tr_rare_avg_fc] =calculate_timelockGrandAverage_and_se(channels,fa_cond_tr_rare_fc,numberSubjects);
%             [fa_cond_tr_freq_avg_fc] =calculate_timelockGrandAverage_and_se(channels,fa_cond_tr_freq_fc,numberSubjects);
%             
%             [fa_cond_tr_long_avg_fc] =calculate_timelockGrandAverage_and_se(channels,fa_cond_tr_long_fc,numberSubjects);
%             [fa_cond_tr_short_avg_fc] =calculate_timelockGrandAverage_and_se(channels,fa_cond_tr_short_fc,numberSubjects);
%             
            
            %%-P--------------------------------------------------------%%%%%
%             channels = [];
%             [fa_cond_tr_rare_avg_p] =calculate_timelockGrandAverage_and_se(channels,fa_cond_tr_rare_p,numberSubjects);
%             [fa_cond_tr_freq_avg_p] =calculate_timelockGrandAverage_and_se(channels,fa_cond_tr_freq_p,numberSubjects);
%             
%             [fa_cond_tr_long_avg_p] =calculate_timelockGrandAverage_and_se(channels,fa_cond_tr_long_p,numberSubjects);
%             [fa_cond_tr_short_avg_p] =calculate_timelockGrandAverage_and_se(channels,fa_cond_tr_short_p,numberSubjects);
%         
%         channels = [];
% [fa_cond_tr_rare_avg_plat] =calculate_timelockGrandAverage_and_se(channels,fa_cond_tr_rare_plat,numberSubjects);
% [fa_cond_tr_freq_avg_plat] =calculate_timelockGrandAverage_and_se(channels,fa_cond_tr_freq_plat,numberSubjects);
% 
% [fa_cond_tr_long_avg_plat] =calculate_timelockGrandAverage_and_se(channels,fa_cond_tr_long_plat,numberSubjects);
% [fa_cond_tr_short_avg_plat] =calculate_timelockGrandAverage_and_se(channels,fa_cond_tr_short_plat,numberSubjects);
% 


        
        end
end


channels = 'eeg';
cond_tr_freq_avg = calculate_timelockGrandAverage_and_se(channels,cond_tr_freq,numberSubjects);
cond_tr_rare_avg = calculate_timelockGrandAverage_and_se(channels,cond_tr_rare,numberSubjects);

cond_tr_short_avg = calculate_timelockGrandAverage_and_se(channels,cond_tr_short,numberSubjects);
cond_tr_long_avg = calculate_timelockGrandAverage_and_se(channels,cond_tr_long,numberSubjects);

%%-CPz--------------------------------------------------------------%%%%%%
channels = [];
[cond_tr_rare_avg_cp] =calculate_timelockGrandAverage_and_se(channels,cond_tr_rare_cp,numberSubjects);
[cond_tr_freq_avg_cp] =calculate_timelockGrandAverage_and_se(channels,cond_tr_freq_cp,numberSubjects);

[cond_tr_long_avg_cp] =calculate_timelockGrandAverage_and_se(channels,cond_tr_long_cp,numberSubjects);
[cond_tr_short_avg_cp] =calculate_timelockGrandAverage_and_se(channels,cond_tr_short_cp,numberSubjects);


% channels =[];
% [coh_30_avg_cp] = calculate_timelockGrandAverage_and_se(channels,coh_03_cp,numberSubjects);
% [coh_40_avg_cp] = calculate_timelockGrandAverage_and_se(channels,coh_04_cp,numberSubjects);
% [coh_50_avg_cp] = calculate_timelockGrandAverage_and_se(channels,coh_05_cp,numberSubjects);
% 



%-Cz--------------------------------------------------------------%%%%%%
channels = [];
[cond_tr_rare_avg_c] =calculate_timelockGrandAverage_and_se(channels,cond_tr_rare_c,numberSubjects);
[cond_tr_freq_avg_c] =calculate_timelockGrandAverage_and_se(channels,cond_tr_freq_c,numberSubjects);

[cond_tr_long_avg_c] =calculate_timelockGrandAverage_and_se(channels,cond_tr_long_c,numberSubjects);
[cond_tr_short_avg_c] =calculate_timelockGrandAverage_and_se(channels,cond_tr_short_c,numberSubjects);


% channels =[];
% [coh_30_avg_c] = calculate_timelockGrandAverage_and_se(channels,coh_03_c,numberSubjects);
% [coh_40_avg_c] = calculate_timelockGrandAverage_and_se(channels,coh_04_c,numberSubjects);
% [coh_50_avg_c] = calculate_timelockGrandAverage_and_se(channels,coh_05_c,numberSubjects);
% 
% 

%%-FCz--------------------------------------------------------------%%%%%%
% channels = [];
% [cond_tr_rare_avg_fc] =calculate_timelockGrandAverage_and_se(channels,cond_tr_rare_fc,numberSubjects);
% [cond_tr_freq_avg_fc] =calculate_timelockGrandAverage_and_se(channels,cond_tr_freq_fc,numberSubjects);
% 
% [cond_tr_long_avg_fc] =calculate_timelockGrandAverage_and_se(channels,cond_tr_long_fc,numberSubjects);
% [cond_tr_short_avg_fc] =calculate_timelockGrandAverage_and_se(channels,cond_tr_short_fc,numberSubjects);


% channels =[];
% [coh_30_avg_fc] = calculate_timelockGrandAverage_and_se(channels,coh_03_fc,numberSubjects);
% [coh_40_avg_fc] = calculate_timelockGrandAverage_and_se(channels,coh_04_fc,numberSubjects);
% [coh_50_avg_fc] = calculate_timelockGrandAverage_and_se(channels,coh_05_fc,numberSubjects);
% 

%%-Pz--------------------------------------------------------------%%%%%%
channels = [];
[cond_tr_rare_avg_p] =calculate_timelockGrandAverage_and_se(channels,cond_tr_rare_p,numberSubjects);
[cond_tr_freq_avg_p] =calculate_timelockGrandAverage_and_se(channels,cond_tr_freq_p,numberSubjects);

[cond_tr_long_avg_p] =calculate_timelockGrandAverage_and_se(channels,cond_tr_long_p,numberSubjects);
[cond_tr_short_avg_p] =calculate_timelockGrandAverage_and_se(channels,cond_tr_short_p,numberSubjects);

% 
% channels =[];
% [coh_30_avg_p] = calculate_timelockGrandAverage_and_se(channels,coh_03_p,numberSubjects);
% [coh_40_avg_p] = calculate_timelockGrandAverage_and_se(channels,coh_04_p,numberSubjects);
% [coh_50_avg_p] = calculate_timelockGrandAverage_and_se(channels,coh_05_p,numberSubjects);
% 

%%--Plateral--------------------------------------------------------------%%
channels = [];
[cond_tr_rare_avg_plat] =calculate_timelockGrandAverage_and_se(channels,cond_tr_rare_plat,numberSubjects);
[cond_tr_freq_avg_plat] =calculate_timelockGrandAverage_and_se(channels,cond_tr_freq_plat,numberSubjects);

[cond_tr_long_avg_plat] =calculate_timelockGrandAverage_and_se(channels,cond_tr_long_plat,numberSubjects);
[cond_tr_short_avg_plat] =calculate_timelockGrandAverage_and_se(channels,cond_tr_short_plat,numberSubjects);

% 
% channels =[];
% [coh_30_avg_plat] = calculate_timelockGrandAverage_and_se(channels,coh_03_plat,numberSubjects);
% [coh_40_avg_plat] = calculate_timelockGrandAverage_and_se(channels,coh_04_plat,numberSubjects);
% [coh_50_avg_plat] = calculate_timelockGrandAverage_and_se(channels,coh_05_plat,numberSubjects);
% 


% permest
% permutation test
% for trials
keyboard;
save_loc = fullfile(options.path.EEG.analysis,'perm_stat_trial_freq_rare.mat');
[stat] = permutation_test(cond_tr_freq_cp, cond_tr_rare_cp, [], 'condition', timelockCondition, save_loc);
[stat_length] = permutation_test(cond_tr_short_cp, cond_tr_long_cp, [], 'condition', timelockCondition, save_loc);


fsample = 100;
pos = plot_stats_permtest(stat, cond_tr_freq_avg, cond_tr_rare_avg,fsample, 0.2, 0, 15);
pos_length = plot_stats_permtest(stat_length, cond_tr_short_avg, cond_tr_long_avg,fsample, 0.2, 0, 15);



%%

plot_ERP_timeseries_stats(pos_length,stat_length,cond_tr_short_avg_cp,cond_tr_long_avg_cp, cond_tr_short_c,cond_tr_long_c)
plot_ERP_timeseries_stats(pos,stat,cond_tr_freq_avg_cp,cond_tr_rare_avg_cp, cond_tr_rare_cp,  cond_tr_freq_cp)


% permutation test false alarm
% for trials
if faFlag
save_loc = fullfile(options.path.EEG.analysis,'fa_perm_stat_trial_freq_rare.mat');
[stat_fa_length] = permutation_test(fa_cond_tr_short_cp, fa_cond_tr_long_cp, [], 'condition', timelockCondition, save_loc);
[stat_fa] = permutation_test(fa_cond_tr_freq_cp, fa_cond_tr_rare_cp, [], 'condition', timelockCondition, save_loc);


% plot stats
fsample = 100;
pos_fa = plot_stats_permtest(stat_fa,fa_cond_tr_freq_avg, fa_cond_tr_rare_avg,fsample, 0.2, 0, 15);
pos_fa_length = plot_stats_permtest(stat_fa_length,fa_cond_tr_short_avg, fa_cond_tr_long_avg,fsample, 0.2, 0, 15);


plot_ERP_timeseries_stats(pos_fa_length,stat_fa_length,fa_cond_tr_short_avg_cp,fa_cond_tr_long_avg_cp,fa_cond_tr_short_cp,fa_cond_tr_long_cp)
plot_ERP_timeseries_stats(pos_fa,stat_fa,fa_cond_tr_freq_avg_cp,fa_cond_tr_rare_avg_cp,fa_cond_tr_freq_cp,fa_cond_tr_rare_cp)


end



%% plotting
figure; 
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [200 100]);    
    set(gcf, 'Position',  [500, 500, 850, 1290])
 %%
 
%      set(gcf, 'PaperUnits', 'inches');
%     set(gcf, 'PaperSize', [200 100]);    
%     set(gcf, 'Position',  [500, 500, 2000, 1290])

 
     set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [200 100]);    
    set(gcf, 'Position',  [100, 100, 400, 245])

cl = cbrewer('div','RdBu',100);
sb_idx = 1;
start_time = -0.5;
%for t = 1:12
    %subplot(2,6,sb_idx)
    sb_idx = sb_idx + 1;
    cfg = [];
    endTime = start_time + 0.5;
    cfg.xlim = [start_time endTime];
    cfg.highlight = 'on'; 
    cfg.highlightchannel = {'P1','P2','Pz'};
    cfg.highlightsymbol = '^';
    %cfg.zlim = [-(2e-4), (2e-4)]; %[-3 3]; %[-(2e-4), (2e-4)]; %[-3 3];
    cfg.layout = 'easycapM1.mat';
    cfg.comment = [num2str(start_time),['s to '],num2str(endTime),'s'];
    cfg.FontName = 'Arial';
    cfg.commentpos = 'title';
    
    
    cfg.colormap = flip(cl);
    ft_topoplotER(cfg,cond_avg_ALL);
    

    tidyfig;
    start_time = start_time + 0.5;
    h = colorbar();
%end


%%     trial startv
figure
 
%      set(gcf, 'PaperUnits', 'inches');
%     set(gcf, 'PaperSize', [200 100]);    
%     set(gcf, 'Position',  [500, 500, 2000, 1290])
    
         set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [200 100]);    
    set(gcf, 'Position',  [100, 100, 400, 245])
    
    
cl = cbrewer('div','RdBu',100);
sb_idx = 1;
start_time = 2.75;
% for t = 1:18
%     subplot(3,6,sb_idx)
    sb_idx = sb_idx + 1;
    cfg = [];
    endTime = start_time + 0.5;
    cfg.xlim = [start_time endTime];
    cfg.highlight = 'on'; 
    cfg.highlightchannel = {'P1','P2','Pz'};
    cfg.highlightsymbol = '^';
    cfg.zlim = [-0.8 0.8];%[-(2e-4), (2e-4)]; %[-0.9 0.9]; %[-1 1]; %[-(2e-4), (2e-4)]; %[-3 3];
    cfg.layout = 'easycapM1.mat';
    cfg.comment = [num2str(start_time),['s to '],num2str(endTime),'s'];
    cfg.commentpos = 'title';
    
    
    cfg.colormap = flip(cl);
    ft_topoplotER(cfg,cond_avg_ALL);
    

    tidyfig;
    start_time = start_time + 0.5;
    colorbar();
% end

%% plot

cl = cbrewer('qual','Set1',3);
subplot(2,2,1)
hold on
h = shadedErrorBar(cond_tr_freq_avg_fc.time, cond_tr_freq_avg_fc.avg, cond_tr_freq_avg_fc.se, 'lineprops', '-k');
h.patch.FaceColor = cl(1,:);
h.mainLine.Color = cl(1,:);
h.mainLine.LineWidth = 0.5;
h.patch.FaceAlpha = 0.3;

h = shadedErrorBar(cond_tr_rare_avg_fc.time, cond_tr_rare_avg_fc.avg, cond_tr_rare_avg_fc.se, 'lineprops', '-k');
h.patch.FaceColor = cl(2,:);
h.mainLine.Color = cl(2,:);
h.mainLine.LineWidth = 0.5;
h.patch.FaceAlpha = 0.3;
%title('difference between freq and rare conditions CPz')
title('FCz')
%legend({'short','long'})
legend({'frequent','rare'})
tidyfig;
hold off





subplot(2,2,2)
hold on
h = shadedErrorBar(cond_tr_freq_avg_c.time, cond_tr_freq_avg_c.avg, cond_tr_freq_avg_c.se, 'lineprops', '-k');
h.patch.FaceColor = cl(1,:);
h.mainLine.Color = cl(1,:);
h.mainLine.LineWidth = 0.5;
h.patch.FaceAlpha = 0.3;

h = shadedErrorBar(cond_tr_rare_avg_c.time, cond_tr_rare_avg_c.avg, cond_tr_rare_avg_c.se, 'lineprops', '-k');
h.patch.FaceColor = cl(2,:);
h.mainLine.Color = cl(2,:);
h.mainLine.LineWidth = 0.5;
h.patch.FaceAlpha = 0.3;
%title('difference between freq and rare conditions CPz')
title(' Cz')
%legend({'short','long'})
legend({'frequent','rare'})
hold off
tidyfig;




subplot(2,2,3)
hold on
h = shadedErrorBar(cond_tr_freq_avg_cp.time, cond_tr_freq_avg_cp.avg, cond_tr_freq_avg_cp.se, 'lineprops', '-k');
h.patch.FaceColor = cl(1,:);
h.mainLine.Color = cl(1,:);
h.mainLine.LineWidth = 0.5;
h.patch.FaceAlpha = 0.3;

h = shadedErrorBar(cond_tr_rare_avg_cp.time, cond_tr_rare_avg_cp.avg, cond_tr_rare_avg_cp.se, 'lineprops', '-k');
h.patch.FaceColor = cl(2,:);
h.mainLine.Color = cl(2,:);
h.mainLine.LineWidth = 0.5;
h.patch.FaceAlpha = 0.3;


title('CP')
xlabel('time(s) 0 start of event')
legend({'frequent','rare'})
hold off
tidyfig;


subplot(2,2,4)
hold on
h = shadedErrorBar(cond_tr_freq_avg_p.time, cond_tr_freq_avg_p.avg, cond_tr_freq_avg_p.se, 'lineprops', '-k');
h.patch.FaceColor = cl(1,:);
h.mainLine.Color = cl(1,:);
h.mainLine.LineWidth = 0.5;
h.patch.FaceAlpha = 0.3;

h = shadedErrorBar(cond_tr_rare_avg_p.time, cond_tr_rare_avg_p.avg, cond_tr_rare_avg_p.se, 'lineprops', '-k');
h.patch.FaceColor = cl(2,:);
h.mainLine.Color = cl(2,:);
h.mainLine.LineWidth = 0.5;
h.patch.FaceAlpha = 0.3;


title('P')
xlabel('time(s) 0 start of event')
legend({'frequent','rare'})
hold off
tidyfig;


%%

figure
hold on
h = shadedErrorBar(cond_avg_ALL.time, cond_tr_freq_avg_plat.avg, cond_avg_ALL.se, 'lineprops', '-k');
h.patch.FaceColor = cl(1,:);
h.mainLine.Color = cl(1,:);
h.mainLine.LineWidth = 0.5;
h.patch.FaceAlpha = 0.3;
% 
% h = shadedErrorBar(cond_tr_rare_avg_plat.time, cond_tr_rare_avg_plat.avg, cond_tr_rare_avg_plat.se, 'lineprops', '-k');
% h.patch.FaceColor = cl(2,:);
% h.mainLine.Color = cl(2,:);
% h.mainLine.LineWidth = 0.5;
% h.patch.FaceAlpha = 0.3;


title('cp average across conditions')
xlabel('time(s) 0 start of event')
%legend({'frequent','rare'})
hold off
tidyfig;
