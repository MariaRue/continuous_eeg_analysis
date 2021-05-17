clear all;
close all;

subj_list = [16, 18:21, 24, 26, 27, 28, 29, 31, 33, 34, 35, 39, 41, 42, 43, 50, 51, 52, 54, 55, 57, 58]; % subject 40 and 47 removed - not working, no idea why, investigate
%subj_list = [40,47];

[EEGdir,EEGdirdata,scriptdir,nSess,nS] = setup_EEG_session(subj_list);
source_density = 0; % flag if source density data is used
EEGpreproc = '/Volumes/LaCie/data_preproc';  % path to behav data all subjs
load_name = fullfile(EEGpreproc,'behav_data_all_subjs_all3');
load(load_name)
nS = length(subj_list);
reference_type = 'LM_RM';
%% get data timelocked to trial start

for  sj = 1:length(subj_list)
    subID = subj_list(sj);
    disp('subject:')
    disp(subID)
    
    EEGdatadir =  fullfile(EEGdir,'subjects',sprintf('sub%03.0f',subID),'eeg');
    clear data_append data
    
    
    switch reference_type
        case 'LM_RM'
            
            if source_density
                eegdat_fname = fullfile(EEGdir,['csd_trial_start_locked_EEG_dat_',sprintf('sub%03.0f',subID),'.mat']);
            else
                eegdat_fname = fullfile(EEGdir,'test_April_2020',['trial_start_locked_EEG_dat_',sprintf('sub%03.0f',subID),'.mat']);
                
            end
        case 'averaged_electrodes'
            if source_density
                eegdat_fname = fullfile(EEGdir,'average_reference',['csd_trial_start_locked_EEG_dat_',sprintf('sub%03.0f',subID),'.mat']);
            else
                eegdat_fname = fullfile(EEGdir,'average_reference',['trial_start_locked_EEG_dat_',sprintf('sub%03.0f',subID),'.mat']);
            end
            
    end
    
    if subID == 55
        nSess = 1:5;
        
    elseif subID == 47
        nSess =[1,2,3,5];
        
            elseif subID == 40
        nSess =[2,3,4,5];
        
    else
        nSess = 1:6;
        
    end
    
    if 1 == 1 %eegdat_fname ~= 2
        
        for sess_count = 1:length(nSess)
            i = nSess(sess_count);
            disp('session:')
            disp(i)
            cfg = [];
            
            switch reference_type
                
                case 'LM_RM'
                    if source_density
                        cfg.dataset = fullfile(EEGdatadir,sprintf('cnanart_cfMdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',subID,i));
                    else
                        cfg.dataset = fullfile(EEGdatadir,sprintf('nanart_fMdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',subID,i));
                        
                        
                    end
                case 'averaged_electrodes'
                    if source_density
                        cfg.dataset = fullfile(EEGdatadir,'average_reference',sprintf('cnanart_cfMdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',subID,i));
                    else
                        
                        cfg.dataset = fullfile(EEGdatadir,'average_reference',sprintf('nanart_fMdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',subID,i));
                    end
            end
            
            cfg.subject_responses = all_responses(all_responses(:,10)==i & all_responses(:,12)==subID,:);
            stream_sj = unique(all_responses(all_responses(:,12)==subID,11));
            disp('stream_sj:')
            disp(stream_sj)
            cfg.mean_stim_streams = mean_stim_streams_org{stream_sj,i};
            cfg.trialdef.prestim = 2;
            cfg.trialdef.poststim = 8;
            cfg.timelock_event = 'trial start';
            
            
            
            cfg.trialfun  = 'ft_trialfun_continuous_eeg';
            cfg = ft_definetrial(cfg);
            
            
            
            data{sess_count} = ft_preprocessing(cfg);
            
            
        end
        
        
        cfg = [];
        
        data_append =  ft_appenddata(cfg,data{:});
        save(eegdat_fname,'data_append');
        
        cd (scriptdir)
        
    end
end

%% put all subjs into one dataframe

for  sj = 1:length(subj_list)
    subID = subj_list(sj);
    clear data_load
    rtFlag = 1;
    switch reference_type
        
        case 'LM_RM'
            if source_density
                eegdat_fname = fullfile(EEGdir,['csd_trial_start_locked_EEG_dat_',sprintf('sub%03.0f',subID),'.mat']);
                save_name_coherence = fullfile(EEGdir,['csd_trial_start_locked_EEG_avg_coherence.mat']);
                save_name_condition = fullfile(EEGdir,['csd_trial_start_locked_EEG_avg_condition.mat']);
            else
                eegdat_fname = fullfile(EEGdir,'test_April_2020',['trial_start_locked_EEG_dat_',sprintf('sub%03.0f',subID),'.mat']);
                save_name_coherence = fullfile(EEGdir,'test_April_2020',['trial_start_locked_EEG_avg_coherenceN.mat']);
                save_name_condition = fullfile(EEGdir,'test_April_2020',['trial_start_locked_EEG_avg_conditionN.mat']);
            end
        case 'averaged_electrodes'
            if source_density
                eegdat_fname = fullfile(EEGdir,'average_reference',['csd_trial_start_locked_EEG_dat_',sprintf('sub%03.0f',subID),'.mat']);
                save_name_coherence = fullfile(EEGdir,'average_reference',['csd_trial_start_locked_EEG_avg_coherence.mat']);
                save_name_condition = fullfile(EEGdir,'average_reference',['csd_trial_start_locked_EEG_avg_condition.mat']);
            else
                eegdat_fname = fullfile(EEGdir,'average_reference',['trial_start_locked_EEG_dat_',sprintf('sub%03.0f',subID),'.mat']);
                save_name_coherence = fullfile(EEGdir,'average_reference',['trial_start_locked_EEG_avg_coherence.mat']);
                save_name_condition = fullfile(EEGdir,'average_reference',['trial_start_locked_EEG_avg_condition.mat']);
            end
            
    end
    
    baseline = 1;
    
    data_load = load(eegdat_fname);
    data{sj} = data_load.data_append;
    % change labels to easycap
    [easy_cap_labels] = change_electrode_labels(data{sj}.label);
    
    data{sj}.label = easy_cap_labels;
    
    % average for each coherence level for each sj
    data_avg =  compute_average_single_subject_level(data{sj},'coherence',rtFlag,[-2 -1],baseline);
    coh_30{sj} = data_avg{1};
    coh_40{sj} = data_avg{2};
    coh_50{sj} = data_avg{3};
    
    
    % average for each condition level for each sj
    data_avg_con =  compute_average_single_subject_level(data{sj},'condition',rtFlag,[-2 -1],baseline);
    cond_1{sj} = data_avg_con{1};
    cond_2{sj} = data_avg_con{2};
    cond_3{sj} = data_avg_con{3};
    cond_4{sj} = data_avg_con{4};
    
    
    % prepare data for permutation test for freq vs rare conditions
    data_avg_con =  compute_average_single_subject_level(data{sj},'perm_test_rare_freq_trial',rtFlag,[-2 -1],baseline);
    cond_tr_freq{sj} = data_avg_con{1};
    cond_tr_rare{sj} = data_avg_con{2};
    
    % prepare data for permutation test for long vs short conditions
    data_avg_con =  compute_average_single_subject_level(data{sj},'perm_test_short_long_trial',rtFlag,[-2 -1],baseline);
    cond_tr_short{sj} = data_avg_con{1};
    cond_tr_long{sj} = data_avg_con{2};
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % condition rare vs freq condition
    
    cfg = [];
    cfg.channel = {'CP1', 'CPz', 'CP2'};
    cfg.avgoverchan = 'yes';
    cfg.nanmean = 'yes';
    [cond_tr_rare_cp{sj}] = ft_selectdata(cfg,  cond_tr_rare{sj});
    [cond_tr_freq_cp{sj}] = ft_selectdata(cfg, cond_tr_freq{sj});
    
    cfg = [];
    cfg.channel = {'FC1', 'FCz', 'FC2'};
    cfg.avgoverchan = 'yes';
    cfg.nanmean = 'yes';
    [cond_tr_rare_fc{sj}] = ft_selectdata(cfg,  cond_tr_rare{sj});
    [cond_tr_freq_fc{sj}] = ft_selectdata(cfg, cond_tr_freq{sj});
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % condition long vs short
    
    cfg = [];
    cfg.channel = {'FC1', 'FCz', 'FC2'};
    cfg.avgoverchan = 'yes';
    cfg.nanmean = 'yes';
    [cond_tr_long_fc{sj}] = ft_selectdata(cfg,  cond_tr_long{sj});
    [cond_tr_short_fc{sj}] = ft_selectdata(cfg, cond_tr_short{sj});
    
    
    
    cfg = [];
    cfg.channel = {'CP1', 'CPz', 'CP2'};
    cfg.avgoverchan = 'yes';
    cfg.nanmean = 'yes';
    [cond_tr_short_cp{sj}] = ft_selectdata(cfg,  cond_tr_short{sj});
    [cond_tr_long_cp{sj}] = ft_selectdata(cfg, cond_tr_long{sj});
    
    
    %%%%%%%%%%%%%%%%%%%
    % coherence
    
    cfg = [];
    cfg.channel = {'FC1', 'FCz', 'FC2'};
    cfg.avgoverchan = 'yes';
    cfg.nanmean = 'yes';
    [coh_03_fc{sj}] = ft_selectdata(cfg, coh_30{sj});
    [coh_04_fc{sj}] = ft_selectdata(cfg, coh_40{sj});
    [coh_05_fc{sj}] = ft_selectdata(cfg, coh_50{sj});
    
    
    cfg = [];
    cfg.channel = {'CP1', 'CPz', 'CP2'};
    cfg.avgoverchan = 'yes';
    cfg.nanmean = 'yes';
    [coh_03_cp{sj}] = ft_selectdata(cfg,  coh_30{sj});
    [coh_04_cp{sj}] = ft_selectdata(cfg, coh_40{sj});
    [coh_05_cp{sj}] = ft_selectdata(cfg, coh_50{sj});
    
    
    %     data_avg_con =  compute_average_single_subject_level(data{sj},'perm_test_rare_freq_fa',0,[-2 -1]);
    %     cond_fa_freq{sj} = data_avg_con{1};
    %     cond_fa_rare{sj} = data_avg_con{2};
    
        cfg = [];
    cfg.channel = {'VEOG'};
    cfg.avgoverchan = 'yes';
    cfg.nanmean = 'yes';
    [cond_tr_rare_VEOG{sj}] = ft_selectdata(cfg,  cond_tr_rare{sj});
    [cond_tr_freq_VEOG{sj}] = ft_selectdata(cfg, cond_tr_freq{sj});
    
    
    
    cfg = [];
    cfg.channel = {'HEOG'};
    cfg.avgoverchan = 'yes';
    cfg.nanmean = 'yes';
    [cond_tr_rare_HEOG{sj}] = ft_selectdata(cfg,  cond_tr_rare{sj});
    [cond_tr_freq_HEOG{sj}] = ft_selectdata(cfg, cond_tr_freq{sj});% 

end
%
% lets calculate grand average for each coherence level
cfg = [];
coh_avg{1} = ft_timelockgrandaverage(cfg,coh_30{:});
coh_avg{2} = ft_timelockgrandaverage(cfg,coh_40{:});
coh_avg{3} = ft_timelockgrandaverage(cfg,coh_50{:});

% and for each condition
cond_avg{1} = ft_timelockgrandaverage(cfg,cond_1{:});


cond_avg{2} = ft_timelockgrandaverage(cfg,cond_2{:});


cond_avg{3} = ft_timelockgrandaverage(cfg,cond_3{:});


cond_avg{4} = ft_timelockgrandaverage(cfg,cond_4{:});
keyboard;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare perm data
cfg = [];
cfg.channels = 'eeg';
cond_tr_freq_avg = ft_timelockgrandaverage(cfg,cond_tr_freq{:});

cfg = [];
cfg.channels = 'eeg';
cond_tr_rare_avg = ft_timelockgrandaverage(cfg,cond_tr_rare{:});

% for short vs long trials

cfg = [];
cfg.channels = 'eeg';
cond_tr_short_avg = ft_timelockgrandaverage(cfg,cond_tr_short{:});

cfg = [];
cfg.channels = 'eeg';
cond_tr_long_avg = ft_timelockgrandaverage(cfg,cond_tr_long{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% conditon freq vs rare average cp and fc

[cond_tr_rare_avg_cp] = ft_timelockgrandaverage(cfg, cond_tr_rare_cp{:});
% calculate se
[se] = calculate_se(cond_tr_rare_avg_cp.var,nS);
cond_tr_rare_avg_cp.se = se;


[cond_tr_freq_avg_cp] = ft_timelockgrandaverage(cfg, cond_tr_freq_cp{:});
% calculate se
[se] = calculate_se(cond_tr_freq_avg_cp.var,nS);
cond_tr_freq_avg_cp.se = se;



[cond_tr_rare_avg_fc] = ft_timelockgrandaverage(cfg, cond_tr_rare_fc{:});
% calculate se
[se] = calculate_se(cond_tr_rare_avg_fc.var,nS);
cond_tr_rare_avg_fc.se = se;


[cond_tr_freq_avg_fc] = ft_timelockgrandaverage(cfg, cond_tr_freq_fc{:});
% calculate se
[se] = calculate_se(cond_tr_freq_avg_fc.var,nS);
cond_tr_freq_avg_fc.se = se;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% long vs short trials cp fc


[cond_tr_long_avg_cp] = ft_timelockgrandaverage(cfg, cond_tr_long_cp{:});
% calculate se
[se] = calculate_se(cond_tr_long_avg_cp.var,nS);
cond_tr_long_avg_cp.se = se;


[cond_tr_short_avg_cp] = ft_timelockgrandaverage(cfg, cond_tr_short_cp{:});
% calculate se
[se] = calculate_se(cond_tr_short_avg_cp.var,nS);
cond_tr_short_avg_cp.se = se;



[cond_tr_long_avg_fc] = ft_timelockgrandaverage(cfg, cond_tr_long_fc{:});
% calculate se
[se] = calculate_se(cond_tr_long_avg_fc.var,nS);
cond_tr_long_avg_fc.se = se;


[cond_tr_short_avg_fc] = ft_timelockgrandaverage(cfg, cond_tr_short_fc{:});
% calculate se
[se] = calculate_se(cond_tr_short_avg_fc.var,nS);
cond_tr_short_avg_fc.se = se;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% coherence

coh_03_avg_fc = ft_timelockgrandaverage(cfg, coh_03_fc{:});
[se] = calculate_se(coh_03_avg_fc.var,nS);
coh_03_avg_fc.se = se;

coh_04_avg_fc = ft_timelockgrandaverage(cfg, coh_04_fc{:});
[se] = calculate_se(coh_04_avg_fc.var,nS);
coh_04_avg_fc.se = se;

coh_05_avg_fc = ft_timelockgrandaverage(cfg, coh_05_fc{:});
[se] = calculate_se(coh_05_avg_fc.var,nS);
coh_05_avg_fc.se = se;



coh_03_avg_cp = ft_timelockgrandaverage(cfg, coh_03_cp{:});
[se] = calculate_se(coh_03_avg_cp.var,nS);
coh_03_avg_cp.se = se;

coh_04_avg_cp = ft_timelockgrandaverage(cfg, coh_04_cp{:});
[se] = calculate_se(coh_04_avg_cp.var,nS);
coh_04_avg_cp.se = se;

coh_05_avg_cp = ft_timelockgrandaverage(cfg, coh_05_cp{:});
[se] = calculate_se(coh_05_avg_cp.var,nS);
coh_05_avg_cp.se = se;


cfg = [];
[cond_tr_freq_avg_VEOG] = ft_timelockgrandaverage(cfg,cond_tr_freq_VEOG{:});
% calculate se
[se] = calculate_se(cond_tr_freq_avg_VEOG.var,nS);
cond_tr_freq_avg_VEOG.se = se;

cfg = [];
[cond_tr_freq_avg_HEOG] = ft_timelockgrandaverage(cfg,cond_tr_freq_HEOG{:});
% calculate se
[se] = calculate_se(cond_tr_freq_avg_HEOG.var,nS);
cond_tr_freq_avg_HEOG.se = se;


cfg = [];
[cond_tr_rare_avg_VEOG] = ft_timelockgrandaverage(cfg,cond_tr_rare_VEOG{:});
% calculate se
[se] = calculate_se(cond_tr_rare_avg_VEOG.var,nS);
cond_tr_rare_avg_VEOG.se = se;

cfg = [];
[cond_tr_rare_avg_HEOG] = ft_timelockgrandaverage(cfg,cond_tr_rare_HEOG{:});
% calculate se
[se] = calculate_se(cond_tr_rare_avg_HEOG.var,nS);
cond_tr_rare_avg_HEOG.se = se;

keyboard;
% cfg = [];
% cond_fa_freq_avg = ft_timelockgrandaverage(cfg,cond_fa_freq{:});
% cond_fa_rare_avg = ft_timelockgrandaverage(cfg,cond_fa_rare{:});

% permutation test
% for trials
save_loc = fullfile(EEGdir,'test_April_2020','perm_stat_trial_freq_rare.mat');
[stat] = permutation_test(cond_tr_freq, cond_tr_rare, [], 'condition', 'trial start', save_loc);

[stat_length] = permutation_test(cond_tr_long, cond_tr_short, [], 'condition', 'trial start', save_loc);
% plot stats
fsample = 100;
pos = plot_stats_permtest(stat, cond_tr_freq_avg, cond_tr_rare_avg,fsample, 0.2, 0, 15);
pos = plot_stats_permtest(stat_length,cond_tr_long_avg, cond_tr_short_avg,fsample, 0.2, 0, 15);
keyboard;

plot_ERP_timeseries_stats(pos,stat,cond_tr_freq_avg_cp,cond_tr_rare_avg_cp,cond_tr_freq_avg_fc, cond_tr_rare_avg_fc)
plot_ERP_timeseries_stats(pos,stat,cond_tr_short_avg_cp,cond_tr_long_avg_cp,cond_tr_short_avg_fc,cond_tr_long_avg_fc)
% delete cfg field - taking up too much memory to save variable...
for l = 1:4
    
    cond_avg{l}.cfg = [];
    
end

for l = 1:3
    coh_avg{l}.cfg = [];
end

save(save_name_condition,'cond_avg' );
save(save_name_coherence,'coh_avg' );
%% average across subjects and conditions for each coherence level - timelocked to trial start


lim = quantile(coh_avg{1}.avg(:),[0.1 0.9]);

cl = cbrewer('seq','Blues',12);
cl =  cl([6 10 12],:);
minlim = lim(1);
maxlim = lim(2);

cfg = [];
cfg.channel = {'CP1','CPz','CP2'};
%cfg.channel = {'CPz'};
cfg.layout = 'easycapM1.mat';
% cfg.ylim = [minlim maxlim];
%cfg.graphcolor = ['b','r','k'];


cfg.graphcolor = cl;
cfg.linewidth = 3;


ft_singleplotER(cfg,coh_avg{:});
%ft_movieplotER(cfg, coh_avg{1});

% YLIM - zero in middle!!! and make legend work!!!
% it *should have zero in the middle.
% and most importantly it should be the same in any related series of plots
% (e.g., all coherence levels)
legend({'30%', '40%', '50%'},'FontSize',14)
title('Averaged ERP across Subjects and conditions for different coherence levels','FontSize',14)
xlabel('time (s) - trial start at 0','FontSize',14)
set(gca,'FontSize',18)

%% 
figure
cl = cbrewer('qual','Set1',3);
figure 
hold on 
h = shadedErrorBar(cond_tr_rare_avg_VEOG.time,cond_tr_rare_avg_VEOG.avg,cond_tr_rare_avg_VEOG.se, 'lineprops', '-k');
h.patch.FaceColor = cl(1,:);
h.mainLine.Color = cl(1,:);

h = shadedErrorBar(cond_tr_freq_avg_VEOG.time,cond_tr_freq_avg_VEOG.avg,cond_tr_freq_avg_VEOG.se, 'lineprops', '-k');
h.patch.FaceColor = cl(2,:);
h.mainLine.Color = cl(2,:);
title('VEOG averaged across participants for freq and rare conditions')
legend({'rare','freq'})
tidyfig; 

%%
% now plot topoplot
figure;
sb_idx = 1;

lim = quantile(coh_avg{1}.avg(:),[0.1 0.9]);

minlim = lim(1);
maxlim = -lim(1);
for i = 1:3
    start_time = -1;
    for t = 1:8
        cfg = [];
        cfg.xlim = [start_time start_time + 1];
        start_time = start_time + 1;
        cfg.zlim = [minlim maxlim];
        cfg.layout = 'easycapM1.mat';
        subplot(3,8,sb_idx)
        sb_idx = sb_idx + 1;
        ft_topoplotER(cfg,coh_avg{i}); colorbar
    end
end

subplot(3,8,1)
title('30% coherence timelocked to trial start')
subplot(3,8,2)
title('0 1')
subplot(3,8,3)
title('1 2sec')
subplot(3,8,4)
title('2 3sec')
subplot(3,8,5)
title('3 4sec')
subplot(3,8,6)
title('4 5sec')
subplot(3,8,7)
title('5 6sec')
subplot(3,8,8)
title('6 7sec')

subplot(3,8,9)
title('40% coherence')

subplot(3,8,17)
title('50% coherence')
%%
for i = 1:24
    subplot(3,8,i)
    tidyfig;
    
end


%% average across subjects and conditions for each coherence level - timelocked to trial start

lim = quantile(coh_avg{1}.avg(:),[0.1 0.9]);

cl = cbrewer('seq','Blues',12);
cl =  cl([6 10 12],:);
minlim = -lim(2);
maxlim = lim(2);

cfg = [];
cfg.channel = {'CP1','CPz','CP2'};
cfg.layout = 'easycapM1.mat';
% cfg.ylim = [minlim maxlim];
%
cfg.graphcolor = ['b','r','k'];
%cfg.graphcolor = cl;
cfg.linewidth = 3;
%%%%%%%%%%%%

figure
ft_multiplotER(cfg,coh_avg{:});


% YLIM - zero in middle!!! and make legend work!!!
% it *should have zero in the middle.
% and most importantly it should be the same in any related series of plots
% (e.g., all coherence levels)
legend({'30%', '40%', '50%'},'FontSize',14)
title('Averaged ERP across Subjects and conditions for different coherence levels','FontSize',14)
xlabel('time (s) - trial start at 0','FontSize',14)
set(gca,'FontSize',18)


%%
% now plot topoplot
figure;
sb_idx = 1;

lim = quantile(coh_avg{1}.avg(:),[0.1 0.9]);



minlim = lim(1);
maxlim = lim(2);
for i = 1:3
    
    start_time = -1;
    for t = 1:8
        cfg = [];
        cfg.xlim = [start_time start_time + 1];
        start_time = start_time + 1;
        cfg.zlim = [minlim maxlim];
        cfg.layout = 'easycapM1.mat';
        subplot(3,8,sb_idx)
        sb_idx = sb_idx + 1;
        ft_topoplotER(cfg,coh_avg{i}); colorbar
    end
end

subplot(3,8,1)
title('30% coherence timelocked to trial start at 0')
subplot(3,8,2)
title('0 1')
subplot(3,8,3)
title('1 2sec')
subplot(3,8,4)
title('2 3sec')
subplot(3,8,5)
title('3 4sec')
subplot(3,8,6)
title('4 5sec')
subplot(3,8,7)
title('5 6sec')
subplot(3,8,8)
title('6 7sec')

subplot(3,8,9)
title('40% coherence')

subplot(3,8,17)
title('50% coherence')
%%
for i = 1:24
    subplot(3,8,i)
    tidyfig;
    
end
%% average across subjects and coherence levels for each condition - timelocked to trial start
figure
lim = quantile(cond_avg{1}.avg(:),[0.1 0.9]);

cl = cbrewer('div','RdBu', 12);
cl = cl([4 1 9 12],:);
minlim = -10e-05;
maxlim = 10e-05;

cfg = [];
cfg.channel = {'CPz','CP1','CP2'};
cfg.layout = 'easycapM1.mat';
% cfg.ylim = [minlim maxlim];
%cfg.graphcolor = ['b','r','k','g'];
cfg.graphcolor = cl;
cfg.linewidth = 3;
%%%%%%%%%%%%

ft_singleplotER(cfg,cond_avg{:});

legend({'ITIS INTS', 'ITIS INTL', 'ITIL INTS','ITIL INTL'},'FontSize',14)
title('Averaged ERP across Subjects and coherences for different conditions','FontSize',14)
xlabel('time (s) - trial start at 0','FontSize',14)
set(gca,'FontSize',18)

%% now plot topoplot
figure;
sb_idx = 1;

lim = quantile(cond_avg{1}.avg(:),[0.1 0.9]);

minlim = lim(1);
maxlim = -lim(1);
for i = 1:4
    start_time = -1;
    for t = 1:8
        cfg = [];
        cfg.xlim = [start_time start_time + 1];
        start_time = start_time + 1;
        cfg.zlim = [minlim maxlim];
        cfg.layout = 'easycapM1.mat';
        subplot(4,8,sb_idx)
        sb_idx = sb_idx + 1;
        ft_topoplotER(cfg,cond_avg{i}); colorbar
    end
end

subplot(4,8,1)
title('ITIS INTS condition timelocked to trial start at 0')
subplot(4,8,2)
title('0 1sec')
subplot(4,8,3)
title('1 2sec')
subplot(4,8,4)
title('2 3sec')
subplot(4,8,5)
title('3 4sec')
subplot(4,8,6)
title('4 5sec')
subplot(4,8,7)
title('5 6sec')

subplot(4,8,9)
title('ITIS INTL')

subplot(4,8,17)
title('ITIL INTS')

subplot(4,8,25)
title('ITIL INTL')

%%
for i = 1:32
    subplot(4,8,i)
    tidyfig;
    
end
