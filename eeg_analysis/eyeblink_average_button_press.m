close all;
clear all;

subj_list = [16, 18:21, 24, 26, 27, 28, 29, 31, 32, 33, 34, 35, 39, 41, 42, 43, 50, 51, 52, 54, 55, 57, 58];
%subj_list = [  50, 51, 52, 54, 55, 57, 58];

[EEGdir,EEGdirdata,scriptdir,nSess,nS] = setup_EEG_session(subj_list);
source_density = 0; % flag if source density data is used
EEGpreproc = '/Volumes/LaCie/data_preproc';  % path to behav data all subjs
load_name = fullfile(EEGpreproc,'behav_data_all_subjs_all3');
load(load_name)
reference_type = 'LM_RM';


%% 
for  sj = 1:length(subj_list)
    subID = subj_list(sj);
    rt_flag = 0;
    baseline = 0;
    clear data_load
    
    switch reference_type
        
        case 'LM_RM'
            if source_density
                eegdat_fname = fullfile(EEGdir,['csd_response_locked_EEG_dat_',sprintf('sub%03.0f',subID),'.mat']);
                save_name_coherence = fullfile(EEGdir,['csd_response_locked_EEG_avg_coherence.mat']);
                save_name_condition = fullfile(EEGdir,['csd_response_locked_EEG_avg_condition.mat']);
            else
                eegdat_fname = fullfile(EEGdir,['response_locked_EEG_dat',sprintf('sub%03.0f',subID),'.mat']);
                save_name_coherence = fullfile(EEGdir,['response_locked_EEG_avg_coherence.mat']);
                save_name_condition = fullfile(EEGdir,['response_locked_EEG_avg_condition.mat']);
            end
            
        case 'averaged_electrodes'
            if source_density
                eegdat_fname = fullfile(EEGdir,'average_reference',['csd_response_locked_EEG_dat_',sprintf('sub%03.0f',subID),'.mat']);
                save_name_coherence = fullfile(EEGdir,'average_reference',['csd_response_locked_EEG_avg_coherence.mat']);
                save_name_condition = fullfile(EEGdir,'average_reference',['csd_response_locked_EEG_avg_condition.mat']);
            else
                eegdat_fname = fullfile(EEGdir,'average_reference',['response_locked_EEG_dat',sprintf('sub%03.0f',subID),'.mat']);
                save_name_coherence = fullfile(EEGdir,'average_reference',['response_locked_EEG_avg_coherence.mat']);
                save_name_condition = fullfile(EEGdir,'average_reference',['response_locked_EEG_avg_condition.mat']);
            end
            
    end
    
    
    
    data_load = load(eegdat_fname);
    data = data_load.data_append;
    
    
    
    conditions_perm = [1 2; 3 4];
    
    for con = 1:2
        if rt_flag
            idx_con = find((data.trialinfo(:,9) == conditions_perm(con,1) | data.trialinfo(:,9) == conditions_perm(con,2))& data.trialinfo(:,7) == 1 & data.trialinfo(:,2) <= 3.5);
        else
            idx_con = find((data.trialinfo(:,9) == conditions_perm(con,1) | data.trialinfo(:,9) == conditions_perm(con,2))& data.trialinfo(:,7) == 1);
        end
        
        cfg = [];
        cfg.trials = idx_con;
        cfg.channel = 'VEOG';
        data_condition{con} = ft_selectdata(cfg,data);
        
        cfg = [];
        data_avg{con} = ft_timelockanalysis(cfg,data_condition{con});
        
        if baseline
            cfg = [];
            cfg.baseline = baseline_window;
            
            data_avg{con} = ft_timelockbaseline(cfg,data_avg{con});
        end
    end
    
    
    VEOG_avg_sj_freq{sj} = data_avg{1};
    VEOG_avg_sj_rare{sj} = data_avg{2};
    
    for con = 1:2
        
        if rt_flag
            idx_con = find((data.trialinfo(:,9) == conditions_perm(con,1) | data.trialinfo(:,9) == conditions_perm(con,2))& data.trialinfo(:,7) == 1 & data.trialinfo(:,2) <= 3.5);
        else
            idx_con = find((data.trialinfo(:,9) == conditions_perm(con,1) | data.trialinfo(:,9) == conditions_perm(con,2))& data.trialinfo(:,7) == 1);
        end
        cfg = [];
        cfg.trials = idx_con;
        cfg.channel = 'HEOG';
        data_condition{con} = ft_selectdata(cfg,data);
        
        cfg = [];
        data_avg{con} = ft_timelockanalysis(cfg,data_condition{con});
        
        
        if baseline
            cfg = [];
            cfg.baseline = baseline_window;
            
            data_avg{con} = ft_timelockbaseline(cfg,data_avg{con});
        end
    end
    
    HEOG_avg_sj_freq{sj} = data_avg{1};
    HEOG_avg_sj_rare{sj} = data_avg{2};
    
end

%% 
cfg = [];
[cond_tr_rare_avg_VEOG] = ft_timelockgrandaverage(cfg,VEOG_avg_sj_rare{:});
% calculate se
[se] = calculate_se(cond_tr_rare_avg_VEOG.var,nS);
cond_tr_rare_avg_VEOG.se = se;

cfg = [];
[cond_tr_freq_avg_VEOG] = ft_timelockgrandaverage(cfg,VEOG_avg_sj_freq{:});
% calculate se
[se] = calculate_se(cond_tr_freq_avg_VEOG.var,nS);
cond_tr_freq_avg_VEOG.se = se;

cfg = [];
[cond_tr_rare_avg_HEOG] = ft_timelockgrandaverage(cfg,HEOG_avg_sj_rare{:});
% calculate se
[se] = calculate_se(cond_tr_rare_avg_HEOG.var,nS);
cond_tr_rare_avg_HEOG.se = se;

cfg = [];
[cond_tr_freq_avg_HEOG] = ft_timelockgrandaverage(cfg,HEOG_avg_sj_freq{:});
% calculate se
[se] = calculate_se(cond_tr_freq_avg_HEOG.var,nS);
cond_tr_freq_avg_HEOG.se = se;



cl = cbrewer('qual','Set1',3);
subplot(1,2,1)
hold on
h = shadedErrorBar(cond_tr_rare_avg_VEOG.time,cond_tr_rare_avg_VEOG.avg,cond_tr_rare_avg_VEOG.se, 'lineprops', '-k');
h.patch.FaceColor = cl(2,:);
h.mainLine.Color = cl(2,:);

h = shadedErrorBar(cond_tr_freq_avg_VEOG.time,cond_tr_freq_avg_VEOG.avg,cond_tr_freq_avg_VEOG.se, 'lineprops', '-k');
h.patch.FaceColor = cl(1,:);
h.mainLine.Color = cl(1,:);
xlabel('time sec, 0 = button press')
ylabel('mV')
title('Vertical eye movements')
legend({'rare','frequent'})
tidyfig;

hold off

subplot(1,2,2)

hold on
h = shadedErrorBar(cond_tr_rare_avg_HEOG.time,cond_tr_rare_avg_HEOG.avg,cond_tr_rare_avg_HEOG.se, 'lineprops', '-k');
h.patch.FaceColor = cl(2,:);
h.mainLine.Color = cl(2,:);

h = shadedErrorBar(cond_tr_freq_avg_HEOG.time,cond_tr_freq_avg_HEOG.avg,cond_tr_freq_avg_HEOG.se, 'lineprops', '-k');
h.patch.FaceColor = cl(1,:);
h.mainLine.Color = cl(1,:);
xlabel('time sec, 0 = button press')
ylabel('mV')
title('horizontal eye movements')
legend({'rare','frequent'})
tidyfig;
hold off

