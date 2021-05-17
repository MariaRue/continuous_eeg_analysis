close all;
clear all;

subj_list = [16,18:21, 24, 26, 27, 28, 29, 31, 33, 34, 35, 39,40,  41, 42, 43, 47,50, 51, 52, 54, 55, 57, 58];
%subj_list = [40, 47];
  baseline = 0;
    rtFlag = 0;
[EEGdir,EEGdirdata,scriptdir,nSess,nS] = setup_EEG_session(subj_list);
source_density = 0; % flag if source density data is used
EEGpreproc = '/Volumes/LaCie/data_preproc';  % path to behav data all subjs
load_name = fullfile(EEGpreproc,'behav_data_all_subjs_all3');
load(load_name)
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
                eegdat_fname = fullfile(EEGdir,['csd_response_locked_EEG_dat_',sprintf('sub%03.0f',subID),'.mat']);
            else
                eegdat_fname = fullfile(EEGdir,['response_locked_EEG_dat',sprintf('sub%03.0f',subID),'.mat']);
            end
        case 'averaged_electrodes'
            if source_density
                eegdat_fname = fullfile(EEGdir,'average_reference',['csd_response_locked_EEG_dat_',sprintf('sub%03.0f',subID),'.mat']);
            else
                eegdat_fname = fullfile(EEGdir,'average_reference',['response_locked_EEG_dat',sprintf('sub%03.0f',subID),'.mat']);
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
            cfg.trialdef.prestim = 7;
            cfg.trialdef.poststim = 4;
            cfg.timelock_event = 'button press';
            
 
            
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
    
    switch reference_type
        
        case 'LM_RM'
            if source_density
                eegdat_fname = fullfile(EEGdir,['csd_response_locked_EEG_dat_',sprintf('sub%03.0f',subID),'.mat']);
                save_name_coherence = fullfile(EEGdir,['csd_response_locked_EEG_avg_coherence.mat']);
                save_name_condition = fullfile(EEGdir,['csd_response_locked_EEG_avg_condition.mat']);
            else
                eegdat_fname = fullfile('/Volumes/crdkData/conventionalEEGAnalysis/LMRM/',['EYEresponse_locked_EEG_dat',sprintf('sub%03.0f',subID),'.mat']);
                %eegdat_fname = fullfile('/Volumes/LaCie/daten_fuer_laptop/',['response_locked_EEG_dat_',sprintf('sub%03.0f',subID),'.mat']);
%                 save_name_coherence = fullfile(EEGdir,['response_locked_EEG_avg_coherence.mat']);
%                 save_name_condition = fullfile(EEGdir,['response_locked_EEG_avg_condition.mat']);
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
    data{sj} = data_load.dataAppend;
    
 
    % change labels to easycap
    [easy_cap_labels] = change_electrode_labels(data{sj}.label);
    
    data{sj}.label = easy_cap_labels;
%  
%  
%     % average for each coherence level for each sj
%     data_avg =  compute_average_single_subject_level(data{sj},'coherence',rtFlag,[4 5],baseline);
%     coh_30{sj} = data_avg{1};
%     coh_40{sj} = data_avg{2};
%     coh_50{sj} = data_avg{3};
%     
%     
%     
% %     average for each condition level for each sj
%     data_avg_con =  compute_average_single_subject_level(data{sj},'condition',rtFlag,[4 5],baseline);
%     cond_1{sj} = data_avg_con{1};
%     cond_2{sj} = data_avg_con{2};
%     cond_3{sj} = data_avg_con{3};
%     cond_4{sj} = data_avg_con{4};
% %     
% %     
%     % average for each condition level for each sj for fa
%     data_avg_con =  compute_average_single_subject_level(data{sj},'false alarm',rtFlag,[4 5],baseline);
%     fa_cond_1{sj} = data_avg_con{1};
%     fa_cond_2{sj} = data_avg_con{2};
%     fa_cond_3{sj} = data_avg_con{3};
%     fa_cond_4{sj} = data_avg_con{4};
% %     
%     
%     % prepare data for permutation test for freq vs rare conditions
%     data_avg_con =  compute_average_single_subject_level(data{sj},'perm_test_rare_freq_trial',rtFlag,[-7 -6],baseline);
%     cond_tr_freq{sj} = data_avg_con{1};
%     cond_tr_rare{sj} = data_avg_con{2};
%     
%     
%     % prepare data for permutation test for freq vs rare conditions false alarms
%     data_avg_con =  compute_average_single_subject_level(data{sj},'perm_test_rare_freq_fa',rtFlag,[-2 -1],baseline);
%     fa_cond_tr_freq{sj} = data_avg_con{1};
%     fa_cond_tr_rare{sj} = data_avg_con{2};
%     
%     
%     % prepare data for permutation test for long vs short conditions
%     data_avg_con =  compute_average_single_subject_level(data{sj},'perm_test_short_long_trial',rtFlag,[-2 -1],baseline);
%     cond_tr_short{sj} = data_avg_con{1};
%     cond_tr_long{sj} = data_avg_con{2};
%     
%     
%     % prepare data for permutation test for long vs short conditions
%     % false alarms
%     data_avg_con =  compute_average_single_subject_level(data{sj},'perm_test_short_long_fa',rtFlag,[-2 -1],baseline);
%     fa_cond_tr_short{sj} = data_avg_con{1};
%     fa_cond_tr_long{sj} = data_avg_con{2};
%     
%     
%     
%     
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % conditions rare vs freq
%     
%     cfg = [];
%     cfg.channel = {'CP1', 'CPz', 'CP2'};
%     cfg.avgoverchan = 'yes';
%     cfg.nanmean = 'yes';
%     [cond_tr_rare_cp{sj}] = ft_selectdata(cfg,  cond_tr_rare{sj});
%     [cond_tr_freq_cp{sj}] = ft_selectdata(cfg, cond_tr_freq{sj});
%     
% %     cfg = [];
% %     cfg.channel = {'FC1', 'FCz', 'FC2'};
% %     cfg.avgoverchan = 'yes';
% %     cfg.nanmean = 'yes';
% %     [cond_tr_rare_fc{sj}] = ft_selectdata(cfg,  cond_tr_rare{sj});
% %     [cond_tr_freq_fc{sj}] = ft_selectdata(cfg, cond_tr_freq{sj});
% %     
% %     cfg = [];
% %     cfg.channel = {'Fp1', 'Fpz', 'Fp2'};
% %     cfg.avgoverchan = 'yes';
% %     cfg.nanmean = 'yes';
% %     [cond_tr_rare_fp{sj}] = ft_selectdata(cfg,  cond_tr_rare{sj});
% %     [cond_tr_freq_fp{sj}] = ft_selectdata(cfg, cond_tr_freq{sj});
%     
% %     cfg = [];
% %     cfg.channel = {'AF3', 'AFz', 'AF4'};
% %     cfg.avgoverchan = 'yes';
% %     cfg.nanmean = 'yes';
% %     [cond_tr_rare_af{sj}] = ft_selectdata(cfg,  cond_tr_rare{sj});
% %     [cond_tr_freq_af{sj}] = ft_selectdata(cfg, cond_tr_freq{sj});
%     
%     % false alarms
%     cfg = [];
%     cfg.channel = {'CP1', 'CPz', 'CP2'};
%     cfg.avgoverchan = 'yes';
%     cfg.nanmean = 'yes';
%     [fa_cond_tr_rare_cp{sj}] = ft_selectdata(cfg,  fa_cond_tr_rare{sj});
%     [fa_cond_tr_freq_cp{sj}] = ft_selectdata(cfg, fa_cond_tr_freq{sj});
%     
% %     cfg = [];
% %     cfg.channel = {'FC1', 'FCz', 'FC2'};
% %     cfg.avgoverchan = 'yes';
% %     cfg.nanmean = 'yes';
% %     [fa_cond_tr_rare_fc{sj}] = ft_selectdata(cfg,  fa_cond_tr_rare{sj});
% %     [fa_cond_tr_freq_fc{sj}] = ft_selectdata(cfg, fa_cond_tr_freq{sj});
%     
% %      cfg = [];
% %     cfg.channel = {'Fp1', 'Fpz', 'Fp2'};
% %     cfg.avgoverchan = 'yes';
% %     cfg.nanmean = 'yes';
% %     [fa_cond_tr_rare_fp{sj}] = ft_selectdata(cfg,  fa_cond_tr_rare{sj});
% %     [fa_cond_tr_freq_fp{sj}] = ft_selectdata(cfg, fa_cond_tr_freq{sj});
% %     
% %     cfg = [];
% %     cfg.channel = {'AF3', 'AFz', 'AF4'};
% %     cfg.avgoverchan = 'yes';
% %     cfg.nanmean = 'yes';
% %     [fa_cond_tr_rare_af{sj}] = ft_selectdata(cfg,  fa_cond_tr_rare{sj});
% %     [fa_cond_tr_freq_af{sj}] = ft_selectdata(cfg, fa_cond_tr_freq{sj});
% %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     %conditions long vs short
%     
%     cfg = [];
%     cfg.channel = {'CP1', 'CPz', 'CP2'};
%     cfg.avgoverchan = 'yes';
%     cfg.nanmean = 'yes';
%     [cond_tr_long_cp{sj}] = ft_selectdata(cfg,  cond_tr_long{sj});
%     [cond_tr_short_cp{sj}] = ft_selectdata(cfg, cond_tr_short{sj});
%     
% %     cfg = [];
% %     cfg.channel = {'FC1', 'FCz', 'FC2'};
% %     cfg.avgoverchan = 'yes';
% %     cfg.nanmean = 'yes';
% %     [cond_tr_long_fc{sj}] = ft_selectdata(cfg,  cond_tr_long{sj});
% %     [cond_tr_short_fc{sj}] = ft_selectdata(cfg, cond_tr_short{sj});
%     
%     % false alarm
%     cfg = [];
%     cfg.channel = {'CP1', 'CPz', 'CP2'};
%     cfg.avgoverchan = 'yes';
%     cfg.nanmean = 'yes';
%     [fa_cond_tr_long_cp{sj}] = ft_selectdata(cfg,  fa_cond_tr_long{sj});
%     [fa_cond_tr_short_cp{sj}] = ft_selectdata(cfg, fa_cond_tr_short{sj});
%     
% %     cfg = [];
% %     cfg.channel = {'FC1', 'FCz', 'FC2'};
% %     cfg.avgoverchan = 'yes';
% %     cfg.nanmean = 'yes';
% %     [fa_cond_tr_long_fc{sj}] = ft_selectdata(cfg,  fa_cond_tr_long{sj});
% %     [fa_cond_tr_short_fc{sj}] = ft_selectdata(cfg, fa_cond_tr_short{sj});
% %     
% %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     % coherences
% %     cfg = [];
% %     cfg.channel = {'FC1', 'FCz', 'FC2'};
% %     cfg.avgoverchan = 'yes';
% %     cfg.nanmean = 'yes';
% %     [coh_03_fc{sj}] = ft_selectdata(cfg, coh_30{sj});
% %     [coh_04_fc{sj}] = ft_selectdata(cfg, coh_40{sj});
% %     [coh_05_fc{sj}] = ft_selectdata(cfg, coh_50{sj});
% %     
%     
%     cfg = [];
%     cfg.channel = {'CP1', 'CPz', 'CP2'};
%     cfg.avgoverchan = 'yes';
%     cfg.nanmean = 'yes';
%     [coh_03_cp{sj}] = ft_selectdata(cfg,  coh_30{sj});
%     [coh_04_cp{sj}] = ft_selectdata(cfg, coh_40{sj});
%     [coh_05_cp{sj}] = ft_selectdata(cfg, coh_50{sj});
%     
%     
%     cfg = [];
%     cfg.channel = {'VEOG'};
%     cfg.avgoverchan = 'yes';
%     cfg.nanmean = 'yes';
%     [cond_tr_rare_VEOG{sj}] = ft_selectdata(cfg,  cond_tr_rare{sj});
%     [cond_tr_freq_VEOG{sj}] = ft_selectdata(cfg, cond_tr_freq{sj});
%     
%     
%     
%     cfg = [];
%     cfg.channel = {'HEOG'};
%     cfg.avgoverchan = 'yes';
%     cfg.nanmean = 'yes';
%     [cond_tr_rare_HEOG{sj}] = ft_selectdata(cfg,  cond_tr_rare{sj});
%     [cond_tr_freq_HEOG{sj}] = ft_selectdata(cfg, cond_tr_freq{sj});% 
% 
%     cfg = [];
%     cfg.channel = {'VEOG'};
%     cfg.avgoverchan = 'yes';
%     cfg.nanmean = 'yes';
%     [fa_cond_tr_rare_VEOG{sj}] = ft_selectdata(cfg,  fa_cond_tr_rare{sj});
%     [fa_cond_tr_freq_VEOG{sj}] = ft_selectdata(cfg, fa_cond_tr_freq{sj});
%     
%     
%     
%     cfg = [];
%     cfg.channel = {'HEOG'};
%     cfg.avgoverchan = 'yes';
%     cfg.nanmean = 'yes';
%     [fa_cond_tr_rare_HEOG{sj}] = ft_selectdata(cfg,  fa_cond_tr_rare{sj});
%     [fa_cond_tr_freq_HEOG{sj}] = ft_selectdata(cfg, fa_cond_tr_freq{sj});
end
keyboard;
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

cfg = [];
fa_cond_avg_button{1} = ft_timelockgrandaverage(cfg,fa_cond_1{:});
fa_cond_avg_button{1}.cfg = [];
fa_cond_avg_button{2} = ft_timelockgrandaverage(cfg,fa_cond_2{:});
fa_cond_avg_button{2}.cfg = [];
fa_cond_avg_button{3} = ft_timelockgrandaverage(cfg,fa_cond_3{:});
fa_cond_avg_button{3}.cfg = [];
fa_cond_avg_button{4} = ft_timelockgrandaverage(cfg,fa_cond_4{:});
cond_avg_button{4}.cfg = [];

keyboard;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare data for permutation plotting topoplot freq vs rare
cfg = [];
cfg.channels = 'eeg';
cond_tr_freq_avg = ft_timelockgrandaverage(cfg,cond_tr_freq{:});
% calculate se
[se] = calculate_se(cond_tr_freq_avg.var,nS);
cond_tr_freq_avg.se = se;

cond_tr_rare_avg = ft_timelockgrandaverage(cfg,cond_tr_rare{:});
% calculate se
[se] = calculate_se(cond_tr_rare_avg.var,nS);
cond_tr_rare_avg.se = se;


cond_multiplot{1} = cond_tr_freq_avg; 
cond_multiplot{2} = cond_tr_rare_avg; 

%%%%%%%%%%% false alarms
cfg = [];
% cfg.channels = 'eeg';
fa_cond_tr_freq_avg = ft_timelockgrandaverage(cfg,fa_cond_tr_freq{:});
% calculate se
[se] = calculate_se(fa_cond_tr_freq_avg.var,nS);
fa_cond_tr_freq_avg.se = se;

fa_cond_tr_rare_avg = ft_timelockgrandaverage(cfg,fa_cond_tr_rare{:});
% calculate se
[se] = calculate_se(fa_cond_tr_rare_avg.var,nS);
fa_cond_tr_rare_avg.se = se;


cfg = [];
% cfg.channels = 'eeg';
fa_cond_tr_short_avg = ft_timelockgrandaverage(cfg,fa_cond_tr_short{:});
% calculate se
[se] = calculate_se(fa_cond_tr_short_avg.var,nS);
fa_cond_tr_short_avg.se = se;

fa_cond_tr_long_avg = ft_timelockgrandaverage(cfg,fa_cond_tr_long{:});
% calculate se
[se] = calculate_se(fa_cond_tr_long_avg.var,nS);
fa_cond_tr_long_avg.se = se;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%prepare data for permutation plotting topoplot long vs short
cfg = [];
cfg.channels = 'eeg';
cond_tr_short_avg = ft_timelockgrandaverage(cfg,cond_tr_short{:});
% calculate se
[se] = calculate_se(cond_tr_short_avg.var,nS);
cond_tr_short_avg.se = se;

cond_tr_long_avg = ft_timelockgrandaverage(cfg,cond_tr_long{:});
% calculate se
[se] = calculate_se(cond_tr_long_avg.var,nS);
cond_tr_long_avg.se = se;



%%%%%%
%false alarm
cfg = [];
cfg.channels = 'eeg';
cond_tr_short_avg = ft_timelockgrandaverage(cfg,fa_cond_tr_short{:});
% calculate se
[se] = calculate_se(cond_tr_short_avg.var,nS);
cond_tr_short_avg.se = se;

cond_tr_long_avg = ft_timelockgrandaverage(cfg,fa_cond_tr_long{:});
% calculate se
[se] = calculate_se(cond_tr_long_avg.var,nS);
cond_tr_long_avg.se = se;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare data for plotting at specific electrodes


%%%%%%% freq vs rare

[cond_tr_rare_avg_cp] = ft_timelockgrandaverage(cfg, cond_tr_rare_cp{:});
% calculate se
[se] = calculate_se(cond_tr_rare_avg_cp.var,nS);
cond_tr_rare_avg_cp.se = se;


[cond_tr_freq_avg_cp] = ft_timelockgrandaverage(cfg, cond_tr_freq_cp{:});
% calculate se
[se] = calculate_se(cond_tr_freq_avg_cp.var,nS);
cond_tr_freq_avg_cp.se = se;

% [cond_tr_rare_avg_fc] = ft_timelockgrandaverage(cfg, cond_tr_rare_fc{:});
% % calculate se
% [se] = calculate_se(cond_tr_rare_avg_fc.var,nS);
% cond_tr_rare_avg_fc.se = se;


% [cond_tr_freq_avg_fc] = ft_timelockgrandaverage(cfg, cond_tr_freq_fc{:});
% % calculate se
% [se] = calculate_se(cond_tr_freq_avg_fc.var,nS);
% cond_tr_freq_avg_fc.se = se;
% 
% 
% 
% [cond_tr_rare_avg_af] = ft_timelockgrandaverage(cfg, cond_tr_rare_af{:});
% % calculate se
% [se] = calculate_se(cond_tr_rare_avg_af.var,nS);
% cond_tr_rare_avg_af.se = se;
% 
% 
% [cond_tr_freq_avg_af] = ft_timelockgrandaverage(cfg, cond_tr_freq_af{:});
% % calculate se
% [se] = calculate_se(cond_tr_freq_avg_af.var,nS);
% cond_tr_freq_avg_af.se = se;
% 
% [cond_tr_rare_avg_fp] = ft_timelockgrandaverage(cfg, cond_tr_rare_fp{:});
% % calculate se
% [se] = calculate_se(cond_tr_rare_avg_fp.var,nS);
% cond_tr_rare_avg_fp.se = se;
% 
% 
% [cond_tr_freq_avg_fp] = ft_timelockgrandaverage(cfg, cond_tr_freq_fp{:});
% % calculate se
% [se] = calculate_se(cond_tr_freq_avg_fp.var,nS);
% cond_tr_freq_avg_fp.se = se;


%%%%%%%%%% short vs long
% 
% [cond_tr_long_avg_fc] = ft_timelockgrandaverage(cfg, cond_tr_long_fc{:});
% % calculate se
% [se] = calculate_se(cond_tr_long_avg_fc.var,nS);
% cond_tr_long_avg_fc.se = se;
% 
% 
% [cond_tr_short_avg_fc] = ft_timelockgrandaverage(cfg, cond_tr_short_fc{:});
% % calculate se
% [se] = calculate_se(cond_tr_short_avg_fc.var,nS);
% cond_tr_short_avg_fc.se = se;

[cond_tr_long_avg_cp] = ft_timelockgrandaverage(cfg, cond_tr_long_cp{:});
% calculate se
[se] = calculate_se(cond_tr_long_avg_cp.var,nS);
cond_tr_long_avg_cp.se = se;


[cond_tr_short_avg_cp] = ft_timelockgrandaverage(cfg, cond_tr_short_cp{:});
% calculate se
[se] = calculate_se(cond_tr_short_avg_cp.var,nS);
cond_tr_short_avg_cp.se = se;



% %%%%%%% false alarm
% 
% % freq vs rare
% 
% [fa_cond_tr_rare_avg_fc] = ft_timelockgrandaverage(cfg, fa_cond_tr_rare_fc{:});
% % calculate se
% [se] = calculate_se(fa_cond_tr_rare_avg_fc.var,nS);
% fa_cond_tr_rare_avg_fc.se = se;
% 
% 
% [fa_cond_tr_freq_avg_fc] = ft_timelockgrandaverage(cfg, fa_cond_tr_freq_fc{:});
% % calculate se
% [se] = calculate_se(fa_cond_tr_freq_avg_fc.var,nS);
% fa_cond_tr_freq_avg_fc.se = se;
% 



[fa_cond_tr_rare_avg_cp] = ft_timelockgrandaverage(cfg, fa_cond_tr_rare_cp{:});
% calculate se
[se] = calculate_se(fa_cond_tr_rare_avg_cp.var,nS);
fa_cond_tr_rare_avg_cp.se = se;


[fa_cond_tr_freq_avg_cp] = ft_timelockgrandaverage(cfg, fa_cond_tr_freq_cp{:});
% calculate se
[se] = calculate_se(fa_cond_tr_freq_avg_cp.var,nS);
fa_cond_tr_freq_avg_cp.se = se;
% 

% [fa_cond_tr_rare_avg_af] = ft_timelockgrandaverage(cfg, fa_cond_tr_rare_af{:});
% % calculate se
% [se] = calculate_se(fa_cond_tr_rare_avg_af.var,nS);
% fa_cond_tr_rare_avg_af.se = se;
% 
% 
% [fa_cond_tr_freq_avg_af] = ft_timelockgrandaverage(cfg,fa_cond_tr_freq_af{:});
% % calculate se
% [se] = calculate_se(fa_cond_tr_freq_avg_af.var,nS);
% fa_cond_tr_freq_avg_af.se = se;
% 
% [fa_cond_tr_rare_avg_fp] = ft_timelockgrandaverage(cfg, fa_cond_tr_rare_fp{:});
% % calculate se
% [se] = calculate_se(fa_cond_tr_rare_avg_fp.var,nS);
% fa_cond_tr_rare_avg_fp.se = se;
% 
% 
% [fa_cond_tr_freq_avg_fp] = ft_timelockgrandaverage(cfg, fa_cond_tr_freq_fp{:});
% % calculate se
% [se] = calculate_se(fa_cond_tr_freq_avg_fp.var,nS);
% fa_cond_tr_freq_avg_fp.se = se;
% 


%%%%% long vs short


% [fa_cond_tr_long_avg_fc] = ft_timelockgrandaverage(cfg, fa_cond_tr_long_fc{:});
% % calculate se
% [se] = calculate_se(fa_cond_tr_long_avg_fc.var,nS);
% fa_cond_tr_long_avg_fc.se = se;
% 
% 
% [fa_cond_tr_short_avg_fc] = ft_timelockgrandaverage(cfg, fa_cond_tr_short_fc{:});
% % calculate se
% [se] = calculate_se(fa_cond_tr_short_avg_fc.var,nS);
% fa_cond_tr_short_avg_fc.se = se;


[fa_cond_tr_long_avg_cp] = ft_timelockgrandaverage(cfg, fa_cond_tr_long_cp{:});
% calculate se
[se] = calculate_se(fa_cond_tr_long_avg_cp.var,nS);
fa_cond_tr_long_avg_cp.se = se;


[fa_cond_tr_short_avg_cp] = ft_timelockgrandaverage(cfg, fa_cond_tr_short_cp{:});
% calculate se
[se] = calculate_se(fa_cond_tr_short_avg_cp.var,nS);
fa_cond_tr_short_avg_cp.se = se;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg = [];
cfg.channels = 'eeg';


[coh_30_avg_cp] = ft_timelockgrandaverage(cfg,coh_03_cp{:});
% calculate se
[se] = calculate_se(coh_30_avg_cp.var,nS);
coh_30_avg_cp.se = se;

[coh_40_avg_cp] = ft_timelockgrandaverage(cfg,coh_04_cp{:});
% calculate se
[se] = calculate_se(coh_40_avg_cp.var,nS);
coh_40_avg_cp.se = se;

[coh_50_avg_cp] = ft_timelockgrandaverage(cfg,coh_05_cp{:});
% calculate se
[se] = calculate_se(coh_50_avg_cp.var,nS);
coh_50_avg_cp.se = se;



% [coh_30_avg_fc] = ft_timelockgrandaverage(cfg,coh_03_fc{:});
% % calculate se
% [se] = calculate_se(coh_30_avg_fc.var,nS);
% coh_30_avg_fc.se = se;
% 
% [coh_40_avg_fc] = ft_timelockgrandaverage(cfg,coh_04_fc{:});
% % calculate se
% [se] = calculate_se(coh_40_avg_fc.var,nS);
% coh_40_avg_fc.se = se;
% 
% [coh_50_avg_fc] = ft_timelockgrandaverage(cfg,coh_05_fc{:});
% % calculate se
% [se] = calculate_se(coh_50_avg_fc.var,nS);
% coh_50_avg_fc.se = se;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cfg = [];
% [cond_tr_rare_avg_VEOG] = ft_timelockgrandaverage(cfg,cond_tr_rare_VEOG{:});
% % calculate se
% [se] = calculate_se(cond_tr_rare_avg_VEOG.var,nS);
% cond_tr_rare_avg_VEOG.se = se;
% 
% cfg = [];
% [cond_tr_freq_avg_VEOG] = ft_timelockgrandaverage(cfg,cond_tr_freq_VEOG{:});
% % calculate se
% [se] = calculate_se(cond_tr_freq_avg_VEOG.var,nS);
% cond_tr_freq_avg_VEOG.se = se;
% 
% cfg = [];
% [cond_tr_rare_avg_HEOG] = ft_timelockgrandaverage(cfg,cond_tr_rare_HEOG{:});
% % calculate se
% [se] = calculate_se(cond_tr_rare_avg_HEOG.var,nS);
% cond_tr_rare_avg_HEOG.se = se;
% 
% cfg = [];
% [cond_tr_freq_avg_HEOG] = ft_timelockgrandaverage(cfg,cond_tr_freq_HEOG{:});
% % calculate se
% [se] = calculate_se(cond_tr_freq_avg_HEOG.var,nS);
% cond_tr_freq_avg_HEOG.se = se;
% 
% 
% cfg = [];
% [fa_cond_tr_rare_avg_VEOG] = ft_timelockgrandaverage(cfg,fa_cond_tr_rare_VEOG{:});
% %calculate se
% [se] = calculate_se(fa_cond_tr_rare_avg_VEOG.var,nS);
% fa_cond_tr_rare_avg_VEOG.se = se;
% 
% cfg = [];
% [fa_cond_tr_freq_avg_VEOG] = ft_timelockgrandaverage(cfg,fa_cond_tr_freq_VEOG{:});
% %calculate se
% [se] = calculate_se(fa_cond_tr_freq_avg_VEOG.var,nS);
% fa_cond_tr_freq_avg_VEOG.se = se;
% 
% cfg = [];
% [fa_cond_tr_rare_avg_HEOG] = ft_timelockgrandaverage(cfg,fa_cond_tr_rare_HEOG{:});
% %calculate se
% [se] = calculate_se(fa_cond_tr_rare_avg_HEOG.var,nS);
% fa_cond_tr_rare_avg_HEOG.se = se;
% 
% cfg = [];
% [fa_cond_tr_freq_avg_HEOG] = ft_timelockgrandaverage(cfg,fa_cond_tr_freq_HEOG{:});
% %calculate se
% [se] = calculate_se(fa_cond_tr_freq_avg_HEOG.var,nS);
% fa_cond_tr_freq_avg_HEOG.se = se;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

keyboard;
% permutation test
% for trials
save_loc = fullfile(EEGdir,'perm_stat_trial_freq_rare.mat');
[stat] = permutation_test(cond_tr_freq_cp, cond_tr_rare_cp, [], 'condition', 'button press', save_loc);
[stat_length] = permutation_test(cond_tr_short, cond_tr_long, [], 'condition', 'button press', save_loc);


% permutation test false alarm
% for trials
save_loc = fullfile(EEGdir,'test_April_2020','fa_perm_stat_trial_freq_rare.mat');
[stat_fa_length] = permutation_test(fa_cond_tr_short, fa_cond_tr_long, [], 'condition', 'button press', save_loc);
[stat_fa] = permutation_test(fa_cond_tr_freq, fa_cond_tr_rare, [], 'condition', 'button press', save_loc);

keyboard;
% plot stats
fsample = 100;
pos_fa = plot_stats_permtest(stat_fa,fa_cond_tr_freq_avg, fa_cond_tr_rare_avg,fsample, 0.2, 0, 15);
pos = plot_stats_permtest(stat, cond_tr_freq_avg, cond_tr_rare_avg,fsample, 0.2, 0, 15);


pos_fa_length = plot_stats_permtest(stat_fa_length,fa_cond_tr_short_avg, fa_cond_tr_long_avg,fsample, 0.2, 0, 15);
pos_length = plot_stats_permtest(stat_length, cond_tr_short_avg, cond_tr_long_avg,fsample, 0.2, 0, 15);

plot_ERP_timeseries_stats(pos_fa_length,stat_fa_length,fa_cond_tr_short_avg_cp,fa_cond_tr_long_avg_cp,fa_cond_tr_short_avg_fc,fa_cond_tr_long_avg_fc)
plot_ERP_timeseries_stats(pos_length,stat_length,cond_tr_short_avg_cp,cond_tr_long_avg_cp, cond_tr_short_avg_fc,cond_tr_long_avg_fc)

plot_ERP_timeseries_stats(pos_fa,stat_fa,fa_cond_tr_freq_avg_cp,fa_cond_tr_rare_avg_cp,fa_cond_tr_freq_avg_fc,fa_cond_tr_rare_avg_fc)
plot_ERP_timeseries_stats(pos,stat,cond_tr_freq_avg_cp,cond_tr_rare_avg_cp, cond_tr_freq_cp,  cond_tr_rare_cp)

% save(save_name_coherence,'coh_avg_button');
% save(save_name_condition,'cond_avg_button');
%% average across subjects and conditions for each coherence level - timelocked to trial start


lim = quantile(coh_avg_button{1}.avg(:),[0.1 0.9]);

cl = cbrewer('seq','Blues',12);
cl =  cl([6 10 12],:);
minlim = -lim(2);
maxlim = lim(2);

cfg = [];
cfg.channel = {'FC1','FCz','FC2'};
cfg.layout = 'easycapM1.mat';
% cfg.ylim = [minlim maxlim];
% cfg.graphcolor = ['b','r','k'];
cfg.graphcolor = [cl];
cfg.linewidth = 3;
%%%%%%%%%%%%


ft_singleplotER(cfg,coh_avg_button{:});


% YLIM - zero in middle!!! and make legend work!!!
% it *should have zero in the middle.
% and most importantly it should be the same in any related series of plots
% (e.g., all coherence levels)
legend({'30%', '40%', '50%'},'FontSize',14)
title('Averaged ERP across Subjects and conditions for different coherence levels','FontSize',14)
xlabel('time (s) - button press at 0','FontSize',14)
set(gca,'FontSize',18)

%%


h = shadedErrorBar(coh_30_avg_cp.time,coh_30_avg_cp.avg,coh_30_avg_cp.se, 'lineprops', '-k');
h.patch.FaceColor = cl(1,:);
h.mainLine.Color = cl(1,:);

h = shadedErrorBar(coh_40_avg_cp.time,coh_40_avg_cp.avg,coh_40_avg_cp.se, 'lineprops', '-k');
h.patch.FaceColor = cl(2,:);
h.mainLine.Color = cl(2,:);

h = shadedErrorBar(coh_50_avg_cp.time,coh_50_avg_cp.avg,coh_50_avg_cp.se, 'lineprops', '-k');
h.patch.FaceColor = cl(3,:);
h.mainLine.Color = cl(3,:);

%%
cl = cbrewer('qual','Set1',3);
figure 
hold on 
h = shadedErrorBar(fa_cond_tr_rare_avg_VEOG.time,fa_cond_tr_rare_avg_VEOG.avg,fa_cond_tr_rare_avg_VEOG.se, 'lineprops', '-k');
h.patch.FaceColor = cl(1,:);
h.mainLine.Color = cl(1,:);

h = shadedErrorBar(fa_cond_tr_freq_avg_VEOG.time,fa_cond_tr_freq_avg_VEOG.avg,fa_cond_tr_freq_avg_VEOG.se, 'lineprops', '-k');
h.patch.FaceColor = cl(2,:);
h.mainLine.Color = cl(2,:);
title('VEOG averaged across participants for freq and rare conditions')
legend({'rare','freq'})
tidyfig; 
%%
% now plot topoplot
figure;
sb_idx = 1;

lim = quantile(coh_avg_button{1}.avg(:),[0.1 0.9]);
% minlim = min(min(coh_avg_button{1}.avg(:)));
% maxlim = max(max(coh_avg_button{1}.avg(:)));
minlim = -5 %lim(1);
maxlim = 5 %-lim(1);
for i = 1:3
    start_time = -3;
    for t = 1:8
        cfg = [];
        cfg.xlim = [start_time start_time + 0.5];
        start_time = start_time + 0.5;
        cfg.zlim = [minlim maxlim];
        cfg.layout = 'easycapM1.mat';
        subplot(3,8,sb_idx)
        sb_idx = sb_idx + 1;
        ft_topoplotER(cfg,coh_avg_button{i}); colorbar
    end
end

subplot(3,8,1)
title('30% coherence timelocked to button press')
subplot(3,8,2)
title('-2.5 -2')
subplot(3,8,3)
title('-2 -1.5sec')
subplot(3,8,4)
title('-1.5 -1sec')
subplot(3,8,5)
title('-1 -0.5sec')
subplot(3,8,6)
title('-0.5 0sec')
subplot(3,8,7)
title('0 0.5sec')
subplot(3,8,8)
title('0.5 1sec')

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

lim = quantile(coh_avg_button{1}.avg(:),[0.1 0.9]);

cl = cbrewer('seq','Blues',12);
cl =  cl([6 10 12],:);
minlim = -lim(2);
maxlim = lim(2);

cfg = [];
cfg.channel = {'CPz'};
cfg.layout = 'easycapM1.mat';
%cfg.ylim = [minlim maxlim];
cfg.xlim = [-2 2];
% cfg.graphcolor = ['b','r','k'];
cfg.graphcolor = cl;
cfg.linewidth = 3;
%%%%%%%%%%%%


ft_singleplotER(cfg,coh_avg_button{:});


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

lim = quantile(coh_avg_button{1}.avg(:),[0.1 0.9]);



minlim = -lim(2);
maxlim = lim(2);

% minlim = -7;
% maxlim = 7;
for i = 1:3
    
    start_time = -2;
    for t = 1:8
        cfg = [];
        cfg.xlim = [start_time start_time + 0.5];
        start_time = start_time + 0.5;
        cfg.zlim = [minlim maxlim];
        cfg.layout = 'easycapM1.mat';
        subplot(3,8,sb_idx)
        sb_idx = sb_idx + 1;
        ft_topoplotER(cfg,coh_avg_button{i}); colorbar
    end
end

subplot(3,8,1)
title('30% coherence timelocked to button press at 0')
subplot(3,8,2)
title('-1.5 -1')
subplot(3,8,3)
title('-1 -0.5sec')
subplot(3,8,4)
title('-0.5 0sec')
subplot(3,8,5)
title('0 0.5sec')
subplot(3,8,6)
title('0.5 1sec')
subplot(3,8,7)
title('1 1.5sec')
subplot(3,8,8)
title('1.5 2sec')

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

%lim = quantile(cond_multiplot{1}.avg(:),[0.1 0.9]);

cl = cbrewer('div','RdBu', 12);
cl = cl([4 1 9 12],:);

% cl = cbrewer('qual','Set1',3);
% 
% 
% minlim = -10e-05;
% maxlim = 10e-05;



cfg = [];
cfg.channel = {'CP1','CPz','CP2'};
cfg.layout = 'easycapM1.mat';
% cfg.ylim = [minlim maxlim];
%cfg.graphcolor = ['b','r','k'];
cfg.graphcolor = cl;
cfg.linewidth = 3;
%cfg.xlim = [-2 2];
%%%%%%%%%%%%

ft_singleplotER(cfg,fa_cond_avg_button{:});
%ft_multiplotER(cfg,cond_multiplot{:});
tidyfig;
legend({'frequent', 'rare'},'FontSize',14)
title('Averaged ERP across Subjects and coherences for different conditions','FontSize',14)
xlabel('time (s) - button press at 0','FontSize',14)
set(gca,'FontSize',18) 
%% plot csd transform data for several subjects to show it is not working 
figure

for sj = 1:length(cond_1)
    subplot(4,7,sj)
   
    plot(cond_1{sj}.time,cond_1{sj}.avg(40,:))
    xlim([-4 4])
    
end 
%legend(num2str(1:26))
ylabel('ERP average per subject at CPz')
xlabel('time in seconds, 0 = button press correct trials')
title('condition: frequent and short trials no csd')
tidyfig;

figure

for sj = 1:length(cond_1)
    subplot(4,7,sj)
   
    plot(cond_3{sj}.time,cond_2{sj}.avg(40,:))
    xlim([-4 4])
    
end 
%legend(num2str(1:26))
ylabel('ERP average per subject at CPz')
xlabel('time in seconds, 0 = button press correct trials')
title('condition: frequent and long trials no csd')
tidyfig;

figure
for sj = 1:length(cond_1)
    subplot(4,7,sj)

    plot(cond_3{sj}.time,cond_3{sj}.avg(40,:))
    xlim([-4 4])
    
end 
%legend(num2str(1:26))
ylabel('ERP average per subject at CPz')
xlabel('time in seconds, 0 = button press correct trials')
title('condition: rare and short trials no csd')
tidyfig;



figure
for sj = 1:length(cond_1)
  subplot(4,7,sj)
  
    plot(cond_3{sj}.time,cond_4{sj}.avg(40,:))
    xlim([-4 4])
    
end 
%legend(num2str(1:26))
ylabel('ERP average per subject at CPz')
xlabel('time in seconds, 0 = button press correct trials')
title('condition: rare and long trials no csd')
tidyfig;
%% now plot topoplot
figure;
sb_idx = 1;

lim = quantile(cond_avg_button{1}.avg(:),[0.1 0.9]);

% minlim = -5 %lim(1);
% maxlim = 5 %-lim(1);

minlim = lim(1);
maxlim = -lim(1);
for i = 1:4
    start_time = -1;
    for t = 1:8
        cfg = [];
        cfg.xlim = [start_time start_time + 0.25];
        start_time = start_time + 0.25;
        cfg.zlim = [minlim maxlim];
        cfg.layout = 'easycapM1.mat';
        subplot(4,8,sb_idx)
        sb_idx = sb_idx + 1;
        ft_topoplotER(cfg,cond_avg_button{i}); colorbar
    end
end

subplot(4,8,1)
title('ITIS INTS condition timelocked to button press at 0')
subplot(4,8,2)
title('-0.75 -0.5sec')
subplot(4,8,3)
title('-0.5 -0.25sec')
subplot(4,8,4)
title('-0.25 0sec')
subplot(4,8,5)
title('0 0.25sec')
subplot(4,8,6)
title('0.25 0.5sec')
subplot(4,8,7)
title('0.5 0.75sec')
subplot(4,8,8)
title('0.75 1sec')
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
