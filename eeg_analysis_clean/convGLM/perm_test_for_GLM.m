function [cond_tr_avg_all_cp] = perm_test_for_GLM(betas_sj,time_idx,chanlabels,EEGdir)

cl = cbrewer('div','RdBu',100);
% average frequent and rare trial conditions
for r = 3%1:length(time_idx)
    % frequent averaged
    freq_avg_sj = mean(betas_sj{r}(:,:,:,1:2),4);
    
    % short averaged
    short_avg_sj = mean(betas_sj{r}(:,:,:,[1,3]),4);
    
    % rare averaged
    rare_avg_sj = mean(betas_sj{r}(:,:,:,3:4),4);
    
    % long averaged
    long_avg_sj = mean(betas_sj{r}(:,:,:,[2,4]),4);
    
    all_avg_sj = mean(betas_sj{r}(:,:,:,:),4);
    
end

% get correct chanlabels
new_labels = change_electrode_labels(chanlabels);
% get electrode information
load('elec_field_for_GLM');

% transform betas into fieldtrip structure


for sj = 1:length(rare_avg_sj(:,1,1))
    
    data_ft_freq{sj}.time = time_idx(r).timeBins/1000; % get timebins in seconds
    data_ft_freq{sj}.label = new_labels;
    data_ft_freq{sj}.avg = squeeze(freq_avg_sj(sj,:,:));
    data_ft_freq{sj}.dimord = 'chan_time';
    data_ft_freq{sj}.elec = elecs;
    
    data_ft_rare{sj}.time = time_idx(r).timeBins/1000; % get timebins in seconds
    data_ft_rare{sj}.label = new_labels;
    data_ft_rare{sj}.avg = squeeze(rare_avg_sj(sj,:,:));
    data_ft_rare{sj}.dimord = 'chan_time';
    data_ft_rare{sj}.elec = elecs;
    
    
    data_ft_short{sj}.time = time_idx(r).timeBins/1000; % get timebins in seconds
    data_ft_short{sj}.label = new_labels;
    data_ft_short{sj}.avg = squeeze(short_avg_sj(sj,:,:));
    data_ft_short{sj}.dimord = 'chan_time';
    data_ft_short{sj}.elec = elecs;
    
    data_ft_long{sj}.time = time_idx(r).timeBins/1000; % get timebins in seconds
    data_ft_long{sj}.label = new_labels;
    data_ft_long{sj}.avg = squeeze(long_avg_sj(sj,:,:));
    data_ft_long{sj}.dimord = 'chan_time';
    data_ft_long{sj}.elec = elecs;
    
    data_ft_all{sj}.time = time_idx(r).timeBins/1000; % get timebins in seconds
    data_ft_all{sj}.label = new_labels;
    data_ft_all{sj}.avg = squeeze(all_avg_sj(sj,:,:));
    data_ft_all{sj}.dimord = 'chan_time';
    data_ft_all{sj}.elec = elecs;
    
    
    
%     cfg = [];
%     cfg.channel = {'FC1', 'FCz', 'FC2'};
%     cfg.avgoverchan = 'yes';
%     cfg.nanmean = 'yes';
%     [cond_tr_rare_fc{sj}] = ft_selectdata(cfg,  data_ft_rare{sj});
%     [cond_tr_freq_fc{sj}] = ft_selectdata(cfg, data_ft_freq{sj});
%     [cond_tr_long_fc{sj}] = ft_selectdata(cfg,  data_ft_long{sj});
%     [cond_tr_short_fc{sj}] = ft_selectdata(cfg, data_ft_short{sj});
    
    
    cfg = [];
    cfg.channel = {'CP1', 'CPz', 'CP2'};
    cfg.avgoverchan = 'yes';
    cfg.nanmean = 'yes';
    [cond_tr_rare_cp{sj}] = ft_selectdata(cfg,  data_ft_rare{sj});
    [cond_tr_freq_cp{sj}] = ft_selectdata(cfg, data_ft_freq{sj});
    [cond_tr_short_cp{sj}] = ft_selectdata(cfg,  data_ft_short{sj});
    [cond_tr_long_cp{sj}] = ft_selectdata(cfg, data_ft_long{sj});
    
    
    cfg = [];
    cfg.channel = {'CP1', 'CPz', 'CP2'};
    cfg.avgoverchan = 'yes';
    cfg.nanmean = 'yes';
   [cond_tr_all_cp{sj}] = ft_selectdata(cfg,  data_ft_all{sj}); 
   
   
%     cfg = [];
%     cfg.channels = {'Pz','P1','P2'};
%     cfg.avgoverchan = 'yes';
%     cfg.nanmean = 'yes';
%     [cond_tr_rare_p{sj}] = ft_selectdata(cfg,  data_ft_rare{sj});
%     [cond_tr_freq_p{sj}] = ft_selectdata(cfg, data_ft_freq{sj});
%     [cond_tr_long_p{sj}] = ft_selectdata(cfg,  data_ft_long{sj});
%     [cond_tr_short_p{sj}] = ft_selectdata(cfg, data_ft_short{sj});
% %     

    cfg = [];
    cfg.channels = {'Pz','P1','P2'};
    cfg.avgoverchan = 'yes';
    cfg.nanmean = 'yes';
   [cond_tr_all_p{sj}] = ft_selectdata(cfg,  data_ft_all{sj}); 
%     

%     
%     cfg = [];
%     cfg.channels = {'CPz','CP1','CP2','C1','C2','Cz'};
%     cfg.avgoverchan = 'yes';
%     cfg.nanmean = 'yes';
%     [cond_tr_rare_cpc{sj}] = ft_selectdata(cfg,  data_ft_rare{sj});
%     [cond_tr_freq_cpc{sj}] = ft_selectdata(cfg, data_ft_freq{sj});
%     [cond_tr_long_cpc{sj}] = ft_selectdata(cfg,  data_ft_long{sj});
%     [cond_tr_short_cpc{sj}] = ft_selectdata(cfg, data_ft_short{sj});
%     
%      cfg = [];
%     cfg.channels = {'C1','C2','Cz'};
%     cfg.avgoverchan = 'yes';
%     cfg.nanmean = 'yes';
%     [cond_tr_rare_c{sj}] = ft_selectdata(cfg,  data_ft_rare{sj});
%     [cond_tr_freq_c{sj}] = ft_selectdata(cfg, data_ft_freq{sj});
%     [cond_tr_long_c{sj}] = ft_selectdata(cfg,  data_ft_long{sj});
%     [cond_tr_short_c{sj}] = ft_selectdata(cfg, data_ft_short{sj});
%     
% 
%      cfg = [];
%     cfg.channels = {'C4'};
%     cfg.avgoverchan = 'yes';
%     cfg.nanmean = 'yes';
%     [cond_tr_rare_Right{sj}] = ft_selectdata(cfg,  data_ft_rare{sj});
%     [cond_tr_freq_Right{sj}] = ft_selectdata(cfg, data_ft_freq{sj});
%     [cond_tr_long_Right{sj}] = ft_selectdata(cfg,  data_ft_long{sj});
%     [cond_tr_short_Right{sj}] = ft_selectdata(cfg, data_ft_short{sj});
%     
%      cfg = [];
%     cfg.channels = {'C3'};
%     cfg.avgoverchan = 'yes';
%     cfg.nanmean = 'yes';
%     [cond_tr_rare_Left{sj}] = ft_selectdata(cfg,  data_ft_rare{sj});
%     [cond_tr_freq_Left{sj}] = ft_selectdata(cfg,  data_ft_freq{sj});
%     [cond_tr_long_Left{sj}] = ft_selectdata(cfg,  data_ft_long{sj});
%     [cond_tr_short_Left{sj}] = ft_selectdata(cfg, data_ft_short{sj});
%     
%     
%          cfg = [];
%     cfg.channels = {'P3','P4','P5','P6'};
%     cfg.avgoverchan = 'yes';
%     cfg.nanmean = 'yes';
%     [cond_tr_rare_Plat{sj}] = ft_selectdata(cfg,  data_ft_rare{sj});
%     [cond_tr_freq_Plat{sj}] = ft_selectdata(cfg, data_ft_freq{sj});
%     [cond_tr_long_Plat{sj}] = ft_selectdata(cfg,  data_ft_long{sj});
%     [cond_tr_short_Plat{sj}] = ft_selectdata(cfg, data_ft_short{sj});
    

    
    
end


% 
% cfg = [];
% %cfg.channels = {'FCZ','FC1','FC2'};
% cond_tr_freq_avg = ft_timelockgrandaverage(cfg,data_ft_freq{:});
% % calculate se
% [se] = calculate_se_convGLM(data_ft_freq);
% cond_tr_freq_avg.se = se;
% 
% 
% 
% cfg = [];
% %cfg.channels = {'FCZ','FC1','FC2'};
% cond_tr_rare_avg = ft_timelockgrandaverage(cfg,data_ft_rare{:});
% % calculate se
% [se] = calculate_se_convGLM(data_ft_rare);
% cond_tr_rare_avg.se = se;

%cfg = [];
%cfg.channels = {'FCZ','FC1','FC2'};
%cond_tr_avg{2} = ft_timelockgrandaverage(cfg,data_ft_short{:});
% calculate se
% [se] = calculate_se_convGLM(data_ft_short);
% cond_tr_avg{1}.se = se;



%cfg = [];
%cfg.channels = {'FCZ','FC1','FC2'};
%cond_tr_avg{1} = ft_timelockgrandaverage(cfg,data_ft_long{:});
% calculate se
% [se] = calculate_se_convGLM(data_ft_long);
% cond_tr_avg{2}.se = se;
cfg = [];
%cfg.channels = {'FCZ','FC1','FC2'};
cond_tr_avg_all = ft_timelockgrandaverage(cfg, data_ft_all{:});
% calculate se
[se] = calculate_se_convGLM( data_ft_all);
cond_tr_avg_all.se = se;
% 



cfg = [];
%cfg.channels = {'FCZ','FC1','FC2'};
cond_tr_avg_all_cp = ft_timelockgrandaverage(cfg,cond_tr_all_cp{:});
% calculate se
[se] = calculate_se_convGLM(cond_tr_all_cp);
cond_tr_avg_all_cp.se = se;
% 
% 

cfg = [];
%cfg.channels = {'FCZ','FC1','FC2'};
cond_tr_avg_all_p = ft_timelockgrandaverage(cfg,cond_tr_all_p{:});
% calculate se
[se] = calculate_se_convGLM(cond_tr_all_p);
cond_tr_avg_all_p.se = se;
% 
% 
% 
% cfg = [];
% cfg.channels = {'FCz','FC1','FC2'};
% cond_tr_freq_avg_fc = ft_timelockgrandaverage(cfg,cond_tr_freq_fc{:});
% % calculate se
% [se] = calculate_se_convGLM(cond_tr_freq_fc);
% cond_tr_freq_avg_fc.se = se;
% 
% 
% 
% cfg = [];
% cfg.channels = {'FCz','FC1','FC2'};
% cond_tr_rare_avg_fc = ft_timelockgrandaverage(cfg,cond_tr_rare_fc{:});
% % calculate se
% [se] = calculate_se_convGLM(cond_tr_rare_fc);
% cond_tr_rare_avg_fc.se = se;
% 
% 
% cfg = [];
% cfg.channels = {'FCz','FC1','FC2'};
% cond_tr_short_avg_fc = ft_timelockgrandaverage(cfg,cond_tr_short_fc{:});
% % calculate se
% [se] = calculate_se_convGLM(cond_tr_short_fc);
% cond_tr_short_avg_fc.se = se;
% 
% 
% 
% cfg = [];
% cfg.channels = {'FCz','FC1','FC2'};
% cond_tr_long_avg_fc = ft_timelockgrandaverage(cfg,cond_tr_long_fc{:});
% % calculate se
% [se] = calculate_se_convGLM(cond_tr_rare_fc);
% cond_tr_long_avg_fc.se = se;
% 
% 
% cfg = [];
% cfg.channels = {'CPz','CP1','CP2'};
% cond_tr_freq_avg_cp = ft_timelockgrandaverage(cfg,cond_tr_freq_cp{:});
% % calculate se
% [se] = calculate_se_convGLM(cond_tr_freq_cp);
% cond_tr_freq_avg_cp.se = se;
% 
% cfg = [];
% cfg.channels = {'CPz','CP1','CP2'};
% cond_tr_rare_avg_cp = ft_timelockgrandaverage(cfg,cond_tr_rare_cp{:});
% % calculate se
% [se] = calculate_se_convGLM(cond_tr_rare_cp);
% cond_tr_rare_avg_cp.se = se;
% 
% 
% cfg = [];
% cfg.channels = {'CPz','CP1','CP2'};
% cond_tr_short_avg_cp = ft_timelockgrandaverage(cfg,cond_tr_short_cp{:});
% % calculate se
% [se] = calculate_se_convGLM(cond_tr_short_cp);
% cond_tr_short_avg_cp.se = se;
% 
% cfg = [];
% cfg.channels = {'CPz','CP1','CP2'};
% cond_tr_long_avg_cp = ft_timelockgrandaverage(cfg,cond_tr_long_cp{:});
% % calculate se
% [se] = calculate_se_convGLM(cond_tr_long_cp);
% cond_tr_long_avg_cp.se = se;
% 
% 
% cfg = [];
% cfg.channels = {'CPz','CP1','CP2'};
% cond_tr_all_avg_cp = ft_timelockgrandaverage(cfg,cond_tr_all_cp{:});
% % calculate se
% [se] = calculate_se_convGLM(cond_tr_all_cp);
% cond_tr_all_avg_cp.se = se;
% 
% cfg = [];
% cfg.channels = {'Pz','P1','P2'};
% cond_tr_freq_avg_p = ft_timelockgrandaverage(cfg,cond_tr_freq_p{:});
% % calculate se
% [se] = calculate_se_convGLM(cond_tr_freq_p);
% cond_tr_freq_avg_p.se = se;
% 
% 
% cfg = [];
% cfg.channels = {'Pz','P1','P2'};
% cond_tr_rare_avg_p = ft_timelockgrandaverage(cfg,cond_tr_rare_p{:});
% % calculate se
% [se] = calculate_se_convGLM(cond_tr_rare_p);
% cond_tr_rare_avg_p.se = se;
% 
% 
% 
% 
% cfg = [];
% cfg.channels = {'Pz','P1','P2'};
% cond_tr_short_avg_p = ft_timelockgrandaverage(cfg,cond_tr_short_p{:});
% % calculate se
% [se] = calculate_se_convGLM(cond_tr_short_p);
% cond_tr_short_avg_p.se = se;
% 
% 
% cfg = [];
% cfg.channels = {'Pz','P1','P2'};
% cond_tr_long_avg_p = ft_timelockgrandaverage(cfg,cond_tr_long_p{:});
% % calculate se
% [se] = calculate_se_convGLM(cond_tr_long_p);
% cond_tr_long_avg_p.se = se;

% 
% % 
% % cfg = [];
% % cfg.channels = {'CPz','CP1','CP2', 'C1','C2','Cz'};
% % cond_tr_freq_avg_cpc = ft_timelockgrandaverage(cfg,cond_tr_freq_cpc{:});
% % % calculate se
% % [se] = calculate_se_convGLM(cond_tr_freq_cpc);
% % cond_tr_freq_avg_cpc.se = se;
% % 
% % cfg = [];
% % cfg.channels = {'CPz','CP1','CP2', 'C1','C2','Cz'};
% % cond_tr_rare_avg_cpc = ft_timelockgrandaverage(cfg,cond_tr_rare_cpc{:});
% % % calculate se
% % [se] = calculate_se_convGLM(cond_tr_rare_cpc);
% % cond_tr_rare_avg_cpc.se = se;
% % 
% % 
% % cfg = [];
% % cfg.channels = {'CPz','CP1','CP2', 'C1','C2','Cz'};
% % cond_tr_short_avg_cpc = ft_timelockgrandaverage(cfg,cond_tr_short_cpc{:});
% % % calculate sec
% % [se] = calculate_se_convGLM(cond_tr_short_cpc);
% % cond_tr_short_avg_cpc.se = se;
% % 
% % cfg = [];
% % cfg.channels = {'CPz','CP1','CP2', 'C1','C2','Cz'};
% % cond_tr_long_avg_cpc = ft_timelockgrandaverage(cfg,cond_tr_long_cpc{:});
% % % calculate se
% % [se] = calculate_se_convGLM(cond_tr_long_cpc);
% % cond_tr_long_avg_cpc.se = se;
% % 
% 
% cfg = [];
% cfg.channels = {'C1','C2','Cz'};
% cond_tr_freq_avg_c = ft_timelockgrandaverage(cfg,cond_tr_freq_c{:});
% % calculate se
% [se] = calculate_se_convGLM(cond_tr_freq_c);
% cond_tr_freq_avg_c.se = se;
% 
% % cfg = [];
% % cfg.channels = { 'C1','C2','Cz'};
% % cond_tr_rare_avg_c = ft_timelockgrandaverage(cfg,cond_tr_rare_c{:});
% % % calculate se
% % [se] = calculate_se_convGLM(cond_tr_rare_c);
% % cond_tr_rare_avg_c.se = se;
% % 
% % 
% % cfg = [];
% % cfg.channels = {'C1','C2','Cz'};
% % cond_tr_short_avg_c = ft_timelockgrandaverage(cfg,cond_tr_short_c{:});
% % % calculate sec
% % [se] = calculate_se_convGLM(cond_tr_short_c);
% % cond_tr_short_avg_c.se = se;
% % 
% % cfg = [];
% % cfg.channels = {'C1','C2','Cz'};
% % cond_tr_long_avg_c = ft_timelockgrandaverage(cfg,cond_tr_long_c{:});
% % % calculate se
% % [se] = calculate_se_convGLM(cond_tr_long_c);
% % cond_tr_long_avg_c.se = se;
% % % % 
% % % 
% % % 
% % % %%%%%%%%%%%%%%%%%%%
% % % 
% % % 
% % % cfg = [];
% % % %LEFT
% % % cond_tr_freq_avg_Left = ft_timelockgrandaverage(cfg,cond_tr_freq_Left{:});
% % % % calculate se
% % % [se] = calculate_se_convGLM(cond_tr_freq_Left);
% % % cond_tr_freq_avg_Left.se = se;
% % % 
% % % cfg = [];
% % % 
% % % cond_tr_rare_avg_Left = ft_timelockgrandaverage(cfg,cond_tr_rare_Left{:});
% % % % calculate se
% % % [se] = calculate_se_convGLM(cond_tr_rare_Left);
% % % cond_tr_rare_avg_Left.se = se;
% % % 
% % % 
% % % cfg = [];
% % % 
% % % cond_tr_short_avg_Left = ft_timelockgrandaverage(cfg,cond_tr_short_Left{:});
% % % % calculate sec
% % % [se] = calculate_se_convGLM(cond_tr_short_Left);
% % % cond_tr_short_avg_Left.se = se;
% % % 
% % % cfg = [];
% % % 
% % % cond_tr_long_avg_Left = ft_timelockgrandaverage(cfg,cond_tr_long_Left{:});
% % % % calculate se
% % % [se] = calculate_se_convGLM(cond_tr_long_Left);
% % % cond_tr_long_avg_Left.se = se;
% % % 
% % % %%%%%%%%%% 
% % % % RIGHT
% % % cfg = [];
% % % 
% % % cond_tr_freq_avg_Right = ft_timelockgrandaverage(cfg,cond_tr_freq_Right{:});
% % % % calculate se
% % % [se] = calculate_se_convGLM(cond_tr_freq_Right);
% % % cond_tr_freq_avg_Right.se = se;
% % % 
% % % cfg = [];
% % % 
% % % cond_tr_rare_avg_Right = ft_timelockgrandaverage(cfg,cond_tr_rare_Right{:});
% % % % calculate se
% % % [se] = calculate_se_convGLM(cond_tr_rare_Right);
% % % cond_tr_rare_avg_Right.se = se;
% % % 
% % % 
% % % cfg = [];
% % % 
% % % cond_tr_short_avg_Right = ft_timelockgrandaverage(cfg,cond_tr_short_Right{:});
% % % % calculate sec
% % % [se] = calculate_se_convGLM(cond_tr_short_Right);
% % % cond_tr_short_avg_Right.se = se;
% % % 
% % % cfg = [];
% % % 
% % % cond_tr_long_avg_Right = ft_timelockgrandaverage(cfg,cond_tr_long_Right{:});
% % % % calculate se
% % % [se] = calculate_se_convGLM(cond_tr_long_Right);
% % % cond_tr_long_avg_Right.se = se;
% % % 
% % % %%%%%%% 
% % % 
% % % cfg = [];
% % % 
% % % cond_tr_freq_avg_Plat = ft_timelockgrandaverage(cfg,cond_tr_freq_Plat{:});
% % % % calculate se
% % % [se] = calculate_se_convGLM(cond_tr_freq_Plat);
% % % cond_tr_freq_avg_Plat.se = se;
% % % 
% % % cfg = [];
% % % 
% % % cond_tr_rare_avg_Plat = ft_timelockgrandaverage(cfg,cond_tr_rare_Plat{:});
% % % % calculate se
% % % [se] = calculate_se_convGLM(cond_tr_rare_Plat);
% % % cond_tr_rare_avg_Plat.se = se;
% % % 
% % % 
% % % cfg = [];
% % % 
% % % cond_tr_short_avg_Plat = ft_timelockgrandaverage(cfg,cond_tr_short_Plat{:});
% % % % calculate sec
% % % [se] = calculate_se_convGLM(cond_tr_short_Plat);
% % % cond_tr_short_avg_Plat.se = se;
% % % 
% % % cfg = [];
% % % 
% % % cond_tr_long_avg_Plat = ft_timelockgrandaverage(cfg,cond_tr_long_Plat{:});
% % % % calculate se
% % % [se] = calculate_se_convGLM(cond_tr_long_Plat);
% % % cond_tr_long_avg_Plat.se = se;
% % 
% % 
%  % keyboard;
% 
% % [stat] = permutation_testGLM(data_ft_freq, data_ft_rare);
% % 
% % 
% % 
% % [stat_len] = permutation_testGLM(data_ft_short, data_ft_long);
% % 
% % 
% % 
 [stat] = permutation_testGLM(cond_tr_freq_cp, cond_tr_rare_cp);
% % 
keyboard;
% % [stat_len] = permutation_testGLM(cond_tr_short_cp, cond_tr_long_cp);
% % keyboard;
% % % 
% % cfg = [];
% % cfg.operation = 'subtract';
% % cfg.parameter = 'avg';
% % dat1_vs_dat2= ft_math(cfg,cond_tr_freq_avg,cond_tr_rare_avg);
% % 
% % %define parameters for plotting
% % timestep =0.2;
% % sampling_rate = 100;
% % sample_count  = length(stat.time);
% % idx_start = 1;
% % idx_end = length(stat.time);
% % j = [idx_start:timestep*sampling_rate:sample_count];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
% % t_real = [stat.time(1):timestep:stat.time(end)]; % real time for whole ERP signal
% % 
% % % start point eeg sample
% % idx_start = 1;
% % idx_end = length(stat.time);
% % m = [idx_start:timestep*sampling_rate:idx_end];  % temporal endpoints in M/EEG samples
% % % get relevant (significant) values
% pos_cluster_pvals = [stat.posclusters(:).prob];
% 
% pos_signif_clust = find(pos_cluster_pvals < stat.cfg.alpha);
% pos = ismember(stat.posclusterslabelmat, pos_signif_clust);
% 
% 
% neg_cluster_pvals = [stat.negclusters(:).prob];
% 
% neg_signif_clust = find(neg_cluster_pvals < stat.cfg.alpha);
% pos = ismember(stat.negclusterslabelmat, neg_signif_clust);
% 
% %pos(neg) = 1;
% % 
% % [i1,i2] = match_str(dat1_vs_dat2.label, stat.label);
% % % plot
% % 
% % lim = quantile(dat1_vs_dat2.avg,[0.01 0.99]);
% % figure
% % minlim = -lim(2);
% % maxlim = lim(2);
% % 
% % for k = 1:10
% %     subplot(2,5,k);
% %     cfg = [];
% %     cfg.xlim=[j(k) j(k+1)];
% %     %cfg.zlim = [minlim maxlim];
% %     
% %     pos_int = zeros(numel(dat1_vs_dat2.label),1);
% %     pos_int(i1) = any(pos(i2, m(k):m(k+1)), 2);
% %     cfg.highlight = 'on';
% %     cfg.highlightchannel = find(pos_int);
% %     %cfg.comment = 'xlim';
% %     cfg.comment = ['time = ',num2str(t_real(k)),[' - '],num2str(t_real(k+1))];
% %     cfg.commentpos = 'title';
% %     cfg.layout = 'easycapM1.mat';
% %     cfg.colormap = flip(cl);
% %     ft_topoplotER(cfg, dat1_vs_dat2);
% %     tidyfig;
% %     if k == 5
% %         colorbar;
% %     end
% % end
% 
% %%
% 
% [dat1_vs_dat2_avg] = calculate_difference_waveform(cond_tr_freq_cp, cond_tr_rare_cp);
% [dat1_vs_dat2_avgLen] = calculate_difference_waveform(cond_tr_short_cp, cond_tr_long_cp);
% 
% figure;
% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperSize', [200 100]);
% set(gcf, 'Position',  [500, 500, 800, 1290])
% 
% 
% cl = cbrewer('qual','Set1',3);
% 
% 
% subplot(2,1,1)
% hold on
% % channel = {'CPz'};
% channel = [39,40, 41]; % 22 - FCZ, 40 = CPZ
% % [t1,t2] = match_str( stat.label, channel);
% hold on
% h = shadedErrorBar(cond_tr_freq_avg_cp.time,cond_tr_freq_avg_cp.avg,dat1_vs_dat2_avg.se, 'lineprops', '-k');
% h.patch.FaceColor = cl(1,:);
% h.mainLine.Color = cl(1,:);
% 
% h = shadedErrorBar(cond_tr_rare_avg_cp.time,cond_tr_rare_avg_cp.avg, dat1_vs_dat2_avg.se, 'lineprops', '-k');
% h.patch.FaceColor = cl(2,:);
% h.mainLine.Color = cl(2,:);
% title('difference between freq and rare conditions correct trial CPZ')
% legend({'frequent','rare'})
% %time points in which either or all of the channes cp1 cp2 cpz are
% %significant
% xlim([-0.1 0.8])
% %ylim([-0.5 1])
% 
% % idx = pos>0;
% % %
% % plot(stat.time(idx),ones(sum(idx),1).*1,'kx')
% xlabel('time(s) 0 start of event')
% % hold off
% tidyfig;
% 
% 
% subplot(2,1,2)
% 
% hold on
% % channel = {'CPz'};
% channel = [21,22, 23]; % 22 - FCZ, 40 = CPZ
% % [t1,t2] = match_str( stat.label, channel);
% hold on
% h = shadedErrorBar(cond_tr_short_avg_cp.time,cond_tr_short_avg_cp.avg,dat1_vs_dat2_avgLen.se, 'lineprops', '-k');
% h.patch.FaceColor = cl(1,:);
% h.mainLine.Color = cl(1,:);
% 
% h = shadedErrorBar(cond_tr_long_avg_cp.time,cond_tr_long_avg_cp.avg, dat1_vs_dat2_avgLen.se, 'lineprops', '-k');
% h.patch.FaceColor = cl(2,:);
% h.mainLine.Color = cl(2,:);
% title('difference between freq and rare conditions correct trial FCZ')
% legend({'short','long'})
% xlim([-0.1 0.8])
% 
% %xlim([0 0.8])
% %ylim([-0.5 1])
% %time points in which either or all of the channes cp1 cp2 cpz are
% %significant
% 
% % idx = pos>0;
% % % %
% % plot(stat_len.time(idx),ones(sum(idx),1).*0.93,'kx')
% % xlabel('time(s) 0 start of event')
% % hold off
% tidyfig;
% 
% 
% %% difference wave 1
% 
% %[dat1_vs_dat2_avgLen] = calculate_difference_waveform(cond_tr_short_cp, cond_tr_long_cp);
% [dat1_vs_dat2_avgLen] = calculate_difference_waveform(cond_tr_long_cp, cond_tr_short_cp);
% figure;
% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperSize', [200 100]);
% set(gcf, 'Position',  [500, 500, 800, 1290])
% 
% 
% cl = cbrewer('qual','Set1',3);
% 
% 
% 
% subplot(2,1,1)
% hold on
% % channel = {'CPz'};
% channel = [39,40, 41]; % 22 - FCZ, 40 = CPZ
% % [t1,t2] = match_str( stat.label, channel);
% hold on
% h = shadedErrorBar(cond_tr_short_avg_cp.time,cond_tr_short_avg_cp.avg,dat1_vs_dat2_avgLen.se, 'lineprops', '-k');
% h.patch.FaceColor = cl(1,:);
% h.mainLine.Color = cl(1,:);
% 
% h = shadedErrorBar(cond_tr_long_avg_cp.time,cond_tr_long_avg_cp.avg, dat1_vs_dat2_avgLen.se, 'lineprops', '-k');
% h.patch.FaceColor = cl(2,:);
% h.mainLine.Color = cl(2,:);
% title('difference between freq and rare conditions correct trial CPZ')
% legend({'frequent','rare'})
% %time points in which either or all of the channes cp1 cp2 cpz are
% %significant
% xlim([-4 0.1])
% %ylim([-3 3])
% 
% % idx = pos>0;
% % 
% % plot(stat_len.time(idx),ones(sum(idx),1).*1,'kx')
% xlabel('time(s) 0 start of event')
% % hold off
% tidyfig;
% 
% 
% subplot(2,1,2)
% 
% hold on
% % channel = {'CPz'};
% channel = [21,22, 23]; % 22 - FCZ, 40 = CPZ
% % [t1,t2] = match_str( stat.label, channel);
% hold on
% h = shadedErrorBar(dat1_vs_dat2_avgLen.time,dat1_vs_dat2_avgLen.avg,dat1_vs_dat2_avgLen.se, 'lineprops', '-k');
% h.patch.FaceColor = cl(1,:);
% h.mainLine.Color = cl(1,:);
% 
% title('difference between freq and rare conditions correct trial FCZ')
% legend({'short','long'})
% xlim([-4 0.1])
% %ylim([-2 2])
% %time points in which either or all of the channes cp1 cp2 cpz are
% %significant
% 
% % idx = pos>0;
% % % %
% % plot(stat_len.time(idx),ones(sum(idx),1).*0.93,'kx')
% % xlabel('time(s) 0 start of event')
% % hold off
% tidyfig;
% 
% %% difference wave 2
% [dat1_vs_dat2_avg] = calculate_difference_waveform(cond_tr_rare_cp, cond_tr_freq_cp);
% figure;
% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperSize', [200 100]);
% set(gcf, 'Position',  [500, 500, 800, 1290])
% 
% 
% cl = cbrewer('qual','Set1',3);
% 
% 
% subplot(2,1,1)
% hold on
% % channel = {'CPz'};
% channel = [39,40, 41]; % 22 - FCZ, 40 = CPZ
% % [t1,t2] = match_str( stat.label, channel);
% hold on
% h = shadedErrorBar(cond_tr_freq_avg_cp.time,cond_tr_freq_avg_cp.avg,dat1_vs_dat2_avg.se, 'lineprops', '-k');
% h.patch.FaceColor = cl(1,:);
% h.mainLine.Color = cl(1,:);
% 
% h = shadedErrorBar(cond_tr_rare_avg_cp.time,cond_tr_rare_avg_cp.avg, dat1_vs_dat2_avg.se, 'lineprops', '-k');
% h.patch.FaceColor = cl(2,:);
% h.mainLine.Color = cl(2,:);
% title('difference between freq and rare conditions correct trial CPZ')
% legend({'frequent','rare'})
% %time points in which either or all of the channes cp1 cp2 cpz are
% %significant
% xlim([-4 0.1])
% %ylim([-3 3])
% %ylim([-0.5 1])
% 
% % idx = pos>0;
% % %
% % plot(stat.time(idx),ones(sum(idx),1).*1,'kx')
% xlabel('time(s) 0 start of event')
% % hold off
% tidyfig;
% 
% subplot(2,1,2)
% 
% hold on
% % channel = {'CPz'};
% channel = [21,22, 23]; % 22 - FCZ, 40 = CPZ
% % [t1,t2] = match_str( stat.label, channel);
% hold on
% h = shadedErrorBar(dat1_vs_dat2_avg.time,dat1_vs_dat2_avg.avg,dat1_vs_dat2_avg.se, 'lineprops', '-k');
% h.patch.FaceColor = cl(1,:);
% h.mainLine.Color = cl(1,:);
% 
% title('difference between freq and rare conditions correct trial FCZ')
% legend({'short','long'})
% xlim([-4 0.1])
% %ylim([-2 2])
% %time points in which either or all of the channes cp1 cp2 cpz are
% %significant
% 
% % idx = pos>0;
% % % %
% % plot(stat_len.time(idx),ones(sum(idx),1).*0.93,'kx')
% % xlabel('time(s) 0 start of event')
% % hold off
% tidyfig;
% 
% 
% %% 
% cl = cbrewer('div','RdBu',100);
% 
%  figure; 
% %      set(gcf, 'PaperUnits', 'inches');
% %     set(gcf, 'PaperSize', [200 100]);    
% %     set(gcf, 'Position',  [500, 500, 2000, 1290])
% cl = cbrewer('div','RdBu',100);
% sb_idx = 1;
% start_time = -0.5;
% %for t = 1:12
%  %   subplot(2,6,sb_idx)
%     sb_idx = sb_idx + 1;
%     cfg = [];
%     endTime = start_time + 0.5;
%     cfg.xlim = [start_time endTime];
%     cfg.highlight = 'on'; 
%     cfg.highlightchannel = {'P1','P2','Pz'};
%     cfg.highlightsymbol = '^';
%     %cfg.zlim =[-0.5 0.5];% [-(2e-4), (2e-4)];%[-0.8 0.8]; %[-(2e-4), (2e-4)]; %[-2 3]; %[-(2e-4), (2e-4)]; %[-3 3];
%     cfg.layout = 'easycapM1.mat';
%     cfg.comment = [num2str(start_time),['s to '],num2str(endTime),'s'];
%     cfg.FontName = 'Arial';
%     cfg.commentpos = 'title';
%     
%     
%     cfg.colormap = flip(cl);
%     ft_topoplotER(cfg,cond_tr_avg_all);
%     
% 
%     tidyfig;
%     start_time = start_time + 0.5;
%     h = colorbar();
% %end
% 
% %%     trial startv
% figure
%  cl = cbrewer('div','RdBu',100);
% %      set(gcf, 'PaperUnits', 'inches');
% %     set(gcf, 'PaperSize', [200 100]);    
% %     set(gcf, 'Position',  [500, 500, 2000, 1290])
% % cl = cbrewer('div','RdBu',100);
% sb_idx = 1;
% start_time = 2.75;
% % for t = 1:18
%     %subplot(3,6,sb_idx)
%     sb_idx = sb_idx + 1;
%     cfg = [];
%     endTime = start_time + 0.5;
%     cfg.xlim = [start_time endTime];
%     cfg.highlight = 'on'; 
%     cfg.highlightchannel = {'P1','P2','Pz'};
%     cfg.highlightsymbol = '^';
%     %cfg.zlim = [-(2e-4), (2e-4)]; %[-0.8 0.8];%[-(2e-4), (2e-4)]; %[-0.9 0.9]; %[-1 1]; %[-(2e-4), (2e-4)]; %[-3 3];
%     cfg.layout = 'easycapM1.mat';
%     cfg.comment = [num2str(start_time),['s to '],num2str(endTime),'s'];
%     cfg.commentpos = 'title';
%     
%     
%     cfg.colormap = flip(cl);
%     ft_topoplotER(cfg,cond_tr_avg_all);
%     
% 
%     tidyfig;
%     start_time = start_time + 0.5;
%     colorbar();
% end





end