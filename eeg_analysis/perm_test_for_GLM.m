function [stat] = perm_test_for_GLM(betas_sj,time_idx,chanlabels,EEGdir)
keyboard;
cl = cbrewer('div','RdBu',100);
% average frequent and rare trial conditions
for r = 2%1:length(time_idx)
    % frequent averaged
    %freq_avg_sj = mean(betas_sj{r}(:,:,:,1:2),4);
    freq_avg_sj = mean(betas_sj{r}(:,:,:,[1,3]),4);
    
    % rare averaged
    %rare_avg_sj = mean(betas_sj{r}(:,:,:,3:4),4);
    rare_avg_sj = mean(betas_sj{r}(:,:,:,[2,4]),4);
end

% get correct chanlabels
new_labels = change_electrode_labels(chanlabels);
% get electrode information
load('elec_field_for_GLM');

% transform betas into fieldtrip structure

for r = 4%1:length(time_idx)
    for sj = 1:length(rare_avg_sj(:,1,1))
        
        data_ft_freq{sj}.time = time_idx(r).timebins/1000; % get timebins in seconds
        data_ft_freq{sj}.label = new_labels;
        data_ft_freq{sj}.avg = squeeze(freq_avg_sj(sj,:,:));
        data_ft_freq{sj}.dimord = 'chan_time';
        data_ft_freq{sj}.elec = elecs;
        
        data_ft_rare{sj}.time = time_idx(r).timebins/1000; % get timebins in seconds
        data_ft_rare{sj}.label = new_labels;
        data_ft_rare{sj}.avg = squeeze(rare_avg_sj(sj,:,:));
        data_ft_rare{sj}.dimord = 'chan_time';
        data_ft_rare{sj}.elec = elecs;
        
        
    cfg = [];
    cfg.channel = {'FC1', 'FCz', 'FC2'};
    cfg.avgoverchan = 'yes';
    cfg.nanmean = 'yes';
    [cond_tr_rare_fc{sj}] = ft_selectdata(cfg,  data_ft_rare{sj});
    [cond_tr_freq_fc{sj}] = ft_selectdata(cfg, data_ft_freq{sj});


    cfg = [];
    cfg.channel = {'FC1', 'FCz', 'FC2'};
    cfg.avgoverchan = 'yes';
    cfg.nanmean = 'yes';
    [cond_tr_rare_cp{sj}] = ft_selectdata(cfg,  data_ft_rare{sj});
    [cond_tr_freq_cp{sj}] = ft_selectdata(cfg, data_ft_freq{sj});
        
        
        
        
        
    end
end


cfg = [];
%cfg.channels = {'FCZ','FC1','FC2'};
cond_tr_freq_avg = ft_timelockgrandaverage(cfg,data_ft_freq{:});
% calculate se
[se] = calculate_se_convGLM(data_ft_freq);
cond_tr_freq_avg.se = se;



cfg = [];
%cfg.channels = {'FCZ','FC1','FC2'};
cond_tr_rare_avg = ft_timelockgrandaverage(cfg,data_ft_rare{:});
% calculate se
[se] = calculate_se_convGLM(data_ft_rare);
cond_tr_rare_avg.se = se;



cfg = [];
cfg.channels = {'FCZ','FC1','FC2'};
cond_tr_freq_avg_fc = ft_timelockgrandaverage(cfg,cond_tr_freq_fc{:});
% calculate se
[se] = calculate_se_convGLM(cond_tr_freq_fc);
cond_tr_freq_avg_fc.se = se;



cfg = [];
cfg.channels = {'FCZ','FC1','FC2'};
cond_tr_rare_avg_fc = ft_timelockgrandaverage(cfg,cond_tr_rare_fc{:});
% calculate se
[se] = calculate_se_convGLM(cond_tr_rare_fc);
cond_tr_rare_avg_fc.se = se;


cfg = [];
cfg.channels = {'CPZ','CP1','CP2'};
cond_tr_freq_avg_cp = ft_timelockgrandaverage(cfg,cond_tr_freq_cp{:});
% calculate se
[se] = calculate_se_convGLM(cond_tr_freq_cp);
cond_tr_freq_avg_cp.se = se;



cfg = [];
cfg.channels = {'CPZ','CP1','CP2'};
cond_tr_rare_avg_cp = ft_timelockgrandaverage(cfg,cond_tr_rare_cp{:});
% calculate se
[se] = calculate_se_convGLM(cond_tr_rare_cp);
cond_tr_rare_avg_cp.se = se;


save_loc = fullfile(EEGdir,'test_April_2020','perm_stat_trial_freq_rare_GLM.mat');
[stat] = permutation_test(data_ft_freq, data_ft_rare,[], 'condition', 'button press', save_loc);

cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
dat1_vs_dat2= ft_math(cfg,cond_tr_freq_avg,cond_tr_rare_avg);

%define parameters for plotting
timestep =0.2;
sampling_rate = 100;
sample_count  = length(stat.time);
idx_start = 1;
idx_end = length(stat.time);
j = [idx_start:timestep*sampling_rate:sample_count];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
t_real = [stat.time(1):timestep:stat.time(end)]; % real time for whole ERP signal

% start point eeg sample
idx_start = 1;
idx_end = length(stat.time);
m = [idx_start:timestep*sampling_rate:idx_end];  % temporal endpoints in M/EEG samples
% get relevant (significant) values
pos_cluster_pvals = [stat.posclusters(:).prob];

pos_signif_clust = find(pos_cluster_pvals < stat.cfg.alpha);
pos = ismember(stat.posclusterslabelmat, pos_signif_clust);


neg_cluster_pvals = [stat.negclusters(:).prob];

neg_signif_clust = find(neg_cluster_pvals < stat.cfg.alpha);
neg = ismember(stat.negclusterslabelmat, neg_signif_clust);

pos(neg) = 1;

[i1,i2] = match_str(dat1_vs_dat2.label, stat.label);
% plot

lim = quantile(dat1_vs_dat2.avg,[0.01 0.99]);

minlim = -1.5; %-lim(2);
maxlim = 1.5; %lim(2);

for k = 1:32;
    subplot(4,8,k);
    cfg = [];
    cfg.xlim=[j(k) j(k+1)];
    cfg.zlim = [minlim maxlim];
    
    pos_int = zeros(numel(dat1_vs_dat2.label),1);
    pos_int(i1) = any(pos(i2, m(k):m(k+1)), 2);
    cfg.highlight = 'on';
    cfg.highlightchannel = find(pos_int);
    %cfg.comment = 'xlim';
    cfg.comment = ['time = ',num2str(t_real(k)),[' - '],num2str(t_real(k+1))];
    cfg.commentpos = 'title';
    cfg.layout = 'easycapM1.mat';
    cfg.colormap = flip(cl);
    ft_topoplotER(cfg, dat1_vs_dat2);
    tidyfig;
    if k == 32
        colorbar;
    end
end

cl = cbrewer('div','RdBu', 12);
cl = cl([4 1 9 12],:);






figure;
subplot(2,1,1)
hold on
% channel = {'CPz'};
channel = [39,40, 41]; % 22 - FCZ, 40 = CPZ
% [t1,t2] = match_str( stat.label, channel);
hold on
h = shadedErrorBar(cond_tr_freq_avg.time,cond_tr_freq_avg_cp.avg,cond_tr_freq_avg_cp.se, 'lineprops', '-k');
h.patch.FaceColor = cl(2,:);
h.mainLine.Color = cl(2,:);

h = shadedErrorBar(cond_tr_freq_avg.time,cond_tr_rare_avg_cp.avg, cond_tr_rare_avg_cp.se, 'lineprops', '-k');
h.patch.FaceColor = cl(4,:);
h.mainLine.Color = cl(4,:);
title('difference between freq and rare conditions correct trial CPZ')
legend({'frequent','rare'})
%time points in which either or all of the channes cp1 cp2 cpz are
%significant

s = pos(channel,:);
idx = sum(s)>0;

plot(stat.time(idx),ones(sum(idx),1).*0.1*10e-3,'kx')
xlabel('time(s) 0 start of event')
hold off
tidyfig;


subplot(2,1,2)

hold on
% channel = {'CPz'};
channel = [21,22, 23]; % 22 - FCZ, 40 = CPZ
% [t1,t2] = match_str( stat.label, channel);
hold on
h = shadedErrorBar(cond_tr_freq_avg.time,cond_tr_freq_avg_fc.avg,cond_tr_freq_avg_fc.se, 'lineprops', '-k');
h.patch.FaceColor = cl(2,:);
h.mainLine.Color = cl(2,:);

h = shadedErrorBar(cond_tr_freq_avg.time,cond_tr_rare_avg_fc.avg, cond_tr_rare_avg_fc.se, 'lineprops', '-k');
h.patch.FaceColor = cl(4,:);
h.mainLine.Color = cl(4,:);
title('difference between freq and rare conditions correct trial FCZ')
legend({'frequent','rare'})
%time points in which either or all of the channes cp1 cp2 cpz are
%significant

s = pos(channel,:);
idx = sum(s)>0;

plot(stat.time(idx),ones(sum(idx),1).*0.1*10e-3,'kx')
xlabel('time(s) 0 start of event')
hold off
tidyfig;








%%%%%% plot difference
% 
% cfg = [];
% cfg.channel = 'VEOG';
% cfg.layout = 'easycapM10';
% ft_singleplotER(cfg,dat1_vs_dat2)
% xlabel('time sec 0 = button press')
% title('freq - rare fa HEOG')
% tidyfig;


keyboard;




end