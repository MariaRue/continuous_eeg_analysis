function [pos] = plot_stats_permtest(stat,data_avg1, data_avg2,fsample,timestep,csd,length_k)

% plot_stats_permtest(stat,data_avg1, data_avg2,fsample,se_data_rare_cp, se_data_freq_cp, se_data_rare_fc, se_data_freq_fc, csd)

cl = cbrewer('div','RdBu',100);    

    % substract conditions from each other
    cfg = [];
    cfg.operation = 'subtract';
    cfg.parameter = 'avg';
    dat1_vs_dat2= ft_math(cfg,data_avg1,data_avg2);
    
    %define parameters for plotting
    %timestep = timestep * fsample;       %(in seconds)
    sampling_rate = fsample;
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
    

   if isfield(stat.negclusters,'prob')
    neg_cluster_pvals = [stat.negclusters(:).prob];
   
    neg_signif_clust = find(neg_cluster_pvals < stat.cfg.alpha);
    neg = ismember(stat.negclusterslabelmat, neg_signif_clust);
     pos(neg) = 1; 
   end 
%    keyboard
%     figure;
%     
% 
%     % First ensure the channels to have the same order in the average and in the statistical output.
%     % This might not be the case, because ft_math might shuffle the order
%     [i1,i2] = match_str(dat1_vs_dat2.label, stat.label);
%     % plot
%    
% 
%     lim = quantile(dat1_vs_dat2.avg(:),[0.1 0.9]);
%     
%     minlim = lim(1);
%     maxlim = -lim(1);
%     for k = 1:25
%         subplot(5,5,k);
%         cfg = [];
%         cfg.xlim=[j(k) j(k+1)];
%         cfg.zlim = [minlim maxlim];
%         
%         pos_int = zeros(numel(dat1_vs_dat2.label(1:61)),1);
%         pos_int(i1) = any(pos(i2, m(k):m(k+1)), 2);
%         cfg.highlight = 'on';
%         cfg.highlightchannel = find(pos_int);
%         %cfg.comment = 'xlim';
%         cfg.comment = ['time = ',num2str(t_real(k)),[' - '],num2str(t_real(k+1))];
%         cfg.commentpos = 'title';
%         cfg.layout = 'easycapM1.mat';
%         cfg.colormap = cl;
%         ft_topoplotER(cfg, dat1_vs_dat2);
%         tidyfig;
%         if k == 40
%             colorbar;
%         end
%     end
%     
%     
%    %%  
%     
% 
% 
% 
% 
% %%% topoplot for csd trial start 
% % % % % % % % % % % % % % % % % %     [i1,i2] = match_str(dat1_vs_dat2.label, stat.label);
% % % % % % % % % % % % % % % % % %     % plot
% % % % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % % % %     lim = quantile(dat1_vs_dat2.avg,[0.01 0.99]);
% % % % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % % % % time_window_start = [1.03 2.28];
% % % % % % % % % % % % % % % % % % time_window_end = [1.14 2.4];
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % m = [4 129];
% % % % % % % % % % % % % % % % % % m_end = [15 141];
% 
% 
%     figure;
%     
% 
%     % First ensure the channels to have the same order in the average and in the statistical output.
%     % This might not be the case, because ft_math might shuffle the order
%     [i1,i2] = match_str(dat1_vs_dat2.label, stat.label);
%     % plot
%    
% 
%     lim = quantile(dat1_vs_dat2.avg(:),[0.1 0.9]);
%     
%     minlim = -lim(2);
%     maxlim = lim(2);
% minlim = -1.5e-4
% maxlim = 1.5e-4
% 
%         cfg = [];
%         cfg.xlim=[151 201];
%         cfg.zlim = [minlim maxlim];
%         
%         pos_int = zeros(numel(dat1_vs_dat2.label(1:61)),1);
%         pos_int(i1) = any(pos(i2, 151:201), 2);
%         cfg.highlight = 'on';
%         cfg.highlightchannel = find(pos_int);
%         %cfg.comment = 'xlim';
%         cfg.comment = ['time = -0.5s to 0s'];
%         cfg.commentpos = 'title';
%         cfg.layout = 'easycapM1.mat';
%         cfg.colormap = cl;
%         ft_topoplotER(cfg, dat1_vs_dat2);
%         tidyfig;
%   
% 
% 
% 
% 
% 
% % 
% %     [i1,i2] = match_str(dat1_vs_dat2.label, stat.label);
% %     % plot
% %     
% %     lim = quantile(dat1_vs_dat2.avg,[0.01 0.99]);
% %     
% % time_window_start = [-1];
% % time_window_end = [0];
% % t_real = stat.time(601:4:701);
% % start_id = 601;
% % end_id = 701; 
% % 
% % m = [start_id : 4 : end_id];
% 
% % m = [601];
% % m_end = [15 141];
% % 
% %     lim = quantile(dat1_vs_dat2.avg,[0.01 0.99]);
% %     
% %     minlim = lim(1);
% %     maxlim = -lim(1);
% %     for k = 1:26
% %         subplot(5,6,k);
% %         cfg = [];
% %         %cfg.xlim=[time_window_start(k) time_window_end(k)];
% %         cfg.zlim = [minlim maxlim];
% %         
% %         pos_int = zeros(numel(dat1_vs_dat2.label),1);
% %         pos_int(i1) = all(pos(i2, m(k):m(k+1)), 2);
% %         cfg.highlight = 'on';
% %         cfg.highlightchannel = find(pos_int);
% %         %cfg.comment = 'xlim';
% %         cfg.comment = ['time = ',num2str(t_real(k)),[' - '],num2str(t_real(k+1))];
% %         cfg.commentpos = 'title';
% %         cfg.layout = 'easycapM1.mat';
% %         ft_topoplotER(cfg, dat1_vs_dat2);
% %         
% %         if k == 2
% %             colorbar;
% %         end
% %     end
%     
%    
% 
%      
%      figure
%      
%      cfg = [];
% cfg.channel = 'CPz'; 
% cfg.layout = 'easycapM10'; 
% ft_singleplotER(cfg,dat1_vs_dat2)
% xlabel('time sec 0 = button press')
% title('freq - rare correct trials FCz')
% tidyfig;


end