clear all;
addpath(genpath(pwd))

glmFlag = 'jumps_absolute';

options = continuous_RDK_set_options('iMac');

% subject list
subjectList = [16 18:21, 24, 26, 27, 28, 29, 31,  33, 34, 35, 42, 43, 47, 50, 51, 52, 54, 55, 57, 58]; %32 taken out
%subjectList = [62:64,66,68,70]; % vertical motion only 

csdFlag = 0; % 1 for csd transformed data
reference = 'LMRM';

for subject = 1:length(subjectList)
    subID = subjectList(subject);
    
    
    disp('subject: ')
    disp(subID)
    [details,paths] =  conrdk_subjects( subID,options,reference,csdFlag);
    
    % load betas for one subject
    betaSubject = load( paths.(reference).subjectLevelGLM.(glmFlag).saveName);
    
    for regressor = 1:length(options.subjectLevelGLM.(glmFlag).regressors)
        
        for sessionCount = 1:length(details.sessionIDs) % loop through sessions
            session = details.sessionIDs(sessionCount);
            
            for condition = 1:4
                
                betas_all_subjects{regressor}(subject,:,:,session,condition) = betaSubject.betas_subject{regressor}(:,:,session,condition);
                
                
            end
        end
    end
end
%%
for r = 1:length(options.subjectLevelGLM.(glmFlag).regressors) % loop over regressors
    betas_all_subjects_sessavg{r} = squeeze(mean(betas_all_subjects{r},4));
end

%[cond_tr_all_avg_cp1] =perm_test_for_GLM(betas_all_subjects_sessavg,options.subjectLevelGLM.(glmFlag).regressors,betaSubject.chanlabels);
%%
[stat] = permutation_testGLM(data_ft_freq, data_ft_rare);



[stat_len] = permutation_testGLM(data_ft_short, data_ft_long);



[stat] = permutation_testGLM(cond_tr_freq_cp, cond_tr_rare_cp);


[stat_len] = permutation_testGLM(cond_tr_short_p, cond_tr_long_p);
keyboard;

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
% pos_cluster_pvals = [stat.posclusters(:).prob];
% 
% pos_signif_clust = find(pos_cluster_pvals < stat.cfg.alpha);
% pos = ismember(stat.posclusterslabelmat, pos_signif_clust);
% 

neg_cluster_pvals = [stat_len.negclusters(:).prob];

neg_signif_clust = find(neg_cluster_pvals < stat_len.cfg.alpha);
pos = ismember(stat_len.negclusterslabelmat, neg_signif_clust);

%pos(neg) = 1;



%%
dm = ones(length(subjectList),1);

for r = 1:length(options.subjectLevelGLM.(glmFlag).regressors) % loop through regressors
    data = betas_all_subjects_sessavg{r};
    [cope,~,tstat] = ols(data(:,:),dm,1); %run GLM across participants
    
    avg_across_subjects{r} = reshape(cope,[size(data,2),size(data,3), size(data,4)]);
    tstat_across_subjects{r} = reshape(tstat,[size(data,2),size(data,3), size(data,4)]);
end


%% plot grand average response at CPz for all regressors

channel_ID = find(strcmp( betaSubject.chanlabels,'CPZ'));
%channel_ID = find(strcmp(chanlabels,'C4_C3_LRP'));

cl = cbrewer('div','RdBu', 12);
cl = cl([2 4 12 10],:);

for r = 1:length(options.subjectLevelGLM.(glmFlag).regressors) % loop through regressors
    figure;
    
    toplot = squeeze(avg_across_subjects{r}(channel_ID,:,:));
    smoothed_toplot = conv2(ones(10,1)/10,1,toplot,'same');
    %plot(time_idx(r).timebins,smoothed_toplot') %each line is now one subject
    for con = 1:4 % loop through conditions and plot in different colours
        plot(options.subjectLevelGLM.(glmFlag).regressors(r).timeBins,smoothed_toplot(:,con), 'LineWidth', 4,'Color',cl(con,:));
        hold on;tidyfig;xlabel('Time relative to event (ms)');
    end
    title(options.subjectLevelGLM.(glmFlag).regressors(r).name);legend(options.conditionLabels)
    tidyfig;
end

%% 
cl = cbrewer('div','RdBu', 12);
cl = cl([2 4 12 10],:);
for sj = 1:27
subplot(4,7,sj)
hold on 
for condition = 1:4
plot(options.subjectLevelGLM.(glmFlag).regressors(6).timeBins, squeeze(betas_all_subjects_sessavg{6}(sj,40,:,condition)),'Color',cl(condition,:))
title(['subject: ', num2str(sj)])
end 
hold off 
end

%% 
% ft_struct is the data structure with fields required by the ft topoplot function
cl = cbrewer('div','RdBu',100);
ft_struct.dimord = 'chan_time';
% transform label names into fieldtrip accepted names:
ft_struct.label = change_electrode_labels(betaSubject.chanlabels(1:61));

figure

fig_id = 1; % subplot index
start_times = [450 ]; %start times of plotting windows, in ms
end_times   = [ 550 ]; %end times of plotting windows, in ms

for r = 4% 1:length(time_idx) % loop through regressors
    ft_struct.avg = mean(avg_across_subjects{r},3); % mean across conditions
    ft_struct.time = options.subjectLevelGLM.(glmFlag).regressors(r).timeBins; % time window
    
    %cax_lim = max(abs(prctile(ft_struct.avg(:),[1 99])));
    
    %first subplot is displaying title of regressor
    subplot(1,length(start_times)+1,fig_id);
    text(0.5,0.5,options.subjectLevelGLM.(glmFlag).regressors(r).name,'FontSize',20);axis off
    fig_id = fig_id + 1;
    
    
    for t = 1:length(start_times) % loop through time windows
        
        
        
        cfg = [];
            cfg.highlight = 'on'; 
    cfg.highlightchannel = {'CP1','CP2','CPz'};
    cfg.highlightsymbol = '^';
        cfg.xlim = [start_times(t) end_times(t)];  % time window for which we create a topo plot
        cfg.zlim = [-cax_lim cax_lim]; %colorbar scale
        cfg.layout = 'easycapM1.mat';
        cfg.colormap = flip(cl);
        subplot(1,length(start_times)+1,fig_id);
        fig_id = fig_id + 1;
        
        ft_topoplotER(cfg,ft_struct); colorbar
        
        
    end
end


