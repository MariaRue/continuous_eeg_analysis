% plotting the results obtained with step 3 - continuous GLM
clear all

subj_list = [ 20,21, 24, 26,  28,  34, 35,   51, 52];
%subj_list = [16,18:19,20,21,24];

%[7 8 9 11 12]




[EEGdir,EEGdirdata,scriptdir,nSess,nS] = setup_EEG_session(subj_list);

design_matrix_type ='jumps_plus_absolute'; % which design matrix was used in step 3?
flag_32_weird = 0; % flag for one weird subjectin LRP signed ;
source_density = 0;
reference_type = 'LM_RM';
%% creating big matrix including all subjects - betas_all_subjects sorted by conditions rather than block order
for sj = 1:length(subj_list)
    subID = subj_list(sj);
    
    %load in output of step 3, GLM:
    switch design_matrix_type
        case 'absolute_coherences'
            save_name = sprintf('sub%03.0f_betas_abs_coh.mat',subID);
        case 'jumps_and_jump_PEs'
            save_name = sprintf('sub%03.0f_betas_all_reg.mat',subID);
        case 'jumps_plus_absolute'
            save_name = sprintf('sub%03.0f_betas_all_reg_plus_abs.mat',subID)%'sub%03.0f_betas_all_reg_plus_abs_NEW_IN3.mat',subID);
        case 'jumps_plus_absolute_CSD'
            save_name = sprintf('sub%03.0f_betas_all_reg_plus_abs_CSD.mat',subID);
        case 'LRP_signed'
            flag_32_weird = 1;
            save_name = sprintf('sub%03.0f_betas_LRP_signed.mat',subID);
        case 'LRP_signed_wo_PE'
            save_name = sprintf('sub%03.0f_betas_LRP_signed_wo_PE.mat',subID);
        case 'LRP_signed_wo_stim_stream'
            save_name = sprintf('sub%03.0f_betas_LRP_signed_wo_stim_stream.mat',subID);
        case 'LRP_signed_wo_PE_wo_stim'
            save_name = sprintf('sub%03.0f_betas_LRP_signed_wo_PE_wo_stim.mat',subID);
        case 'split_PE'
            save_name = sprintf('sub%03.0f_betas_split_PE_IN.mat',subID);
        case 'jumps_plus_absolute_vertical'
            save_name = sprintf('sub%03.0f_betas_all_reg_plus_abs_plus_vertical.mat',subID);
        case 'all_regs_demeaned'
            save_name = sprintf('sub%03.0f_betas_all_reg_plus_abs_demean.mat',subID);
        case 'trial_start_response'
            save_name = sprintf('sub%03.0f_betas_trial_start_response_N.mat',subID);
        case 'trial_start_response_coherence'
            save_name = sprintf('sub%03.0f_betas_trial_start_response_coherence.mat',subID);
        case 'trial_start_response_coherence_in_one_reg'
            save_name = sprintf('sub%03.0f_betas_trial_start_response_coherence_in_one_reg_N.mat',subID);
        case 'trial_start_button_response_left_right'
            save_name = sprintf('sub%03.0f_betas_trial_start_button_response_left_right_N.mat',subID);
        case 'trial_start_button_response_left_minus_right'
            save_name = sprintf('sub%03.0f_betas_trial_start_button_response_left_minus_right_N.mat',subID);
        case 'trial_start_button_response_left_minus_right_trial'
            save_name = sprintf('sub%03.0f_betas_trial_start_button_response_left_minus_right_trial_N.mat',subID);
        case 'trial_start_button_response_left_minus_right_plus_mean_response'
            save_name = sprintf('sub%03.0f_betas_trial_start_button_response_left_minus_right_plus_mean_response.mat',subID);
        case 'trial_start_button_response_left_minus_right_plus_mean_response_coherence'
            save_name = sprintf('sub%03.0f_betas_trial_start_button_response_left_minus_right_plus_mean_coherence_response.mat',subID);
    end
    
    if source_density %append 'csd' to save_name if working with CSD transformed data
        save_name(end-3:end+4) = '_csd.mat';
    end
    
    switch reference_type
        case 'LM_RM'
            %fname_betas = fullfile(EEGdir, 'preprocessed_EEG_dat_new',save_name);
            fname_betas = fullfile(EEGdir, 'preprocessed_EEG_dat',save_name);
            load(fname_betas,'betas','time_idx','chanlabels')
        case 'average_electrodes'
            fname_betas = fullfile(EEGdir, 'preprocessed_EEG_dat_new','average_reference',save_name);
            load(fname_betas,'betas','time_idx','chanlabels')
    end
   
    %load in output of step 2, telling us which block corresponds to which condition:
    fname_BS = fullfile(EEGdir,'preprocessed_EEG_dat',[sprintf('sub%03.0f',subID),'_block_sess_id']);
    load(fname_BS,'block_Session_ID');
    condition_labels{1} = 'ITIS_INTS';
    condition_labels{2} = 'ITIS_INTL';
    condition_labels{3} = 'ITIL_INTS';
    condition_labels{4} = 'ITIL_INTL';
    
    for r = 1:length(time_idx) % loop through regressors
        for i = 1:4 % loop over conditions
            
            for s = 1:6 % loop through sessions
                
                if (subID == 70 && s == 4) || (subID == 29 && s == 2) || (subID == 42 && s == 5) || (subID == 47 && s == 4)|| (subID == 55 && s == 1)|| (subID == 55 && s == 6)|| ((subID == 32 || subID == 33 || subID == 16 || subID == 24 || subID == 25) && flag_32_weird == 1)
                    % something is weird about this session
                    
                    
                else
                    
                    block_idx = find(block_Session_ID(s,:,2) == i); %which block corresponds to condition i?
                    betas_all_subjects{r}(sj,:,:,s,i) = betas{r}(:,:,s,block_idx);
                end
                
            end
        end
    end
end

%% first, average betas_all_subjects across all 6 sessions

for r = 1:length(time_idx) % loop over regressors
    betas_all_subjects_sessavg{r} = squeeze(mean(betas_all_subjects{r},4));
end


% plot PE for each subject subplot 
figure
for sj = 1:11
    subplot(3,3,sj)
    data = squeeze(betas_all_subjects_sessavg{3}(sj,40,:,1));
     smoothed_toplot = conv2(ones(10,1)/10,1,data,'same');
    plot(time_idx(3).timebins, smoothed_toplot)
    title(num2str(subj_list(sj)))
    tidyfig;
end 


%[stat] = perm_test_for_GLM(betas_all_subjects_sessavg,time_idx,chanlabels,EEGdir);

%% run a group GLM across all participants

% make a design matrix across participants - just calculates the mean
dm = ones(nS,1);

for r = 1:length(time_idx) % loop through regressors
    data = betas_all_subjects_sessavg{r};
    [cope,~,tstat] = ols(data(:,:),dm,1); %run GLM across participants
    
    avg_across_subjects{r} = reshape(cope,[size(data,2),size(data,3), size(data,4)]);
    tstat_across_subjects{r} = reshape(tstat,[size(data,2),size(data,3), size(data,4)]);
end

%% plot grand average response at CPz for all regressors

channel_ID = find(strcmp(chanlabels,'CPZ'));
%channel_ID = find(strcmp(chanlabels,'C4_C3_LRP'));

cl = cbrewer('div','RdBu', 12);
cl = cl([2 4 12 10],:);

for r = 1:length(time_idx) % loop through regressors
    figure;
    
    toplot = squeeze(avg_across_subjects{r}(channel_ID,:,:));
    smoothed_toplot = conv2(ones(10,1)/10,1,toplot,'same');
    %plot(time_idx(r).timebins,smoothed_toplot') %each line is now one subject
    for con = 1:4 % loop through conditions and plot in different colours
        plot(time_idx(r).timebins,smoothed_toplot(:,con), 'LineWidth', 4,'Color',cl(con,:));
        hold on;tidyfig;xlabel('Time relative to event (ms)');
    end
    title(time_idx(r).name);legend(condition_labels)
    tidyfig;
end

%% if trial responses with coherence plot regressors for different coherences together for each condition
channel_ID = find(strcmp(chanlabels,'FCZ'));
%channel_ID = find(strcmp(chanlabels,'C4_C3_LRP'));
coherence_labels = {'30', '40', '50'};
cl = cbrewer('seq','Blues', 9);
cl = cl([3 6 9],:);
col = 1;
for r = 1:length(time_idx) % loop through regressors
    
    
    toplot = squeeze(avg_across_subjects{r}(channel_ID,:,:));
    smoothed_toplot = conv2(ones(10,1)/10,1,toplot,'same');
    %plot(time_idx(r).timebins,smoothed_toplot') %each line is now one subject
    for con = 1:4 % loop through conditions and plot in different colours
        if r <=3
            figure(con)
        else
            figure(con + 4)
        end
        hold on
        plot(time_idx(r).timebins,smoothed_toplot(:,con), 'LineWidth', 4,'Color',cl(col,:));
        hold off;tidyfig;xlabel('Time relative to event (ms)');
        title(condition_labels(con));legend(coherence_labels)
    end
    
    tidyfig;
    
    col = col + 1;
    
    if col == 4
        col = 1;
    end
end



%% plot topoplots for each regressor

% ft_struct is the data structure with fields required by the ft topoplot function

ft_struct.dimord = 'chan_time';
% transform label names into fieldtrip accepted names:
ft_struct.label = change_electrode_labels(chanlabels(1:61));

figure

fig_id = 1; % subplot index
start_times = [100 260 450]; %start times of plotting windows, in ms
end_times   = [240 350 520]; %end times of plotting windows, in ms

for r = 4% 1:length(time_idx) % loop through regressors
    ft_struct.avg = mean(avg_across_subjects{r},3); % mean across conditions
    ft_struct.time = time_idx(r).timebins; % time window
    
    cax_lim = max(abs(prctile(ft_struct.avg(:),[5 95])));
    
    %first subplot is displaying title of regressor
    subplot(1,length(start_times)+1,fig_id);
    text(0.5,0.5,time_idx(r).name,'FontSize',20);axis off
    fig_id = fig_id + 1;
    
    
    for t = 1:length(start_times) % loop through time windows
        
        
        
        cfg = [];
        cfg.xlim = [start_times(t) end_times(t)];  % time window for which we create a topo plot
        cfg.zlim = [-cax_lim cax_lim]; %colorbar scale
        cfg.layout = 'easycapM1.mat';
        
        subplot(1,length(start_times)+1,fig_id);
        fig_id = fig_id + 1;
        
        ft_topoplotER(cfg,ft_struct); colorbar
        
        
    end
end

%% plot topoplots for each regressor and condition

% ft_struct is the data structure with fields required by the ft topoplot function

ft_struct.dimord = 'chan_time';
% transform label names into fieldtrip accepted names:
ft_struct.label = change_electrode_labels(chanlabels(1:61));

figure

fig_id = 1; % subplot index
start_times = [-250 ]; %start times of plotting windows, in ms
end_times   = [-40  ]; %end times of plotting windows, in ms

for r = 3% 1:length(time_idx) % loop through regressors
    
    ft_struct.avg = mean(avg_across_subjects{r},3); % mean across conditions
    ft_struct.time = time_idx(r).timebins; % time window
    
    cax_lim = max(abs(prctile(ft_struct.avg(:),[5 95])));
    
    %first subplot is displaying title of regressor
    %     subplot(1,length(start_times)+1,fig_id);
    %     text(0.5,0.5,time_idx(r).name,'FontSize',20);axis off
    %     fig_id = fig_id + 1;
    
    for con = 1:4
        
        ft_struct.avg = avg_across_subjects{r}(:,:,con); % mean across conditions
        subplot(2,2,con)
        if con == 1
            title(time_idx(r).name,'FontSize',20);
        else
            
            title(condition_labels{con})
            
        end
        
        for t = 1:length(start_times) % loop through time windows
            
            
            
            cfg = [];
            cfg.xlim = [start_times(t) end_times(t)];  % time window for which we create a topo plot
            cfg.zlim = [-cax_lim cax_lim]; %colorbar scale
            cfg.layout = 'easycapM1.mat';
            %
            %         subplot(1,length(start_times)+1,fig_id);
            %         fig_id = fig_id + 1;
            
            ft_topoplotER(cfg,ft_struct); colorbar
            
            
        end
    end
end

%% multiplotER - looking at waveform of abs coherence at every single electrode

% ft_struct is the data structure with fields required by the ft topoplot function

ft_struct.dimord = 'chan_time';
% transform label names into fieldtrip accepted names:
ft_struct.label = change_electrode_labels(chanlabels(1:61));



for r = 1:length(time_idx)
    ft_struct.avg = mean(avg_across_subjects{r},3); % mean across conditions
    ft_struct.time = time_idx(r).timebins; % time window
    
    
    cfg = [];
    cfg.layout = 'easycapM1.mat';
    
    figure
    ft_multiplotER(cfg,ft_struct);
    
end

%% consistency across subjects?

r = 1; %regressor of interest
toplot = squeeze(mean(betas_all_subjects_sessavg{r}(:,40,:,:),4));
smoothed_toplot = conv2(1,ones(10,1)/10,toplot,'same');
plot(time_idx(r).timebins,smoothed_toplot') %each line is now one subject

%% difference waveform for different conditions

% difference waveform for ITI short - ITI long

for r = 1:length(time_idx)
    
    ITI_difference_waveform{r} = mean(betas_all_subjects_sessavg{r}(:,:,:,[1 2]),4) - mean(betas_all_subjects_sessavg{r}(:,:,:,[3 4]),4);
    
end
r = 3;
toplot = squeeze(ITI_difference_waveform{r}(:,40,:,:));
smoothed_toplot = conv2(1,ones(10,1)/10,toplot,'same');
plot(time_idx(r).timebins,smoothed_toplot') %each line is now one subject
plotmse(smoothed_toplot,1,time_idx(r).timebins)
%% topoplots for difference waveform

% ft_struct is the data structure with fields required by the ft topoplot function

ft_struct.dimord = 'chan_time';
% transform label names into fieldtrip accepted names:
ft_struct.label = change_electrode_labels(chanlabels(1:61));

for r = 1:length(time_idx)
    ft_struct.avg = squeeze(mean(ITI_difference_waveform{r},1)); % mean across conditions
    ft_struct.time = time_idx(r).timebins; % time window
    
    
    cfg = [];
    cfg.layout = 'easycapM1.mat';
    cfg.comment = time_idx(r).name;
    
    figure
    ft_multiplotER(cfg,ft_struct);
    
end