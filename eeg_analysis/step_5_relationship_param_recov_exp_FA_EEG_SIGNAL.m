% plotting the results obtained with step 3 - continuous GLM and
% behavioural analysis obtained with calculate_mean_coh_leading_toFA.m

subj_list = [16,18:21,24,26,32, 52, 41, 34, 35, 51, 40, 28]; %2 missing
%subj_list = [62:64,66,68,70];
%subj_list(1) = []; %problem with CSD transformation in this subject?

[EEGdir,EEGdirdata,scriptdir,nSess,nS] = setup_EEG_session(subj_list);
design_matrix_type = 'jumps_plus_absolute'; % which design matrix was used in step 3?


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
            save_name = sprintf('sub%03.0f_betas_all_reg_plus_abs.mat',subID);
        case 'jumps_plus_absolute_CSD'
            save_name = sprintf('sub%03.0f_betas_all_reg_plus_abs_CSD.mat',subID);
        case 'LRP_signed'
            save_name = sprintf('sub%03.0f_betas_LRP_signed.mat',subID);
        case 'LRP_signed_wo_PE'
            save_name = sprintf('sub%03.0f_betas_LRP_signed_wo_PE.mat',subID);
        case 'LRP_signed_wo_stim_stream'
            save_name = sprintf('sub%03.0f_betas_LRP_signed_wo_stim_stream.mat',subID);
        case 'LRP_signed_wo_PE_wo_stim'
            save_name = sprintf('sub%03.0f_betas_LRP_signed_wo_PE_wo_stim.mat',subID);
        case 'split_PE'
            save_name = sprintf('sub%03.0f_betas_split_PE.mat',subID);
        case 'jumps_plus_absolute_vertical'
            save_name = sprintf('sub%03.0f_betas_all_reg_plus_abs_plus_vertical.mat',subID);        
    end
    fname_betas = fullfile(EEGdir, 'preprocessed_EEG_dat',save_name);
    load(fname_betas,'betas','time_idx','chanlabels')
    
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
                
                if subID == 70 && s == 4
                    % something is weird about this session
                else
                 
                block_idx = find(block_Session_ID(s,:,2) == i); %which block corresponds to condition i?
                betas_all_subjects{r}(sj,:,:,s,i) = betas{r}(:,:,s,block_idx);
                end

            end
        end
    end
end

%% %% first, average betas_all_subjects across all 6 sessions

for r = 1:length(time_idx) % loop over regressors
    betas_all_subjects_sessavg{r} = squeeze(mean(betas_all_subjects{r},4));
end

%% prepare the behavioural param recov results 

tau = squeeze(pnew_sj(:,2,:));
tau = 1./tau; 
median_tau = median(tau(:)); 

idx_high = tau >= median_tau; 
idx_low  =~ idx_high; 

%% extract beta weights from regrossors 

for r = 1:length(time_idx) % loop through regressors
    
    for sj = 1:nS
    
        high_tau_all{r}(sj,:,:) = mean(squeeze(betas_all_subjects_sessavg{r}(sj,:,:,idx_high(:,sj))),3);
        
        low_tau_all{r}(sj,:,:) = mean(squeeze(betas_all_subjects_sessavg{r}(sj,:,:,idx_low(:,sj))),3);
        
    end 
    
end 
%% plot low-high taus across conditions for different subjects for prediction error and abs 

subplot(1,2,1)
plot(time_idx(3).timebins, squeeze(high_tau_all{3}(:,40,:)),'b')
hold on
plot(time_idx(3).timebins, nanmean(squeeze(high_tau_all{3}(:,40,:))),'b','LineWidth',4);


plot(time_idx(3).timebins, squeeze(low_tau_all{3}(:,40,:)),'k')

plot(time_idx(3).timebins, nanmean(squeeze(low_tau_all{3}(:,40,:))),'k','LineWidth',4);

title('prediction error')
tidyfig

subplot(1,2,2)
plot(time_idx(7).timebins, squeeze(high_tau_all{7}(:,40,:)),'b')
hold on
plot(time_idx(7).timebins, nanmean(squeeze(high_tau_all{7}(:,40,:))),'b','LineWidth',4);


plot(time_idx(7).timebins, squeeze(low_tau_all{7}(:,40,:)),'k')

plot(time_idx(7).timebins, nanmean(squeeze(low_tau_all{7}(:,40,:))),'k','LineWidth',4);
tidyfig
title('abs stim')
%% for each condition
for r = 1:length(time_idx)
    for con = 1:4
        reshp = permute(betas_all_subjects_sessavg{r},[2 3 4 1]);
        high_tau{r,con}(:,:) = squeeze(reshp(40,:,con,idx_high(con,:)));
        
        low_tau{r,con}(:,:) = squeeze(reshp(40,:,con,idx_low(con,:)));
    end
    
end
%% plot prediction error for low and high tau
condition = {'ITIS INTS', 'ITIS INTL','ITIL INTS', 'ITIL INTL'};
figure
for con = 1:4
   
    subplot(2,2,con)
     hold on
    plot(time_idx(3).timebins, high_tau{3,con},'b')
    plot(time_idx(3).timebins,low_tau{3,con},'k')
    plot(time_idx(3).timebins,mean(high_tau{3,con},2),'Color','b','LineWidth',4)
    plot(time_idx(3).timebins,mean(low_tau{3,con},2),'Color','k','LineWidth',4)
    hold off
    tidyfig
    title(condition{con})
    
end
%% same as above for abs stim 
condition = {'ITIS INTS', 'ITIS INTL','ITIL INTS', 'ITIL INTL'};
figure
for con = 1:4
   
    subplot(2,2,con)
     hold on
    plot(time_idx(7).timebins, high_tau{7,con},'b')
    plot(time_idx(7).timebins,low_tau{7,con},'k')
    plot(time_idx(7).timebins,mean(high_tau{7,con},2),'Color','b','LineWidth',4)
    plot(time_idx(7).timebins,mean(low_tau{7,con},2),'Color','k','LineWidth',4)
    hold off
    tidyfig
    title(condition{con})
    
end

%% 
cl = cbrewer('div','RdBu', 12);
cl = cl([2 4 12 10],:);


figure
subplot(1,2,1)
for con = 1:4
    hold on
    plot(time_idx(7).timebins,mean(high_tau{3,con},2),'Color',cl(con,:),'LineWidth',4)
    
end 
hold off 


subplot(1,2,2)
for con = 1:4
    hold on
    plot(time_idx(7).timebins,mean(low_tau{3,con},2),'Color',cl(con,:),'LineWidth',4)
    
end 
hold off 