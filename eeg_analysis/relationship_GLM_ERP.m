subj_list = [16, 18:21, 24, 26, 27, 28, 29, 31, 33, 34, 35, 39, 41, 42, 43, 50, 51, 52, 54, 55, 57, 58]; % subject 40 and 47 removed - not working, no idea why, investigate
% subj_list = [ 50, 51, 52, 54, 55, 57, 58];

[EEGdir,EEGdirdata,scriptdir,nSess,nS] = setup_EEG_session(subj_list);
source_density = 0; % flag if source density data is used
EEGpreproc = '/Volumes/LaCie/data_preproc';  % path to behav data all subjs
load_name = fullfile(EEGpreproc,'behav_data_all_subjs_all3');
load(load_name)
condition = {'ITIS INTS', 'ITIS INTL','ITIL INTS', 'ITIL INTL'};
%% load ERP data

for  sj = 1:length(subj_list)
    subID = subj_list(sj);
    clear data_load
    
    if source_density
        eegdat_fname = fullfile(EEGdir,['csd_trial_start_locked_EEG_dat_',sprintf('sub%03.0f',subID),'.mat']);
        
    else
        eegdat_fname = fullfile(EEGdir,['trial_start_locked_EEG_dat_',sprintf('sub%03.0f',subID),'.mat']);
        
    end
    
    
    
    data_load = load(eegdat_fname);
    data{sj} = data_load.data_append;
    % change labels to easycap
    [easy_cap_labels] = change_electrode_labels(data{sj}.label);
    data{sj}.label = easy_cap_labels;
    
   
    % average for each condition level for each sj
    data_avg_con =  compute_average_for_single_participant(data{sj},0, [-2 -1],0);
    cond_1{sj} = data_avg_con{1};
    cond_2{sj} = data_avg_con{2};
    cond_3{sj} = data_avg_con{3};
    cond_4{sj} = data_avg_con{4};
    
end
%% load GLM DATA for 'jumps_plus_absolute' model
for sj = 1:length(subj_list)
    subID = subj_list(sj);
    save_name = sprintf('sub%03.0f_betas_all_reg_plus_abs.mat',subID);
    
    
    if source_density %append 'csd' to save_name if working with CSD transformed data
        save_name(end-3:end+4) = '_csd.mat';
    end
    
    fname_betas = fullfile(EEGdir, 'preprocessed_EEG_dat_new',save_name);
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
                
                if (subID == 29 && s == 2) || (subID == 42 && s == 5) || (subID == 47 && s == 4)|| (subID == 55 && s == 1)|| (subID == 55 && s == 6)
                    
                else
                    
                    block_idx = find(block_Session_ID(s,:,2) == i); %which block corresponds to condition i?
                    betas_all_subjects{r}(sj,:,:,s,i) = betas{r}(:,:,s,block_idx);
                end
            end
        end
    end
end



for r = 1:length(time_idx) % loop over regressors
    betas_all_subjects_sessavg{r} = squeeze(mean(betas_all_subjects{r},4));
end

%% extract data from subjects from ERP data 
% re-organise data
for sj = 1:length(cond_1)

data_erp(sj,:,1) = conv(cond_1{sj}.avg(40,:),ones(10,1),'same'); 
data_erp(sj,:,2) = conv(cond_2{sj}.avg(40,:),ones(10,1),'same'); 
data_erp(sj,:,3) = conv(cond_3{sj}.avg(40,:),ones(10,1),'same'); 
data_erp(sj,:,4) = conv(cond_4{sj}.avg(40,:),ones(10,1),'same');
end 
time_1 = find(cond_1{1}.time == 1); 
time_5 = find(cond_1{1}.time == 5); 

for con = 1:4
[~,idx] = max(squeeze(data_erp(:,time_1:time_5,con)),[],2);

indices = [time_1 + idx-5 time_1 + idx + 5];
 
sj_mean_erp(:,con) = mean(squeeze(data_erp(:,indices(:,1):indices(:,2),con)),2);
end 
%% extract data from GLM data 
% select data from coherence jump and CPz 
for con = 1:4
data_glm(:,:,con) = squeeze( betas_all_subjects_sessavg{1}(:,40,:,1));

% smooth the data for all subjects 
for sj = 1:length(data_glm(:,1,1))
    
    smoothed_data(sj,:,con) = conv(squeeze(data_glm(sj,:,con)),ones(20,1),'same');
    
   
    
    
end 

idx_230ms = find(time_idx(1).timebins == 230);
idx_630ms = find(time_idx(1).timebins == 630);

 % select peak between 230 and 630ms and cut 100 frame around it 
[~,idx] = max(squeeze(smoothed_data(:,idx_230ms:idx_630ms,con)),[],2);

indices = [idx_230ms + idx-5 idx_230ms+idx + 5];
 
sj_mean_glm(:,con) = mean(squeeze(smoothed_data(:,indices(:,1):indices(:,2),con)),2);


end
%% calculate correlation for each condition 
cl = cbrewer('div','RdBu', 12);
cl = cl([2 4 12 10],:);

for con = 1:4
    [R,P] = corrcoef(sj_mean_glm(:,con),sj_mean_erp(:,con));
    subplot(2,2,con)
    plot(sj_mean_glm(:,con),sj_mean_erp(:,con),'x','Color',cl(con,:),'LineWidth',3,'MarkerSize',10)
    title([condition{con},' ','R= ',num2str(round(R(1,2),2)),' ','P= ', num2str(round(P(1,2),2))])
    xlabel('glm')
    ylabel('erp')
    tidyfig; 

end 
%% extract data from ERP when ramp occurs - between 1.5 and 2.5 seconds 
for sj = 1:length(cond_1)

data_erp(sj,:,1) = conv(cond_1{sj}.avg(40,:),ones(10,1),'same'); 
data_erp(sj,:,2) = conv(cond_2{sj}.avg(40,:),ones(10,1),'same'); 
data_erp(sj,:,3) = conv(cond_3{sj}.avg(40,:),ones(10,1),'same'); 
data_erp(sj,:,4) = conv(cond_4{sj}.avg(40,:),ones(10,1),'same');
end 
time_15 = find(cond_1{1}.time == 1.5); 
time_25 = find(cond_1{1}.time == 2.5); 

for con = 1:4

indices = [time_15 time_25];
 
sj_mean_erp(:,con) = mean(squeeze(data_erp(:,indices(1):indices(2),con)),2);
end 


