function nan_perc_Matrix = check_badsamples(csd)
subj_list = [16,18:21, 24, 26, 27, 28, 29, 31, 33, 34, 35, 39,40,  41, 42, 43, 47,50, 51, 52, 54, 55, 57, 58];
condition = {'frequent and short', 'frequent and long', 'rare and short', 'rare and long'};
con_count = 1;
for sj = 1:length(subj_list)
    subID = subj_list(sj);
    if csd
        eegdat_fname = fullfile(    '/Users/maria/Documents/data/data_preproc/',['csd_response_locked_EEG_dat_',sprintf('sub%03.0f',subID),'.mat']);
    else
        eegdat_fname = fullfile(    '/Users/maria/Documents/data/data_preproc/',['response_locked_EEG_dat',sprintf('sub%03.0f',subID),'.mat']);
    end
    
    data_load = load(eegdat_fname);
    data = data_load.data_append;
 
    


    
    figure
    for con = 1:4 
    idx_trials =   find(data.trialinfo(:,9) == con);    
    num_trials = length(idx_trials); 
    
    trialMatrix = [];
    for trials = 1:num_trials
        trialMatrix(trials,:) = data.trial{idx_trials(trials)}(40,:);
    end 
    
  nan_perc_Matrix = sum(isnan(trialMatrix));
  
        
    subplot(2,2,con)
    con_count = con_count + 1;
    plot(data.time{1},nan_perc_Matrix,'-x')
    if con == 1 
        title([condition{1},' total num trials= ',num2str(num_trials),' subject: ',num2str(sj)]);
    else 
         title([condition{1},' total num trials= ',num2str(num_trials)]);
    end 
    
    end
savefig(num2str(sj))
end



end