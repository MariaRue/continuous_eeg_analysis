%% add paths
addpath(genpath('/Users/maria/Documents/Matlab/cbrewer'))
addpath('/Users/maria/Documents/MATLAB/raacampbell-shadedErrorBar-9b19a7b')
EEGpreproc = '/Volumes/LaCie 1/data_preproc';  % path to behav data all subjs

condition = {'Tr frequent TR short', 'Tr frequent Tr short','Tr rare Tr short', 'Tr rare Tr long'};
% plot each subject separately for conditions
cl = cbrewer('div','RdBu', 12);
cl = cl([2 4 12 10],:);
%% load data
load_name = fullfile(EEGpreproc,'behav_data_all_subjs_pilot');
load(load_name)

lags = 400; % frames back in time that are leading up to FA

which_responses = 'false alarms';  % calculating integration kernel for either FA button presses or button presses during trials options: 'false alarms' or 'trials', 'behviour only FA', behaviour only TR' - the last two are needed for poorly behavioural trianing sessions etc where no triggers are generated and we have to go with the all - responses matrix

nS = max(all_responses(:,11)); % number of subjects
%%
% loop through subjects and find button presses
mean_coherences = [];
for sj = 1 : 1
    
    
    
    
    combined_long = [];
    combined_short = [];
    for se = 1:6
        
        condition_ID = unique(all_responses(all_responses(:,11) == sj & all_responses(:,10) == se, 9));
        for bl = 1:3
            % select only stimstreams from all sessions that belong to specific
            % block
            
            % select all stim streams that belong to one subject
            stim_streams_sj = [];
            stim_streams_sj = stim_streams{sj,se}(:,bl);
            

            
            
            % find all triggers that lead to a button press
            switch which_responses
                
                
                case 'false alarms'
                    
                    % find triggers right and left button press (202 and 206)
                    mean_stim = mean_stim_streams{sj,se}(:,bl);
                    start_trial_idx = mean_stim(1:end-1) == 0 & mean_stim(2:end) ~= 0;
                    start_trial_frames = find(start_trial_idx)+1; % frames of trials at start
                    end_trial_idx = mean_stim(1:end-1) ~= 0 & mean_stim(2:end) == 0;
                    end_trial_frames = find(end_trial_idx); % last frames of trial
                    
                    FA_frames_responses = all_responses(all_responses(:,7)==2 & all_responses(:,13) == bl & all_responses(:,11) == sj & all_responses(:,10) == se,[3,6]);
                    
                    % in this case flex feedback = 50 frames
                    % trial length  = 100 frames - but save this with response function
                    
                    fb = 50;
                    
                    
                    
                    if condition_ID == 1
                        tr_len = 100;
                    else
                        tr_len = 300;
                        
                    end
                    % check whether FA happened during trial period - which would be super bad
                    t = zeros(length(start_trial_frames),1);
                    for tr = 1:length(start_trial_frames)
                        
                        t(tr) = any(FA_frames_responses(:,2) <= end_trial_frames(tr) & FA_frames_responses(:,2) >= start_trial_frames(tr));
                        
                    end
                    
                    if sum(t) ~= 0
                        keyboard;
                    end
                    
                    % check which FAs are within a trial period - but check that the trial actuall was 1sec long and there wasn't a trial response yet
                    FA_TR_ID = zeros(length(end_trial_frames),1);
                    coherence_TR = zeros(length(end_trial_frames),1);
                    for tr = 1:length(end_trial_frames)
                        
                        
                        length_tr = end_trial_frames(tr) - start_trial_frames(tr);
                        
                        
                        if length_tr == tr_len
                            
                            
                            ID = find(FA_frames_responses(:,2) <= end_trial_frames(tr) + fb & FA_frames_responses(:,2) >= end_trial_frames(tr));
                            
                        else
                            ID = [];
                        end
                        
                        
                        if ID
                            FA_TR_ID(tr) = ID;
                            coherence_TR(tr) = mean_stim(start_trial_frames(tr));
                            
                        end
                    end
                    
                    
                    FA_TR_ID(FA_TR_ID  == 0) = [];
                    coherence_TR(coherence_TR == 0) = [];
                    
                    FA_TR_FR = FA_frames_responses(FA_TR_ID,:);
                    FA_TR_FR = [FA_TR_FR,coherence_TR];
                    
                    % updated FA frames
                    FA_frames = FA_frames_responses;
                    FA_frames(FA_TR_ID,:) = [];
                    
                    triggers_left = FA_frames(FA_frames(:,1)==0,2);
                    triggers_right = FA_frames(FA_frames(:,1)==1,2);
                    
                    
                    
                    
                    
                case 'trials'
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
            end
            % only choose triggers that are bigger than the lags we go
            % back
            triggers_right(triggers_right(:,1)<=lags,:) = [];
            triggers_left(triggers_left(:,1)<=lags,:) = [];
            
            % loop through triggers for right and left button presses and
            % select coherences from stim_streams_bl
            matrix_right = [];
            for i = 1:length(triggers_right(:,1))
                
                matrix_right(:,i) = stim_streams_sj(triggers_right(i,1) - lags : triggers_right(i,1));
            end
            
            matrix_left = [];
            for i = 1:length(triggers_left(:,1))
                
                matrix_left(:,i) = stim_streams_sj(triggers_left(i,1) - lags : triggers_left(i,1)) .* -1 ; % to be able to combine with right button presses * -1
            end
            
            if condition_ID == 2
                combined_long = [combined_long,matrix_right, matrix_left];
            elseif condition_ID == 1
                
                combined_short = [combined_short,matrix_right, matrix_left];
                
            else
                
                keyboard;
            end
            
            
            
        end
        
    end
    
    mean_coherences(:,1,sj) = mean(combined_short,2);
    mean_coherences(:,2,sj) = mean(combined_long,2);
    sem_coherence(:,1,sj) = std(combined_short')/sqrt(size(mean_coherences,1));
    sem_coherence(:,2,sj) = std(combined_long')/sqrt(size(mean_coherences,1));
    % num_false_alarms_short(sj,condition_ID) = size(combined_short,2) + size(combined_long,2);
    
    
end



sem_across_subjects = squeeze(std(permute(mean_coherences,[3,1,2]))/sqrt(size(mean_coherences,3)));
mean_across_subjects = mean(mean_coherences,3);
