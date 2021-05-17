function step_3_continuous_GLM_diff_intertrial_trial_trial(design_matrix_type,source_density,reference_type,ITI_flag)
% this function runs a conv GLM. Regressors used are different for
% different design matrices. GLM can be run on either csd transformed data
% nor normal EEG data.


% Input
% design_matrix_type - string defining model
% source_density - flag, 1 if csd transformed data is used
% reference_type - during preprocessing step data was either referenced to
%                  left and right mastoids or average of all electrodes
%                  (LM_RM or average_electrodes)
% ITI_flag - flag, 1 to exclude events happening during trial periods



% subjects to evaluate
subj_list = [16, 18:21, 24, 26, 27, 28, 29, 31,  33, 34, 35, 39, 40, 41, 42, 43, 47, 50, 51, 52, 54, 55, 57, 58]; %32 taken out
% subj_list = [62:64,66,68,70]; % for vertical motion experiment



[EEGdir,EEGdirdata,scriptdir,nSess,nS] = setup_EEG_session(subj_list); % set up fieldtrip and some paths

% load behavioural data for button press and trial start regressors
EEGpreproc = '/Volumes/LaCie/data_preproc';  % path to behav data all subjs
load_name = fullfile('/Volumes/LaCie/data_preproc/behav_data_all_subjs_all3'); % load behavioural data
load(load_name)



% loop through subjects
for sj = 1:length(subj_list)
    
    
    subID = subj_list(sj);
    disp('subject:')
    disp(subID);
    
    if subID == 55 % comment here
        nSess = 5;
    else
        nSess = 6;
        
    end
    
    %% load in aligned data depending on reference type
    
    switch reference_type %
        case 'LM_RM'
            if source_density %%% more comments here for each line
                eegdat_fname = fullfile(EEGdir,'preprocessed_EEG_dat_new','eeg_matched_data',[sprintf('sub%03.0f',subID),'_csdEEGdat.mat']);
                load(fullfile(EEGdir,'preprocessed_EEG_dat_new','eeg_matched_data',[sprintf('sub%03.0f',subID),'_block_sess_id_csd'])); % load into structure not into workspace
                load(fullfile(EEGdir,'preprocessed_EEG_dat_new','eeg_matched_data',[sprintf('sub%03.0f',subID),'_trigger_vals_eegmatch_csd']));
            else
                keyboard;
                eegdat_fname = fullfile(EEGdir,'preprocessed_EEG_dat_new','eeg_matched_data',[sprintf('sub%03.0f',subID),'_EEGdat.mat']);
                load(fullfile(EEGdir,'preprocessed_EEG_dat_new','eeg_matched_data',[sprintf('sub%03.0f',subID),'_block_sess_id']));
                load(fullfile(EEGdir,'preprocessed_EEG_dat_new','eeg_matched_data',[sprintf('sub%03.0f',subID),'_trigger_vals_eegmatch']));
            end
            load(eegdat_fname);
            
            
            
            
            EEGdatadir =  fullfile(EEGdir,sprintf('sub%03.0f',subID),'eeg');
            STdatadir = fullfile(EEGdir,sprintf('sub%03.0f',subID),'stim');
            BHVdatadir = fullfile(EEGdir,sprintf('sub%03.0f',subID),'behaviour');
        case 'average_electrodes'
            
            if source_density
                eegdat_fname = fullfile(EEGdir,'preprocessed_EEG_dat','average_reference',[sprintf('sub%03.0f',subID),'_csdEEGdat.mat']);
                load(fullfile(EEGdir,'preprocessed_EEG_dat','average_reference',[sprintf('sub%03.0f',subID),'_block_sess_id_csd']));
                load(fullfile(EEGdir,'preprocessed_EEG_dat','average_reference',[sprintf('sub%03.0f',subID),'_trigger_vals_eegmatch_csd']));
            else
                eegdat_fname = fullfile(EEGdir,'preprocessed_EEG_dat','average_reference',[sprintf('sub%03.0f',subID),'_EEGdat.mat']);
                load(fullfile(EEGdir,'preprocessed_EEG_dat','average_reference',[sprintf('sub%03.0f',subID),'_block_sess_id']));
                load(fullfile(EEGdir,'preprocessed_EEG_dat','average_reference',[sprintf('sub%03.0f',subID),'_trigger_vals_eegmatch']));
            end
            load(eegdat_fname);
            
            
            EEGdatadir =  fullfile(EEGdir,sprintf('sub%03.0f',subID),'eeg','average_reference');
            STdatadir = fullfile(EEGdir,sprintf('sub%03.0f',subID),'stim');
            BHVdatadir = fullfile(EEGdir,sprintf('sub%03.0f',subID),'behaviour');
            
            
    end
    
    
    nBlocks = 4;
    nChannels = 63;
    
    % specify file names for conv GLM results depending on design matrix
    switch design_matrix_type
        case 'absolute_coherences'
            save_name = sprintf('sub%03.0f_betas_abs_coh_IN.mat',subID);
        case 'jumps_and_jump_PEs'
            save_name = sprintf('sub%03.0f_betas_all_reg_IN.mat',subID);
        case 'jumps_plus_absolute'
            save_name = sprintf('sub%03.0f_betas_all_reg_plus_abs_NEW_IN4.mat',subID);
        case 'LRP_signed'
            save_name = sprintf('sub%03.0f_betas_LRP_signed_IN.mat',subID);
        case 'LRP_signed_wo_PE'
            save_name = sprintf('sub%03.0f_betas_LRP_signed_wo_PE_IN.mat',subID);
        case 'LRP_signed_wo_stim_stream'
            save_name = sprintf('sub%03.0f_betas_LRP_signed_wo_stim_stream_IN.mat',subID);
        case 'LRP_signed_wo_PE_wo_stim'
            save_name = sprintf('sub%03.0f_betas_LRP_signed_wo_PE_wo_stim_IN.mat',subID);
        case 'split_PE'
            save_name = sprintf('sub%03.0f_betas_split_PE_IN.mat',subID);
        case 'jumps_plus_absolute_vertical'
            save_name = sprintf('sub%03.0f_betas_all_reg_plus_abs_plus_vertical_IN.mat',subID);
        case 'trial_start_response'
            save_name = sprintf('sub%03.0f_betas_trial_start_response_N.mat',subID);
        case 'trial_start_response_coherence'
            save_name = sprintf('sub%03.0f_betas_trial_start_response_coherence_N.mat',subID);
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
    
    
    
    
    for i = 1:nSess % loop through sessions
        
        % load stim data and aligned EEG data
        fname_behav = fullfile(BHVdatadir,sprintf('sub%03.0f_sess%03.0f_behav.mat',subID,i));
        bhv{i} = load(fname_behav);
        fname_stim = fullfile(STdatadir,sprintf('sub%03.0f_sess%03.0f_stim.mat',subID,i));
        stim{i} = load(fname_stim);
        
        disp('session:')
        disp(i);
        fname_target = fullfile(EEGdatadir,sprintf('fMdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',subID,i));
        D{i} = spm_eeg_load(fname_target);
        
        
        keyboard; 
        
        %
        for b = 1:nBlocks % loop through blocks
            
            
            
            
            
            
            
            if block_Session_ID(i,b,1) % for some blocks the alignment of EEG data and triggers did not work, this is to exclude those blocks from the analysis
                
                blockID(i,b) = str2num(stim{i}.S.block_ID_cells{b});
                
                
                %
                %                 loading in the stim data from the orignal behavioural
                %                 files
                coherence = bhv{i}.B.coherence_frame{b}; %vector of coherence levels for this block
                coherence(coherence>1) = 1; coherence(coherence<-1) = -1; % in presentation code, if abs(coherence) is >1
                % then *all* dots move in same direction, i.e. coherence = 1
                % sj id to get right mean stim stream
                
                % this is loading in the stimulus stream from a response file that also has trigges for responses
                stream_sj = unique(all_responses(all_responses(:,12)==subID,11)); % I need to make sure I identify the correct subject
                mean_stim = mean_stim_streams_org{stream_sj,i}(:,blockID(i,b)); % indicating trials
                coherencet = stim_streams{stream_sj,i}(:,blockID(i,b));%vector of coherence levels for this block
                coherencet(coherencet>1) = 1; coherencet(coherencet<-1) = -1; % in presentation code, if abs(coherence) is >1
                %                 % then *all* dots move in same direction, i.e. coherence = 1
                
                if ~all(coherence == coherencet)
                    keyboard;
                end
                
                %               %%%%%%%%%%%%%%%%%%%%%%%%%%%% start of defining regressors
                %               for each block
                
                coherence_jump = abs([0; diff(coherence)])>0; %vector of coherence 'jumps'
                coherence_jump_level = coherence_jump.*abs(coherence); %vector of coherences of those 'jumps'
                
                
                
                % demean coherence jump level - but only non-zero values -
                mean_coh_jump_lev = mean(coherence_jump_level(coherence_jump_level ~= 0));
                coherence_jump_level(coherence_jump_level ~= 0) = coherence_jump_level(coherence_jump_level ~= 0) - mean_coh_jump_lev;
                
                
                % signed coherence jump level for LRP analysis
                coherence_jump_level_signed = coherence_jump.*coherence; %vector of coherence 'jumps'
                
                mean_coherence_jump_level_signed = mean(coherence_jump_level_signed(coherence_jump_level_signed~=0));
                coherence_jump_level_signed(coherence_jump_level_signed~=0) = coherence_jump_level_signed(coherence_jump_level_signed~=0) - mean_coherence_jump_level_signed;
                
                
                
                
                
                
                                correct_responses = trigger_vals_eegmatch{i}{b} == 201 |... % vector of button presses during trial periods
                                    trigger_vals_eegmatch{i}{b} == 205;
                                %                        trigger_vals_eegmatch{i}{b} == 205 |...
                                %                        trigger_vals_eegmatch{i}{b} == 206;
                
                                false_alarm_responses = trigger_vals_eegmatch{i}{b} == 202 |...
                                    trigger_vals_eegmatch{i}{b} == 206; % vector of button presses during intertrial periods
                
%                                 button_press_left = trigger_vals_eegmatch{i}{b} == 205 |... % vector of button presses left
%                                     trigger_vals_eegmatch{i}{b} == 206;
%                 
%                                 button_press_right = trigger_vals_eegmatch{i}{b} == 201 |... % vector of button presses right
%                                     trigger_vals_eegmatch{i}{b} == 202;
                
                               correct_trial_starts = trigger_vals_eegmatch{i}{b} == 30 |...  % get start of each trial for all coherence levels
                                    trigger_vals_eegmatch{i}{b} == 40 |...
                                    trigger_vals_eegmatch{i}{b} == 50 |...
                                    trigger_vals_eegmatch{i}{b} == 130 |...
                                    trigger_vals_eegmatch{i}{b} == 140 |...
                                    trigger_vals_eegmatch{i}{b} == 150;
                
                
                % regressor for prediciton error - absoluted difference in
                % coherences of current and previous jump
                diff_coherences = diff(coherence(coherence_jump));
                diff_coherences = [coherence(1); diff_coherences]; % differnce to prev cohernce for first coherence is that coherence itself
                jump_idx = find(coherence_jump);
                coherence_level_difference = zeros(size(coherence,1),1);
                coherence_level_difference(jump_idx) = abs(diff_coherences);
                
                % de-mean coherence Level difference (pred error)
                coherence_level_difference_mean = mean(coherence_level_difference(coherence_level_difference ~= 0));
                coherence_level_difference(coherence_level_difference ~= 0) =  coherence_level_difference(coherence_level_difference ~= 0) - coherence_level_difference_mean;
                
                
                % PE for LRP (same as above but signed)
                coherence_level_difference_signed = zeros(size(coherence,1),1);
                coherence_level_difference_signed(jump_idx) = diff_coherences;
                %demean
                coherence_level_difference_signed_mean = mean(coherence_level_difference_signed(coherence_level_difference_signed ~= 0));
                coherence_level_difference_signed(coherence_level_difference_signed ~= 0) =  coherence_level_difference_signed(coherence_level_difference_signed ~= 0) - coherence_level_difference_signed_mean;
                
                
                
                %%%%% different regressors for prediction errors towards
                %%%%% and away from zero
                
                jump_coherences = coherence(coherence_jump);
                % PE jump towards zero coherence
                % find jumps towards 0
                jumps_towards_zero = [0;sign(diff_coherences(2:end)) ~= sign(jump_coherences(1:end-1))];
                jumps_towards_zero_difference = zeros(size(coherence,1),1);
                jumps_towards_zero_difference (coherence_jump) = jumps_towards_zero .* diff_coherences;
                
                %demeaned
                jumps_towards_zero_difference_mean = mean(jumps_towards_zero_difference(jumps_towards_zero_difference ~= 0));
                jumps_towards_zero_difference(jumps_towards_zero_difference ~= 0) =  jumps_towards_zero_difference(jumps_towards_zero_difference ~= 0) - jumps_towards_zero_difference_mean;
                
                % PE jump away from zero
                jumps_away_from_zero = [0;sign(diff_coherences(2:end)) == sign(jump_coherences(1:end-1))];
                jumps_away_from_zero_difference = zeros(size(coherence,1),1);
                jumps_away_from_zero_difference (coherence_jump) = jumps_away_from_zero .* diff_coherences;
                
                %demeaned
                jumps_away_from_zero_difference_mean = mean(jumps_away_from_zero_difference(jumps_away_from_zero_difference ~= 0));
                jumps_away_from_zero_difference( jumps_away_from_zero_difference ~= 0) =   jumps_away_from_zero_difference( jumps_away_from_zero_difference ~= 0) -  jumps_away_from_zero_difference_mean;
                
                
                
                %%%%%%%%%%%%%%%%%% button press plus trial start regressors
                
                % get the responses for this subject, block and session
                id = all_responses(:,10)==i & all_responses(:,12)==subID & all_responses(:,9) == blockID(i,b);
                responses = all_responses(id,:);
                
                % trial responses - we need these later to match them with start frame
                % of trial
                id_trials = (responses(:,7) == 0 | responses(:,7) == 1 | responses(:,7) == 3);
                
                trial_frame = find(mean_stim(2:end)~= 0 & mean_stim(1:end-1) == 0);
                
                responses(id_trials,13) = trial_frame;
                
                %%%%% define regressors for correct responses
%                 correct_responses = zeros(length(mean_stim),1);
%                 correct_responses_diff = zeros(length(mean_stim),1);
%                 
%                 frames_correct_responses_right = responses(responses(:,7) == 1 & responses(:,3) == 1,6);
%                 frames_correct_responses_left = responses(responses(:,7) == 1 & responses(:,3) == 0,6);
%                 
%                 %%% regressor with 1 for right and -1 for left
%                 %%% responses
%                 correct_responses_diff(frames_correct_responses_right) = 1;
%                 correct_responses_diff(frames_correct_responses_left) = -1;
%                 
%                 correct_responses(frames_correct_responses_right) = 1;
%                 correct_responses(frames_correct_responses_left) =  1;
%                 
%                 
%                 %%% repeat for falmse alarms
%                 false_alarm_responses = zeros(length(mean_stim),1);
%                 false_alarm_responses_diff = zeros(length(mean_stim),1);
%                 
%                 frames_false_alarms_right = responses(responses(:,7) == 2 & responses(:,3) == 1,6);
%                 frames_false_alarms_left = responses(responses(:,7) == 2 & responses(:,3) == 0,6);
%                 
%                 false_alarm_responses_diff(frames_false_alarms_right) = 1;
%                 false_alarm_responses_diff(frames_false_alarms_left) = -1;
%                 
%                 false_alarm_responses(frames_false_alarms_right) = 1;
%                 false_alarm_responses(frames_false_alarms_left) = 1;
%                 
%                 % correct trial starts
%                 correct_trial_starts = zeros(length(mean_stim),1);
%                 frames_correct_trial_starts = responses(responses(:,7) == 1,13);
%                 correct_trial_starts(frames_correct_trial_starts) = 1;
%                 
%                 
%                 
%                 %%%% if we want to include a regressor for different coherences
%                 coherence_responses = zeros(length(mean_stim),1);
%                 
%                 frames_correct_responses_coh = responses(responses(:,7) == 1 ,[4 6]);
%                 coh_tr_id = abs(frames_correct_responses_coh(:,1)) == 0.3;
%                 coherence_responses(frames_correct_responses_coh(coh_tr_id,2),1) = -1;
%                 coh_tr_id = abs(frames_correct_responses_coh(:,1)) == 0.5;
%                 coherence_responses(frames_correct_responses_coh(coh_tr_id,2),1) = 1;
%                 
%                 
%                 correct_trial_starts_coh = zeros(length(mean_stim),1);
%                 
%                 frames_correct_trial_starts_coh = responses(responses(:,7) == 1,[4 13]);
%                 coh_tr_id = abs(frames_correct_trial_starts_coh(:,1)) == 0.3;
%                 correct_trial_starts_coh(frames_correct_trial_starts_coh(coh_tr_id,2),1) = -1;
%                 coh_tr_id = abs(frames_correct_trial_starts_coh(:,1)) == 0.5;
%                 correct_trial_starts_coh(frames_correct_trial_starts_coh(coh_tr_id,2),1) = 1;
                
                
                %%%
                if ITI_flag % flag to remove events during trial periods
                    mean_stim(mean_stim ~= 0) = 1;
                    mean_coherence = logical(mean_stim); % logical pointing to trial periods
                    
                    % exclude trial periods for all regressors
                    coherence_jump(mean_coherence) = 0;
                    coherence_jump_level(mean_coherence) = 0;
                    coherence_level_difference(mean_coherence) = 0;
                    jumps_towards_zero_difference(mean_coherence) = 0;
                    jumps_away_from_zero_difference(mean_coherence) = 0;
                    coherence_level_difference_signed(mean_coherence) = 0;
                    coherence_jump_level_signed(mean_coherence) = 0;
                end
                
                %%%%%%%%%%%%%%%%%%%%%
                % vertical motion regressors
                
                switch design_matrix_type
                    case 'jumps_plus_absolute_vertical'
                        coherence_vertical = stim{i}.S.coherence_frame_v{b}; %vector of coherence levels for this block
                        coherence_vertical(coherence_vertical>1) = 1; coherence_vertical(coherence_vertical<-1) = -1; % in presentation code, if abs(coherence) is >1
                        
                        
                        coherence_jump_vertical = abs([0; diff(coherence_vertical)])>0; %vector of coherence 'jumps'
                        coherence_jump_level_vertical = coherence_jump_vertical .* abs(coherence_vertical);
                        
                        diff_coherences_vertical = diff(coherence_vertical(coherence_jump_vertical));
                        diff_coherences_vertical = [coherence_vertical(1); diff_coherences_vertical]; % differnce to prev cohernce for first coherence is that coherence itself
                        jump_idx_v = find(coherence_jump_vertical);
                        coherence_level_difference_vertical = zeros(size(coherence_vertical,1),1);
                        coherence_level_difference_vertical(jump_idx_v) = abs(diff_coherences_vertical);
                        
                        % demean coherence jump level - but only non-zero values -
                        % not the zeros!!!!
                        mean_coh_jump_lev_vertical = mean(coherence_jump_level_vertical(coherence_jump_level_vertical ~= 0));
                        coherence_jump_level_vertical(coherence_jump_level_vertical ~= 0) = coherence_jump_level_vertical(coherence_jump_level_vertical ~= 0) - mean_coh_jump_lev_vertical;
                        
                        mean_coherence_level_difference_vertical = mean(coherence_level_difference_vertical(coherence_level_difference_vertical ~= 0));
                        coherence_level_difference_vertical(coherence_level_difference_vertical ~= 0) = coherence_level_difference_vertical(coherence_level_difference_vertical ~= 0) - mean_coherence_level_difference_vertical;
                    otherwise
                        % do nothing
                end
                
                nF = length(coherence);
                %
                switch design_matrix_type
                    case 'absolute_coherences'
                        regressor_list(1).value = abs(coherence);
                        regressor_list(1).nLagsBack = 150;
                        regressor_list(1).nLagsForward = 150;
                        regressor_list(1).name = 'abs stimulus';
                        
                        regressor_list(2).value = correct_responses;
                        regressor_list(2).nLagsBack = 150;
                        regressor_list(2).nLagsForward = 150;
                        regressor_list(2).name = 'button_press';
                        
                        regressor_list(3).value = false_alarm_responses;
                        regressor_list(3).nLagsBack = 150;
                        regressor_list(3).nLagsForward = 150;
                        regressor_list(3).name = 'iti button press';
                        
                        regressor_list(4).value = correct_trial_starts;
                        regressor_list(4).nLagsBack = 50;
                        regressor_list(4).nLagsForward = 700;
                        regressor_list(4).name = 'trial start';
                        
                    case 'jumps_and_jump_PEs'
                        
                        regressor_list(1).value = coherence_jump;
                        regressor_list(1).nLagsBack = 100; %unit
                        regressor_list(1).nLagsForward = 150;
                        regressor_list(1).name = 'coherence_jump';
                        
                        regressor_list(2).value = coherence_jump_level;
                        regressor_list(2).nLagsBack = 100;
                        regressor_list(2).nLagsForward = 150;
                        regressor_list(2).name = 'coherence_jump_level';
                        
                        regressor_list(3).value = coherence_level_difference;
                        regressor_list(3).nLagsBack = 150;
                        regressor_list(3).nLagsForward = 150;
                        regressor_list(3).name = 'prediction error';
                        
                        %                         regressor_list(4).value = button_press;
                        %                         regressor_list(4).nLagsBack = 150;
                        %                         regressor_list(4).nLagsForward = 150;
                        %                         regressor_list(4).name = 'button_press';
                        %
                        %                         regressor_list(5).value = button_press_incoh_motion;
                        %                         regressor_list(5).nLagsBack = 150;
                        %                         regressor_list(5).nLagsForward = 150;
                        %                         regressor_list(5).name = 'iti button press';
                        %
                        %                         regressor_list(6).value = trial_start;
                        %                         regressor_list(6).nLagsBack = 50;
                        %                         regressor_list(6).nLagsForward = 800;
                        %                         regressor_list(6).name = 'trial start';
                        
                    case 'jumps_plus_absolute'
                        
                        regressor_list(1).value = coherence_jump;
                        regressor_list(1).nLagsBack = 100;
                        regressor_list(1).nLagsForward = 150;
                        regressor_list(1).name = 'coherence_jump';
                        
                        regressor_list(2).value = coherence_jump_level;
                        regressor_list(2).nLagsBack = 100;
                        regressor_list(2).nLagsForward = 150;
                        regressor_list(2).name = 'coherence_jump_level';
                        
                        regressor_list(3).value = coherence_level_difference;
                        regressor_list(3).nLagsBack = 150;
                        regressor_list(3).nLagsForward = 150;
                        regressor_list(3).name = 'prediction error';
                        
                        regressor_list(4).value = abs(coherence);
                        regressor_list(4).nLagsBack = 150;
                        regressor_list(4).nLagsForward = 150;
                        regressor_list(4).name = 'abs stimulus';
                        
                        regressor_list(5).value = correct_responses;
                        regressor_list(5).nLagsBack = 150;
                        regressor_list(5).nLagsForward = 150;
                        regressor_list(5).name = 'button_press';
                        
                        regressor_list(6).value = false_alarm_responses;
                        regressor_list(6).nLagsBack = 150;
                        regressor_list(6).nLagsForward = 150;
                        regressor_list(6).name = 'iti button press';
                        
                        regressor_list(7).value = correct_trial_starts;
                        regressor_list(7).nLagsBack = 50;
                        regressor_list(7).nLagsForward = 800;
                        regressor_list(7).name = 'trial start';
                        %
                        
                        %
                        % %                         if we want to include vertical and horizontal EOG as
                        % %                         confounds:
                        % %
                        %                         regressor_list(8).value = EEGdat{i}{b}(63,:,:)';
                        %                         regressor_list(8).nLagsBack = 150;
                        %                         regressor_list(8).nLagsForward = 150;
                        %                         regressor_list(8).name = 'confound_EOG_reg_ver';
                        %
                        %                         regressor_list(9).value = EEGdat{i}{b}(64,:,:)';
                        %                         regressor_list(9).nLagsBack = 150;
                        %                         regressor_list(9).nLagsForward = 150;
                        %                         regressor_list(9).name = 'confound_EOG_reg_hor';
                        
                        
                        %}
                        
                    case 'LRP_signed'
                        
                        regressor_list(1).value = coherence_jump;
                        regressor_list(1).nLagsBack = 100;
                        regressor_list(1).nLagsForward = 150;
                        regressor_list(1).name = 'coherence_jump';
                        
                        regressor_list(2).value = coherence_jump_level_signed;
                        regressor_list(2).nLagsBack = 100;
                        regressor_list(2).nLagsForward = 150;
                        regressor_list(2).name = 'coherence_jump_level_signed';
                        
                        %                         regressor_list(3).value = button_press_right;
                        %                         regressor_list(3).nLagsBack = 150;
                        %                         regressor_list(3).nLagsForward = 150;
                        %                         regressor_list(3).name = 'button press right';
                        %
                        %                         regressor_list(4).value = button_press_left;
                        %                         regressor_list(4).nLagsBack = 150;
                        %                         regressor_list(4).nLagsForward = 150;
                        %                         regressor_list(4).name = 'button press left';
                        
                        regressor_list(3).value = coherence;
                        regressor_list(3).nLagsBack = 150;
                        regressor_list(3).nLagsForward = 150;
                        regressor_list(3).name = 'stimulus stream';
                        
                        regressor_list(4).value = coherence_level_difference_signed;
                        regressor_list(4).nLagsBack = 150;
                        regressor_list(4).nLagsForward = 150;
                        regressor_list(4).name = 'prediction error';
                        
                        
                    case 'LRP_signed_wo_PE'
                        regressor_list(1).value = coherence_jump;
                        regressor_list(1).nLagsBack = 100;
                        regressor_list(1).nLagsForward = 150;
                        regressor_list(1).name = 'coherence_jump';
                        
                        regressor_list(2).value = coherence_jump_level_signed;
                        regressor_list(2).nLagsBack = 100;
                        regressor_list(2).nLagsForward = 150;
                        regressor_list(2).name = 'coherence_jump_level_signed';
                        
                        %                         regressor_list(3).value = button_press_right;
                        %                         regressor_list(3).nLagsBack = 150;
                        %                         regressor_list(3).nLagsForward = 150;
                        %                         regressor_list(3).name = 'button press right';
                        %
                        %                         regressor_list(4).value = button_press_left;
                        %                         regressor_list(4).nLagsBack = 150;
                        %                         regressor_list(4).nLagsForward = 150;
                        %                         regressor_list(4).name = 'button press left';
                        
                        regressor_list(3).value = coherence;
                        regressor_list(3).nLagsBack = 150;
                        regressor_list(3).nLagsForward = 150;
                        regressor_list(3).name = 'stimulus stream';
                        
                    case 'LRP_signed_wo_stim_stream'
                        regressor_list(1).value = coherence_jump;
                        regressor_list(1).nLagsBack = 100;
                        regressor_list(1).nLagsForward = 150;
                        regressor_list(1).name = 'coherence_jump';
                        
                        regressor_list(2).value = coherence_jump_level_signed;
                        regressor_list(2).nLagsBack = 100;
                        regressor_list(2).nLagsForward = 150;
                        regressor_list(2).name = 'coherence_jump_level_signed';
                        
                        %                         regressor_list(3).value = button_press_right;
                        %                         regressor_list(3).nLagsBack = 150;
                        %                         regressor_list(3).nLagsForward = 150;
                        %                         regressor_list(3).name = 'button press right';
                        %
                        %                         regressor_list(4).value = button_press_left;
                        %                         regressor_list(4).nLagsBack = 150;
                        %                         regressor_list(4).nLagsForward = 150;
                        %                         regressor_list(4).name = 'button press left';
                        
                        regressor_list(3).value = coherence_level_difference_signed;
                        regressor_list(3).nLagsBack = 150;
                        regressor_list(3).nLagsForward = 150;
                        regressor_list(3).name = 'prediction error';
                        
                    case 'LRP_signed_wo_PE_wo_stim'
                        
                        regressor_list(1).value = coherence_jump;
                        regressor_list(1).nLagsBack = 100;
                        regressor_list(1).nLagsForward = 150;
                        regressor_list(1).name = 'coherence_jump';
                        
                        regressor_list(2).value = coherence_jump_level_signed;
                        regressor_list(2).nLagsBack = 100;
                        regressor_list(2).nLagsForward = 150;
                        regressor_list(2).name = 'coherence_jump_level_signed';
                        
                        %                         regressor_list(3).value = button_press_right;
                        %                         regressor_list(3).nLagsBack = 150;
                        %                         regressor_list(3).nLagsForward = 150;
                        %                         regressor_list(3).name = 'button press right';
                        %
                        %                         regressor_list(4).value = button_press_left;
                        %                         regressor_list(4).nLagsBack = 150;
                        %                         regressor_list(4).nLagsForward = 150;
                        %                         regressor_list(4).name = 'button press left';
                        
                    case 'split_PE'
                        
                        regressor_list(1).value = coherence_jump;
                        regressor_list(1).nLagsBack = 100;
                        regressor_list(1).nLagsForward = 150;
                        regressor_list(1).name = 'coherence_jump';
                        
                        regressor_list(2).value = coherence_jump_level;
                        regressor_list(2).nLagsBack = 100;
                        regressor_list(2).nLagsForward = 150;
                        regressor_list(2).name = 'coherence_jump_level';
                        
                        regressor_list(3).value = jumps_towards_zero_difference;
                        regressor_list(3).nLagsBack = 150;
                        regressor_list(3).nLagsForward = 150;
                        regressor_list(3).name = 'prediction error towards zero';
                        
                        regressor_list(4).value = jumps_away_from_zero_difference;
                        regressor_list(4).nLagsBack = 150;
                        regressor_list(4).nLagsForward = 150;
                        regressor_list(4).name = 'prediction error away from zero';
                        
                        %                         regressor_list(5).value = button_press;
                        %                         regressor_list(5).nLagsBack = 150;
                        %                         regressor_list(5).nLagsForward = 150;
                        %                         regressor_list(5).name = 'button_press';
                        %
                        %                         regressor_list(6).value = button_press_incoh_motion;
                        %                         regressor_list(6).nLagsBack = 150;
                        %                         regressor_list(6).nLagsForward = 150;
                        %                         regressor_list(6).name = 'iti button press';
                        %
                        %                         regressor_list(7).value = trial_start;
                        %                         regressor_list(7).nLagsBack = 50;
                        %                         regressor_list(7).nLagsForward = 800;
                        %                         regressor_list(7).name = 'trial start';
                        
                        
                    case 'jumps_plus_absolute_vertical'
                        
                        
                        regressor_list(1).value = coherence_level_difference;
                        regressor_list(1).nLagsBack = 150;
                        regressor_list(1).nLagsForward = 150;
                        regressor_list(1).name = 'prediction error';
                        
                        regressor_list(2).value = correct_responses;
                        regressor_list(2).nLagsBack = 150;
                        regressor_list(2).nLagsForward = 150;
                        regressor_list(2).name = 'button_press';
                        
                        regressor_list(3).value = abs(coherence);
                        regressor_list(3).nLagsBack = 150;
                        regressor_list(3).nLagsForward = 150;
                        regressor_list(3).name = 'abs stimulus';
                        
                        regressor_list(4).value = coherence_level_difference_vertical;
                        regressor_list(4).nLagsBack = 150;
                        regressor_list(4).nLagsForward = 150;
                        regressor_list(4).name = 'prediction error vertical';
                        
                        regressor_list(5).value = abs(coherence_vertical);
                        regressor_list(5).nLagsBack = 150;
                        regressor_list(5).nLagsForward = 150;
                        regressor_list(5).name = 'abs stimulus vertical';
                        
                        regressor_list(6).value = coherence_jump;
                        regressor_list(6).nLagsBack = 100;
                        regressor_list(6).nLagsForward = 150;
                        regressor_list(6).name = 'coherence_jump';
                        
                        regressor_list(7).value = coherence_jump_level;
                        regressor_list(7).nLagsBack = 100;
                        regressor_list(7).nLagsForward = 150;
                        regressor_list(7).name = 'coherence_jump_level';
                        
                        regressor_list(8).value = coherence_jump_vertical;
                        regressor_list(8).nLagsBack = 100;
                        regressor_list(8).nLagsForward = 150;
                        regressor_list(8).name = 'coherence_jump_vertical';
                        
                        regressor_list(9).value = coherence_jump_level_vertical;
                        regressor_list(9).nLagsBack = 100;
                        regressor_list(9).nLagsForward = 150;
                        regressor_list(9).name = 'coherence jump level vertical';
                        
                    case 'trial_start_response'
                        
                        regressor_list(1).value = correct_trial_starts;
                        regressor_list(1).nLagsBack = 50;
                        regressor_list(1).nLagsForward = 800;
                        regressor_list(1).name = 'triall start';
                        
                        regressor_list(2).value = correct_responses;
                        regressor_list(2).nLagsBack = 150;
                        regressor_list(2).nLagsForward = 150;
                        regressor_list(2).name = 'correct responses';
                        
                        regressor_list(3).value = false_alarm_responses;
                        regressor_list(3).nLagsBack = 150;
                        regressor_list(3).nLagsForward = 150;
                        regressor_list(3).name = 'false alarm';
                        
                    case 'trial_start_response_coherence'
                        
                        
                        
                        
                        
                        
                        correct_responses = zeros(length(mean_stim),3);
                        mean_correct_responses= zeros(length(mean_stim),1);
                        
                        frames_correct_responses = responses(responses(:,7) == 1, [4 6]);
                        mean_correct_responses(frames_correct_responses(:,2)) = 1;
                        
                        
                        correct_trial_starts = zeros(length(mean_stim),3);
                        mean_correct_trial_starts = zeros(length(mean_stim),1);
                        
                        frames_correct_trial_starts = responses(responses(:,7) == 1,[4 13]);
                        mean_correct_trial_starts(frames_correct_trial_starts(:,2)) = 1;
                        
                        id_col = 1;
                        for coh = [0.3 0.4 0.5]
                            
                            coh_tr_id = abs(frames_correct_trial_starts(:,1)) == coh;
                            coh_rp_id = abs(frames_correct_responses(:,1)) == coh;
                            correct_trial_starts(frames_correct_trial_starts(coh_tr_id,2), id_col) = 1;
                            correct_responses(frames_correct_responses(coh_rp_id,2), id_col) = 1;
                            id_col = id_col + 1;
                        end
                        
                        
                        regressor_list(1).value = correct_trial_starts(:,1);
                        regressor_list(1).nLagsBack = 50;
                        regressor_list(1).nLagsForward = 800;
                        regressor_list(1).name = 'triall start 0.3 coh';
                        
                        regressor_list(2).value = correct_trial_starts(:,2);
                        regressor_list(2).nLagsBack = 50;
                        regressor_list(2).nLagsForward = 800;
                        regressor_list(2).name = 'triall start 0.4 coh';
                        
                        regressor_list(3).value = correct_trial_starts(:,3);
                        regressor_list(3).nLagsBack = 50;
                        regressor_list(3).nLagsForward = 800;
                        regressor_list(3).name = 'triall start 0.5 coh';
                        
                        regressor_list(4).value = correct_responses(:,1);
                        regressor_list(4).nLagsBack = 150;
                        regressor_list(4).nLagsForward = 150;
                        regressor_list(4).name = 'correct responses coh 0.3';
                        
                        regressor_list(5).value = correct_responses(:,2);
                        regressor_list(5).nLagsBack = 150;
                        regressor_list(5).nLagsForward = 150;
                        regressor_list(5).name = 'correct responses coh 0.4';
                        
                        regressor_list(6).value = correct_responses(:,3);
                        regressor_list(6).nLagsBack = 150;
                        regressor_list(6).nLagsForward = 150;
                        regressor_list(6).name = 'correct responses coh 0.5';
                        
                        
                        %                                              regressor_list(7).value = mean_correct_responses;
                        %                         regressor_list(7).nLagsBack = 150;
                        %                         regressor_list(7).nLagsForward = 150;
                        %                         regressor_list(7).name = 'correct responses coh 0.5';
                        %
                        %                                                 regressor_list(8).value = mean_correct_trial_starts;
                        %                         regressor_list(8).nLagsBack = 50;
                        %                         regressor_list(8).nLagsForward = 800;
                        %                         regressor_list(8).name = 'trial start 0.5 coh';
                        
                    case 'trial_start_response_coherence_in_one_reg'
                        
                        
                        
                        regressor_list(1).value = correct_trial_starts;
                        regressor_list(1).nLagsBack = 50;
                        regressor_list(1).nLagsForward = 800;
                        regressor_list(1).name = 'triall start mean';
                        
                        regressor_list(2).value = correct_trial_starts_coh;
                        regressor_list(2).nLagsBack = 50;
                        regressor_list(2).nLagsForward = 800;
                        regressor_list(2).name = 'triall start coherence';
                        
                        regressor_list(3).value = correct_responses;
                        regressor_list(3).nLagsBack = 700;
                        regressor_list(3).nLagsForward = 150;
                        regressor_list(3).name = 'correct responses';
                        
                        regressor_list(4).value = correct_responses_diff;
                        regressor_list(4).nLagsBack = 700;
                        regressor_list(4).nLagsForward = 150;
                        regressor_list(4).name = 'correct responses right - left';
                        
                        
                        
                        
                    case 'trial_start_button_response_left_right'
                        
                        
                        regressor_list(1).value = correct_trial_starts;
                        regressor_list(1).nLagsBack = 50;
                        regressor_list(1).nLagsForward = 800;
                        regressor_list(1).name = 'triall start';
                        
                        regressor_list(2).value = correct_responses_right;
                        regressor_list(2).nLagsBack = 150;
                        regressor_list(2).nLagsForward = 150;
                        regressor_list(2).name = 'correct responses right';
                        
                        regressor_list(3).value = correct_responses_left;
                        regressor_list(3).nLagsBack = 150;
                        regressor_list(3).nLagsForward = 150;
                        regressor_list(3).name = 'correct responses left';
                        
                        regressor_list(4).value = false_alarm_responses_right;
                        regressor_list(4).nLagsBack = 150;
                        regressor_list(4).nLagsForward = 150;
                        regressor_list(4).name = 'FA responses right';
                        
                        regressor_list(5).value = false_alarm_responses_left;
                        regressor_list(5).nLagsBack = 150;
                        regressor_list(5).nLagsForward = 150;
                        regressor_list(5).name = 'FA responses left';
                        
                    case 'trial_start_button_response_left_minus_right'
                        
                        
                        regressor_list(1).value = correct_trial_starts;
                        regressor_list(1).nLagsBack = 50;
                        regressor_list(1).nLagsForward = 800;
                        regressor_list(1).name = 'triall start';
                        
                        regressor_list(2).value = correct_responses;
                        regressor_list(2).nLagsBack = 150;
                        regressor_list(2).nLagsForward = 150;
                        regressor_list(2).name = 'correct responses difference left right';
                        
                        regressor_list(3).value = false_alarm_responses;
                        regressor_list(3).nLagsBack = 150;
                        regressor_list(3).nLagsForward = 150;
                        regressor_list(3).name = 'false alarm responsed difference left right';
                        
                    case 'trial_start_button_response_left_minus_right_trial'
                        
                        
                        regressor_list(1).value = correct_trial_starts;
                        regressor_list(1).nLagsBack = 50;
                        regressor_list(1).nLagsForward = 800;
                        regressor_list(1).name = 'triall start left right';
                        
                        regressor_list(2).value = correct_responses;
                        regressor_list(2).nLagsBack = 150;
                        regressor_list(2).nLagsForward = 150;
                        regressor_list(2).name = 'correct responses difference left right';
                        
                        regressor_list(3).value = false_alarm_responses;
                        regressor_list(3).nLagsBack = 150;
                        regressor_list(3).nLagsForward = 150;
                        regressor_list(3).name = 'false alarm responsed difference left right';
                        
                    case 'trial_start_button_response_left_minus_right_plus_mean_response'
                        
                        regressor_list(1).value = correct_trial_starts;
                        regressor_list(1).nLagsBack = 50;
                        regressor_list(1).nLagsForward = 800;
                        regressor_list(1).name = 'triall start';
                        
                        regressor_list(2).value = correct_responses;
                        regressor_list(2).nLagsBack = 700;
                        regressor_list(2).nLagsForward = 150;
                        regressor_list(2).name = 'main effect correct responses';
                        
                        regressor_list(3).value = correct_responses_diff;
                        regressor_list(3).nLagsBack = 700;
                        regressor_list(3).nLagsForward = 150;
                        regressor_list(3).name = 'correct responses difference left right';
                        
                        regressor_list(4).value = false_alarm_responses;
                        regressor_list(4).nLagsBack = 700;
                        regressor_list(4).nLagsForward = 150;
                        regressor_list(4).name = 'main effect false alarm responses';
                        
                        regressor_list(5).value = false_alarm_responses_diff;
                        regressor_list(5).nLagsBack = 700;
                        regressor_list(5).nLagsForward = 150;
                        regressor_list(5).name = 'false alarm responsed difference left right';
                        
                        
                    case 'trial_start_button_response_left_minus_right_plus_mean_response_coherence'
                        
                        
                        
                        
                        
                        regressor_list(1).value = correct_trial_starts;
                        regressor_list(1).nLagsBack = 50;
                        regressor_list(1).nLagsForward = 800;
                        regressor_list(1).name = 'triall start';
                        
                        regressor_list(2).value = correct_responses;
                        regressor_list(2).nLagsBack = 700;
                        regressor_list(2).nLagsForward = 150;
                        regressor_list(2).name = 'main effect correct responses';
                        
                        regressor_list(3).value = correct_responses_diff;
                        regressor_list(3).nLagsBack = 700;
                        regressor_list(3).nLagsForward = 150;
                        regressor_list(3).name = 'correct responses difference left right';
                        
                        regressor_list(4).value = coherence_responses;
                        regressor_list(4).nLagsBack = 700;
                        regressor_list(4).nLagsForward = 150;
                        regressor_list(4).name = 'correct responses coherences';
                        
                        
                end
                %
                
                
                %                 if length(mean_stim_streams{stream_sj,i}(:,blockID(i,b))) ~= length(coherence)
                %                     faulty(fault_id,:) = [subID, i, blockID(i,b)]
                %
                %                     fault_id = fault_id + 1;
                %                 end
                
                Fs = D{i}.fsample;
                [lagged_design_matrix, time_idx] = create_lagged_design_matrix(regressor_list, Fs);
                
                VEOG_indx = selectchannels(D{i},'VEOG');
                pDM = geninv(lagged_design_matrix'); %pseudoinvesrse of design matrix
                for ch = 1:length(D{1}.chanlabels)
                    
                    % indices of badsamples for this channel and eyeblinks
                    removed_samples = badsamples{i}{b}(ch,:) | badsamples{i}{b}(VEOG_indx,:);
                    
                    tmp(ch,:) = pDM(:,~removed_samples)*EEGdat{i}{b}(ch,~removed_samples)';
                end
                for r = 1:length(regressor_list)
                    betas{r}(:,:,i,b) = tmp(:,time_idx(r).dm_row_idx); %betas (indexed by regressors): channels * lags * sessions * blocks
                end
                
            end
        end
    end
    
    
    channel_ind = 40; %channel of interest (CPz = 40);
    chanlbCPZ = D{1}.chanlabels(channel_ind); chanlbCPZ = chanlbCPZ{1};
    chanlabels = D{1}.chanlabels;
    
    switch reference_type
        case 'average_electrodes'
            save(fullfile(EEGdir,'preprocessed_EEG_dat_new','average_reference',save_name),'betas', 'time_idx','chanlbCPZ','channel_ind','chanlabels');
        case 'LM_RM'
            save(fullfile(EEGdir,'preprocessed_EEG_dat_new',save_name),'betas', 'time_idx','chanlbCPZ','channel_ind','chanlabels');
    end
    
end




end