%%%%
% save matrix for all subjects in one file
% add column with subject ID, and block ID as well as trial num and session
% num

% save different matrix with coherence traces, subject ID and block ID and
% session num

which_data = 'main experiment';

%% define subjects and and pathways

switch which_data
    
    case 'main experiment'
        
        % EEGdir = '/Volumes/LaCie 1/data/EEG';
        EEGpreproc = '/Volumes/LaCie/data/EEG';
        EEGdir = '/Volumes/LaCie/data/EEG';
%         subj_list = [16,18:21,24,26,32, 52, 41, 34, 35, 51, 40, 28];% old - best subjects
subj_list = [16, 18:21, 24, 26, 27, 28, 29, 31, 32, 33, 34, 35, 39, 40, 41, 42, 43, 47, 50, 51, 52, 54, 55, 57, 58];

        num_blocks = 4;
    case 'vertical control'
        
        
        
    case 'short trials pilot'
        EEGpreproc = '/Volumes/LaCie 1/data_preproc';
        EEGdir = '/Volumes/LaCie 1/pilot_short_trials';
        subj_list = [201,202];
        num_blocks = 3;
end

%% 

idx_block_start = 1; % counter for indexing into all_responses
all_responses = nan(10000,11); % big matrix with info across subjects

sample_rate = nan(length(subj_list),6);

%stim_streams = nan(31000,4,6,length(subj_list));s s
%trigger_streams = nan(31000,4,6,length(subj_list));


for sj = 1:length(subj_list) % loop through subjects
    subID = subj_list(sj);
    disp('subject:')
    disp(subID)
    
    BehaveDir = fullfile(EEGdir,sprintf('sub%03.0f',subID),'behaviour'); % pathway to subject behaviour folder
    StimDir = fullfile(EEGdir,sprintf('sub%03.0f',subID),'stim'); % pathway to subject stim folder
    
    se_total = 6; 
    
    if subID == 55 
        
        se_total = 5; 
        
    end 
    
    for se = 1:se_total % loop through sessions
        
        disp('session:')
        disp(se)
        % define filename and load relevant stimulus and behavioural variables for each subject
        file_ID_behaviour = fullfile(BehaveDir,sprintf('sub%03.0f_sess%03.0f_behav.mat',subID,se));
        file_ID_stim = fullfile(StimDir,sprintf('sub%03.0f_sess%03.0f_stim.mat',subID,se));
        
        
        response_mat= load(file_ID_behaviour,'respMat','B'); % responses
        responses = response_mat.respMat;
        behav_stim_trig = response_mat.B;
        
        
        bhv_stim_struct = load(file_ID_behaviour,'B'); % other stim info that changed during session
        coherence_stream = bhv_stim_struct.B.coherence_frame; % this is the coherence stream particiaptns have seen
        
        % load stim info
        
        stim = load(file_ID_stim,'S','tconst'); % original info - for example contains sequence of conditions
        sample_rate(sj,se) = stim.tconst.framerate;
        
        stim = load(file_ID_stim,'S'); % original info - for example contains sequence of conditions
        
        
        
        for bl = 1:num_blocks
            
            num_trials_per_block = max(max(stim.S.blocks_shuffled{bl}));
            
           
          
            [new_response_matrix,nResponsesDiscarded(sj,se,bl)] = return_new_response_matrix(behav_stim_trig.mean_coherence{bl},responses{bl});
 
            
            
             [new_trigger_list] = relabel_EEG_button_press_triggers( behav_stim_trig.trigger_vals{bl}, new_response_matrix);
            
            block = str2double(stim.S.block_ID_cells{bl}); % get condition ID
            
            num_responses = length(new_response_matrix(:,1)); % get num of butten presses and missed trials
            
            % all_responses is matrix button press/missed, respMat from
            % responses file, trials per block, condition ID, session,
            % subject

            
            %                       1: points won on current trial or for current
%                           response (check paramte.csv and create stimuli for
%                           more details on exact points that can be won or
%                           lost for a response)
%                       2: reaction time in secs 
%                       3: choice, 0 left, 1 right 
%                       4: coherence of dots 
%                       5: choice correct 1, incorrect 0  
%                       6: frame on which response
%                          occured 
%                       7: flag for type of response,
%                          0 incorrect response during coherent motion, 
%                          1 for correct response during coherent motion, 
%                          2 response during incoherent motion, 
%                          3 missed response to coherent motion 
%                          4 correct suppression of response during single
%                          trial version
            all_responses(idx_block_start:idx_block_start + num_responses-1,1:7) = new_response_matrix(:,1:7);
            all_responses(idx_block_start:idx_block_start + num_responses-1,8) = ones(num_responses,1) * num_trials_per_block;
            all_responses(idx_block_start:idx_block_start + num_responses-1,9) = ones(num_responses,1) * block;
            all_responses(idx_block_start:idx_block_start + num_responses-1,10) = ones(num_responses,1) * se;
            all_responses(idx_block_start:idx_block_start + num_responses-1,11) = ones(num_responses,1) * sj;
            all_responses(idx_block_start:idx_block_start + num_responses-1,12) = ones(num_responses,1) * subID;
            
                        % figure out 

            
            
            switch which_data
                case 'short trials pilot' % to know in pilot when one block ends and other begins - cause they are all the same condition
                    
                    all_responses(idx_block_start:idx_block_start + num_responses-1,13) = ones(num_responses,1) * bl;
                    
            end
            
            

            
            
            
            idx_block_start = idx_block_start + num_responses;
            
            
            % stim streams = coherence frames displayed, trigger streams =
            % vector with triggers for events such as button presses
            
            switch which_data
                
                case 'main experiment'
                    stim_streams{sj,se}(:,block) = behav_stim_trig.coherence_frame{bl};
                    
                    trigger_streams{sj,se}(:,block) = new_trigger_list;
                    
                    mean_stim_streams{sj, se}(:,block) = behav_stim_trig.mean_coherence{bl};
                    
                    mean_stim_streams_org{sj,se}(:,block) = stim.S.mean_coherence_org{bl};
                    
                    stim_streams_org{sj,se}(:,block) = stim.S.coherence_frame_org{bl};
                    
                    noise_streams{sj,se}(:,block) = stim.S.coherence_frame_noise{bl};
                    
                case 'short trials pilot'
                    
                    stim_streams{sj,se}(:,bl) = behav_stim_trig.coherence_frame{bl};
                    
                    trigger_streams{sj,se}(:,bl) = new_trigger_list;
                    
                    mean_stim_streams{sj, se}(:,bl) = behav_stim_trig.mean_coherence{bl};
                    
                    
            end
        end
        
        
        
    end
    
end

% remove nan values
all_responses = all_responses(~isnan(all_responses(:,1)),:);



% at some point save these matrices at a place that makes sense

switch which_data
    
    case 'short trials pilot'
        save_name = fullfile(EEGpreproc,'behav_data_all_subjs_pilot');
    case 'main experiment'
        save_name = fullfile(EEGpreproc,'behav_data_all_subjs_all');
end

save(save_name,'trigger_streams','stim_streams','all_responses','sample_rate','mean_stim_streams','mean_stim_streams_org','stim_streams_org','noise_streams');
