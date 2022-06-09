%%%%
% save matrix for all subjects in one file
% add column with subject ID, and block ID as well as trial num and session
% num

% save different matrix with coherence traces, subject ID and block ID and
% session num



%% define subjects and and pathways
options = continuous_RDK_set_options('LTHiMac');
subjectList = [62:64,66,68,70];

reference = 'LMRM';
csdFlag = 0; 

idx_block_start = 1; % counter for indexing into all_responses
all_responses = nan(10000,11); % big matrix with info across subjects

sample_rate = nan(length(subjectList),6);

%stim_streams = nan(31000,4,6,length(subj_list));s s
%trigger_streams = nan(31000,4,6,length(subj_list));


for sj = 1:length(subjectList) % loop through subjects
    subID = subjectList(sj);
    disp('subject:')
    disp(subID)
    [details,paths] =  conrdk_subjects(subID,options,reference,csdFlag); 

    for se = 1:length(details.sessionIDs) % loop through sessions
       
        sessionID = details.sessionIDs(se);
        disp('session:')
        disp(se)
        % define filename and load relevant stimulus and behavioural variables for each subject
        file_ID_behaviour = paths.behaviour(se).sessionList;
        file_ID_stim = paths.stimulus(se).sessionList;
        
        
        response_mat= load(file_ID_behaviour,'respMat','B'); % responses
        responses = response_mat.respMat;
        behav_stim_trig = response_mat.B;
        
        
        bhv_stim_struct = load(file_ID_behaviour,'B'); % other stim info that changed during session
        coherence_stream = bhv_stim_struct.B.coherence_frame; % this is the coherence stream particiaptns have seen
        
        % load stim info
        
        stim = load(file_ID_stim,'S','tconst'); % original info - for example contains sequence of conditions
        sample_rate(sj,sessionID) = stim.tconst.framerate;
        
        stim = load(file_ID_stim,'S'); % original info - for example contains sequence of conditions
        
        
        
        for bl = 1:4
            
            num_trials_per_block = max(max(stim.S.blocks_shuffled{bl}));
            
           
          
            [new_response_matrix,nResponsesDiscarded(sj,sessionID,bl)] = return_new_response_matrix(behav_stim_trig.mean_coherence{bl},responses{bl});
 
            
            
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
            all_responses(idx_block_start:idx_block_start + num_responses-1,10) = ones(num_responses,1) * sessionID;
            all_responses(idx_block_start:idx_block_start + num_responses-1,11) = ones(num_responses,1) * sj;
            all_responses(idx_block_start:idx_block_start + num_responses-1,12) = ones(num_responses,1) * subID;
            
                        % figure out 

  
            
            idx_block_start = idx_block_start + num_responses;
            
            
            % stim streams = coherence frames displayed, trigger streams =
            % vector with triggers for events such as button presses
            
         
                    stim_streams{sj,sessionID}(:,block) = behav_stim_trig.coherence_frame{bl};
                    
                    trigger_streams{sj,sessionID}(:,block) = new_trigger_list;
                    
                    mean_stim_streams{sj, sessionID}(:,block) = behav_stim_trig.mean_coherence{bl};
                    
                    mean_stim_streams_org{sj,sessionID}(:,block) = stim.S.mean_coherence_org{bl};
                    
                    stim_streams_org{sj,sessionID}(:,block) = stim.S.coherence_frame_org{bl};
                    
                    noise_streams{sj,sessionID}(:,block) = stim.S.coherence_frame_noise{bl};
                    
                   
                    if isfield(stim.S,'coherence_frame_v')
                        
                        vertical_stim_streams{sj,sessionID}(:,block) = stim.S.coherence_frame_v{bl};
                    end 

           
        end
        
        
        
    end
    
end

% remove nan values
all_responses = all_responses(~isnan(all_responses(:,1)),:);



% at some point save these matrices at a place that makes sense


        save_name = 'behav_data_all_subjs_allVertical';

if isfield(stim.S,'coherence_frame_v')
    save(save_name,'trigger_streams','stim_streams','all_responses','sample_rate','mean_stim_streams','mean_stim_streams_org','stim_streams_org','noise_streams','vertical_stim_streams');

else 
save(save_name,'trigger_streams','stim_streams','all_responses','sample_rate','mean_stim_streams','mean_stim_streams_org','stim_streams_org','noise_streams');
end