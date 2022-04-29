function all_regressors = select_and_prepare_regressor_data_for_subject_level_glm(subID,session,condition,all_responses,stim_streams, mean_stim_streams,vertical,vertical_stim_streams)
%PREPARE_REGRESSOR_DATA_FOR_SUBJECT_LEVEL_GLM
% Collects all regressors for the convolutional GLM from behavioural and stimulus data
%
%IN: subID = subject identifier
%    session = session id
%    condition = 1 frequent/short
%                2 frequent/long
%                3 rare/short
%                4 rare/long
%    all_responses = response matrix with information to button presses and
%    trial information such as coherence, left/right choice
%    stim_stream = stimulus participant has seen
%    mean_stim_stream = indicates trial periods and average coherence
%
%OUT: all_regressors = struct with all regressors that might be needed by
%GLM with


%--select data for given subject, session and condition-------------------%
stream_sj = unique(all_responses(all_responses(:,12)==subID,11)); % I need to make sure I identify the correct subject
mean_stim = mean_stim_streams{stream_sj,session}(:,condition); % indicating trials
coherence = stim_streams{stream_sj,session}(:,condition);%vector of coherence levels for this block
coherence(coherence>1) = 1;
coherence(coherence<-1) = -1; % in presentation code, if abs(coherence) is >1
%                 % then *all* dots move in same direction, i.e. coherence = 1

% get the responses for this subject, session and condition
id = all_responses(:,10)==session & all_responses(:,12)==subID & all_responses(:,9) == condition;
responses = all_responses(id,:);

%%%% define regressors for correct responses and false alarms
correct_responses = zeros(length(mean_stim),1); % allocate space for regressor for an average impulse response function for a correct button press during trial (including left and right)
correct_responses_diff = zeros(length(mean_stim),1); % allocate space for a regressor looking at the difference waveform between left and right correct responses

% select frames of correct responses
frames_correct_responses_right = responses(responses(:,7) == 1 & responses(:,3) == 1,6);
frames_correct_responses_left = responses(responses(:,7) == 1 & responses(:,3) == 0,6);

%%% regressor with 1 for right and -1 for left
%%% responses to look at difference waveform between left and right
%%% responses
correct_responses_diff(frames_correct_responses_right) = 1;
correct_responses_diff(frames_correct_responses_left) = -1;

% average impulse response function for a correct response
correct_responses(frames_correct_responses_right) = 1;
correct_responses(frames_correct_responses_left) =  1;


%%% repeat for falmse alarms
false_alarm_responses = zeros(length(mean_stim),1); % allocate space for regressor for an average impulse response function for a false alarm button press during intertrial time (including left and right)
false_alarm_responses_diff = zeros(length(mean_stim),1);% allocate space for a regressor looking at the difference waveform between left and right false alarms

frames_false_alarms_right = responses(responses(:,7) == 2 & responses(:,3) == 1,6);
frames_false_alarms_left = responses(responses(:,7) == 2 & responses(:,3) == 0,6);

false_alarm_responses_diff(frames_false_alarms_right) = 1;
false_alarm_responses_diff(frames_false_alarms_left) = -1;

false_alarm_responses(frames_false_alarms_right) = 1;
false_alarm_responses(frames_false_alarms_left) = 1;

%%%%%%%%%%%%%%%%%%%% coherence correct responses

coherence_responses = zeros(length(mean_stim),1);
frames_correct_responses_coh = responses(responses(:,7) == 1 ,[4 6]);
coh_tr_id = abs(frames_correct_responses_coh(:,1)) == 0.3;
coherence_responses(frames_correct_responses_coh(coh_tr_id,2),1) = -1;
coh_tr_id = abs(frames_correct_responses_coh(:,1)) == 0.5;
coherence_responses(frames_correct_responses_coh(coh_tr_id,2),1) = 1;



%%%%%%%%%%%%%%%%%%% correct trial starts


% find indices of correct trial responses
id_trials = (responses(:,7) == 0 | responses(:,7) == 1 | responses(:,7) == 3);

% find start frame of trials
trial_frame = find(mean_stim(2:end)~= 0 & mean_stim(1:end-1) == 0);

% create a additional column to save the start frame of trials
responses(id_trials,13) = trial_frame;

% regressor for trial start
correct_trial_starts = zeros(length(mean_stim),1);
frames_correct_trial_starts = responses(responses(:,7) == 1,13);
correct_trial_starts(frames_correct_trial_starts) = 1;

%%%%%%%%%%% intertrial regressors for noise

coherence_jump = abs([0; diff(coherence)])>0; %regressor values for coherence 'jumps'

coherence_jump_level = coherence_jump.*abs(coherence); %regressor values for coherence levels of those 'jumps'

coherence_jump_level_signed = coherence_jump.*coherence; %regressor values for coherence levels of those 'jumps', signed by left/right motion

% demean coherence jump level - but only non-zero values -
mean_coh_jump_lev = mean(coherence_jump_level(coherence_jump_level ~= 0));
coherence_jump_level(coherence_jump_level ~= 0) = coherence_jump_level(coherence_jump_level ~= 0) - mean_coh_jump_lev;


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





%%----------VERTICAL----------------------------------------------%%%%%%%%%
if vertical
coherence_vertical = vertical_stim_streams{stream_sj,session}(:,condition); %vector of coherence levels for this block
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



all_regressors.coherence_jump_vertical           = coherence_jump_vertical;
all_regressors.coherence_jump_level_vertical     = coherence_jump_level_vertical;
all_regressors.prediction_error_vertical         = coherence_level_difference_vertical;
all_regressors.absoluted_stimulus_vertical       = abs(coherence_vertical);


end 
% 

%%--------------------------------------------------------------%%%%%%%%%%%







% remove events trom trial periods for all initertrial regressors
mean_stim(mean_stim ~= 0) = 1;
mean_coherence = logical(mean_stim); % logical pointing to trial periods
coherence_jump(mean_coherence) = 0;
coherence_jump_level(mean_coherence) = 0;
coherence_jump_level_signed(mean_coherence) = 0;
coherence_level_difference(mean_coherence) = 0;
signed_stimulus = coherence; coherence(mean_coherence) = 0;

% build structure will all regressors from which we later select the
% regressors we need for a specific GLM model
all_regressors.coherence_jump           = coherence_jump;
all_regressors.coherence_jump_level     = coherence_jump_level;
all_regressors.coherence_jump_level_signed= coherence_jump_level_signed;
all_regressors.prediction_error         = coherence_level_difference;
all_regressors.trial_start              = correct_trial_starts;
all_regressors.correct_trial_response   = correct_responses;
all_regressors.false_alarm              = false_alarm_responses;
all_regressors.difference_waveform_correct_trial_response = correct_responses_diff;
all_regressors.difference_waveform_false_alarm = false_alarm_responses_diff;
all_regressors.absoluted_stimulus       = abs(coherence);
all_regressors.signed_stimulus          = signed_stimulus;
all_regressors.coherence_responses      = coherence_responses; 
end