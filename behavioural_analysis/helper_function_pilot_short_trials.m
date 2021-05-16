% some fals alarms were in feedback time of stim - find those 

% D = response output matrix, lines indicating trial or response (in continuous motin version)  
%               column: 
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
%                           trial version

EEGpreproc = '/Volumes/LaCie 1/data_preproc';  % path to behav data all subjs
load_name = fullfile(EEGpreproc,'behav_data_all_subjs_pilot');
load(load_name)

bl = 1; 
se = 1; 
sj = 1; 
%%
mean_stim = mean_stim_streams{sj,se}(:,bl);
start_trial_idx = mean_stim(1:end-1) == 0 & mean_stim(2:end) ~= 0;
start_trial_frames = find(start_trial_idx)+1; % frames of trials at start
end_trial_idx = mean_stim(1:end-1) ~= 0 & mean_stim(2:end) == 0; 
end_trial_frames = find(end_trial_idx); % last frames of trial 

FA_frames_responses = all_responses(all_responses(:,7)==2 & all_responses(:,13) == bl & all_responses(:,11) == sj & all_responses(:,10) == se,[3,6]);

% in this case flex feedback = 50 frames
% trial length  = 100 frames - but save this with response function 

fb = 50; 



if unique(all_responses(all_responses(:,11) == sj & all_responses(:,10) == se, 9)) == 1
tr_len = 100; 
else 
    tr_len = 300;
    
end 
%% check whether FA happened during trial period - which would be super bad 
t = zeros(length(start_trial_frames),1);
for tr = 1:length(start_trial_frames)
    
    t(tr) = any(FA_frames_responses(:,2) <= end_trial_frames(tr) & FA_frames_responses(:,2) >= start_trial_frames(tr)); 
    
end

if sum(t) ~= 0
    keyboard; 
end 

%% check which FAs are within a trial period - but check that the trial actuall was 1sec long and there wasn't a trial response yet 
FA_TR_ID = zeros(length(end_trial_frames),1);
coherence_TR = zeros(length(end_trial_frames),1);
for tr = 1:length(end_trial_frames)
    
    
    length_tr = end_trial_frames(tr) - start_trial_frames(tr); 

    
    if length_tr == tr_len
    
        
  ID = find(FA_frames_responses(:,2) <= end_trial_frames(tr) + fb & FA_frames_responses(:,2) >= end_trial_frames(tr));  
  
  
  if ID 
      FA_TR_ID(tr) = ID; 
      coherence_TR(tr) = mean_stim(start_trial_frames(tr)); 
      
  end 
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

