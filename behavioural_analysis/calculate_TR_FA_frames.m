function [triggers_left, triggers_right] = calculate_TR_FA_frames(mean_stim, all_responses,which_case, fb, bl, sj, se,condition_ID)



start_trial_idx = mean_stim(1:end-1) == 0 & mean_stim(2:end) ~= 0;
start_trial_frames = find(start_trial_idx)+1; % frames of trials at start
end_trial_idx = mean_stim(1:end-1) ~= 0 & mean_stim(2:end) == 0;
end_trial_frames = find(end_trial_idx); % last frames of trial

FA_frames_responses = all_responses(all_responses(:,7)==2 & all_responses(:,13) == bl & all_responses(:,11) == sj & all_responses(:,10) == se,[3,6]);
TR_frames_responses = all_responses(all_responses(:,7)==1 & all_responses(:,13) == bl & all_responses(:,11) == sj & all_responses(:,10) == se,[3,6]);
% in this case flex feedback = 50 frames
% trial length  = 100 frames - but save this with response function

if condition_ID == 1
    tr_len = 100;
elseif condition_ID == 2
    tr_len = 300;
    
else 
    
    keyboard; 
    
end
%check whether FA happened during trial period - which would be super bad
t = zeros(length(start_trial_frames),1);
for tr = 1:length(start_trial_frames)
    
    t(tr) = any(FA_frames_responses(:,2) <= end_trial_frames(tr) & FA_frames_responses(:,2) >= start_trial_frames(tr));
    
end

if sum(t) ~= 0
    keyboard;
end

%check which FAs are within a trial period - but check that the trial period actually is full length and there wasn't a trial response yet
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
    

    

    
    
    
    switch which_case
        
        case 'FA'
            
                % updated FA frames
    FA_frames = FA_frames_responses;
    FA_frames(FA_TR_ID,:) = [];
    
            
            triggers_left = FA_frames(FA_frames(:,1)==0,2);
            triggers_right = FA_frames(FA_frames(:,1)==1,2);
            
        case 'Trial'
            
            % first check whether all wrong FAs are correct responses to
            % trial period
            coherence_TR(coherence_TR == 0) = [];
            
            FA_TR_FR = FA_frames_responses(FA_TR_ID,:);
            
            try
            FA_TR_FR = [FA_TR_FR,coherence_TR];
            catch 
                keyboard; 
            end 
            
            % for left responses
            idx_wrong_fa_press = (FA_TR_FR(:,1)==0 & FA_TR_FR(:,3) > 0);
            
            if any(idx_wrong_fa_press)
                
                FA_TR_FR(idx_wrong_fa_press) = [];
                
            end
            
            % for right responses
            
            idx_wrong_fa_press = (FA_TR_FR(:,1)==1 & FA_TR_FR(:,3) < 0);
            
            if any(idx_wrong_fa_press)
                
                FA_TR_FR(idx_wrong_fa_press) = [];
                
            end
            
            
            triggers_left = sort([TR_frames_responses(TR_frames_responses(:,1)==0,2);FA_TR_FR(FA_TR_FR(:,1)==0,2)]);
            triggers_right = sort([TR_frames_responses(TR_frames_responses(:,1)==1,2);FA_TR_FR(FA_TR_FR(:,1)==1,2)]);
            
    end
    
end