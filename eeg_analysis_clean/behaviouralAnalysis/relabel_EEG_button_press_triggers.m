function [new_trigger_list] = relabel_EEG_button_press_triggers(recorded_trigger_list, new_response_matrix)



% make a copy of the old trigger list
new_trigger_list = recorded_trigger_list;

% trigger labels

fa_left = 206;
fa_right = 202;
trial_left = 205;
trial_right = 201;

% select false alarm and trial frames 
behav_button_responses = new_response_matrix(new_response_matrix(:,7)==0 |...
    new_response_matrix(:,7)==1 | new_response_matrix(:,7)==2,:);



% make a list of button press frame number from the trigger list
recorded_response_EEG_triggers = find(recorded_trigger_list == 201 | recorded_trigger_list == 202 | recorded_trigger_list == 205 | recorded_trigger_list == 206);

% if there are more EEG triggers than behav triggers I would be quite
% worried and this function should throw an error to investigate the issue
if length(recorded_response_EEG_triggers) > length(behav_button_responses(:,6))
    
    warning('There are more EEG triggers for button presses than behavourial triggers')
    
end

% loop through behav_button_responses and compare with EEG responses and
% label correctly

EEG_tr = 1; % tracks idx of EEG triggers in case we have more behav triggers than eeg triggers 
for tr = 1:length(behav_button_responses(:,1))
    
    if behav_button_responses(tr,6) == recorded_response_EEG_triggers(EEG_tr)
        
        
        if behav_button_responses(tr,7) == 2
            
            if behav_button_responses(tr,3) ==0 % is a left response
                
                new_trigger_list(recorded_response_EEG_triggers(EEG_tr)) = fa_left;
                
            else  % it is a right choice
                new_trigger_list(recorded_response_EEG_triggers(EEG_tr)) = fa_right;
                
            end
            
        else % it is a trial
            
            if behav_button_responses(tr,3) ==0 % is a left response
                
                new_trigger_list(recorded_response_EEG_triggers(EEG_tr)) = trial_left;
                
            else  % it is a right choice
                new_trigger_list(recorded_response_EEG_triggers(EEG_tr)) = trial_right;
                
            end
            
            
        end
        
        EEG_tr = EEG_tr + 1; 
    else % loop through nearby frames and see whether we can find trigger there, if not then insert a new trigger at that time point
        
%         % check whether EEG trigger occurs in -+10 frames away from behave
%         % trigger - set starting frame f to behav trigger frame - 10
%         warning('eeg and behave triggers do not match exactly')
        
        f = behav_button_responses(tr,8) - 10;
        
        EEG_trigger_ID = 0;
        
        while f >= behav_button_responses(tr,6) - 10 && f <= behav_button_responses(tr,6) + 10
            
            if f == recorded_response_EEG_triggers(EEG_tr)
                
                EEG_trigger_ID = 1;
                
            else
                
                f = f+1;
                
            end
            
            
            
        end
        
        
        if EEG_trigger_ID % if triggers matched assign correct label
            
            EEG_tr = EEG_tr + 1; 
            
            if behav_button_responses(tr,7) == 2 % if false alarm
                
                if behav_button_responses(tr,3) ==0 % is a left response
                    
                    new_trigger_list(recorded_response_EEG_triggers(EEG_tr)) = fa_left;
                    
                else  % it is a right choice
                    new_trigger_list(recorded_response_EEG_triggers(EEG_tr)) = fa_right;
                    
                end
                
            else % it is a trial
                
                if behav_button_responses(tr,3) ==0 % is a left response
                    
                    new_trigger_list(recorded_response_EEG_triggers(EEG_tr)) = trial_left;
                    
                else  % it is a right choice
                    new_trigger_list(recorded_response_EEG_triggers(EEG_tr)) = trial_right;
                    
                end
                
                
            end
            
            
        else % if trigger is missing in EEG 
            
            if behav_button_responses(tr,7) == 2 % if false alarm
                
                if behav_button_responses(tr,3) ==0 % is a left response
                    
                    new_trigger_list(behav_button_responses(tr,6)) = fa_left;
                    
                else  % it is a right choice
                    new_trigger_list(behav_button_responses(tr,6)) = fa_right;
                    
                end
                
            else % it is a trial
                
                if behav_button_responses(tr,3) ==0 % is a left response
                    
                    new_trigger_list(behav_button_responses(tr,6)) = trial_left;
                    
                else  % it is a right choice
                    new_trigger_list(behav_button_responses(tr,6)) = trial_right;
                    
                end
                
                
            end
            
            
            
            warning(' we could not match all the behav triggers with the EEG triggers. There are some behavioural triggers that are not in the EEG data. These triggers have been inserted now')
            
        end
        
        
        
    end
    
    
    
end







end