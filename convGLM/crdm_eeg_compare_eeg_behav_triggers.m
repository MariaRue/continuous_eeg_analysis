function [eegTriggerVals, flagBlockMatched] = ...
    crdm_eeg_compare_eeg_behav_triggers( eegTriggers, bhvTriggers )
%COMPARE_EEG_AND_BEHAVIOURAL_TRIGGERS A number of triggers in S.trigger_vals 
%haven't made it into the EEG data. We now try to correct this problem, making
%'trigger_vals_eegmatch' - a version of S.trigger_vals that only contains the 
%triggers that are also in the EEG data

idxNonzeroTriggers = find(bhvTriggers); 
bhvTriggerList = bhvTriggers(idxNonzeroTriggers); % a list of all the non-zero
% trigger values in S.trigger_vals

try 
    eegTriggerList = [eegTriggers.value]'; % a list of all trigger values that are in the EEG data

    if length(bhvTriggerList) <= length(eegTriggerList)
        % quick sanity check - can't have fewer behav triggers for this to
        % work
        error('Fewer behavioural triggers detected than EEG triggers.')
        
    elseif bhvTriggerList(end) ~= 210 || eegTriggerList(end) ~= 210 
        % second sanity check - final trigger must be 210 in both lists
        error('Final trigger in EEG and/or behav trigger list is not 210');
    end
catch 
    % Why are we using a try-catch here?
end 

% now loop backwards from end of EEG trigger list and try to match with
% behavioural list
currentEegIndex = length(eegTriggerList);
currentBhvIndex = length(bhvTriggerList);

keep = []; % this will be a list of all entries in behav triggers that we keep
while currentBhvIndex > 0 && currentEegIndex > 0
    if bhvTriggerList(currentBhvIndex) == eegTriggerList(currentEegIndex)
        % we have a match, so we keep this trigger
        keep(currentBhvIndex) = 1;
        % we have matched this EEG trigger, move to the next one
        currentEegIndex = currentEegIndex - 1;
    else
        % behav trigger doesn't correspond to EEG trigger - stay with this
        % EEG trigger and try again with the next behav trigger
        keep(currentBhvIndex) = 0;
    end
    % move to the next behav trigger - either to try again on the current
    % EEG trigger, or to compare with the next EEG trigger if the previous
    % has been matched
    currentBhvIndex = currentBhvIndex-1;
end

% check whether we matched all the triggers 
if currentEegIndex > 0
    warning('Unmatched EEG triggers remaining');
end

if length(bhvTriggerList(keep==1)) == length(eegTriggerList) ...
        && all(bhvTriggerList(keep==1) == eegTriggerList) % we've matched them all!
    
    eegTriggerVals = bhvTriggers;
    eegTriggerVals(idxNonzeroTriggers(find(~keep))) = 0; % set the other ones to 0

    flagBlockMatched = 1; 
  
else 
    flagBlockMatched = 0; 
	eegTriggerVals = [];
end

end