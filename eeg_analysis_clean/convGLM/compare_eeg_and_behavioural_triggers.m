function [eegTriggerVals,flagBlockMatched] = compare_eeg_and_behavioural_triggers(blockEventsEEG,behaviouralTriggers)
%COMPARE_EEG_AND_ BEHAVIOURAL_TRIGGERS



%a number of triggers in S.trigger_vals haven't made it into
%the EEG data. We now try to correct this problem, making
%'trigger_vals_eegmatch' - a version of S.trigger_vals that
%only contains the triggers that are also in the EEG data

triggerValuesInd = find(behaviouralTriggers); triggerValuesList = behaviouralTriggers(triggerValuesInd); % a list of all the trigger values in S.trigger vals

try 
triggerValuesListEEG = [blockEventsEEG.value]'; % a list of all trigger values that are in the EEG data

if length(triggerValuesList)<=length(triggerValuesListEEG)
    % keyboard; % quick sanity check - are there more really triggers in S.tvlist? (answer = yes for pilot dataset)
    error;
elseif triggerValuesList(end)~=210||triggerValuesListEEG(end)~=210 % second sanity check - is final trigger 210 in both lists? (answer is yes, i.e. 210 always found
    error;
end
catch 
   
end 

%now loop backwards from end of tvlist_eeg and try to match with tvlist
teegind = length(triggerValuesListEEG);

tind = length(triggerValuesList);

keep = []; %this will be a list of all entries in tvlist_eeg that we keep
while tind>0 & teegind>0
    if triggerValuesList(tind)==triggerValuesListEEG(teegind)
        keep(tind) = 1;
        teegind = teegind - 1;
    else
        keep(tind) = 0;
    end
    tind = tind-1;
    
end

% check whether we matched all the triggers 

if length(triggerValuesList(find(keep)))==length(triggerValuesListEEG) ...
        && all(triggerValuesList(find(keep))==triggerValuesListEEG) %we've matched them all!
    
    eegTriggerVals = behaviouralTriggers;
    eegTriggerVals(triggerValuesInd(find(~keep))) = 0; %set the other ones to 0

    flagBlockMatched = 1; 
  
else 
    
        flagBlockMatched = 0; 
        eegTriggerVals = [];
       
 end
end