function [blockStartEEGIdx,blockEndEEGIdx] = match_eeg_data_with_behavioural_triggers(blockEventsEEG,behaviouralTriggers,D)
%MATCH_EEG_DATA_WITH_BEHAVIOURAL_TRIGGERS

blocklength = blockEventsEEG(end).time-blockEventsEEG(1).time; %length of block, in seconds
nSamplesEEG = round(blocklength*D.fsample)+1; %number of samples in block, at 100 Hz (fsample)



%%%% this part here should be eegTriggerVals? Therefore, we don't need the
%%%% part above

%make a vector from the eeg trigger channel that corresponds to S.trigger_vals (with 0s where no trigger occurs)
trigger_vals_eeg = zeros(nSamplesEEG,1);
for s = 1:length(blockEventsEEG)
    trigger_vals_eeg(round((blockEventsEEG(s).time-blockEventsEEG(1).time)*D.fsample)+1) = ...
        blockEventsEEG(s).value;
    
end




%trigger values from behaviour
nSamplesBehav = length(behaviouralTriggers);

%now 'zero-pad' the EEG triggers, as the behavioural triggers may run over in length
zp_size = 500; %number of samples to zero pad by
trigger_vals_eeg_zp = [zeros(zp_size,1); trigger_vals_eeg; zeros(zp_size,1)]; %zeropadded eeg triggervalues
nSamplesEEG_zp = length(trigger_vals_eeg_zp);



for c = 1:(nSamplesEEG_zp + 1 - nSamplesBehav)
    nMatch(c) = sum(behaviouralTriggers==trigger_vals_eeg_zp(c:c+nSamplesBehav-1));
    
end



%         plot(nMatch); disp(b); pause; % this reveals a clear 'spike' in every session -
% where the triggers in the EEG data match the behavioural triggers
%but a bit strangely, it doesn't always seem
%to be at 500 - it is sometimes up to half a
%second earlier - Maria to investigate?

[~,best_match] = max(nMatch);

%we may be out here by one sample - I can't quite work the indexing out, but
%it won't matter in the grand scheme of things...
blockStartEEGIdx = findc(D.time,blockEventsEEG(1).time)-zp_size+best_match;
%

blockEndEEGIdx   = blockStartEEGIdx + nSamplesBehav - 1;


end