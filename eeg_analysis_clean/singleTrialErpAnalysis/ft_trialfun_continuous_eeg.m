function [trl] = ft_trialfun_continuous_eeg(cfg)

timelock_event = cfg.timelock_event; 

% read the header information and the events from the data
hdr   = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset);

% search for "trigger" events
value  = [event(find(strcmp('trigger', {event.type}))).value]';
sample = [event(find(strcmp('trigger', {event.type}))).sample]';

condition_id = unique(cfg.subject_responses(:,9),'stable');

% determine the number of samples before and after the trigger
pretrig  = -round(cfg.trialdef.prestim  * hdr.Fs);
posttrig =  round(cfg.trialdef.poststim * hdr.Fs);


% find block starts and ends

sob = sample(value == 11);

eob = sample(value == 210);

blocks = length(sob);

if blocks ~= 4
    
    block_id = 2;
    
else
    
    block_id = 1;
    
end

trl = [];
for c = 1 :blocks
   
    % find all response events
    id = cfg.subject_responses(:,9) == condition_id(block_id);
    responses = cfg.subject_responses(id,:);
    
    % trial responses - we need these later to match them with start frame
    % of trial
    id_trials = (responses(:,7) == 0 | responses(:,7) == 1 | responses(:,7) == 3);
    
    % find trial starts for all trials
    trial_start_frames = find(cfg.mean_stim_streams(2:end,condition_id(block_id)) ~=0 & cfg.mean_stim_streams(1:end-1,condition_id(block_id)) ==0);
    
    % convert frames into actual EEG sample id for fieldtrip
    responses(:,6) = sob(c) + responses(:,6);
    trial_start_samples = sob(c) + trial_start_frames;
    
    responses(id_trials,13) = trial_start_samples;
    
    switch timelock_event
        
        case 'trialStart'
            
            trlbegin = trial_start_samples + pretrig;
            trlend   = trial_start_samples + posttrig;
            offset   = pretrig;
            newtrl   = [trlbegin trlend ones(length(trlend),1) * offset responses(id_trials,:)];
            

        case 'buttonPress'
            
            trlbegin = responses(:,6) + pretrig;
            trlend   = responses(:,6) + posttrig;
            offset   = pretrig;
            newtrl   = [trlbegin trlend ones(length(trlend),1) * offset responses];
            
    end
    
    trl      = [trl; newtrl];

   
    
    block_id = block_id + 1;
end



end