% This is a analysis pipeline for EEG data that is in SPM format but not
% preprocessed yet. This scirpt should automate the process of reading in
% the spm data and running several GLMs on the data, plotting the betas and
% topoplots

%% specify the datasets to read in by identifying the subid
% in terminal go to the data folder containing sufolders for each
% participant --> ls > txt (make sure that list contains only direct
% subject directories)
%scriptdir = fullfile('/Users/maria/Documents/Matlab/continuous_eeg_analysis/eeg_analysis');
%EEGdir= fullfile('/Users/maria/Documents/data/data.continuous_rdk','data','EEG');
%subj_list = [16,18:21,24,26,32, 52, 41, 34, 35, 51, 40, 28]; %2 missing
%subj_list = [16, 18:21, 24, 26, 27, 28, 29, 31, 33, 34, 35, 39, 40, 41, 42, 43, 47, 50, 51, 52, 54, 55, 57, 58];
clear all 
subj_list = [42];
[EEGdir,EEGdirdata,scriptdir,nSess,nS] = setup_EEG_session(subj_list);
source_density = 0; %use current source density transformed data?
 reference_type = 'LM_RM';

%% loop through each subject and load the data and match the behavioural triggers with the EEG triggers and save the resulting eegdatasets

for sj = 1:length(subj_list)
    subID = subj_list(sj);
    
    if subID == 55
        nSess = 5;
        
    else
        nSess = 6;
    end
    for i = 1:nSess % load in behavioural and stim data for each session
        fprintf('Subject %0.0f - run %0.0f\n',subID,i);
        
        BHVdatadir = fullfile(EEGdir,'subjects',sprintf('sub%03.0f',subID),'behaviour');
        STdatadir = fullfile(EEGdir,'subjects',sprintf('sub%03.0f',subID),'stim');
        EEGdatadir =  fullfile(EEGdir,'subjects',sprintf('sub%03.0f',subID),'eeg');
        
        switch reference_type
            case 'averaged_reference'
                if source_density
                    fname_target = fullfile(EEGdatadir,'average_reference',sprintf('cfMdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',subID,i));
                else
                    fname_target = fullfile(EEGdatadir,'average_reference',sprintf('fMdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',subID,i));
                end
                
            case 'LM_RM'
                if source_density
                    fname_target = fullfile(EEGdatadir,sprintf('cfMdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',subID,i));
                else
                    fname_target = fullfile(EEGdatadir,sprintf('fMdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',subID,i));
                end
                
        end
        
        D{i} = spm_eeg_load(fname_target);
        fname_behav = fullfile(BHVdatadir,sprintf('sub%03.0f_sess%03.0f_behav.mat',subID,i));
        bhv{i} = load(fname_behav);
        fname_stim = fullfile(STdatadir,sprintf('sub%03.0f_sess%03.0f_stim.mat',subID,i));
        stim{i} = load(fname_stim);
    
    end
    
 
    %
    %     cd (scriptdir)
    %
    %
    
    % match eeg triggers and behavioural triggers if EEGdat doesn't exist.
    % EEGdat is EEG data that matches the trigger vals in the behaviour
    switch reference_type
        case 'averaged_reference'
            if source_density
                eegdat_fname = fullfile(EEGdir,'preprocessed_EEG_dat','average_reference',[sprintf('sub%03.0f',subID),'_csdEEGdat.mat']);
            else
                eegdat_fname = fullfile(EEGdir,'preprocessed_EEG_dat','average_reference',[sprintf('sub%03.0f',subID),'_EEGdat.mat']);
            end
        case 'LM_RM'
            
            if source_density
                eegdat_fname = fullfile(EEGdir,'preprocessed_EEG_dat_new',[sprintf('sub%03.0f',subID),'_csdEEGdat.mat']);
            else
                eegdat_fname = fullfile(EEGdir,'preprocessed_EEG_dat_new',[sprintf('sub%03.0f',subID),'_EEGdat.mat']);
            end
            
            
    end
    
    % if exist(eegdat_fname) ~= 2
    
    
    nBlocks = 4;
    
    % initiliaise matrix with logicals for which blocks could be matched
    % for a session and blockID (1st - 4th Column =
    % block logical, sheet 2 = 1-4th column
    % BlockID corresponding to block in first sheet)
    % rows indicate session
    
    block_Session_ID = zeros(6,4,2);
    
    
    for i = 1:nSess %loop over sessions
        
        eeg_events = D{i}.events; %events in EEG data
        eeg_events = eeg_events(strmatch('trigger', {eeg_events.type})); %triggers in EEG data
        
        eob = find([eeg_events.value]==210); %end of block trigger
        sob = find([eeg_events.value]== 11);
        
        if i == 5
            keyboard; 
        end 
        if length(eob)~=nBlocks
            sprintf('didn''t find 4 end of blocks, subID %d, session %d',subID, i);
        elseif length(sob) ~= nBlocks
            sprintf('didn''t find 4 starts of blocks, subID %d, session %d',subID, i);
        else
            
            for b = 1:nBlocks % loop over blocks
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
%%--------------------matche_eeg_and_behaviour_triggers-------------------    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
                
                
                %a number of triggers in S.trigger_vals haven't made it into
                %the EEG data. We now try to correct this problem, making
                %'trigger_vals_eegmatch' - a version of S.trigger_vals that
                %only contains the triggers that are also in the EEG data
                trigger_vals_behav = bhv{i}.B.trigger_vals{b};
                tvind = find(trigger_vals_behav); tvlist = trigger_vals_behav(tvind); % a list of all the trigger values in S.trigger vals
                
                % eeg_events_block = eeg_events(eob(b)+1:eob(b+1)); %eeg_events just corresponding to this block
                eeg_events_block = eeg_events(sob(b):eob(b)); %eeg_events just corresponding to this block
                tvlist_eeg = [eeg_events_block.value]'; % a list of all trigger values that are in the EEG data
                
                if length(tvlist)<=length(tvlist_eeg)
                    % keyboard; % quick sanity check - are there more really triggers in S.tvlist? (answer = yes for pilot dataset)
                    error;
                elseif tvlist(end)~=210||tvlist_eeg(end)~=210 % second sanity check - is final trigger 210 in both lists? (answer is yes, i.e. 210 always found
                    error;
                end
                
                %now loop backwards from end of tvlist_eeg and try to match with tvlist
                teegind = length(tvlist_eeg);
                
                tind = length(tvlist);
                keep = []; %this will be a list of all entries in tvlist_eeg that we keep
                while tind>0 & teegind>0
                    if tvlist(tind)==tvlist_eeg(teegind)
                        keep(tind) = 1;
                        teegind = teegind - 1;
                    else
                        keep(tind) = 0;
                    end
                    tind = tind-1;
                    
                    %                     %these are the three places where there are oth bugs - might be worth further investigation by Maria
                    %                     if i==1&b==4&teegind==693 % at this time point a trigger has been send to the EEG recorder that does not exist or is not defined, which was number 2, correct trigger before would have been 26
                    %                         %   keyboard;
                    %                     elseif i==6&b==2&teegind==877 % very weird trigger 19 only occurd once in the tvlist vector at 886 - so in the past from 877 in the eeg list - maybe all triggers are sort of delayed in the eeg trig list?
                    %                         %  keyboard
                    %                     elseif i==6&b==3&teegind==345 % trigger 8 has been recorded in the eeg recording file but that trigger doesn't exist!
                    %                         %   keyboard
                    %                     elseif i == 1 & b == 1 && teegind == 2
                    %                         % keyboard;
                    %                     end
                end % while loop  % This might also explain the shifts we find in the lag between EEG and behav data? It also seems that over time the the lag between eeg triggers and behav triggers increases from 1 to 2 frames or more
                
                if length(tvlist(find(keep)))==length(tvlist_eeg) ...
                        && all(tvlist(find(keep))==tvlist_eeg) %we've matched them all!
                    trigger_vals_eegmatch{i}{b} = trigger_vals_behav;
                    trigger_vals_eegmatch{i}{b}(tvind(find(~keep))) = 0; %set the other ones to 0
                    
                    % block_Session_ID = session x block, first layer
                    % indicating blocks where trigger match worked,
                    % second layer indicates condition ID
                    
                    block_Session_ID(i,b,1) = 1;
                else
                    fprintf('Haven''t matched session %0.0f, block %0.0f\n',i,b);
                    
                end
                
                block_Session_ID(i,b,2) = str2double(stim{i}.S.block_ID_cells{b});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
%%--------------------matche_eeg_and_behaviour_triggers-------------------    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
            end
        end
       
        if i == 5
            keyboard; 
        end 
        
    end
    
    switch reference_type
        case 'averaged_reference'
            if source_density 
            save(fullfile(EEGdir,'preprocessed_EEG_dat','average_reference',[sprintf('sub%03.0f',subID),'_block_sess_id_csd']),'block_Session_ID','trigger_vals_eegmatch');
            save(fullfile(EEGdir,'preprocessed_EEG_dat','average_reference',[sprintf('sub%03.0f',subID),'_trigger_vals_eegmatch_csd']),'trigger_vals_eegmatch');
            else 
            save(fullfile(EEGdir,'preprocessed_EEG_dat','average_reference',[sprintf('sub%03.0f',subID),'_block_sess_id']),'block_Session_ID','trigger_vals_eegmatch');
            save(fullfile(EEGdir,'preprocessed_EEG_dat','average_reference',[sprintf('sub%03.0f',subID),'_trigger_vals_eegmatch']),'trigger_vals_eegmatch');
            end 
        case 'LM_RM'
            if source_density
            save(fullfile(EEGdir,'preprocessed_EEG_dat_new',[sprintf('sub%03.0f',subID),'_block_sess_id_csd']),'block_Session_ID','trigger_vals_eegmatch');
            save(fullfile(EEGdir,'preprocessed_EEG_dat_new',[sprintf('sub%03.0f',subID),'_trigger_vals_eegmatch_csd']),'trigger_vals_eegmatch');            
                
            else
            save(fullfile(EEGdir,'preprocessed_EEG_dat_new',[sprintf('sub%03.0f',subID),'_block_sess_id']),'block_Session_ID','trigger_vals_eegmatch');
            save(fullfile(EEGdir,'preprocessed_EEG_dat_new',[sprintf('sub%03.0f',subID),'_trigger_vals_eegmatch']),'trigger_vals_eegmatch');
            end 
    end
    % now we find the eeg data corresponding to the relevant time-periods in S, check that trigger channel is well aligned, and snip out this eeg data
    
    
    
    for i = 1:nSess
        eeg_events = D{i}.events; %events in EEG data
        eeg_events = eeg_events(strmatch('trigger', {eeg_events.type})); %triggers in EEG data
        eob = find([eeg_events.value]==210); %end of block trigger
        sob = find([eeg_events.value]==11);
        
        
        for b = 1:nBlocks
            
            if block_Session_ID(i,b,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
%%--------------------Match_eeg_data_with_stimulus_stream-------------------    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
                
                eeg_events_block = eeg_events(sob(b):eob(b)); %eeg_events just corresponding to this block
                blocklength = eeg_events_block(end).time-eeg_events_block(1).time; %length of block, in seconds
                nSamplesEEG = round(blocklength*D{1}.fsample)+1; %number of samples in block, at 100 Hz (fsample)
                
                
                %make a vector from the eeg trigger channel that corresponds to S.trigger_vals (with 0s where no trigger occurs)
                trigger_vals_eeg = zeros(nSamplesEEG,1);
                for s = 1:length(eeg_events_block)
                    trigger_vals_eeg(round((eeg_events_block(s).time-eeg_events_block(1).time)*D{i}.fsample)+1) = ...
                        eeg_events_block(s).value;
                    
                end
                
                
                
                trigger_vals_behav = bhv{i}.B.trigger_vals{b}; %trigger values from behaviour
                nSamplesBehav = length(trigger_vals_behav);
                
                %now 'zero-pad' the EEG triggers, as the behavioural triggers may run over in length
                zp_size = 500; %number of samples to zero pad by
                trigger_vals_eeg_zp = [zeros(zp_size,1); trigger_vals_eeg; zeros(zp_size,1)]; %zeropadded eeg triggervalues
                nSamplesEEG_zp = length(trigger_vals_eeg_zp);
                
                
                
                for c = 1:(nSamplesEEG_zp + 1 - nSamplesBehav)
                    nMatch(c) = sum(trigger_vals_behav==trigger_vals_eeg_zp(c:c+nSamplesBehav-1));
                    
                end
                
                
                
                %         plot(nMatch); disp(b); pause; % this reveals a clear 'spike' in every session -
                % where the triggers in the EEG data match the behavioural triggers
                %but a bit strangely, it doesn't always seem
                %to be at 500 - it is sometimes up to half a
                %second earlier - Maria to investigate?
                
                [~,best_match(i,b)] = max(nMatch);
                
                %we may be out here by one sample - I can't quite work the indexing out, but
                %it won't matter in the grand scheme of things...
                block_start_eeg_idx(i,b) = findc(D{1}.time,eeg_events_block(1).time)-zp_size+best_match(i,b);
                %
                
                block_end_eeg_idx(i,b)   = block_start_eeg_idx(i,b) + nSamplesBehav - 1;
                
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
%%--------------------Match_eeg_data_with_stimulus_stream-------------------    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
                
                %
                EEGdat{i}{b} = D{i}(:,block_start_eeg_idx(i,b):block_end_eeg_idx(i,b),1);
                badsamples{i}{b} = D{i}.badsamples(:,block_start_eeg_idx(i,b):block_end_eeg_idx(i,b),1);
                
                keyboard;
                
            end
        end
    end
    save(eegdat_fname,'EEGdat','badsamples');
    
    %end
    
end




