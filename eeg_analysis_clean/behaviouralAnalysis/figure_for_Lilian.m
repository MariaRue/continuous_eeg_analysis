addpath(genpath('/Users/maria/Documents/Matlab/cbrewer'))
addpath('/Users/maria/Documents/MATLAB/raacampbell-shadedErrorBar-9b19a7b')
EEGpreproc = '/Volumes/crdkData/preprocessedData/behaviour';  % path to behav data all subjs


condition = {'Tr frequent TR short', 'Tr frequent Tr short','Tr rare Tr short', 'Tr rare Tr long'};

% plot each subject separately for conditions
cl = cbrewer('div','RdBu', 12);
cl = cl([2 4 12 10],:);
%% load data
load_name = fullfile(EEGpreproc,'behav_data_all_subjs_all3');
load(load_name)

lags = 500; % frames back in time that are leading up to FA


which_responses = 'false alarms';  % calculating integration kernel for either FA button presses or button presses during trials options: 'false alarms' or 'trials', '3sec_rts',


with_coherence = 'without coherence level'; % if we want to know the coherence levels for the trials version
nS = max(all_responses(:,11)); % number of subjects

coherence = [0.3 0.4 0.5];

counter = 0;

%% loop through subjects and find button presses
mean_coherences = [];
for sj = 1 : nS
    
    
    for bl = 1 :4
        combined = [];
        combined_coh = [];
        
        se_total = 6;
        if sj == 26
            
            se_total = 5 ;
            
        end
        for se = 1:se_total
            % select only stimstreams from all sessions that belong to specific
            % block
            
            % select all stim streams that belong to one subject
            stim_streams_sj = [];
            stim_streams_sj = stim_streams{sj,se}(:,bl);
            
            
            % select trigger streams that belong to one subject
            trigger_streams_sj = [];
            trigger_streams_sj = trigger_streams{sj,se}(:,bl);
            
            mean_streams = [];
            mean_streams = mean_stim_streams{sj,se}(:,bl);
            
            
            responses = all_responses((all_responses(:,9)== bl & all_responses(:,10) == se & all_responses(:,11) == sj),:);
            
            
            
            % find all triggers that lead to a button press
            switch which_responses
                
                
                case 'false alarms'
                    
                    % find triggers right and left button press (202 and 206)
                    
                    triggers_right = [];
                    triggers_left = [];
                    
                    
                    % this is with eeg triggers - don't use it
                    %                     triggers_right = find(trigger_streams_sj == 202);
                    %                     triggers_left = find(trigger_streams_sj == 206);
                    
                    
                    % this is with triggers from the response matrix
                    %
                    triggers_right = responses((responses(:,7) == 2 & responses(:,3) == 1),6);
                    triggers_left = responses((responses(:,7) == 2 & responses(:,3) == 0),6);
                    
                    
                    
                case 'trials'
                    % this is with EEG triggers - don't use it
                    
                    % this is with triggers from the EEG
                    %                     triggers_right = find(trigger_streams_sj == 202);
                    %                     triggers_left = find(trigger_streams_sj == 206);
                    
                    % this is with triggers from the response matrix
                    %
                    triggers_right = responses((responses(:,7) == 1 & responses(:,3) == 1),6);
                    triggers_left = responses((responses(:,7) == 1 & responses(:,3) == 0),6);
                    
                    
                    
                case '3sec_rts'
                    
                    % this is for eeg triggers - don't use it
                    
                    
                    % this is with triggers from the EEG
                    
                    %                     triggers_right = find(trigger_streams_sj == 201);
                    %                     triggers_left = find(trigger_streams_sj == 205);
                    % find all
                    
                    
                    
                    % this is for frames taken from response matrix
                    triggers_right = responses((responses(:,7) == 1 & responses(:,3) == 1),6);
                    triggers_left = responses((responses(:,7) == 1 & responses(:,3) == 0),6);
                    
                    % right trials
                    rts_rigth = zeros(length(triggers_right(:,1)),1);
                    for i = 1:length(triggers_right(:,1))
                        
                        t_org = triggers_right(i,1);
                        t = triggers_right(i,1);
                        
                        while ~(mean_streams(t) ~= 0 && mean_streams(t-1) == 0)
                            t = t-1;
                            
                        end
                        
                        rts_rigth(i) = t_org - t;
                        
                        if rts_rigth(i) > 300
                            
                            triggers_right(i,:) = nan;
                            counter = counter + 1;
                        end
                        
                    end
                    
                    
                    % left trials
                    rts_left = zeros(length(triggers_left(:,1)),1);
                    for i = 1:length(triggers_left(:,1))
                        
                        t_org = triggers_left(i,1);
                        t = triggers_left(i,1);
                        
                        while ~(mean_streams(t) ~= 0 && mean_streams(t-1) == 0)
                            t = t-1;
                            
                        end
                        
                        rts_left(i) = t_org - t;
                        
                        if rts_left(i) > 300
                            
                            triggers_left(i,:) = nan;
                            
                        end
                        
                    end
                    
                    
                    % remove nan trials - trials with rt > 300
                    triggers_right(isnan(triggers_right(:,1)),:) = [];
                    triggers_left(isnan(triggers_left(:,1)),:) = [];
            end
            
            
            
            
            % only choose triggers that are bigger than the lags we go
            % back
            
            
            triggers_right(triggers_right(:,1)<=lags,:) = [];
            triggers_left(triggers_left(:,1)<=lags,:) = [];
            
            
            switch with_coherence
                
                case 'with coherence levels'
                    
                    coh_val_right = [];
                    for i = 1:length(triggers_right(:,1))
                        
                        % this for EEG triggers
                        if mean_streams(triggers_right(i,1)) ~= 0
                            coh_val_right(i) = mean_streams(triggers_right(i,1));
                        else
                            t = 0; coh_val = 0;
                            while t <= 51 && coh_val == 0
                                t = t+1;
                                coh_val = mean_streams(triggers_right(i)-t,1);
                                
                                
                            end
                            
                            coh_val_right(i) = coh_val;
                        end
                        
                        
                        
                        
                    end
                    
                    coh_val_left = [];
                    for i = 1:length(triggers_left(:,1))
                        
                        if mean_streams(triggers_left(i,1)) ~= 0
                            coh_val_left(i) = mean_streams(triggers_left(i,1));
                        else
                            t = 0; coh_val = 0;
                            while t <= 51 && coh_val == 0
                                t = t+1;
                                coh_val = mean_streams(triggers_left(i) - t,1);
                                
                                
                            end
                            coh_val_left(i) = coh_val;
                        end
                        
                    end
                    
                    
                    
                    
                    
                    % get rid of incorrect trials
                    
                    triggers_right(coh_val_right < 0) = [];
                    triggers_left(coh_val_left > 0) = [];
                    coh_val_right(coh_val_right < 0) = [];
                    coh_val_left(coh_val_left > 0) = [];
                    
                    
                    
                    
                    
                    matrix_right = [];
                    
                    if any(triggers_right)
                        for i = 1:length(triggers_right(:,1))
                            
                            matrix_right(:,i) = stim_streams_sj(triggers_right(i,1) - lags : triggers_right(i,1));
                        end
                    end
                    
                    matrix_left = [];
                    
                    if any(triggers_left)
                        for i = 1:length(triggers_left(:,1))
                            
                            matrix_left(:,i) = stim_streams_sj(triggers_left(i,1) - lags : triggers_left(i,1)) .* -1 ; % to be able to combine with right button presses * -1
                        end
                    end
                    
                    combined = [combined,matrix_right, matrix_left];
                    combined_coh = [combined_coh,coh_val_right,abs(coh_val_left)];
                case 'without coherence level'
                    % loop through triggers for right and left button presses and
                    % select coherences from stim_streams_bl
                    matrix_right = [];
                    for i = 1:length(triggers_right(:,1))
                        
                        matrix_right(:,i) = stim_streams_sj(triggers_right(i,1) - lags : triggers_right(i,1));
                    end
                    
                    matrix_left = [];
                    for i = 1:length(triggers_left(:,1))
                        
                        matrix_left(:,i) = stim_streams_sj(triggers_left(i,1) - lags : triggers_left(i,1)) .* -1 ; % to be able to combine with right button presses * -1
                    end
                    
                    combined = [combined,matrix_right, matrix_left];
                    
            end
            
            
            
        end
        
        switch with_coherence
            
            case 'without coherence level'
                mean_coherences(:,bl,sj) = nanmean(combined,2);
                sem_coherence(:,bl,sj) = nanstd(combined')/sqrt(size(mean_coherences,1));
                num_false_alarms(sj,bl) = size(combined,2);
                
            case 'with coherence levels'
                
                for coh = 1:3
                    
                    if any(combined_coh == coherence(coh))
                        
                        idx_comb =combined_coh == coherence(coh);
                        
                        
                        mean_coherences(:,bl,sj,coh) = nanmean(combined(:,idx_comb),2);
                        sem_coherence(:,bl,sj,coh) = nanstd(combined(:,idx_comb)')/sqrt(size(mean_coherences,1));
                        num_false_alarms(sj,coh,bl) = size(combined,2);
                        
                        
                        
                    end
                    
                end
        end
        
        
    end
    
    
end





switch with_coherence
    case 'without coherence level'
        sem_across_subjects = squeeze(std(permute(mean_coherences,[3,1,2]))/sqrt(size(mean_coherences,3)));
        mean_across_subjects = nanmean(mean_coherences,3);
        
    case 'with coherence levels'
        
        for coh = 1:3
            sem_across_subjects(:,:,coh) = squeeze(std(permute(squeeze(mean_coherences(:,:,:,coh)),[3,1,2]))/sqrt(size(mean_coherences,3)));
            mean_across_subjects(:,:,coh) = nanmean(squeeze(mean_coherences(:,:,:,coh)),3);
        end
end


%%
% for each subject calculate mean rare and mean freq

for subject = 1:28
    meanFreqSubj(:,subject) = nanmean(squeeze(mean_coherences(:,1:2,subject)),2);
    
    meanRareSubj(:,subject) = nanmean(squeeze(mean_coherences(:,3:4,subject)),2);
    
end

meanFreq = nanmean(meanFreqSubj,2);
meanRare = nanmean(meanRareSubj,2);

seFreq = std(meanFreqSubj')/sqrt(28);

seRare = std(meanRareSubj')/sqrt(28);

%% plot
figure 
    %subplot(2,2,i)
    hold on
    plot(meanFreq, 'LineWidth', 3,'Color',cl(1,:));
    hold on
    h = shadedErrorBar(1:lags+1, meanFreq, seFreq, 'lineprops', '-k');
    h.patch.FaceColor = cl(1,:);
    h.mainLine.Color = cl(1,:);
    
        plot(meanRare, 'LineWidth', 3,'Color',cl(3,:));
    
       i = shadedErrorBar(1:lags+1, meanRare, seRare, 'lineprops', '-k');
    i.patch.FaceColor = cl(3,:);
    i.mainLine.Color = cl(3,:);
    
    
    legend({'frequent','','rare'})
    xlim([0 501])
    
    tidyfig
    xticks([0:100:500])
    xticklabels([-5 -4 -3 -2 -1 0])
   
        xlabel('time to button press [s]')
        ylabel('mean coherence')
   

plot(1:501,zeros(501,1),'k--','LineWidth',3)



