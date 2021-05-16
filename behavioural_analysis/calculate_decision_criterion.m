% this script takes in streams and triggers from the matrixes genarated
% with 'read_in_behav_data'

%% add paths to plotting packages/data

addpath(genpath('/Users/maria/Documents/Matlab/cbrewer'))
addpath('/Users/maria/Documents/MATLAB/raacampbell-shadedErrorBar-9b19a7b')
EEGpreproc = '/Volumes/LaCie 1/data_preproc';  % path to behav data all subjs

cl = cbrewer('div','RdBu', 12);
cl = cl([2 4 12 10],:);
condition = {'ITIS INTS', 'ITIS INTL','ITIL INTS', 'ITIL INTL'};
%% load data
load_name = fullfile(EEGpreproc,'behav_data_all_subjs');
load(load_name)

lags = 100; % frames back in time that are leading up to peak of FA over which we calculate the mean coherence for FA
% with respect to plots of mean curves for all 4 conditions leading up to
% FA, peak coherence is always at around 35-40 frames before the button press
% (from calculate_mean_coh_leading_to_FA script) - removing 33 frames
% before button press)
nS = size(stim_streams,4); % num subjects
% subract first 33 frames which participants don't use any more for
% making a decision
% with respect to plots of mean curves for all 4 conditions leading up to
% FA, peak coherence is always at around 35-40 frames before the button press
% (from calculate_mean_coh_leading_to_FA script) - removing 33 frames
% before button press)
button_resp = 33;
%% calculate mean coherence and round to nearest 100ms for FAs
for sj = 1 : nS
    
    all_streams = [];
    all_streams = permute(squeeze(stim_streams(:,:,:,sj)),[1 3 2]); % permute order of entries so that we can cat across sessions for one subject
    all_streams = reshape(all_streams, [], size(stim_streams,2),1); % reshape is actually doing the concatenating
    
    % repeat the same for the triggers
    all_triggers = [];
    all_triggers = permute(squeeze(trigger_streams(:,:,:,sj)),[1 3 2]);
    all_triggers = reshape(all_triggers, [], size(trigger_streams,2),1);
    
    
    for bl  = 1:4 % loop through the conditions
        
        
        % find all triggers for a left and right button press FA
        % triggers for FA = 202 right button press, 206 = left button press
        triggers_right =  find(all_triggers(:,bl) == 202);
        triggers_left =  find(all_triggers(:,bl) == 206);
        
        
        
        triggers_right(triggers_right <= lags + button_resp) = [];
        triggers_left(triggers_left <= lags + button_resp) = [];
        
        
        triggers_right = triggers_right - button_resp;
        triggers_left = triggers_left - button_resp;
        
        matrix_coherences_right = [];
        for i = 1:length(triggers_right)
            
            matrix_coherences_right(i) = round(mean(all_streams(triggers_right(i) - lags : triggers_right(i),bl)),1);
        end
        
        matrix_coherences_left = [];
        for ii = 1:length(triggers_left)
            matrix_coherences_left(ii) = round(mean(all_streams(triggers_left(ii) - lags : triggers_left(ii),bl) .* -1),1);     % multiply by -1 to cmopare to button press to the riht
        end
        
        FA_coherences{sj,bl} = [matrix_coherences_right,matrix_coherences_left];
        sec_iti_time = sum(sum(squeeze(mean_stim_streams(:,bl,:,sj)==0)))/100;
        out_FA{sj,bl} = [unique(FA_coherences{sj,bl})',histc(FA_coherences{sj,bl},unique(FA_coherences{sj,bl}))'./sec_iti_time];
    end
    
    
    
    
end

%% repeat the same for button presses during trials
for sj = 1 : nS
    
    all_streams = [];
    all_streams = permute(squeeze(stim_streams(:,:,:,sj)),[1 3 2]); % permute order of entries so that we can cat across sessions for one subject
    all_streams = reshape(all_streams, [], size(stim_streams,2),1); % reshape is actually doing the concatenating
    
    
    
    % repeat the same for the triggers
    all_triggers = [];
    all_triggers = permute(squeeze(trigger_streams(:,:,:,sj)),[1 3 2]);
    all_triggers = reshape(all_triggers, [], size(trigger_streams,2),1);
    
    
    for bl  = 1:4 % loop through the conditions
        
        
        % find all triggers for a left and right button press FA
        % triggers for FA = 202 right button press, 206 = left button press
        triggers_right =  find(all_triggers(:,bl) == 201);
        triggers_left =  find(all_triggers(:,bl) == 205);
        
        
        
        triggers_right(triggers_right <= lags + button_resp) = [];
        triggers_left(triggers_left <= lags + button_resp) = [];
        
        
        triggers_right = triggers_right - button_resp;
        triggers_left = triggers_left - button_resp;
        
        matrix_coherences_right = [];
        for i = 1:length(triggers_right)
            
            matrix_coherences_right(i) = round(mean(all_streams(triggers_right(i) - lags : triggers_right(i),bl)),1);
        end
        
        matrix_coherences_left = [];
        for ii = 1:length(triggers_left)
            matrix_coherences_left(ii) = round(mean(all_streams(triggers_left(ii) - lags : triggers_left(ii),bl) .* -1),1);     % multiply by -1 to cmopare to button press to the riht
        end
        
        TR_coherences{sj,bl} = [matrix_coherences_right,matrix_coherences_left];
        num_trials = sum(all_triggers(:,bl) == 30 | all_triggers(:,bl) == 130 | all_triggers(:,bl) == 40 | all_triggers(:,bl) == 140 | all_triggers(:,bl) == 50 | all_triggers(:,bl) == 150);
        out_tr{sj,bl} = [unique(TR_coherences{sj,bl})',histc(TR_coherences{sj,bl},unique(TR_coherences{sj,bl}))'./num_trials];
    end
    
    
    
    
end
%% for each subject and each coherence level make a plot across conditions
list_coh = 0.1:0.1:1;
for sj = 1:nS
    
    hit_rate = zeros(10,4);
    fa_rate = zeros(10,4);
    
    for bl = 1:4
        [~,idx_hit,~] = intersect( round(out_tr{sj,bl}(:,1),2), round(list_coh',2) );
        [~,idx_fa,~] = intersect( round(out_FA{sj,bl}(:,1),2), round(list_coh',2) );
        hit_rate(idx_hit,bl) = out_tr{sj,bl}(idx_hit,2);
        fa_rate(idx_fa,bl) = out_FA{sj,bl}(idx_fa,2);
        
    end
    
    fig_idx = 1; 
    figure 
    title(['subject: ',num2str(sj)])
    for r = 3:8
        subplot(3,2,fig_idx)
        
        for con = 1:4
            hold on
        plot(fa_rate(r,con), hit_rate(r,con), 'x', 'Color',cl(con,:), 'MarkerSize',10, 'LineWidth',5)
        hold off
        end
        ylabel('hit rate')
        xlabel('FA rate')
        if fig_idx == 1
        legend(condition)
        end
        tidyfig
        title(['coherence level ', num2str(list_coh(r))])
        fig_idx = fig_idx + 1; 
    end
    
end