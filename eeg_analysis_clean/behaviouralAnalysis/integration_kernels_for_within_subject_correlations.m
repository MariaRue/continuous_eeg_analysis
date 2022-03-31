
%% add paths
%addpath(genpath('/Users/maria/Documents/Matlab/cbrewer'))
%addpath('/Users/maria/Documents/MATLAB/raacampbell-shadedErrorBar-9b19a7b')
options = continuous_RDK_set_options('LTHiMac');

EEGpreproc = options.path.preproc.behaviour;  % path to behav data all subjs


condition = {'Tr frequent TR short', 'Tr frequent Tr short','Tr rare Tr short', 'Tr rare Tr long'};

% plot each subject separately for conditions
cl = cbrewer('div','RdBu', 12);
cl = cl([2 4 12 10],:);
%% load data
load_name = fullfile(EEGpreproc,'behav_data_all_subjs_all3');
load(load_name)

%% 

subj_list_behav = unique(all_responses(:,12));

%%

lags = 500; % frames back in time that are leading up to FA


which_responses = 'false alarms';  % calculating integration kernel for either FA button presses or button presses during trials options: 'false alarms' or 'trials', '3sec_rts',

with_coherence = 'without coherence level'; % if we want to know the coherence levels for the trials version
nS = length(unique(all_responses(:,12))); % number of subjects

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

kernels(:,1,:) = squeeze(mean(mean_coherences(:,[1,3],:),2)); % short kernels
kernels(:,2,:) = squeeze(mean(mean_coherences(:,[2,4],:),2)); % long kernels
kernels(:,3,:) = squeeze(mean(mean_coherences(:,[1,2],:),2)); % frequent kernels
kernels(:,4,:) = squeeze(mean(mean_coherences(:,[3,4],:),2)); % rare kernels 





%% fit curves do invidual subject data
options = optimset('MaxFunEvals',1000000,'MaxIter',100000);
for con = 1 : 4
    for sj = 1 : size(mean_coherences,3)
        data = kernels(:,con,sj); % data to fit exp model to
        
        % find the peak of the data as starting point
        [val,idx_max] = max(data);
        data_new{con,sj} = data(1:idx_max);
        % time steps
        dt = 0.01;
        num_steps = length(data_new{con,sj});
        
        t_sj{con,sj} = dt : dt : num_steps * dt;
        
        
        % initial param guesses for exp model
        pstart(1) = 1; % Amplitude
        pstart(2) = 1; % 1/tau
        % pstart(3) = 1; % offset
        
        fun = @(p)exp_residual(p,data_new{con,sj},t_sj{con,sj}); % this is the correct cost function that works
        [pnew_sj(con,:,sj),~,exitflag(con,sj)] = fminsearch(fun,pstart,options);
        
    end
end

% find subject with less than 20 responses to exclude from further analysis
[r,c] = find(num_false_alarms<=20);
r = unique(r);

%%
for con = 1:4
    figure;
    for sj = 1:size(mean_coherences,3)
        model{con,sj} = eval_exp(pnew_sj(con,:,sj),t_sj{con,sj});
        subplot(7,4,sj)
        title(['subject', ' ', num2str(subj_list_behav(sj))])
        if sj == 1
            title(condition{con})
            ylabel('mean coherence')
            xlabel('t to FA [s]')
        end
        hold on
        plot(data_new{con,sj},'Color',cl(con,:),'LineWidth',3);
        plot(model{con,sj},'k', 'LineWidth', 3);
        
        txt = ['FA: ', num2str(num_false_alarms(sj,con))];
        text(20,0.5,txt,'FontSize',14)
        xticks([0:100:500])
        xticklabels([ 5000 4000 3000 2000 1000 0])
        set(gca, 'XDir','reverse')
        tidyfig
        
        hold off
    end
end
%% 
%maria was originally removing subjects with less than 20 responses and
%>2SD away from mean. I will keep them in (and do rank correlations where
%needed)
tau_short = squeeze(pnew_sj(1,2,:));
%tau_short(r) = []; %removes subjects with <20 responses
tau_long = squeeze(pnew_sj(2,2,:));
%tau_long(r) = []; %removes subjects with <20 responses
tau_freq = squeeze(pnew_sj(3,2,:));
%tau_freq(r) = []; %removes subjects with <20 responses
tau_rare = squeeze(pnew_sj(4,2,:));
%tau_rare(r) = []; %removes subjects with <20 responses

% remove taus that lie more than 2 SD away from mean (probably not necessary?)
allTau = [tau_short,tau_long,tau_freq,tau_rare];
mean_tau = mean(allTau(:));
std_tau = std(allTau(:));

std_dist = mean_tau + 2*std_tau; 
std_dist2 = mean_tau - 2*std_tau;

[sj_1_tau,l1] = find(allTau >= std_dist);
[sj_2,l] = find(allTau <= std_dist2);

%keep these subjects in for now? LH edit
%tau_short(unique(sj_1_tau)) = [];
%tau_long(unique(sj_1_tau)) = []; 
%tau_freq(unique(sj_1_tau)) = []; 
%tau_rare(unique(sj_1_tau)) = []; 

amp_short = squeeze(pnew_sj(1,1,:));
amp_short(r) = [];
amp_long = squeeze(pnew_sj(2,1,:));
amp_long(r) = [];
amp_freq = squeeze(pnew_sj(3,1,:));
amp_freq(r) = [];
amp_rare = squeeze(pnew_sj(4,1,:));
amp_rare(r) = []; 




allAmp = [amp_short,amp_long,amp_freq,amp_rare];
mean_amp = mean(allAmp(:));
std_amp = std(allAmp(:));

std_dist = mean_amp + 2*std_amp; 
std_dist2 = mean_amp - 2*std_amp;

[sj_1_amp,l1] = find(allAmp >= std_dist);
[sj_2_amp,l] = find(allTau <= std_dist2);


tau_short(unique(sj_1_amp)) = [];
tau_long(unique(sj_1_amp)) = []; 
tau_freq(unique(sj_1_amp)) = []; 
tau_rare(unique(sj_1_amp)) = []; 



[R_tau_len, P_tau_len, RLO_tau_len, RUP_tau_len] = corrcoef(tau_short, tau_long); 
[R_amp_len, P_amp_len, RLO_amp_len, RUP_amp_len] = corrcoef(amp_short, amp_long);

[R_tau_freq, P_tau_freq, RLO_tau_freq, RUP_tau_freq] = corrcoef(tau_freq, tau_rare); 
[R_amp_freq, P_amp_freq, RLO_amp_freq, RUP_amp_freq] = corrcoef(amp_freq, amp_rare);

figure

subplot(2,2,1)
plot(tau_short,tau_long,'kd','LineWidth',1,'MarkerSize',8)
title([condition{con},' ','R= ',num2str(round(R_tau_len(1,2),2)),' ','P= ', num2str(P_tau_len(1,2),2)])
xlabel('short')
ylabel('long')
tidyfig;

subplot(2,2,2)
plot(amp_short,amp_long,'kd','LineWidth',1,'MarkerSize',8)
title([condition{con},' ','R= ',num2str(round(R_amp_len(1,2),2)),' ','P= ', num2str(P_amp_len(1,2),2)])
xlabel('short')
ylabel('long')
tidyfig;

subplot(2,2,3)
plot(tau_freq,tau_rare,'kd','LineWidth',1,'MarkerSize',8)
title([condition{con},' ','R= ',num2str(round(R_tau_freq(1,2),2)),' ','P= ', num2str(P_tau_freq(1,2),2)])
xlabel('frequent')
ylabel('rare')
tidyfig;

subplot(2,2,4)
plot(amp_freq,amp_rare,'kd','LineWidth',1,'MarkerSize',8)
title([condition{con},' ','R= ',num2str(round(R_amp_freq(1,2),2)),' ','P= ', num2str(P_amp_freq(1,2),2)])
xlabel('frequent')
ylabel('rare')
tidyfig;

% to save for subsequent neural analysis:
% save('stored_integration_kernels.mat','allTau','subj_list_behav')


%% helper function for exp fit

function [residual_reg] = exp_residual(params, data, t)
% this function computes the residuals for the exponential fit which is the
% sum of errers data - exp_model
% params is a vector with
%  amplitude, tau, offset


Amp = params(1);
tau = params(2);


model = Amp .* (exp(t/tau)); % compute the model

model(isinf(model)) = 0;
residual = sum((data - model').^2); % compute the error
residual_reg = residual + sum(params.^2) .* 0.01;

end




%% calculate model for estimated paramters

function model = eval_exp(params,t)
Amp = params(1);
tau = params(2);


model = Amp .* (exp(t/tau)); % compute the model
end




