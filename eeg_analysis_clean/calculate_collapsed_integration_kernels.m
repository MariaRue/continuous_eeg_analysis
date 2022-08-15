function [ExpParameters,Correlation] = calculate_collapsed_integration_kernels(all_responses,SubjectList,nS, mean_stim_streams, stim_streams, trigger_streams,lags)


%% loop through subjects and find button presses
mean_coherences = [];
for subject = 1 : nS
    subjectID = SubjectList(subject);
    
    for condition = 1 :4
        combined = [];
        
        
        numSessions = 6;
        if subjectID == 26
            
            numSessions = 5 ;
            
        end
        for session = 1:numSessions
            % select only stimstreams from all sessions that belong to specific
            % block
            
            % select all stim streams that belong to one subject
            stim_streams_sj = [];
            stim_streams_sj = stim_streams{subjectID,session}(:,condition);
            
            
            % select trigger streams that belong to one subject
            trigger_streams_sj = [];
            trigger_streams_sj = trigger_streams{subjectID,session}(:,condition);
            
            mean_streams = [];
            mean_streams = mean_stim_streams{subjectID,session}(:,condition);
            
            
            responses = all_responses((all_responses(:,9)== condition & all_responses(:,10) == session & all_responses(:,11) == subjectID),:);
            
            
            
            % find all triggers that lead to a button press
            
            
            % find triggers right and left button press (202 and 206)
            
            triggers_right = [];
            triggers_left = [];
            
            
            
            
            % this is for eeg triggers - don't use it
            
            
            % this is with triggers from the EEG
            
            %                     triggers_right = find(trigger_streams_sj == 201);
            %                     triggers_left = find(trigger_streams_sj == 205);
            % find all
            
            
            
            % this is for frames taken from response matrix
            triggers_right = responses((responses(:,7) == 2 & responses(:,3) == 1),6);
            triggers_left = responses((responses(:,7) == 2 & responses(:,3) == 0),6);
            
            
         
            
            
            % only choose triggers that are bigger than the lags we go
            % back
            
            
            triggers_right(triggers_right(:,1)<=lags,:) = [];
            triggers_left(triggers_left(:,1)<=lags,:) = [];
            
            
            
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
        
        
        mean_coherences(:,condition,subject) = nanmean(combined,2);
        sem_coherence(:,condition,subject) = nanstd(combined')/sqrt(size(mean_coherences,1));
        num_false_alarms(subject,condition) = size(combined,2);
        
        
        
        
    end
    
    
end






sem_across_subjects = squeeze(std(permute(mean_coherences,[3,1,2]))/sqrt(size(mean_coherences,3)));
mean_across_subjects = nanmean(mean_coherences,3);



subjectIntegrationKernels(:,1,:) = squeeze(mean(mean_coherences(:,[1,3],:),2)); % short kernels
subjectIntegrationKernels(:,2,:) = squeeze(mean(mean_coherences(:,[2,4],:),2)); % long kernels
subjectIntegrationKernels(:,3,:) = squeeze(mean(mean_coherences(:,[1,2],:),2)); % frequent kernels
subjectIntegrationKernels(:,4,:) = squeeze(mean(mean_coherences(:,[3,4],:),2)); % rare kernels


%% fit curves do invidual subject data
options = optimset('MaxFunEvals',1000000,'MaxIter',100000);
for condition = 1 : 4
    for subject = 1 : size(mean_coherences,3)
        data = subjectIntegrationKernels(:,condition,subject); % data to fit exp model to
        
        % find the peak of the data as starting point
        [val,idx_max] = max(data);
        data_new{condition,subject} = data(1:idx_max);
        % time steps
        dt = 0.01;
        num_steps = length(data_new{condition,subject});
        
        t_sj{condition,subject} = -(num_steps-1)*dt:dt:0; %LH edit 280422 to make amplitude interpretable 
        
        
        % initial param guesses for exp model
        pstart(1) = 1; % Amplitude
        pstart(2) = 1; % 1/tau
        % pstart(3) = 1; % offset
        keyboard;
        fun = @(p)calculate_residuals_for_exponential_fit(p,data_new{condition,subject},t_sj{condition,subject}); % this is the correct cost function that works
        [ExpParameters.parameters(condition,:,subject),~,exitflag(condition,subject)] = fminsearch(fun,pstart,options);

    end
end

SubjectIntegrationKernels.DataForModelFit = data_new;
ExpParameters.time = t_sj;

%maria was originally removing subjects with less than 20 responses and
%>2SD away from mean. I will keep them in (and do rank correlations where
%needed)
ExpParameters.sorted.short.tau = squeeze(ExpParameters.parameters(1,2,:));
%tau_short(r) = []; %removes subjects with <20 responses
ExpParameters.sorted.long.tau = squeeze(ExpParameters.parameters(2,2,:));
%tau_long(r) = []; %removes subjects with <20 responses
ExpParameters.sorted.frequent.tau = squeeze(ExpParameters.parameters(3,2,:));
%tau_freq(r) = []; %removes subjects with <20 responses
ExpParameters.sorted.rare.tau = squeeze(ExpParameters.parameters(4,2,:));
%tau_rare(r) = []; %removes subjects with <20 responses

% remove taus that lie more than 2 SD away from mean (probably not necessary?)
% allTau = [tau_short,tau_long,tau_freq,tau_rare];
% mean_tau = mean(allTau(:));
% std_tau = std(allTau(:));
% 
% std_dist = mean_tau + 2*std_tau; 
% std_dist2 = mean_tau - 2*std_tau;
% 
% [sj_1_tau,l1] = find(allTau >= std_dist);
% [sj_2,l] = find(allTau <= std_dist2);

%keep these subjects in for now? LH edit
%tau_short(unique(sj_1_tau)) = [];
%tau_long(unique(sj_1_tau)) = []; 
%tau_freq(unique(sj_1_tau)) = []; 
%tau_rare(unique(sj_1_tau)) = []; 

ExpParameters.sorted.short.amplitude = squeeze(ExpParameters.parameters(1,1,:));
% amp_short(r) = [];
ExpParameters.sorted.long.amplitude = squeeze(ExpParameters.parameters(2,1,:));
% amp_long(r) = [];
ExpParameters.sorted.frequent.amplitude = squeeze(ExpParameters.parameters(3,1,:));
% amp_freq(r) = [];
ExpParameters.sorted.rare.amplitude = squeeze(ExpParameters.parameters(4,1,:));
% amp_rare(r) = []; 




% allAmp = [amp_short,amp_long,amp_freq,amp_rare];
% mean_amp = mean(allAmp(:));
% std_amp = std(allAmp(:));
% 
% std_dist = mean_amp + 2*std_amp; 
% std_dist2 = mean_amp - 2*std_amp;
% 
% [sj_1_amp,l1] = find(allAmp >= std_dist);
% [sj_2_amp,l] = find(allTau <= std_dist2);
% 
% 
% tau_short(unique(sj_1_amp)) = [];
% tau_long(unique(sj_1_amp)) = []; 
% tau_freq(unique(sj_1_amp)) = []; 
% tau_rare(unique(sj_1_amp)) = []; 

confidenceInterval = 0.95; 
[Correlation.length.coef.tau, Correlation.length.pvalue.tau] = corr(ExpParameters.sorted.short.tau, ExpParameters.sorted.long.tau, 'type', 'Spearman'); 
[Correlation.length.upperBound.tau, Correlation.length.lowerBound.tau] = calculate_confidence_interval_for_rank_correlation(Correlation.length.coef.tau, nS, confidenceInterval);

[Correlation.length.coef.amplitude, Correlation.length.pvalue.amplitude] = corr(ExpParameters.sorted.short.amplitude, ExpParameters.sorted.long.amplitude, 'type', 'Spearman');
[Correlation.length.upperBound.amplitude, Correlation.length.lowerBound.amplitude] = calculate_confidence_interval_for_rank_correlation(Correlation.length.coef.amplitude, nS, confidenceInterval);



[Correlation.frequency.coef.tau, Correlation.frequency.pvalue.tau] = corr(ExpParameters.sorted.frequent.tau, ExpParameters.sorted.rare.tau, 'type', 'Spearman'); 
[Correlation.frequency.upperBound.tau, Correlation.frequency.lowerBound.tau] = calculate_confidence_interval_for_rank_correlation(Correlation.frequency.coef.tau, nS, confidenceInterval);

[Correlation.frequency.coef.amplitude, Correlation.frequency.pvalue.amplitude] = corr(ExpParameters.sorted.frequent.amplitude, ExpParameters.sorted.rare.amplitude, 'type', 'Spearman');
[Correlation.frequency.upperBound.amplitude, Correlation.frequency.lowerBound.amplitude] = calculate_confidence_interval_for_rank_correlation(Correlation.frequency.coef.amplitude, nS, confidenceInterval);


% 
%RHO = corr(a.',b.','Type','Spearman');
% n = numel(a);
% STE = 1/sqrt(n-3);
% % here the input is 95% confidence interval, for 99% use 0.99:
% CI = norminv(0.95); 
% upper_bound = tanh(atanh(RHO)+CI*STE);
% lower_bound = tanh(atanh(RHO)-CI*STE);


end