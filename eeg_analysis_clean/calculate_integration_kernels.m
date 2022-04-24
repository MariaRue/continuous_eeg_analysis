function [GroupIntegrationKernels, SubjectIntegrationKernels, SignificantTimePoints, ExpParameters] = calculate_integration_kernels(all_responses,SubjectList,nS, mean_stim_streams, stim_streams, trigger_streams,lags)



ExpParameters = [];
%% loop through subjects and find button presses
mean_coherences = [];
for subject = 1 : nS
    subjectID = SubjectList(subject);
    
    for condition = 1 :4
        combined = [];
       
        
        se_total = 6;
        if subjectID == 26
            
            se_total = 5 ;
            
        end
        for se = 1:se_total
            % select only stimstreams from all sessions that belong to specific
            % block
            
            % select all stim streams that belong to one subject
            stim_streams_sj = [];
            stim_streams_sj = stim_streams{subjectID,se}(:,condition);
            
            
            % select trigger streams that belong to one subject
            trigger_streams_sj = [];
            trigger_streams_sj = trigger_streams{subjectID,se}(:,condition);
            
            mean_streams = [];
            mean_streams = mean_stim_streams{subjectID,se}(:,condition);
            
            
            responses = all_responses((all_responses(:,9)== condition & all_responses(:,10) == se & all_responses(:,11) == subjectID),:);
            
            
            
            % find all triggers that lead to a button press
            
            
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
        
        
        SubjectIntegrationKernels.mean(:,condition,subject) = nanmean(combined,2);
        SubjectIntegrationKernels.sem(:,condition,subject) = nanstd(combined')/sqrt(size(SubjectIntegrationKernels.mean,1));
        num_false_alarms(subject,condition) = size(combined,2);
        
        
        
        
    end
    
    
end



GroupIntegrationKernels.sem = squeeze(std(permute(SubjectIntegrationKernels.mean,[3,1,2]))/sqrt(size(SubjectIntegrationKernels.mean,3)));
GroupIntegrationKernels.mean = nanmean(SubjectIntegrationKernels.mean,3);

%% multiple comparison 

repetitions = 100;
thres = round(.05 * repetitions);
anova_samples = nS;
plt = 1;

%[SignificantTimePoints,p,tbl] = shuffled_permutation_test_Fscore(SubjectIntegrationKernels.mean, lags, thres, repetitions, anova_samples, plt,nS);

SignificantTimePoints = [];
%% fitting of exponentials

options = optimset('MaxFunEvals',1000000,'MaxIter',100000);
for condition = 1 : 4
    for subject = 1 : size(SubjectIntegrationKernels.mean,3)
        data = SubjectIntegrationKernels.mean(:,condition,subject); % data to fit exp model to
        
        % find the peak of the data as starting point
        [val,idx_max] = max(data);
      
        data_new{condition,subject} = data(1:idx_max);
        % time steps
        dt = 0.01;
        num_steps = length(data_new{condition,subject});
        
        t_sj{condition,subject} = dt : dt : num_steps * dt;
        
        
        % initial param guesses for exp model
        pstart(1) = 1; % Amplitude
        pstart(2) = 1; % 1/tau
        % pstart(3) = 1; % offset
        
      fun = @(p)calculate_residuals_for_exponential_fit(p,data_new{condition,subject},t_sj{condition,subject}); % this is the correct cost function that works
      [ExpParameters.parameters(condition,:,subject),~,exitflag(condition,subject)] = fminsearch(fun,pstart,options);
        
    end
end
SubjectIntegrationKernels.DataForModelFit = data_new;
ExpParameters.time = t_sj;
end