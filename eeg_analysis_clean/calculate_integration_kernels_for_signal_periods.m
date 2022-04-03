function [GroupIntegrationKernels, SubjectIntegrationKernels, SignificantTimePoints] = calculate_integration_kernels_for_signal_periods(all_responses,SubjectList,nS, mean_stim_streams, stim_streams, trigger_streams,lags)

ExpParameters = 0;
%% loop through subjects and find button presses

for subject = 1 : nS
    subjectID = SubjectList(subject);
    
    for condition = 1 :4
        combined = [];
        combined_coh = [];
        
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




repetitions = 100;
thres = round(.05 * repetitions);
anova_samples = nS;
plt = 1;

[SignificantTimePoints,p,tbl] = shuffled_permutation_test_Fscore(SubjectIntegrationKernels.mean, lags, thres, repetitions, anova_samples, plt,nS);




end
