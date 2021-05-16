function D = rereference_EEG(D,reference)
% rereference the electrodes and remove artefacts including
        % eyeblinks
        
        %%%% 1) reference - choose between LM-RM to reference to average
        % of left and right mastoids or average_reference which
        % rereferences to average of electrodes across scalp
        
            S = [];
            S.D = D;
            
            S.mode = 'write';
            
            
            S.keepothers = 1; % to keep the EOG channels
            switch reference
                
                
                case 'LMRM' % rereference the 61 EEG channels to L+R mastoid average
                    
                    
                    %in our montage we have 61 EEG channels, plus the right mastoid:
                    S.montage.labelorg = D.chanlabels(1:62);
                    %we keep all EEG channels, but throw out the right mastoid:
                    S.montage.labelnew = D.chanlabels(1:61);
                    
                    %build our M*N matrix for montaging:
                    S.montage.tra = eye(62);
                    %subtract 0.5 of the right mastoid in the re-reference - see
                    %p.108 of Luck book!
                    S.montage.tra(:,62) = -0.5;
                    %get rid of the right mastoid in the new channels:
                    S.montage.tra(62,:) = [];
                    
                    
                    %add in extra electrode for lateralised readiness potential (C4-C3)
                    S.montage.tra(62,:) = S.montage.tra(selectchannels(D,'C4'),:) - S.montage.tra(selectchannels(D,'C3'),:);
                    S.montage.labelnew{62} = 'C4_C3_LRP';
                    
                case 'average_reference'
                   
                    %in our montage we have 61 EEG channels. We can ignore
                    %the right mastoid, as it isn't in our equation.
                    S.montage.labelorg = D.chanlabels(1:61);
                    %we keep all EEG channels
                    S.montage.labelnew = D.chanlabels(1:61);
                    
                    %build our M*N matrix for montaging:
                    S.montage.tra = eye(61)-(1/61);
                    
                    %add in extra electrode for lateralised readiness potential (C4-C3)
                    S.montage.tra(62,:) = S.montage.tra(selectchannels(D,'C4'),:) - S.montage.tra(selectchannels(D,'C3'),:);
                    S.montage.labelnew{62} = 'C4_C3_LRP';
                    
                otherwise
                    error('unrecognised re-referencing type')
            end
            
            D = spm_eeg_montage(S);
            %set the LRP electrode to be EEG chantype
            D = chantype(D,selectchannels(D,'C4_C3_LRP'),'EEG');
            D.save();
            
end 




