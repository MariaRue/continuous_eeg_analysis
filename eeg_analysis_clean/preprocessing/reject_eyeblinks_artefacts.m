function  D = reject_eyeblinks_artefacts(D,id,session,options,threshold)
  %1. grab the EOG data, and set anything exceeding 1000 microvolts
            %   to 0
%             VEOG_data = D(selectchannels(D,'VEOG'),:,:);
%             D(selectchannels(D,'VEOG'),find(abs(VEOG_data)>1000),:) = 0;
%             D.save;
%             
%             %2. use spm_eeg_artefact_eyeblink to identify eyeblinks and add
%             %them as events.
%             
%             % get threshold
%             
%             S = [];
%             S.D = D;
%             S.chanind = selectchannels(D,'VEOG'); %VEOG channel in the
%             S.threshold = threshold; %standard deviations away from the me an
%             
%             S.mode = 'mark';
%             S.append = 1; %this adds the eyeblinks as events
%             S.excwin = options.preproc.eyeblink.excwin; %remove 500ms around each event
%             fprintf('Subject %0.0f - run %0.0f\n',id,session);
%             
%             D = spm_eeg_artefact_eyeblink(S);
%             
%             
%             D.save; %stores the eyeblink events
            %
        
            S = [];
            S.D = D;
            S.chanind = selectchannels(D,'EEG');
            S.threshold = options.preproc.artefact.threshold;
            S.badchanthresh = options.preproc.artefact.badchanthresh; % possibly number of artefacts/length of samples that lead to exclusion of whole channel?
            S.mode = 'mark';
            S.append = 1;
            S.excwin = options.preproc.artefact.excwin;
            
            D = spm_eeg_artefact_threshchan(S);
            bs = sum(D.badsamples([selectchannels(D,'CPZ') selectchannels(D,'VEOG')],:,1))>=1; %find artefacts in any channels
            mean(bs)
%             figure (3)
%             if sum(bs > 0)
%             plot(D.time,D(40,:,1)); hold on; plot(D.time(find(bs==1)), D(40,find(bs==1),1),'r.');
%             end 
            D.save; %stores the events above 100 µm Volt
            %close all

end 