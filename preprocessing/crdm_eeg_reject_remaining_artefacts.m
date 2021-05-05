function  D = crdm_eeg_reject_remaining_artefacts( D, options )
%CRDM_EEG_REJECT_REMAINING_ARTEFACTS Applies a threshold to all EEG channels to 
%detect remaining artefacts (after steps like filtering, eyeblink
%correction, etc. have been performed). Marks these events as bad.
%   IN:     D
%           options
%   OUT:    D

S = [];
S.D = D;

S.chanind = selectchannels(D,'EEG'); % use all channels
S.threshold = options.preproc.artefact.threshold; % threshold in muV
S.badchanthresh = options.preproc.artefact.badchanthresh; % percent of trials 
% marked as bad that leads to exclusion of a whole channel

S.mode = 'mark';
S.append = 1;
S.excwin = options.preproc.artefact.excwin;

D = spm_eeg_artefact_threshchan(S);
%bs = sum(D.badsamples([selectchannels(D,'CPZ') selectchannels(D,'VEOG')],:,1))>=1; %find artefacts in any channels
%mean(bs)
% figure (3)
% if sum(bs > 0)
% 	plot(D.time,D(40,:,1)); hold on; plot(D.time(find(bs==1)), D(40,find(bs==1),1),'r.');
% end 
D.save; % events above the threshold are marked as bad in D

%close all

end 