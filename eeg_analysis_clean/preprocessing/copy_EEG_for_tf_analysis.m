function D = copy_EEG_for_tf_analysis(D,options)

S = [];
S.D = D;
S.outfile = ['for_tf_analysis' D.fname]; %append 'for_tf_analysis' at start of trial name to mark it for tf analysis in contrast to orignal analysis
D = spm_eeg_copy(S);


% S = [];
% S.D = D;
% S.band = 'bandpass';
% S.freq = options.preproc.bandpass_for_tf_analysis;
% D = spm_eeg_filter(S);

end