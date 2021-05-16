options.preproc.eyeblinktreatment   = 'ssp'; % 'reject', 'ssp'
options.preproc.eyeblinkchannels    = {'VEOG'};
options.preproc.eyeblinkthreshold   = 5; % for SD thresholding: in standard deviations, for amp in uV
options.preproc.windowForEyeblinkdetection = 3; % first event of interest (and optionally last)

options.preproc.eyeblinkmode        = 'eventbased'; % uses EEG triggers for trial onsets
options.preproc.eyeblinkwindow      = 0.5; % in s around blink events
options.preproc.eyeblinktrialoffset = 0.1; % in s: EBs won't hurt <100ms after tone onset
options.preproc.eyeblinkEOGchannel  = 'VEOG'; % EOG channel (name/idx) to plot
options.preproc.eyebadchanthresh    = 0.4; % prop of bad trials due to EBs

options.preproc.eyeconfoundcomps    = 1;
options.preproc.eyecorrectionchans  = {'Fp1', 'Fz', 'AF8', 'T7', 'Oz'};
options.preproc.preclean.doFilter           = true;
options.preproc.preclean.lowPassFilterFreq  = 10;
options.preproc.preclean.doBadChannels      = false;
options.preproc.preclean.doBadTrials        = true;
options.preproc.preclean.doRejection        = true;
options.preproc.preclean.badtrialthresh     = 500;
options.preproc.preclean.badchanthresh      = 0.5;
options.preproc.preclean.rejectPrefix       = 'cleaned_';
