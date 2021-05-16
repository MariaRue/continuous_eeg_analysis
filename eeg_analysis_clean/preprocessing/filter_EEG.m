function D = filter_EEG(D,options)


            S = [];
            S.D = D;
            S.band = 'bandpass';
            S.freq = options.preproc.bandpass;
            D = spm_eeg_filter(S);


end 