function D =  downsample_EEG(D,options) 

            S = [];
            S.D = D;
            S.fsample_new = options.preproc.fsample;
            D = spm_eeg_downsample(S);
            
end 