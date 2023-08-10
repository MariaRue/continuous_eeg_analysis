function D = beta_power_EEG(D)

            S = [];
            S.D = D;
            S.method = 'morlet';
            S.frequencies = 20;
            D = spm_eeg_tf(S);


end 