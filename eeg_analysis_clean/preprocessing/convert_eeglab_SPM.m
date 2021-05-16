function D = convert_eeglab_SPM(pathSetFile,pathPreprocessed,SetFile,scriptdir)
cd (pathSetFile) 

 S = [];
            
            S.dataset = SetFile;
            
            
            S.mode = 'continuous';
            D = spm_eeg_convert(S);
            D = move(D,pathPreprocessed);

%cd (scriptdir)
end