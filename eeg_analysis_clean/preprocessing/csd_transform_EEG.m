function D = csd_transform_EEG(D)
 S = [];
            S.D = D;
            S.outfile = ['csd' D.fname]; %append 'csd' on start of filename for CSD transformation
            
            D = spm_eeg_copy(S);
            
            %1) transform data into fieldtrip format
            ft_data = D.ftraw();
            
            %2 source density analysis
            
            % remove non-EEG electrodes as ft_scalpcurrentdensity doesn't
            % use these
            cfg = [];

            cfg.channel = {'all','-RM', '-VEOG', '-HEOG', '-LM', '-C4_C3_LRP'};
            [data_pre] = ft_preprocessing(cfg, ft_data);
            
            cfg = [];
            cfg.elec = ft_data.elec;
            cfg.degree = 14;
            source_density_data =  ft_scalpcurrentdensity(cfg, data_pre);
            
            D(1:61,:,1) = source_density_data.trial{1}; %put the source density transformed data into cfMdspmeeg...
            D.save;



end 