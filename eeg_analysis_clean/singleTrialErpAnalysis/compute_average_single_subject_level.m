function [data_avg] =  compute_average_single_subject_level(data,average_type, rt_flag, baseline_window,baseline)%


coherence = [0.3,0.4,0.5];


switch average_type
    
    case 'coherence'
        for coh = 1:3
            if rt_flag
                
                idx_coh = find(abs(data.trialinfo(:,4)) == coherence(coh) & data.trialinfo(:,7) == 1 & data.trialinfo(:,2) <= 3.5);
                
            else
                idx_coh = find(abs(data.trialinfo(:,4)) == coherence(coh) & data.trialinfo(:,7) == 1);
                
            end
            
            
            cfg = [];
            cfg.trials = idx_coh;
            cfg.channel = 'eeg';
            data_coherence{coh} = ft_selectdata(cfg,data);
            
            try
                cfg = [];
                data_avg{coh} = ft_timelockanalysis(cfg,data_coherence{coh});
            catch
                
                keyboard;
            end
            
            if baseline
                        cfg = [];
                        cfg.baseline = baseline_window;
            
                        data_avg{coh} = ft_timelockbaseline(cfg, data_avg{coh});
            end
        end
        
    case 'condition'
        
        for con = 1:4
            
            if rt_flag
                idx_con = find(data.trialinfo(:,9) == con & data.trialinfo(:,7) == 1 & data.trialinfo(:,2) <= 3.5);
                
            else
                
                idx_con = find(data.trialinfo(:,9) == con & data.trialinfo(:,7) == 1);
            end
            cfg = [];
            cfg.trials = idx_con;
            % cfg.channel = 'eeg';
            data_condition{con} = ft_selectdata(cfg,data);
            
            cfg = [];
            data_avg{con} = ft_timelockanalysis(cfg,data_condition{con});
            
            if baseline
                        cfg = [];
                        cfg.baseline = baseline_window;
            
                        data_avg{con} = ft_timelockbaseline(cfg,data_avg{con});
            end
            
            
        end
        
        
        
    case 'false alarm'
        
        for con = 1:4
            
            idx_con = find(data.trialinfo(:,9) == con & data.trialinfo(:,7) == 2);
            
            
            cfg = [];
            cfg.trials = idx_con;
            %cfg.channel = 'eeg';
            data_condition{con} = ft_selectdata(cfg,data);
            
            cfg = [];
            data_avg{con} = ft_timelockanalysis(cfg,data_condition{con});
            
            if baseline
                        cfg = [];
                        cfg.baseline = baseline_window;
                         cfg.channel = 'eeg';
                     data_avg{con} = ft_timelockbaseline(cfg,data_avg{con});
            end 
            end
        
    case 'perm_test_rare_freq_fa'
        
        
        for con = 1:2
            
            conditions_perm = [1 2; 3 4];
            idx_con = find((data.trialinfo(:,9) == conditions_perm(con,1) | data.trialinfo(:,9) == conditions_perm(con,2))& data.trialinfo(:,7) == 2);
            
            
            cfg = [];
            cfg.trials = idx_con;
            %cfg.channel = 'eeg';
            data_condition{con} = ft_selectdata(cfg,data);
            
            cfg = [];
            data_avg{con} = ft_timelockanalysis(cfg,data_condition{con});
            if baseline
                cfg = [];
                cfg.baseline = baseline_window;
                
                data_avg{con} = ft_timelockbaseline(cfg, data_avg{con});
            end
        end
        
        
    case 'perm_test_rare_freq_trial'
        
        conditions_perm = [1 2; 3 4];
        
        for con = 1:2
            if rt_flag
                idx_con = find((data.trialinfo(:,9) == conditions_perm(con,1) | data.trialinfo(:,9) == conditions_perm(con,2))& data.trialinfo(:,7) == 1 & data.trialinfo(:,2) <= 3.5);
            else
                idx_con = find((data.trialinfo(:,9) == conditions_perm(con,1) | data.trialinfo(:,9) == conditions_perm(con,2))& data.trialinfo(:,7) == 1);
            end
            
            cfg = [];
            cfg.trials = idx_con;
            % cfg.channel = 'eeg';
            data_condition{con} = ft_selectdata(cfg,data);
            
            cfg = [];
            data_avg{con} = ft_timelockanalysis(cfg,data_condition{con});
            
            if baseline
                cfg = [];
                cfg.baseline = baseline_window;
                
                data_avg{con} = ft_timelockbaseline(cfg,data_avg{con});
            end
        end
        
        
    case 'perm_test_short_long_trial'
        
        conditions_perm = [1 3; 2 4];
        
        for con = 1:2
            if rt_flag
                idx_con = find((data.trialinfo(:,9) == conditions_perm(con,1) | data.trialinfo(:,9) == conditions_perm(con,2))& data.trialinfo(:,7) == 1 & data.trialinfo(:,2) <= 3.5);
            else
                idx_con = find((data.trialinfo(:,9) == conditions_perm(con,1) | data.trialinfo(:,9) == conditions_perm(con,2))& data.trialinfo(:,7) == 1);
            end
            
            cfg = [];
            cfg.trials = idx_con;
            cfg.channel = 'eeg';
            data_condition{con} = ft_selectdata(cfg,data);
            
            cfg = [];
            data_avg{con} = ft_timelockanalysis(cfg,data_condition{con});
            
            if baseline
                
                cfg = [];
                cfg.baseline = baseline_window;
                
                data_avg{con} = ft_timelockbaseline(cfg,data_avg{con});
                
            end
            
        end
        
        
    case 'perm_test_short_long_fa'
        
        
        for con = 1:2
            
            conditions_perm = [1 3; 2 4];
            idx_con = find((data.trialinfo(:,9) == conditions_perm(con,1) | data.trialinfo(:,9) == conditions_perm(con,2))& data.trialinfo(:,7) == 2);
            
            
            cfg = [];
            cfg.trials = idx_con;
            cfg.channel = 'eeg';
            data_condition{con} = ft_selectdata(cfg,data);
            
            cfg = [];
            data_avg{con} = ft_timelockanalysis(cfg,data_condition{con});
            if baseline
                cfg = [];
                cfg.baseline = baseline_window;
                
                data_avg{con} = ft_timelockbaseline(cfg, data_avg{con});
            end
        end
        
        
end


end






