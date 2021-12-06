function [diffWaveLeftAvg diffWaveRightAvg] = calculate_diffWaves_tf_data_between_conditions(dataLeftBaseCorr, dataRightBaseCorr)



% average over frequencies for each subject and condition 

for subject = 1:19 
    
    
    for condition = 1:4
            cfg = [];
        
        
        % for left
        cfg.channel = {'C4' 'CP4'};
        cfg.avgoverchan = 'yes';
        cfg.frequency = [13 30];%[10  15];
        cfg.avgoverfreq = 'yes';
        cfg.nanmean = 'yes';
        
        [selectDataLeft{subject,condition}] = ft_selectdata(cfg, dataLeftBaseCorr{subject}{condition});
        
        
        
        % for right
        cfg.channel = {'C3' 'CP3'};
        cfg.avgoverchan = 'yes';
        cfg.frequency = [13 30];%[28  33];
        cfg.avgoverfreq = 'yes';
        cfg.nanmean = 'yes';
        [selectDataRight{subject, condition}] = ft_selectdata(cfg, dataRightBaseCorr{subject}{condition});
    
    end
    
end 

% average across conditions within subjects 

for subject = 1:19
    
    
        cfg = [];
        dataPerConditionLeft{subject, 1} = ft_freqgrandaverage(cfg,selectDataLeft{subject, 1:2});  
        dataPerConditionLeft{subject, 2} = ft_freqgrandaverage(cfg,selectDataLeft{subject, 3:4});
        dataPerConditionLeft{subject, 3} = ft_freqgrandaverage(cfg,selectDataLeft{subject, [1 3]});
        dataPerConditionLeft{subject, 4} = ft_freqgrandaverage(cfg,selectDataLeft{subject, [2 4]});
        
        cfg = [];
        dataPerConditionRight{subject, 1} = ft_freqgrandaverage(cfg,selectDataRight{subject, 1:2});  
        dataPerConditionRight{subject, 2} = ft_freqgrandaverage(cfg,selectDataRight{subject, 3:4});
        dataPerConditionRight{subject, 3} = ft_freqgrandaverage(cfg,selectDataRight{subject, [1 3]});
        dataPerConditionRight{subject, 4} = ft_freqgrandaverage(cfg,selectDataRight{subject, [2 4]});


% 
        cfg = []; 
        cfg.operation = 'x1-x2';
        cfg.parameter = 'powspctrm';
        diffWaveLeft{subject,1} = ft_math(cfg, dataPerConditionLeft{subject,1}, dataPerConditionLeft{subject,2}); 
        
        diffWaveRight{subject,1} = ft_math(cfg, dataPerConditionRight{subject, 1}, dataPerConditionRight{subject, 2}); 
        
        diffWaveLeft{subject,2} = ft_math(cfg, dataPerConditionLeft{subject,3}, dataPerConditionLeft{subject,4}); 
        
        diffWaveRight{subject,2} = ft_math(cfg, dataPerConditionRight{subject, 3}, dataPerConditionRight{subject, 4}); 
        



end 

% grandaverage across subjects 
cfg = [];
diffWaveLeftAvg{1} = ft_freqgrandaverage(cfg, diffWaveLeft{:,1});
diffWaveLeftAvg{2} = ft_freqgrandaverage(cfg, diffWaveLeft{:,2});
diffWaveRightAvg{1} = ft_freqgrandaverage(cfg, diffWaveRight{:,1});
diffWaveRightAvg{2} = ft_freqgrandaverage(cfg, diffWaveRight{:,2});
end 


