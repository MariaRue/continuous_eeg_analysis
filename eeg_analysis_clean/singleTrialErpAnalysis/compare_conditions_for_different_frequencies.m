function [dataPerConditionLeft, dataPerConditionRight] = compare_conditions_for_different_frequencies(dataLeft, dataRight, ifDiffWave, conditions4)

if ifDiffWave
    
    
else
    
    for condition = 1:4
        cfg = [];
        
        
        % for left
        cfg.channel = {'C4' 'CP4'};
        cfg.avgoverchan = 'yes';
        cfg.frequency = [13 30];%[10  15];
        cfg.avgoverfreq = 'yes';
        cfg.nanmean = 'yes';
        
        [selectDataLeft{condition}] = ft_selectdata(cfg, dataLeft{condition});
        
        
        
        % for right
        cfg.channel = {'C3' 'CP3'};
        cfg.avgoverchan = 'yes';
        cfg.frequency = [13 30];%[28  33];
        cfg.avgoverfreq = 'yes';
        cfg.nanmean = 'yes';
        [selectDataRight{condition}] = ft_selectdata(cfg, dataRight{condition});
        
        
        
    end
    
    
    if conditions4
        
        dataPerConditionRight = selectDataRight;
        dataPerConditionLeft = selectDataLeft;
        
    else
        
        cfg = [];
        dataPerConditionLeft{1} = ft_freqgrandaverage(cfg,selectDataLeft{1:2});  
        dataPerConditionLeft{2} = ft_freqgrandaverage(cfg,selectDataLeft{3:4});
        dataPerConditionLeft{3} = ft_freqgrandaverage(cfg,selectDataLeft{[1 3]});
        dataPerConditionLeft{4} = ft_freqgrandaverage(cfg,selectDataLeft{[2 4]});
        
        cfg = [];
        dataPerConditionRight{1} = ft_freqgrandaverage(cfg,selectDataRight{1:2});  
        dataPerConditionRight{2} = ft_freqgrandaverage(cfg,selectDataRight{3:4});
        dataPerConditionRight{3} = ft_freqgrandaverage(cfg,selectDataRight{[1 3]});
        dataPerConditionRight{4} = ft_freqgrandaverage(cfg,selectDataRight{[2 4]});
        
    end
    
    
    
end


end