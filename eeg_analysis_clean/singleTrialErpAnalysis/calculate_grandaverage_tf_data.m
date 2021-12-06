function [grandavgLeft, grandavgRight, grandAvgDiffWave] = calculate_grandaverage_tf_data(dataLeftBaseCorr, dataRightBaseCorr, diffWaveLeftRight, per_condition, for_diffWaves)


grandavgLeft = [];
grandavgRight = [];
grandAvgDiffWave =[];

if per_condition
    
    for condition = 1:4
        
        cfg = [];
        if for_diffWaves
            
              [grandAvgDiffWave{condition}] = ft_freqgrandaverage(cfg,  diffWaveLeftRight{:,condition});

            
        else
            
            for subject = 1:19
                
                dataLeftPerCondition{subject} = dataLeftBaseCorr{subject}{condition};
                dataRightPerCondition{subject} = dataRightBaseCorr{subject}{condition};
            end 
            
            
            
            [grandavgLeft{condition}] = ft_freqgrandaverage(cfg, dataLeftPerCondition{:});

            [grandavgRight{condition}] = ft_freqgrandaverage(cfg, dataRightPerCondition{:});
            
        end
        
        
    end
    
    
else
    cfg = [];
    
    if for_diffWaves
        
      [grandAvgDiffWave] = ft_freqgrandaverage(cfg,  diffWaveLeftRight{:});
   
        
    else
        
    [grandavgLeft] = ft_freqgrandaverage(cfg, dataLeftBaseCorr{:});

    [grandavgRight] = ft_freqgrandaverage(cfg, dataRightBaseCorr{:});
    
    end
    
end


end