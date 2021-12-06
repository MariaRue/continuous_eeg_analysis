function [diffWaveLeftRight] = calculate_diff_wave_tf_data(dataLeftBaseCorr,dataRightBaseCorr,for_different_conditions)

if for_different_conditions
    for condition = 1:4
        
        for subject = 1:19
            
            
            dataLeftPerCondition = dataLeftBaseCorr{subject}{condition};
            
            dataRightPerCondition = dataRightBaseCorr{subject}{condition};
            
             % difference wave left vs right for each condition
        cfg = []; 
        cfg.operation = '(x1-x2)/(x1+x2)';
        cfg.parameter = 'powspctrm';
        diffWaveLeftRight{subject,condition} = ft_math(cfg, dataLeftPerCondition, dataRightPerCondition); 
        
            
        end

        
    end
    
    
else
    for subject = 1:19
      cfg = []; 
        cfg.operation = '(x1-x2)/(x1+x2)';
        cfg.parameter = 'powspctrm';
        diffWaveLeftRight{subject} = ft_math(cfg, dataLeftBaseCorr{subject}, dataRightBaseCorr{subject}); 
        
    end 
    
    
    
end





end