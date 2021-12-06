function [dataLeftBaseCorr, dataRightBaseCorr] = baseline_correct_tf_data(dataLeft, dataRight, per_condition) 

if per_condition

     for condition = 1:4
  cfg = []; 
cfg.baseline     = [-4 -3];
cfg.baselinetype = 'db';

 [dataLeftBaseCorr{condition}] = ft_freqbaseline(cfg, dataLeft{condition});
 
 [dataRightBaseCorr{condition}] = ft_freqbaseline(cfg, dataRight{condition});
     end 
    
else 
    
    cfg = []; 
cfg.baseline     = [-4 -3];
cfg.baselinetype = 'db';

 [dataLeftBaseCorr] = ft_freqbaseline(cfg, dataLeft);
 
 [dataRightBaseCorr] = ft_freqbaseline(cfg, dataRight);

end 
end 