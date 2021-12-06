function [dataLeft, dataRight] = select_trials_for_tf_analysis(data, per_condition)
% left trials

if per_condition
    
    
        % left trials + condition
for condition = 1:4
    
   cfg = [];

cfg.trials =  data.trialinfo(:,5) == 1 & sign(data.trialinfo(:,4)) == -1 & data.trialinfo(:,9) == condition; % correct responses to the left
cfg.avgoverrpt = 'yes';
cfg.nanmean = 'yes';

dataLeft{condition} = ft_selectdata(cfg,data); 



% right trials + condition

cfg.trials = data.trialinfo(:,5) == 1 & sign(data.trialinfo(:,4)) == 1 & data.trialinfo(:,9) == condition; % correct responses to the right 

dataRight{condition} = ft_selectdata(cfg,data);  
    
end 




else 

cfg = [];

cfg.trials =  data.trialinfo(:,5) == 1 & sign(data.trialinfo(:,4)) == -1; % correct responses to the left
cfg.avgoverrpt = 'yes';
cfg.nanmean = 'yes';

dataLeft = ft_selectdata(cfg,data); 



% right trials

cfg.trials = data.trialinfo(:,5) == 1 & sign(data.trialinfo(:,4)) == 1; % correct responses to the right 

dataRight = ft_selectdata(cfg,data); 

end 

end