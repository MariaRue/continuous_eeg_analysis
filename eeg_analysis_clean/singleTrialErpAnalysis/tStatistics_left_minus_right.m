function [stat] =  tStatistics_left_minus_right(dataLeft, dataRight, subjectNumber, per_condition)

if per_condition 
    
   for condition = 1:4 
       
       for subject = 1:subjectNumber 
           
            dataLeftPerCondition{subject} = dataLeft{subject}{condition};
            
            dataRightPerCondition{subject} = dataRight{subject}{condition};
           
       end 
       
       cfg = [];
  cfg.keepindividual = 'yes';
 [grandavgLeft] = ft_freqgrandaverage(cfg, dataLeftPerCondition{:});

 [grandavgRight] = ft_freqgrandaverage(cfg, dataRightPerCondition{:});
 
 
cfg = [];
cfg.method = 'analytic'; 
cfg.statistic = 'ft_statfun_depsamplesT';

 subj = subjectNumber;
design = zeros(2,2*subj);
for i = 1:subj
  design(1,i) = i;
end
for i = 1:subj
  design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;

cfg.design   = design;
cfg.uvar     = 1;
cfg.ivar     = 2;





 [stat{condition}] = ft_freqstatistics(cfg, grandavgLeft, grandavgRight);
       
   end 
    
else 

% transform data into grandavg structure for stats analysis
cfg = [];
  cfg.keepindividual = 'yes';
 [grandavgLeft] = ft_freqgrandaverage(cfg, dataLeft{:});

 [grandavgRight] = ft_freqgrandaverage(cfg, dataRight{:});
 
 
cfg = [];
cfg.method = 'analytic'; 
cfg.statistic = 'ft_statfun_depsamplesT';

 subj = subjectNumber;
design = zeros(2,2*subj);
for i = 1:subj
  design(1,i) = i;
end
for i = 1:subj
  design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;

cfg.design   = design;
cfg.uvar     = 1;
cfg.ivar     = 2;





 [stat] = ft_freqstatistics(cfg, grandavgLeft, grandavgRight);
 

end




end 