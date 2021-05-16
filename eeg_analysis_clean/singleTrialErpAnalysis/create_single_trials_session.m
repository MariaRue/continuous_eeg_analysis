function [data] = create_single_trials_session(responses,meanStimSession,options,timelockCondition,dataName)


            cfg = [];
            cfg.dataset = dataName;
            cfg.subject_responses = responses; 
            cfg.mean_stim_streams = meanStimSession; %mean_stim_streams_org{stream_sj,i};
            cfg.trialdef.prestim = options.singleTrial.(timelockCondition).preStim; %7;
            cfg.trialdef.poststim = options.singleTrial.(timelockCondition).postStim; %4;
            cfg.timelock_event = options.singleTrial.(timelockCondition).name; %'button press';
              
            cfg.trialfun  = 'ft_trialfun_continuous_eeg';
            cfg = ft_definetrial(cfg);
            
            
            data = ft_preprocessing(cfg);
            
end 