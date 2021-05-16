for sj = 1:26
    
    data_sj = data{sj};
    % find frequent trials
    frequentTrials = []; 
    rareTrials = [];
    idx_con = data_sj.trialinfo(:,9) == 1 | data_sj.trialinfo(:,9) == 2;
    trialCountRare = 1;
    trialCountFreq = 1;
    for trial = 1:length(data_sj.trial)
        if idx_con(trial)
            
            frequentTrials(trialCountFreq,:) = zscore(data_sj.trial{trial}(63,:));
            trialCountFreq = trialCountFreq +1;
        else
            
            rareTrials(trialCountRare,:) = zscore(data_sj.trial{trial}(63,:));
            trialCountRare = trialCountRare + 1;
        end
        
        
        
    end
    
    meanfrequentTrials(sj,:) = nanmean(frequentTrials,1);
    meanrareTrials(sj,:) = nanmean(rareTrials,1);
end


    keyboard;
    figure
    subplot(1,2,1)
plot(data_sj.time{1},meanfrequentTrials)
title(['freq VEOG'])
   xlabel('time(s) 0 = button press')
 subplot(1,2,2)
plot(data_sj.time{1},meanrareTrials)
    title(['rare VEOG'])
      xlabel('time(s) 0 = button press')

    
    figure
    subplot(2,1,1)
    plot(data_sj.time{1},frequentTrials)
    xlabel('time(s) 0 = button press')
    title('VEOG for frequent trials')
    
    
   
    subplot(2,1,2)
    plot(data_sj.time{1},rareTrials)
    xlabel('time(s) 0 = button press')
    title('VEOG for rare trials')
    
    
    
keyboard;
badSRare = sum(isnan(rareTrials));
badSFreq = sum(isnan( frequentTrials));   

figure 
subplot(1,2,1)
plot(data_sj.time{1},badSRare,'-x')
title(['rare ', num2str(trialCountRare-1)])

subplot(1,2,2)
plot(data_sj.time{1},badSFreq,'-x')
title(['freq ',num2str(trialCountFreq-1)])
%     stdFreq = nanstd(frequentTrials(:));
%     stdRare = nanstd(rareTrials(:));
%     
%     FreqThresIdx = find(frequentTrials >= 2*stdFreq);
%     RareThresIdx = find(rareTrials >= 2*stdRare);
%     
%    figure (sj)
%    subplot(2,1,1)
%    vectorRare = [1:length(RareThresIdx(:,1))]';
%    plot(data_sj.time{1},vectorRare .* RareThresIdx,'b.')
%     subplot(2,1,2)
%    vectorFreq = [1:length(FreqThresIdx(:,1))]';
%    plot(data_sj.time{1},vectorFreq.*FreqThresIdx,'r.')
%     
