for sj = 1:26
    
%     data_sj = data{sj};
%     % find frequent trials
%     
%     idx_con = data_sj.trialinfo(:,9) == 1 | data_sj.trialinfo(:,9) == 2;
%     trialCountRare = 1;
%     trialCountFreq = 1;
%     for trial = 1:length(data_sj.trial)
%         if idx_con(trial)
%             
%             frequentTrials(trialCountFreq,:) = zscore(data_sj.trial{trial}(40,:));
%             trialCountFreq = trialCountFreq +1;
%         else
%             
%             rareTrials(trialCountRare,:) = zscore(data_sj.trial{trial}(40,:));
%             trialCountRare = trialCountRare + 1;
%         end
%         
%         
%         
%     end
%     
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
    
end