function [statsGroupLevelSummary] = regression_analysis_for_median_rts(options, all_responses)

for subject = 1:options.totalNumberofSubjects
    
    % get rts for regression < 3.5s
    idxRt = all_responses(:,7) == 1 & all_responses(:,11) == subject & all_responses(:,2)<= 3.5;
    rts  = all_responses(idxRt,2);
    
    % define Regressors for coherence, Length, frequency and interaction
    % lengthXFrequency 
    
    idxRtLong = all_responses(idxRt,9) == 2 |  all_responses(idxRt,9) == 4;
    idxRtShort = all_responses(idxRt,9) == 1 |  all_responses(idxRt,9) == 3;
    
    idxRtFreq = all_responses(idxRt,9) == 1 |  all_responses(idxRt,9) == 2;
    idxRtRare = all_responses(idxRt,9) == 3 |  all_responses(idxRt,9) == 4;
    
    LenghtRegressor = zeros(length(rts),1);
    LenghtRegressor(idxRtLong) = -1;
    LenghtRegressor(idxRtShort) = 1;
    
    frequencyRegressor = zeros(length(rts),1);
    frequencyRegressor(idxRtRare) = -1;
    frequencyRegressor(idxRtFreq) = 1;
    
    interactionFrequencyXLengthRegressor = frequencyRegressor .* LenghtRegressor;
    
    coherenceRegressor = abs(all_responses(idxRt,4));
    
    
    % remove nan Rts
    nanRt = isnan(rts);
    
    rts(nanRt) = [];
    coherenceRegressor(nanRt) = [];
    LenghtRegressor(nanRt) = [];
    frequencyRegressor(nanRt) = [];
    interactionFrequencyXLengthRegressor(nanRt) = [];
    
    % demean coherence regressor
    coherenceRegressor = coherenceRegressor - mean(coherenceRegressor);
    
    Y = rts; % dependent
    X = [coherenceRegressor LenghtRegressor  frequencyRegressor interactionFrequencyXLengthRegressor];
    
    % regression subject Level
    [betaRT(:,subject),~,~] = glmfit(X,Y,'normal','constant','off');
    
end


% ttest for group level 
for regressor = 1:size(X,2)
    [~,PValue(regressor),~,statsGroupLevel(regressor)] = ttest(betaRT(regressor,:));
end


% PValues and tStats structure for each regressor

statsGroupLevelSummary.coherence.pValue = PValue(1);
statsGroupLevelSummary.coherence.tStatistic = statsGroupLevel(1).tstat;

statsGroupLevelSummary.length.pValue = PValue(2);
statsGroupLevelSummary.length.tStatistic = statsGroupLevel(2).tstat;

statsGroupLevelSummary.frequency.pValue = PValue(3);
statsGroupLevelSummary.frequency.tStatistic = statsGroupLevel(3).tstat;

statsGroupLevelSummary.lengthXFrequency.pValue = PValue(4);
statsGroupLevelSummary.lengthXFrequency.tStatistic = statsGroupLevel(4).tstat;




end