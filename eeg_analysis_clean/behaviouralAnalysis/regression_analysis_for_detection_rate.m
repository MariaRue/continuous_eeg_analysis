function [statsGroupLevelSummary] = regression_analysis_for_detection_rate(options,detectRateSubject)


for subject = 1:options.totalNumberofSubjects
    
   
Y = detectRateSubject(:,:,subject);
Y = Y(:); % dependent variable

% Regressors 

coherenceRegressor = zeros(length(Y),1);
coherenceRegressor(1:3:end) = 0.3; % mean ... coherence 
coherenceRegressor(2:3:end) = 0.4; 
coherenceRegressor(3:3:end) = 0.5; 

lengthRegressor = zeros(length(Y),1);
lengthRegressor(1:6:end) = 1; % short
lengthRegressor(2:6:end) = 1; 
lengthRegressor(3:6:end) = 1; 
lengthRegressor(4:6:end) = -1; % long
lengthRegressor(5:6:end) = -1; 
lengthRegressor(6:6:end) = -1; 

frequencyRegressor = zeros(length(Y),1);
frequencyRegressor(1:12:end) = 1; % frequent 
frequencyRegressor(2:12:end) = 1; 
frequencyRegressor(3:12:end) = 1; 
frequencyRegressor(4:12:end) = 1; 
frequencyRegressor(5:12:end) = 1; 
frequencyRegressor(6:12:end) = 1; 

frequencyRegressor(7:12:end) = -1; % rare 
frequencyRegressor(8:12:end) = -1; 
frequencyRegressor(9:12:end) = -1; 
frequencyRegressor(10:12:end) = -1; 
frequencyRegressor(11:12:end) = -1; 
frequencyRegressor(12:12:end) = -1; 



LenghtFrequencyInteractionRegressor = frequencyRegressor .* lengthRegressor; 


X = [coherenceRegressor, frequencyRegressor, lengthRegressor, LenghtFrequencyInteractionRegressor];
% single subject level GLM
[beta(:,subject),~,statsGLM(subject)] = glmfit(X,Y,'normal');
end 

% Group Level ttest across betas
for regressor = 1:5
    [~,PValue(regressor),~,statsGroupLevel(regressor)] = ttest(beta(regressor,:));
end 

% PValues and tStats structure for each regressor
statsGroupLevelSummary.constant.pValue = PValue(1);
statsGroupLevelSummary.constant.tStatistic = statsGroupLevel(1).tstat;

statsGroupLevelSummary.coherence.pValue = PValue(2);
statsGroupLevelSummary.coherence.tStatistic = statsGroupLevel(2).tstat;

statsGroupLevelSummary.frequency.pValue = PValue(3);
statsGroupLevelSummary.frequency.tStatistic = statsGroupLevel(3).tstat;

statsGroupLevelSummary.length.pValue = PValue(4);
statsGroupLevelSummary.length.tStatistic = statsGroupLevel(4).tstat;

statsGroupLevelSummary.lengthXFrequency.pValue = PValue(5);
statsGroupLevelSummary.lengthXFrequency.tStatistic = statsGroupLevel(5).tstat;


end 
