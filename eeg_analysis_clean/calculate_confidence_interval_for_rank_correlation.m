function [upperBound, lowerBound] = calculate_confidence_interval_for_rank_correlation(RHO, n, confidenceInterval)

%RHO = corr(a.',b.','Type','Spearman');
%n = numel(a);
STE = 1/sqrt(n-3);
% here the input is 95% confidence interval, for 99% use 0.99:
CI = norminv(confidenceInterval); 
upperBound = tanh(atanh(RHO)+CI*STE);
lowerBound = tanh(atanh(RHO)-CI*STE);






end 