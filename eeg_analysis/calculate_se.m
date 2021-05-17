function [se] = calculate_se(variance,nS)

standard_deviation = sqrt(variance);
se = standard_deviation./sqrt(nS);



end 