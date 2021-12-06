function stats = rm_Anova2_for_false_alarm_rate(FalseAlarmRatePerSubject)

% rm Analysis 
% rm_anvoa2 function code taken from https://de.mathworks.com/matlabcentral/fileexchange/6874-two-way-repeated-measures-anova?s_tid=mwa_osa_a

Y = FalseAlarmRatePerSubject'; 
Y = Y(:);

count = 1; 
for subject = 1:size(FalseAlarmRatePerSubject,1)
    
    S(count : count+3,1) = subject; 
    count = count + 4; 
end 

F1 = zeros(size(FalseAlarmRatePerSubject,1)*4,1);
F2 = zeros(size(FalseAlarmRatePerSubject,1)*4,1);

F1(1:4:end) = 1; % frequent
F1(2:4:end) = 1; 
F1(3:4:end) = 2; % rare
F1(4:4:end) = 2; 

F2(1:4:end) = 1; % short
F2(2:4:end) = 2; % long
F2(3:4:end) = 1; 
F2(4:4:end) = 2; 
FACTNAMES = {'frequency', 'length'};
stats = rm_anova2(Y,S,F1,F2,FACTNAMES);

end 