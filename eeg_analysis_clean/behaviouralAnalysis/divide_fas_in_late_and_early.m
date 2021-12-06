function [falseAlarmsSorted] = divide_fas_in_late_and_early(triggers,MeanLengthOfInterval,StartOfIntertrialIntervals,EndOfIntertrialIntervals)


for falseAlarm = 1:length(triggers)
                
                % find interval
                
                IntervalStartAlarm = StartOfIntertrialIntervals(triggers(falseAlarm) < EndOfIntertrialIntervals & triggers(falseAlarm) > StartOfIntertrialIntervals);
                
                % early or late?
       
                
                if triggers(falseAlarm) > (IntervalStartAlarm + MeanLengthOfInterval)
                    
                    falseAlarmsSorted(falseAlarm,1) = 1;
                    falseAlarmsSorted(falseAlarm,2) = triggers(falseAlarm);

                else
                    
                    falseAlarmsSorted(falseAlarm,1) = 0;
                    falseAlarmsSorted(falseAlarm,2) = triggers(falseAlarm);
                end
  
                
end
             



end