function [kernelsLate, kernelsEarly] = get_kernels(options,stimStream,triggersSorted) 

kernelsLate = nan(options.lags+1,1); 
kernelsEarly = nan(options.lags+1,1); 
countEarly = 1; 
countLate = 1; 

for triggers = 1:length(triggersSorted(:,1))
    
    if triggersSorted(triggers,1) == 1
        
        kernelsLate(:,countLate) = stimStream(triggersSorted(triggers,2) - options.lags : triggersSorted(triggers,2)); % to be able to combine with right button presses * -1
        
        countLate = countLate+1;
    else
        kernelsEarly(:,countEarly) = stimStream(triggersSorted(triggers,2) - options.lags : triggersSorted(triggers,2)); % to be able to combine with right button presses * -1
        countEarly = countEarly+1;
    end
    
end

end 