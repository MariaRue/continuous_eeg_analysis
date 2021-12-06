function [meanKernelLate, meanKernelEarly] = extract_kernels_from_stimStreams(stimStream,options,triggersSortedLeft,triggersSortedRight)

kernelsEarlyLeft = []; 
kernelsEarlyRight = []; 
kernelsLateLeft = []; 
kernelsLateRight = []; 
if ~isnan(triggersSortedRight) 
[kernelsLateRight, kernelsEarlyRight] = get_kernels(options,stimStream,triggersSortedRight); 
end 

if ~isnan(triggersSortedLeft)
[kernelsLateLeft, kernelsEarlyLeft] = get_kernels(options,stimStream,triggersSortedLeft); 



    
    kernelsEarlyLeft = kernelsEarlyLeft .* -1;
    kernelsLateLeft = kernelsLateLeft .* -1;
    
    end 

kernelsEarly = [kernelsEarlyLeft, kernelsEarlyRight];
kernelsLate = [kernelsLateLeft, kernelsLateRight];


meanKernelLate = nanmean(kernelsLate,2); 
meanKernelEarly = nanmean(kernelsEarly,2); 



end