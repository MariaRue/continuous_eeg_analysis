function [MeanLengthOfInterval,StartOfIntertrialIntervals,EndOfIntertrialIntervals] = calculate_mean_length_of_interval(mean_stim_streams,details)

LengthOfIntertrialInterval = cell(4,1);

    for session = 1:length(details.sessionIDs)
        
        for condition = 1:4 % loop conditions
            
       
            StartOfIntertrialIntervals{details.sessionIDs(session),condition} = find( (mean_stim_streams{details.sessionIDs(session)}(1:end-1,condition) ~= 0 & mean_stim_streams{details.sessionIDs(session)}(2:end,condition) == 0));
            EndOfIntertrialIntervals{details.sessionIDs(session),condition} = find( (mean_stim_streams{details.sessionIDs(session)}(1:end-1,condition) == 0 & mean_stim_streams{details.sessionIDs(session)}(2:end,condition) ~= 0));
            
            StartOfIntertrialIntervals{details.sessionIDs(session),condition} = [1; StartOfIntertrialIntervals{details.sessionIDs(session),condition}];
            EndOfIntertrialIntervals{details.sessionIDs(session),condition} = [EndOfIntertrialIntervals{details.sessionIDs(session),condition}; length(mean_stim_streams{details.sessionIDs(session)}(:,condition))];
            
            LengthOfIntertrialInterval{condition} = [LengthOfIntertrialInterval{condition}; (EndOfIntertrialIntervals{details.sessionIDs(session),condition} - StartOfIntertrialIntervals{details.sessionIDs(session),condition})];
            
        end
        
    end
    
    for condition = 1:4
        MeanLengthOfInterval(condition) = mean(LengthOfIntertrialInterval{condition}).*0.5;
    end
    
end