function [threshold] = threshold_fun(subID,sess)
% This function returns the correct threshold to detect artefacts and
% eyeblinks in the data. Threshold is a 2x1 vector with the first entry
% being the threshold for eyeblink detection and the second for artefact
% detection

switch subID
    
    case 16
        
        if  sess == 1
            threshold(1,1) = 4;
            threshold(2,1) = 2;
        elseif sess == 2
            
            % huge artefact at start - won't allow me to set an appropriate
            % eyeblink threshold - this needs to be removed before the
            % eyeblink detection is run again
            
            % for now
            threshold(1,1) = 0.9;
            threshold(2,1) = 2;
            
        elseif sess == 3
            threshold(1,1) = 3;
            threshold(2,1) = 1.5;
            
        elseif sess == 4
            threshold(1,1) = 2;
            threshold(2,1) = 1.5;
            
        elseif sess == 5 || sess == 6
            threshold(1,1) = 1.5; % not sure
            threshold(2,1) = 2; % data towards end quite noise - not sure what to make of this

        end
        
    case 18
        
        
        if  sess == 1 || sess == 2 || sess == 3 || sess == 5 || sess == 6
            
            
            % sess 2 good example of more data thrown out during artefact
            % rejection - is that correct?
            
            threshold(1,1) = 1.5;
            threshold(2,1) = 1.5;
            
            
        elseif sess == 4
            threshold(1,1) = 1.5;
            threshold(2,1) = 2;

        end
        

    case 19
        % session 2/3 eyeblinks weird? are these other artefacts?
        if sess == 1
            threshold(1,1) = 3;
            threshold(2,1) = 1.5;
            
        elseif sess == 2
            threshold(1,1) = 5;
            threshold(2,1) = 1.5;
            
        elseif sess == 3 || sess == 4 || sess == 5 || sess == 6
            threshold(1,1) = 10;
            threshold(2,1) = 1.5;
            
        end
        
    case 20
        
        if sess == 1 || sess == 2
            
            threshold(1,1) = 2;
            threshold(2,1) = 3;
            
            
        elseif sess == 3 || sess == 4 % not sure about eyeblink threshold probably ok?
            threshold(1,1) = 2;
            threshold(2,1) = 1.5;
            
            
            
        elseif sess == 5
            threshold(1,1) = 2;
            threshold(2,1) = 5; %CPZ super noisy? discard? - we have results from this participant though from the GLM that look alright...
            
        elseif sess == 6
            
            threshold(1,1) = 2;
            threshold(2,1) = 5; % same as above  - not sure about threshold
            
        end
        
    case 21
        
        if sess == 1 || sess == 2 || sess == 3
            
            threshold(1,1) = 2;
            threshold(2,1) = 3;
            
        elseif sess >= 4
            
            threshold(1,1) = 2;
            threshold(2,1) = 2;
            
        end
        
    case 24
        
        if sess == 1 || sess == 2 || sess == 3 || sess == 4
            
            threshold(1,1) = 2;
            threshold(2,1) = 2;
            
            % with around 10 eyeblinks per minute if reduced we
            % catch each little bump. Right now we might miss one
            % every now and then
            
            
        elseif sess == 5 || sess == 6
            
            threshold(1,1) = 3;
            threshold(2,1) = 2;
            
        end
        
        
    case 26
        
        threshold(1,1) = 3;
        threshold(2,1) = 2;
        
        % sess 2 super weird for threshold for eyeblinks -
        % missing big peaks - why?
        % same for 3 and 4, 5 threshold for normal artefacts for 5
        % ok?
        
        if sess == 6
            threshold(2,1) = 3.5;
        end

    case 28
        threshold(1,1) = 3;
        threshold(2,1) = 3.5;
        
    case 32
        
        threshold(1,1) = 3;
        threshold(2,1) = 3.5;
        
        
        if sess >= 2
            threshold(1,1) = 4;
            
        end
        
        
    case 34
        threshold(1,1) = 2;
        threshold(2,1) = 3.5;
    case 35
        threshold(1,1) = 2;
        threshold(2,1) = 3.5;
        % good example where I think we don't catch artefacts in sessions 4
        % or 5
    case 40
        % eyeblink detection - is that correct - maybe the signal is just
        % distorted because of big peaks? - after inspecting channel, I
        % believe that there are much less eyeblinks - but not sure
        threshold(1,1) = 2;
        threshold(2,1) = 3.5;

    case 41
        threshold(1,1) = 2;
        threshold(2,1) = 3.5;
        
        % allowed of spikes specifically in last sessions - maybe check it
        % out?
    case 51
        threshold(1,1) = 2;
        threshold(2,1) = 3.5;        
    case 52
        
        threshold(1,1) = 2;
        threshold(2,1) = 3.5;    
        
    case 55 
        if sess == 3 || sess == 5 || sess == 4
        threshold(1,1) = 4;
         threshold(2,1) = 3.5; 
        else 
                    threshold(1,1) = 2;
        threshold(2,1) = 3.5; 
        end 
        % high num of eyeblinks per minute... like  over 25 per minute
    otherwise 
        
        threshold(1,1) = 2; 
        threshold(1,1) = 3.5; 
end
end
