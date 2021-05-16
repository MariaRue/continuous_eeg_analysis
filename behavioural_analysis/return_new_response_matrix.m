function [new_response_matrix,nResponsesDiscarded] = return_new_response_matrix(observed_stimulus_stream,orig_response_matrix)

%make a copy of the original response matrix
response_matrix = orig_response_matrix;

% Remove all rows of response matrix where choice (3rd column)==2 (get rid of missed trials)
response_matrix(response_matrix(:,3)==2,:) = [];

% Remove rows where the choice is a nan (these were pre-allocated for memory issues)
response_matrix(isnan(response_matrix(:,3)),:) = [];

%
% Set current_trial_coherence to 0;
current_trial_coherence = 0;

% Set time_since_last_response to 0;
time_since_last_response = 0;

% Set k to 0; %this is a counter for where we are in the new response matrix
k = 0;

%set up variables that will change during loop
current_trial_start_frame = nan;
in_trial = 0;
recently_in_trial = 0;
recently_responded_in_trial = 0; 

% counter for number of responses thrown away by our new fancy algorithm
nResponsesDiscarded = 0;


% Loop i from 1 to end of stimulus stream

for i = 1:length(observed_stimulus_stream)
    
    if abs(observed_stimulus_stream(i))>0
        if in_trial==0
            current_trial_start_frame = i;
        end
        in_trial = 1;
        recently_in_trial = 0;
        current_trial_coherence =  observed_stimulus_stream(i);
    elseif any(abs(observed_stimulus_stream(max(1,i-50):(i-1)))>0) % was there coherence in last 50 frames...
        in_trial = 0;
        recently_in_trial = 1;
    else
        if recently_in_trial == 1 && ...
                recently_responded_in_trial == 0 %we have just transitioned out of a 'recently_in_trial' period but didn't make a response - this means we missed a trial?
            k = k + 1; %counter for new response matrix
            
            if current_trial_coherence==0
                error('Current trial coherence is not known - BUG!!');
            end
            new_response_matrix(k,1) = -1.5; %cost for missed trial
            new_response_matrix(k,2) = (i-current_trial_start_frame-1)*10/1000; % "reaction time" for missed trial - should be 5.5 or 3.5
            new_response_matrix(k,3) = 2; % response type for missed trial
            new_response_matrix(k,4) = current_trial_coherence; %coherence
            new_response_matrix(k,5) = 0; % incorrect
            new_response_matrix(k,6) = i;
            new_response_matrix(k,7) = 3; % label this as a missed trial
        end
        current_trial_start_frame = nan;
        in_trial = 0;
        recently_in_trial = 0;
        current_trial_coherence = 0;
    end
    
    if any(response_matrix(:,6)==i) % we made a response!

            if recently_responded_in_trial %We throw away these events
     			nResponsesDiscarded = nResponsesDiscarded + 1;
            else
                k = k + 1; %index of current row for new_response_matrix
                
                % obtain j, which element it is that matches ? i.e. which response number
                j = find(response_matrix(:,6)==i);
                
                %number of points given in task (obtained from experimental code)
                new_response_matrix(k,1) = response_matrix(j,1);
                
                %check reaction time in the old response matrix (because it
                %is more accurrate than we can compute here)
                newRT = ((i - current_trial_start_frame)*10)/1000;
                oldRT = response_matrix(j,2);
                if ~isnan(newRT)&&abs((newRT-oldRT)<0.2)
                    %we can use the oldRT because it is more accurate
                    new_response_matrix(k,2) = oldRT;
                elseif ~isnan(newRT) & (isnan(oldRT)|abs((newRT-oldRT)>=0.2))
                    warning('Using new RT! May be less precise or inaccurate')
                    new_response_matrix(k,2) = newRT;
                elseif isnan(newRT)
                    new_response_matrix(k,2) = nan;
                end
                
                %copy across choice
                new_response_matrix(k,3) = response_matrix(j,3);
                
                %set coherence
                new_response_matrix(k,4) = current_trial_coherence;
                
                %set current frame number
                new_response_matrix(k,6) = i;
                
                if in_trial || recently_in_trial %are we currently in a trial period?
                    if current_trial_coherence==0
                        error('Current trial coherence is not known - BUG!!');
                    end
                    % now evaluate whether this is a correct response
                    if new_response_matrix(k,3)==1 && current_trial_coherence > 0 || ... %right response, rightward coherence
                            new_response_matrix(k,3)==0 && current_trial_coherence < 0 %left response, leftward coherence
                        %set response as correct:
                        new_response_matrix(k,5) = 1;
                        new_response_matrix(k,7) = 1;
                    else
                        %set response as incorrect:
                        new_response_matrix(k,5) = 0;
                        new_response_matrix(k,7) = 0;
                    end
                    
                    recently_responded_in_trial = 1;
                else % we made a false alarm OH DEAR
                    if observed_stimulus_stream(i)~=0
                        error('Current trial coherence is not 0 and we think we made a false alarm!');
                    end
                    new_response_matrix(k,5) = 0;
                    new_response_matrix(k,7) = 2;
                    recently_responded_in_trial = 0;
                end
            end
            time_since_last_response = 0;
    else %we didn?t make a response!
        time_since_last_response = time_since_last_response + 1;
    	if recently_responded_in_trial == 1 && time_since_last_response>50
            recently_responded_in_trial = 0;
        end
    end % of response
end % end of i loop
