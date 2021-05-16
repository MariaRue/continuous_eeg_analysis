% this script simulates the integration kernel: Every time the stim stream
% hits a coherence of abs(0.6) there is some probability that a FA is
% assigned - take these and calculate the integration kernel for these to
% see whether they look the same as for the real participants



%% add paths
addpath(genpath('/Users/maria/Documents/Matlab/cbrewer'))
addpath('/Users/maria/Documents/MATLAB/raacampbell-shadedErrorBar-9b19a7b')
EEGpreproc = '/Volumes/crdkData/preprocessedData/behaviour'; ;  % path to behav data all subjs

condition = {'Tr frequent TR short', '','Tr frequent Tr short','','Tr rare Tr short', '','Tr rare Tr long'};
% plot each subject separately for conditions
cl = cbrewer('div','RdBu', 12);
cl = cl([2 4 12 10],:);

%% load data
load_name = fullfile(EEGpreproc,'behav_data_all_subjs_all3');
load(load_name)

p = 0.2; % frames back in time that are leading up to FA

lags = 500; % amount of lags used to calculate kernel

nS = 25%max(all_responses(:,11)); % number of subjects

%% go through stim streams (only ITIs) and assign FAs

for sj = 1:nS
    
    for bl = 1:4
        combined = [];
        for se = 1:6
            
            
            stim = stim_streams{sj,se}(:,bl);
            mean_stim = mean_stim_streams{sj,se}(:,bl);
            % find frames with coherences that are 0.6 or higher?
            
            idx_coh = abs(stim(:)) > 0.5;
            
            % idx_frame = find(idx_coh(2:end) == 1 & idx_coh(1:end-1) == 0); % find the first frame of each of those steps
            
        idx_frame = find(idx_coh);
            idx_frame(mean_stim(idx_frame)~=0) = []; % only select frames that are actually in an iti period
            
            % assign probability
            
            FA_idx = rand(size(idx_frame)) <= p;
            
         
            
            FA_frames = idx_frame(FA_idx);
            
            signed_coh = sign(stim(FA_frames)); % left or right coherences?
            
            
            
            % calculate kernel
            
            triggers_right = [];
            triggers_left = [];
            
            
            triggers_right = FA_frames(signed_coh > 0);
            triggers_left = FA_frames(signed_coh < 0);
            
            
            triggers_right(triggers_right(:,1)<=lags,:) = [];
            triggers_left(triggers_left(:,1)<=lags,:) = [];
            
            
            % loop through triggers for right and left button presses and
            % select coherences from stim
            matrix_right = [];
            for i = 1:length(triggers_right(:,1))
                
                matrix_right(:,i) = stim(triggers_right(i,1) - lags : triggers_right(i,1));
            end
            
            matrix_left = [];
            for i = 1:length(triggers_left(:,1))
                
                matrix_left(:,i) = stim(triggers_left(i,1) - lags : triggers_left(i,1)) .* -1 ; % to be able to combine with right button presses * -1
            end
            
            combined = [combined,matrix_right, matrix_left];
            
            
            
        end
        
        
        mean_coherences(:,bl,sj) = mean(combined,2);
        sem_coherence(:,bl,sj) = std(combined')/sqrt(size(mean_coherences,1));
        num_false_alarms(sj,bl) = size(combined,2);
        
    end
    
end
sem_across_subjects = squeeze(std(permute(mean_coherences,[3,1,2]))/sqrt(size(mean_coherences,3)));
mean_across_subjects = mean(mean_coherences,3);

%% 

figure (1)
for i = 1:4
    %subplot(2,2,i)
    hold on
    plot(mean_across_subjects(:,i), 'LineWidth', 3,'Color',cl(i,:));
    hold on
    h = shadedErrorBar(1:lags+1,mean_across_subjects(:,i),(sem_across_subjects(:,i)), 'lineprops', '-k');
    h.patch.FaceColor = cl(i,:);
    h.mainLine.Color = cl(i,:);
  
    title('mean coherence leading to a button press')
    xlim([0 501])
    
    tidyfig
    xticks([0:100:500])
    xticklabels([-5 -4 -3 -2 -1 0])
    if i == 1
        xlabel('time to button press [s]')
        ylabel('mean coherence')
    end
end
  plot(1:501,zeros(501,1),'k--','LineWidth',3)
hold off

legend(condition)
%% 





