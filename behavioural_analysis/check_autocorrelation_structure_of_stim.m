% check autocorrelation structure of each stimulus stream for each subject 

%% add paths
addpath(genpath('/Users/maria/Documents/Matlab/cbrewer'))
addpath('/Users/maria/Documents/MATLAB/raacampbell-shadedErrorBar-9b19a7b')
EEGpreproc = '/Volumes/LaCie 1/data_preproc';  % path to behav data all subjs

condition = {'Tr frequent TR short', 'Tr frequent Tr short','Tr rare Tr short', 'Tr rare Tr long'};
% plot each subject separately for conditions
cl = cbrewer('div','RdBu', 12);
cl = cl([2 4 12 10],:);
%% load data
load_name = fullfile(EEGpreproc,'behav_data_all_subjs');
load(load_name)
nS = max(max(all_responses(:,11)));
%% 
for sj = 1:nS 
    
    subplot(3,5,sj)
   
    hold on
    for se = 1:6
        for bl = 1:4
    plot(xcorr(stim_streams{sj,se}(:,bl)),'Color',cl(bl,:))
        end 
    end
    title(num2str(sj))
    hold off 
end 

%% check for each subject whether on shorter time scales there is some structure in the data 

for sj = 1:nS
    
    for se = 1:6
       
        figure
        
        for bl = 1:4
         subplot(2,2,bl)
         
         for split = 1:5
             
         end 
            
            
            
        end 
        
    end
    
    
end 