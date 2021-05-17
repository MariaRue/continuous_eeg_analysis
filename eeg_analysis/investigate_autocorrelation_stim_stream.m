%% add paths
addpath(genpath('/Users/maria/Documents/Matlab/cbrewer'))
addpath('/Users/maria/Documents/MATLAB/raacampbell-shadedErrorBar-9b19a7b')
EEGpreproc = '/Volumes/LaCie 1/data_preproc';  % path to behav data all subjs

condition = {'Tr frequent TR short', 'Tr frequent Tr short','Tr rare Tr short', 'Tr rare Tr long'};
% plot each subject separately for conditions
cl = cbrewer('div','RdBu', 12);
cl = cl([2 4 12 10],:);
load_name = fullfile('/Volumes/crdkData/preprocessedData/behaviour','behav_data_all_subjs_all3');
load(load_name)
%% 
clear C
C(:,1) = xcorr(stim_streams{4,2}(:,1));
C(:,2) = xcorr(stim_streams{4,2}(:,2));
C(:,3) = xcorr(stim_streams{4,2}(:,3));
C(:,4) = xcorr(stim_streams{4,2}(:,4));

num = length(stim_streams{2,2}(:,1));
for topl = 1:4 
    subplot(2,2,topl)
plot(C(:,topl),'Color',cl(topl,:),'LineWidth',3)
xlim([num num+10000])
%ylim([-1000,1000])
legend(condition(topl))
end


