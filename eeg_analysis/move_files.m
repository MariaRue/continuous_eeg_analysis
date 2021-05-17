%subj_list = [16,18:21, 24, 26, 27, 28, 29, 31, 33, 34, 35, 39,40,  41, 42, 43, 47,50, 51, 52, 54, 55, 57, 58];
subj_list = [50, 51, 52, 54, 55, 57, 58];
%subj_list = [62:64,66,68,70];
% for sj = 1:length(subj_list)
%     subID = subj_list(sj);
% 
%     sourcePath = fullfile('/Users/maria/Documents/data/data.continuous_rdk/data/EEG/',sprintf('sub%03.0f',subID));
%     destinationPath = fullfile('/Volumes/crdkData/rawData/experiment',sprintf('sub%03.0f',subID));
%     
%     if exist(destinationPath) ~= 7
%         mkdir(destinationPath)
%     end 
%     copyfile(fullfile(sourcePath,'behaviour'),fullfile(destinationPath,'behaviour'));
%     copyfile(fullfile(sourcePath,'eye'),fullfile(destinationPath,'eye'));
%     copyfile(fullfile(sourcePath,'stim'),fullfile(destinationPath,'stim'));
%     
%     mkdir (fullfile(destinationPath,'eeg'))
%     destinationEegPath = fullfile(destinationPath,'eeg');
%     sourceEegPath = fullfile(sourcePath,'eeg');
%     
%     copyfile(fullfile(sourceEegPath,'*.cdt'),destinationEegPath);
%     copyfile(fullfile(sourceEegPath,'*.ceo'),destinationEegPath);
%     copyfile(fullfile(sourceEegPath,'*.dpa'),destinationEegPath);
%     
% end 

    sourcePath = '/Users/maria/Documents/data/data_preproc';
    destinationPath = '/Volumes/crdkData/conventionalEEGAnalysis/LMRM';
for sj = 1:length(subj_list)
    subID = subj_list(sj);
    
copyfile(fullfile(sourcePath,'test_April_2020',['trial_start_locked_EEG_dat_',sprintf('sub%03.0f',subID),'.mat']),destinationPath);
copyfile(fullfile(sourcePath,['csd_trial_start_locked_EEG_dat_',sprintf('sub%03.0f',subID),'.mat']),destinationPath);    

copyfile(fullfile(sourcePath,['csd_response_locked_EEG_dat_',sprintf('sub%03.0f',subID),'.mat']),destinationPath);    
copyfile(fullfile(sourcePath,['response_locked_EEG_dat',sprintf('sub%03.0f',subID),'.mat']),destinationPath);     
    
    
    
end