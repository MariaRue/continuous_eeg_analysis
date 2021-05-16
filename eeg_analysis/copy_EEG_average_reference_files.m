% this script copys files that were generated with averaged references into
% a specific folder to separate them from the LM_RM references

subj_list = [16, 18:21, 24, 26, 27, 28, 29, 31, 32, 33, 34, 35, 39, 40, 41, 42, 43, 47, 50, 51, 52, 54, 55, 57, 58];

user = 'MR_iMac';
switch user
    case 'MR'
        EEGdir = '/Volumes/LaCie/data_preproc';
        
        
    case 'MR_iMac'
        
        
        EEGdir = '/Users/maria/Documents/data/data_preproc';
        
        
end


%% loop through subjects

for sj = 1:length(subj_list)
    
    subID = subj_list(sj);
    
    if subID == 55
        nSess = 5;
    else
        nSess = 6;
    end
    
 
    for i = 1:nSess
        org_EEGdatadir =  fullfile(EEGdir,sprintf('sub%03.0f',subID),'eeg');
        dest_EEGdatadir = fullfile(EEGdir,sprintf('sub%03.0f',subID),'eeg','average_reference');
        
        cd (org_EEGdatadir)
        if exist('average_reference','dir') ~= 7
            
            mkdir('average_reference')
            
        end
        
        fname = fullfile(org_EEGdatadir,sprintf('Mdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',subID,i));
        [SUCCESS] = copyfile(fname,dest_EEGdatadir)
        delete(fname)
        fname = fullfile(org_EEGdatadir,sprintf('fMdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',subID,i));
        [SUCCESS] = copyfile(fname,dest_EEGdatadir)
        delete(fname)
        fname = fullfile(org_EEGdatadir,sprintf('nanart_fMdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',subID,i));
        [SUCCESS] = copyfile(fname,dest_EEGdatadir)
        delete(fname)
        fname = fullfile(org_EEGdatadir,sprintf('cfMdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',subID,i));
        [SUCCESS] = copyfile(fname,dest_EEGdatadir)
        delete(fname)
        
        
        
        fname = fullfile(org_EEGdatadir,sprintf('Mdspmeeg_sub%03.0f_sess%03.0f_eeg.dat',subID,i));
        [SUCCESS] = copyfile(fname,dest_EEGdatadir)
        delete(fname)
        fname = fullfile(org_EEGdatadir,sprintf('fMdspmeeg_sub%03.0f_sess%03.0f_eeg.dat',subID,i));
        [SUCCESS] = copyfile(fname,dest_EEGdatadir)
        delete(fname)
        fname = fullfile(org_EEGdatadir,sprintf('nanart_fMdspmeeg_sub%03.0f_sess%03.0f_eeg.dat',subID,i));
        [SUCCESS] = copyfile(fname,dest_EEGdatadir)
        delete(fname)
        fname = fullfile(org_EEGdatadir,sprintf('cfMdspmeeg_sub%03.0f_sess%03.0f_eeg.dat',subID,i));
        [SUCCESS] = copyfile(fname,dest_EEGdatadir)
        delete(fname)
    end
    
end