%% specify the datasets to read in by identifying the subid
% in terminal go to the data folder containing sufolders for each
% participant --> ls > txt (make sure that list contains only direct
% subject directories)
% EEGdir= fullfile('/Users/maria/Documents/data/data.continuous_rdk','data','EEG');

% subj_list = [16,18:21,24,26,32, 52, 41, 34, 35, 51, 40, 28]; %2 missing
%subj_list = [16, 18:21, 24, 26, 27, 28, 29, 31, 32, 33, 34, 35, 39, 40, 41, 42, 43, 47, 50, 51, 52, 54, 55, 57, 58];
subj_list = [40,47];


[EEGdir,EEGdirdata,scriptdir,nSess,nS] = setup_EEG_session(subj_list);

reference_type = 'LM_RM'; %average of left and right mastoid
source_density = 1; %compute current source density transformation?

force_downsample = 0; %force downsampling to be (re)run
force_montage = 0; %force montaging to be (re)run
force_filter = 0; %force filtering and artifact detection to be (re)run
force_CSD_transform = 1; %force CSD transform to be (re)run

user = 'MR';

%subj_list = [16];
%%
for sj = 1:nS
    
    
    
    subID = subj_list(sj);
    BHVdatadir = fullfile(EEGdir,sprintf('sub%03.0f',subID),'behaviour');
    STdatadir = fullfile(EEGdir,sprintf('sub%03.0f',subID),'stim');
    EEGdatadir =  fullfile(EEGdir,sprintf('sub%03.0f',subID),'eeg');
    
    subID
    
    
    cd (EEGdatadir)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%--------------- convert_raw_data_to_eeglab_format -------------------
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    eeg_files = dir('*set');
    
    % if raw data hasn't been transformed into eeglab set files then run
    % the following code and save raw data as set files
    if isempty({eeg_files.name})
        switch user
            case 'LH'
                %do nothing - EEGlab is already in my path
            case 'MR'
                addpath('/Users/maria/Documents/MATLAB/eeglab14_1_2b/');
        end
        eeglab
        
        eeg_raw = dir('*dpa');
        nSess = 6; %length({eeg_raw.name});
        if subID == 55 
            nSess = 5; 
        end 
    
        for l = 1:nSess
            
            
            fname_load = fullfile(EEGdatadir,...
                sprintf('sub%03.0f_sess%03.0f_eeg.cdt',subID,l));
            EEG = loadcurry(fname_load, 'CurryLocations', 'False');
            
            pop_saveset(EEG,'filename',sprintf('sub%03.0f_sess%03.0f_eeg.set',subID,l));
        end
        
    else
        
        nSess = length({eeg_files.name});
        
        if subID == 55 
            nSess = 5; 
            
        end 
        
    end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%--------------- convert_raw_data_to_eeglab_format -------------------
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % load in data from one subject or created SPM file if that hasn't
    % happened yet.
    for i = 1:nSess
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%--------------- convert_downsample_EEG -------------------
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % name of converted eeg data set from eeglab format
        fname_target = fullfile(EEGdatadir,sprintf('spmeeg_sub%03.0f_sess%03.0f_eeg.mat',subID,i));
        
        
        % if data is already converted to spm, just load the file,
        % otherwise convert it
        if exist(fname_target,'file')
            D{i} = spm_eeg_load(fname_target);
            
        else
            S = [];
            
            S.dataset = fullfile(EEGdatadir,sprintf('sub%03.0f_sess%03.0f_eeg.set',subID,i));
            
            
            S.mode = 'continuous';
            D{i} = spm_eeg_convert(S);
        end
        
        % this is the file name of downsampled file
        fname_target = fullfile(EEGdatadir,sprintf('dspmeeg_sub%03.0f_sess%03.0f_eeg.mat',subID,i));
        
        % check whether this file already exists and load it in otherwise
        % downsmaple data
        if exist(fname_target,'file')&~force_downsample
            D{i} = spm_eeg_load(fname_target);
        else
            S = [];
            S.D = D{i};
            S.fsample_new = 100;
            D{i} = spm_eeg_downsample(S);
            
        end
                
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%--------------- convert_downsample_EEG -------------------
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   
   
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%--------------- rereference_EEG_data -------------------
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % rereference the electrodes and remove artefacts including
        % eyeblinks
        
        %%%% 1) reference - choose between LM-RM to reference to average
        % of left and right mastoids or average_reference which
        % rereferences to average of electrodes across scalp
        
        
        fname_target = fullfile(EEGdatadir,sprintf('Mdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',subID,i));
        
        keyboard;    
        
        if exist(fname_target,'file')&~force_montage
            D{i} = spm_eeg_load(fname_target);
        else
            
            S = [];
            S.D = D{i};
            
            S.mode = 'write';
            
            
            S.keepothers = 1; % to keep the EOG channels
            switch reference_type
                
                
                case 'LM_RM' % rereference the 61 EEG channels to L+R mastoid average
                    
                    
                    %in our montage we have 61 EEG channels, plus the right mastoid:
                    S.montage.labelorg = D{i}.chanlabels(1:62);
                    %we keep all EEG channels, but throw out the right mastoid:
                    S.montage.labelnew = D{i}.chanlabels(1:61);
                    
                    %build our M*N matrix for montaging:
                    S.montage.tra = eye(62);
                    %subtract 0.5 of the right mastoid in the re-reference - see
                    %p.108 of Luck book!
                    S.montage.tra(:,62) = -0.5;
                    %get rid of the right mastoid in the new channels:
                    S.montage.tra(62,:) = [];
                    
                    
                    %add in extra electrode for lateralised readiness potential (C4-C3)
                    S.montage.tra(62,:) = S.montage.tra(selectchannels(D{i},'C4'),:) - S.montage.tra(selectchannels(D{i},'C3'),:);
                    S.montage.labelnew{62} = 'C4_C3_LRP';
                    
                case 'average_reference'
                   
                    %in our montage we have 61 EEG channels. We can ignore
                    %the right mastoid, as it isn't in our equation.
                    S.montage.labelorg = D{i}.chanlabels(1:61);
                    %we keep all EEG channels
                    S.montage.labelnew = D{i}.chanlabels(1:61);
                    
                    %build our M*N matrix for montaging:
                    S.montage.tra = eye(61)-(1/61);
                    
                    %add in extra electrode for lateralised readiness potential (C4-C3)
                    S.montage.tra(62,:) = S.montage.tra(selectchannels(D{i},'C4'),:) - S.montage.tra(selectchannels(D{i},'C3'),:);
                    S.montage.labelnew{62} = 'C4_C3_LRP';
                    
                otherwise
                    error('unrecognised re-referencing type')
            end
            
            D{i} = spm_eeg_montage(S);
            %set the LRP electrode to be EEG chantype
            D{i} = chantype(D{i},selectchannels(D{i},'C4_C3_LRP'),'EEG')
            D{i}.save();
        end % referencing
        
           
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%--------------- rereference_EEG_data -------------------
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   
   
      
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%--------------- artefact_and_eyeblink_labelling-------------------
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%% 2) filter the data bandpass between 0.1 and 30Hz
        
        fname_target = fullfile(EEGdatadir,sprintf('fMdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',subID,i));
        % this is filename of filtered data with artefacts removed
        
        if exist(fname_target,'file')&~force_filter
            D{i} = spm_eeg_load(fname_target);
        else
            
            S = [];
            S.D = D{i};
            S.band = 'bandpass';
            S.freq = [0.1 30];
            D{i} = spm_eeg_filter(S);
            
            
            
            %%%%% 3) % detect eyeblinks and add these to the filtered data file
            
            %1. grab the EOG data, and set anything exceeding 1000 microvolts
            %   to 0
            VEOG_data = D{i}(selectchannels(D{i},'VEOG'),:,:);
            D{i}(selectchannels(D{i},'VEOG'),find(abs(VEOG_data)>1000),:) = 0;
            D{i}.save;
            
            %2. use spm_eeg_artefact_eyeblink to identify eyeblinks and add
            %them as events.
            
            % get threshold
            [threshold] = threshold_fun(subID,i);
            S = [];
            S.D = D{i};
            S.chanind = selectchannels(D{i},'VEOG'); %VEOG channel in the
            S.threshold = threshold(1,1); %standard deviations away from the me an
            
            S.mode = 'mark';
            S.append = 1; %this adds the eyeblinks as events
            S.excwin = 500; %remove 500ms around each event
            fprintf('Subject %0.0f - run %0.0f\n',sj,i);
            
            D{i} = spm_eeg_artefact_eyeblink(S);
            
            
            D{i}.save; %stores the eyeblink events
            %
            %
            S = [];
            S.D = D{i};
            S.chanind = selectchannels(D{i},'EEG');
            S.threshold = 100;
            S.badchanthresh = 1000; % possibly number of artefacts/length of samples that lead to exclusion of whole channel?
            S.mode = 'mark';
            S.append = 1;
            S.excwin = 500;
            
            D{i} = spm_eeg_artefact_threshchan(S);
            bs = sum(D{i}.badsamples([selectchannels(D{i},'CPZ') selectchannels(D{i},'VEOG')],:,1))>=1; %find artefacts in any channels
            mean(bs)
            figure (3)
            plot(D{i}.time,D{i}(40,:,1)); hold on; plot(D{i}.time(find(bs==1)), D{i}(40,find(bs==1),1),'r.');
            D{i}.save; %stores the events above 100 µm Volt
            close all
        end
        
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%--------------- artefact_and_eyeblink_labelling-------------------
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%--------------- make_eeg_copy_with_ artefacts_for_fieldtrip_timelock_analysis-------------------
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
        
        % make a copy of the artefact labelled file where artefacts are
        % replaced with NaNs, for reading into fieldtrip
        fname_target = fullfile(EEGdatadir,sprintf('nanart_fMdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',subID,i));
        if exist(fname_target,'file')
            %do nothing - we don't pass this down to the next stage
        else
            S = [];
            S.D = D{i};
            S.outfile = ['nanart_' D{i}.fname]; %append 'c' on start of filename for CSD transformation
            
            tmpD = spm_eeg_copy(S);
            
            bs = tmpD.badsamples(:,:,:);
            for ch = 1:size(bs,1)
               tmpD(ch,find(bs(ch,:)),1) = nan;
            end
            tmpD.save; clear tmpD
        end
        
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%--------------- make_eeg_copy_with_ artefacts_for_fieldtrip_timelock_analysis-------------------
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
   
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%---------------transform_csd_eeg-------------------
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%% source_density ANALYSIS 
        fname_target = fullfile(EEGdatadir,sprintf('cfMdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',subID,i));
        if exist(fname_target,'file')&~force_CSD_transform  % source density transformation of data
            D{i} = spm_eeg_load(fname_target);
        else
            S = [];
            S.D = D{i};
            S.outfile = ['c' D{i}.fname]; %append 'c' on start of filename for CSD transformation
            
            D{i} = spm_eeg_copy(S);
            
            %1) transform data into fieldtrip format
            ft_data = D{i}.ftraw();
            
            %2 source density analysis
            
            % remove non-EEG electrodes as ft_scalpcurrentdensity doesn't
            % use these
            cfg = [];

            cfg.channel = {'all','-RM', '-VEOG', '-HEOG', '-LM', '-C4_C3_LRP'};
            [data_pre] = ft_preprocessing(cfg, ft_data);
            
            cfg = [];
            cfg.elec = ft_data.elec;
            cfg.degree = 14;
            source_density_data =  ft_scalpcurrentdensity(cfg, data_pre);
            
            D{i}(1:61,:,1) = source_density_data.trial{1}; %put the source density transformed data into cfMdspmeeg...
            D{i}.save;
           
        end
        
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%---------------transform_csd_eeg-------------------
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%---------------make_eeg_copy_with_ artefacts_for_fieldtrip_timelock_analysis-------------------
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% Maria edit - insert nan in csd data 
            fname_target = fullfile(EEGdatadir,sprintf('cnanart_fMdspmeeg_sub%03.0f_sess%03.0f_eeg.mat',subID,i));
        if exist(fname_target,'file')&~force_CSD_transform  % source density transformation of data
           % do nothing - we don't pass this down to the next stage 
        else
            
            
            S = [];
            S.D = D{i};
            S.outfile = ['cnanart_' D{i}.fname]; %append 'c' on start of filename for CSD transformation
            
            tmpD = spm_eeg_copy(S);
        
            bs = tmpD.badsamples(:,:,:);
            for ch = 1:size(bs,1)
               tmpD(ch,find(bs(ch,:)),1) = nan;
            end
            tmpD.save; clear tmpD
       
        end
    end
end
cd(scriptdir)
