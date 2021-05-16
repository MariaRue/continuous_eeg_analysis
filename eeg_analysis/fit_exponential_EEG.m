% fit an exponential function to the eeg data
%% fit to CPz for now

%% init paths, toolboxes, plotting legends
addpath(genpath('/Users/maria/Documents/Matlab/cbrewer'))
addpath('/Users/maria/Documents/MATLAB/raacampbell-shadedErrorBar-9b19a7b')


EEGpreproc = '/Volumes/LaCie 1/data_preproc';  % path to behav data all subjs
EEGdir = '/Volumes/LaCie 1/data/';

condition = {'Tr frequent TR short', 'Tr frequent Tr short','Tr rare Tr short', 'Tr rare Tr long'};
% plot each subject separately for conditions
cl = cbrewer('div','RdBu', 12);
cl = cl([2 4 12 10],:);


%num of lags used for calculating the kernel
    lags = 100;
%% load behav data (all subjects and sessions, as well as EEG data for a specific subject and session)

% behavioural data
load_name = fullfile(EEGpreproc,'behav_data_all_subjs');
load(load_name)
nS = max(max(all_responses(:,11)));
subj_list = unique(all_responses(:,12),'stable');

for sj = 1:nS
    % select Subject with correct ID that matches raw data
    
    
    subID = subj_list(sj);
    
    disp(subID);
    % load the EEG data for that subject - all sessions
    eegdat_fname = fullfile(EEGdir,'preprocessed_EEG_dat',[sprintf('sub%03.0f',subID),'_EEGdat.mat']);
    load(eegdat_fname);
    
    
   
    
  
    
    % define condition (bl) and session (se) 
    for se = 1:6
        % contains order of conditions to select correct block from stim
        % stream data - EEGdat is not ordered for conditions but sequence
        % of conditions as they appeared in session 
        block_IDs = unique(all_responses(all_responses(:,10)==se & all_responses(:,12) == subID,9),'stable');
      
        disp(se)
        for bl = 1:4
            
            clear stim_sum
            % determine the condition 
            condition_ID =  block_IDs(bl);
            
            % select the EEG and stim data
            try
                data = EEGdat{se}{bl};
            catch
                data = 0;
            end
            if any(data)%for some blocks matching to stim didn't work - excluding those from fit
                data = data(40,:);
                
                stim = stim_streams{sj,se}(:,condition_ID);
                
               
                
                
                % fit exponential to EEG data
                
                options = optimset('MaxFunEvals',1000000,'MaxIter',100000); % otherwise no convergence
                
                % initial param guesses for exp model
                pstart(1) = 1; % Amplitude
                pstart(2) = 1; % 1/tau
                
                
                fun = @(params)eval_exp_model(params,data,stim,lags); %
                pnew(:,condition_ID,se,sj) = fminsearch(fun,pstart,options);
                
                
            end
        end
    end
end

%% plot taus

% calculate mean taus for each condition and subject

taus = squeeze(pnew(2,:,:,:));
taus(taus==0) = nan(1,1); 
mean_taus = squeeze(nanmean(taus,2));
mean_taus = mean_taus./mean_taus(2,:);

cl1 = cbrewer('qual','Set1', 15);

figure
for sj = 1:nS
    
    hold on
    
    plot([1:4],mean_taus(:,sj),'x-','Color',cl1(sj,:),'LineWidth',3)
    
end
xlim([0 5])
%ylim([-0.5 2.4])
xticks([0:1:5])
xticklabels({'','ITIS INTS','ITIS INTL','ITIL INTS','ITIL INTL',''})
tidyfig;
ylabel('Amplitude')
legend(num2str([1:nS]'));
title('Amplitude across conditions for each subject')

%% 

%% cost function
function [cost] = eval_exp_model(params,data,stim,lags)

A = params(1); % amplitude exponential
tau = params(2); % decay exponential


model = exp(-[0:lags]/tau);
model = model/sum(model);
model = A*model;

model = A * exp(-[0:lags]/tau); % calculate exp across number of lags (number of samples we sum over)

model_vals = conv(stim,model); % multiply exp with stim 

keyboard; 

model_vals = abs(model_vals(lags:end));
cost = sum((data - model_vals').^2) +  sum(params.^2) .* 0.01; % calculate squared error


end
