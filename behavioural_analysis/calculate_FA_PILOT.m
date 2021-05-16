%% integration kernels example task

addpath(genpath('/Users/maria/Documents/Matlab/cbrewer'))
addpath('/Users/maria/Documents/MATLAB/raacampbell-shadedErrorBar-9b19a7b')
EEGpreproc = '/Volumes/LaCie 1/data_preproc';  % path to behav data all subjs

condition = {'Tr frequent TR short', 'Tr frequent Tr short','Tr rare Tr short', 'Tr rare Tr long'};
% plot each subject separately for conditions
cl = cbrewer('div','RdBu', 12);
cl = cl([2 4 12 10],:);


EEGpreproc = '/Volumes/LaCie 1/data_preproc';  % path to behav data all subjs
load_name = fullfile(EEGpreproc,'behav_data_all_subjs_pilot');
load(load_name)

fb = 50;
lags = 400;

which_case = 'Trial'; % or 'Trial' or 'FA'
%%

for sj = 1: 2
session_combined_short = [];
session_combined_long = [];
sess = [1 2 3 4 5 6];
for se = 1:length(sess)
    
    % find condition of session
    condition_ID = unique(all_responses(all_responses(:,11) == sj & all_responses(:,10) == sess(se), 9));

    
    % loop through each block of a session and get FAs and coherences leading
    % up to it
    
    combined_short = [];
    combined_long = [];
    
    for bl = 1:3
        
        
                    % select all stim streams that belong to one subject
            stim_streams_sj = [];
            stim_streams_sj = stim_streams{sj,se}(:,bl);
            
            
            mean_stim = [];
            mean_stim = mean_stim_streams{sj,se}(:,bl);
        
       [triggers_left, triggers_right] = calculate_TR_FA_frames(mean_stim, all_responses,which_case, fb, bl, sj, se,condition_ID);
        
        triggers_left(triggers_left <= lags) = [];
        triggers_right(triggers_right <= lags) = [];
         matrix_right = [];
            for i = 1:length(triggers_right(:,1))
                
                matrix_right(:,i) = stim_streams_sj(triggers_right(i,1) - lags : triggers_right(i,1));
            end
            
            matrix_left = [];
            for i = 1:length(triggers_left(:,1))
                
                matrix_left(:,i) = stim_streams_sj(triggers_left(i,1) - lags : triggers_left(i,1)) .* -1 ; % to be able to combine with right button presses * -1
            end
            
        
            if condition_ID == 1
        combined_short = [combined_short, matrix_left, matrix_right];
            elseif condition_ID == 2
        combined_long = [combined_long, matrix_left, matrix_right]; 
            else 
                keyboard; 
            end
    end
    
    if condition_ID == 1
    
    session_combined_short = [session_combined_short, combined_short];  
    elseif condition_ID == 2
    session_combined_long = [session_combined_long, combined_long];
    end
end
    mean_coherences(:,1,sj) = mean(session_combined_short,2); 
    mean_coherences(:,2,sj) = mean(session_combined_long,2);
    sem_coherence(:,1,sj) = std(session_combined_short')/sqrt(size(mean_coherences,1));
    sem_coherence(:,2,sj) = std(session_combined_long')/sqrt(size(mean_coherences,1));
end
%% plot mean across sj 
figure
for bl = 1:2
    hold on 
    plot(mean(squeeze(mean_coherences(:,bl,:)),2),'Color',cl(bl,:),'LineWidth',3)
    
    
    
    
end 

plot([1:lags+1],zeros(lags+1),'k','LineWidth',3)
 
        tidyfig
        xticks([0:100:400])
        xticklabels([-4 -3 -2 -1 0])
    
            xlabel('time to FA [s]')
            ylabel('mean coherence')
        
%% for different subjects 
figure

for sj = 1:2
subplot(1,2,sj)
for bl = 1:2
    hold on 
    plot(mean_coherences(:,bl,sj),'Color',cl(bl,:),'LineWidth',3)
    
    
    
    
end 

plot([1:lags+1],zeros(lags+1),'k','LineWidth',3)
if sj == 1
 title('Maria')
else 
    title('Laurence')
end 
        tidyfig
        xticks([0:100:400])
        xticklabels([-4 -3 -2 -1 0])
    
            xlabel('time to FA [s]')
            ylabel('mean coherence')
   legend('short','long')         
            
end

%% fit curves do invidual subject data
options = optimset('MaxFunEvals',1000000,'MaxIter',100000);
for con = 1 : 2
    for sj = 1 : 2
        data = mean_coherences(:,con,sj); % data to fit exp model to
        
        % find the peak of the data as starting point
        [val,idx_max] = max(data);
        data_new{con,sj} = data(1:idx_max);
        % time steps
        dt = 0.01;
        num_steps = length(data_new{con,sj});
        
        t_sj{con,sj} = dt : dt : num_steps * dt;
        
        
        % initial param guesses for exp model
        pstart(1) = 1; % Amplitude
        pstart(2) = 1; % 1/tau
        % pstart(3) = 1; % offset
        
        fun = @(p)exp_residual(p,data_new{con,sj},t_sj{con,sj}); % this is the correct cost function that works
        [pnew_sj(con,:,sj),~,exitflag(con,sj)] = fminsearch(fun,pstart,options);
        
    end
end


%% plot single subject fits
for con = 1:2
    figure;
    for sj = 1:2
        model{con,sj} = eval_exp(pnew_sj(con,:,sj),t_sj{con,sj});
        subplot(1,2,sj)
        title(['Laurence'])
        if sj == 1
            title(condition{con})
            ylabel('mean coherence')
            xlabel('t to FA [s]')
        end
        hold on
        plot(data_new{con,sj},'Color',cl(con,:),'LineWidth',3);
        plot(model{con,sj},'k', 'LineWidth', 3);

        xticks([0:100:500])
        xticklabels([-5 -4 -3 -2 -1 0])
        
        tidyfig
        
        hold off
    end
end
%%
figure 
for sj = 1:2
    subplot(1,2,sj)
    for con = 1:2
        hold on
        plot([nan(500-length(model{con,sj}),1); model{con,sj}'],'Color',cl(con,:),'LineWidth',3)
    end
    hold off
    xticks([0:100:500])
    xticklabels([-5 -4 -3 -2 -1 0])
    if sj == 1 
    title(['Maria'])
    else
         title(['Laurence'])
    end
    if sj == 1
        ylabel('mean coherence')
        xlabel('t to FA [s]')
        legend(condition)
    end
    tidyfig;
end


%%
function [residual_reg] = exp_residual(params, data, t)
% this function computes the residuals for the exponential fit which is the
% sum of errers data - exp_model
% params is a vector with
%  amplitude, tau, offset


Amp = params(1);
tau = params(2);


model = Amp .* (exp(t/tau)); % compute the model

model(isinf(model)) = 0;
residual = sum((data - model').^2); % compute the error
residual_reg = residual + sum(params.^2) .* 0.01;

end




%% calculate model for estimated paramters

function model = eval_exp(params,t)
Amp = params(1);
tau = params(2);


model = Amp .* (exp(t/tau)); % compute the model
end
