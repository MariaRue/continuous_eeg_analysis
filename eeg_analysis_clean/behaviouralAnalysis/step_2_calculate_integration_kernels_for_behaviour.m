% this script takes behavioural data created with step_1_read_in_behav_data
% and calculates integration kernels leading up to a button press and
% fitting an exponential to it

%% add paths
addpath(genpath('/Users/maria/Documents/Matlab/cbrewer'))
addpath('/Users/maria/Documents/MATLAB/raacampbell-shadedErrorBar-9b19a7b')
EEGpreproc = '/Volumes/crdkData/preprocessedData/behaviour';  % path to behav data all subjs


condition = {'Tr frequent TR short', 'Tr frequent Tr short','Tr rare Tr short', 'Tr rare Tr long'};

% plot each subject separately for conditions
cl = cbrewer('div','RdBu', 12);
cl = cl([2 4 12 10],:);
%% load data
load_name = fullfile(EEGpreproc,'behav_data_all_subjs_all3');
load(load_name)

lags = 500; % frames back in time that are leading up to FA


which_responses = 'false alarms';  % calculating integration kernel for either FA button presses or button presses during trials options: 'false alarms' or 'trials', '3sec_rts',

with_coherence = 'without coherence level'; % if we want to know the coherence levels for the trials version
nS = max(all_responses(:,11)); % number of subjects

coherence = [0.3 0.4 0.5];

counter = 0; 

%% loop through subjects and find button presses
mean_coherences = [];
for sj = 1 : nS
    
    
    for bl = 1 :4
        combined = [];
        combined_coh = [];
        
        se_total = 6; 
        if sj == 26 
            
            se_total = 5 ;
            
        end 
        for se = 1:se_total
            % select only stimstreams from all sessions that belong to specific
            % block
            
            % select all stim streams that belong to one subject
            stim_streams_sj = [];
            stim_streams_sj = stim_streams{sj,se}(:,bl);
            
            
            % select trigger streams that belong to one subject
            trigger_streams_sj = [];
            trigger_streams_sj = trigger_streams{sj,se}(:,bl);
            
            mean_streams = [];
            mean_streams = mean_stim_streams{sj,se}(:,bl);
            

            responses = all_responses((all_responses(:,9)== bl & all_responses(:,10) == se & all_responses(:,11) == sj),:);


            
            % find all triggers that lead to a button press
            switch which_responses
                
                
                case 'false alarms'
                    
                    % find triggers right and left button press (202 and 206)
                    
                    triggers_right = [];
                    triggers_left = [];
                    

                    % this is with eeg triggers - don't use it
%                     triggers_right = find(trigger_streams_sj == 202);
%                     triggers_left = find(trigger_streams_sj == 206);
                    
                    
                     % this is with triggers from the response matrix 
%                     
                    triggers_right = responses((responses(:,7) == 2 & responses(:,3) == 1),6);
                    triggers_left = responses((responses(:,7) == 2 & responses(:,3) == 0),6);


                    
                case 'trials'
                    % this is with EEG triggers - don't use it 

                    % this is with triggers from the EEG 
%                     triggers_right = find(trigger_streams_sj == 202);
%                     triggers_left = find(trigger_streams_sj == 206);

                    % this is with triggers from the response matrix 
%                     
                    triggers_right = responses((responses(:,7) == 1 & responses(:,3) == 1),6);
                    triggers_left = responses((responses(:,7) == 1 & responses(:,3) == 0),6);


                    
                case '3sec_rts'

                    % this is for eeg triggers - don't use it 
                    

                               % this is with triggers from the EEG

%                     triggers_right = find(trigger_streams_sj == 201);
%                     triggers_left = find(trigger_streams_sj == 205);
                    % find all
   

                    
                    % this is for frames taken from response matrix
                         triggers_right = responses((responses(:,7) == 1 & responses(:,3) == 1),6);
                    triggers_left = responses((responses(:,7) == 1 & responses(:,3) == 0),6);
                    
                    % right trials
                    rts_rigth = zeros(length(triggers_right(:,1)),1);
                    for i = 1:length(triggers_right(:,1))
                        
                        t_org = triggers_right(i,1);
                        t = triggers_right(i,1);
                        
                        while ~(mean_streams(t) ~= 0 && mean_streams(t-1) == 0)
                            t = t-1;
                            
                        end
                        
                        rts_rigth(i) = t_org - t;
                        
                        if rts_rigth(i) > 300
                            
                            triggers_right(i,:) = nan;
                            counter = counter + 1; 
                        end
                        
                    end
                    
                    
                    % left trials
                    rts_left = zeros(length(triggers_left(:,1)),1);
                    for i = 1:length(triggers_left(:,1))
                        
                        t_org = triggers_left(i,1);
                        t = triggers_left(i,1);
                        
                        while ~(mean_streams(t) ~= 0 && mean_streams(t-1) == 0)
                            t = t-1;
                            
                        end
                        
                        rts_left(i) = t_org - t;
                        
                        if rts_left(i) > 300
                            
                            triggers_left(i,:) = nan;
                            
                        end
                        
                    end
                    
                    
                    % remove nan trials - trials with rt > 300
                    triggers_right(isnan(triggers_right(:,1)),:) = [];
                    triggers_left(isnan(triggers_left(:,1)),:) = [];
            end
            
            
            
            
            % only choose triggers that are bigger than the lags we go
            % back
            
            
            triggers_right(triggers_right(:,1)<=lags,:) = [];
            triggers_left(triggers_left(:,1)<=lags,:) = [];
            
            
            switch with_coherence
                
                case 'with coherence levels'
                    
                    coh_val_right = [];
                    for i = 1:length(triggers_right(:,1))
                        
                        % this for EEG triggers 
                        if mean_streams(triggers_right(i,1)) ~= 0
                            coh_val_right(i) = mean_streams(triggers_right(i,1));
                        else
                            t = 0; coh_val = 0;
                            while t <= 51 && coh_val == 0
                                t = t+1;
                                coh_val = mean_streams(triggers_right(i)-t,1);
                                
                                
                            end
                            
                            coh_val_right(i) = coh_val;
                        end
                        
                        
                        

                    end
                    
                    coh_val_left = [];
                    for i = 1:length(triggers_left(:,1))
                        
                        if mean_streams(triggers_left(i,1)) ~= 0
                            coh_val_left(i) = mean_streams(triggers_left(i,1));
                        else
                            t = 0; coh_val = 0;
                            while t <= 51 && coh_val == 0
                                t = t+1;
                                coh_val = mean_streams(triggers_left(i) - t,1);
                                
                                
                            end
                            coh_val_left(i) = coh_val;
                        end
                        
                    end
                    
                    
                    
                    
                    
                    % get rid of incorrect trials
                    
                    triggers_right(coh_val_right < 0) = [];
                    triggers_left(coh_val_left > 0) = [];
                    coh_val_right(coh_val_right < 0) = [];
                    coh_val_left(coh_val_left > 0) = [];
                    
                    
                    
                    
                    
                    matrix_right = [];
                    
                    if any(triggers_right)
                    for i = 1:length(triggers_right(:,1))
                        
                        matrix_right(:,i) = stim_streams_sj(triggers_right(i,1) - lags : triggers_right(i,1));
                    end
                    end
                    
                    matrix_left = [];
                    
                     if any(triggers_left)
                    for i = 1:length(triggers_left(:,1))
                        
                        matrix_left(:,i) = stim_streams_sj(triggers_left(i,1) - lags : triggers_left(i,1)) .* -1 ; % to be able to combine with right button presses * -1
                    end
                     end
                    
                    combined = [combined,matrix_right, matrix_left];
                    combined_coh = [combined_coh,coh_val_right,abs(coh_val_left)];
                case 'without coherence level'
                    % loop through triggers for right and left button presses and
                    % select coherences from stim_streams_bl
                    matrix_right = [];
                    for i = 1:length(triggers_right(:,1))
                        
                        matrix_right(:,i) = stim_streams_sj(triggers_right(i,1) - lags : triggers_right(i,1));
                    end
                    
                    matrix_left = [];
                    for i = 1:length(triggers_left(:,1))
                        
                        matrix_left(:,i) = stim_streams_sj(triggers_left(i,1) - lags : triggers_left(i,1)) .* -1 ; % to be able to combine with right button presses * -1
                    end
                    
                    combined = [combined,matrix_right, matrix_left];
                    
            end
            
            
            
        end
        
        switch with_coherence
            
            case 'without coherence level'
                mean_coherences(:,bl,sj) = nanmean(combined,2);
                sem_coherence(:,bl,sj) = nanstd(combined')/sqrt(size(mean_coherences,1));
                num_false_alarms(sj,bl) = size(combined,2);
                
            case 'with coherence levels'
                
                for coh = 1:3
                    
                    if any(combined_coh == coherence(coh))
                     
                     idx_comb =combined_coh == coherence(coh);

                  
                    mean_coherences(:,bl,sj,coh) = nanmean(combined(:,idx_comb),2);
                    sem_coherence(:,bl,sj,coh) = nanstd(combined(:,idx_comb)')/sqrt(size(mean_coherences,1));
                    num_false_alarms(sj,coh,bl) = size(combined,2);
                
                        
                    
                    end 
                    
                end
        end
        
        
    end
    
    
end





switch with_coherence
    case 'without coherence level'
        sem_across_subjects = squeeze(std(permute(mean_coherences,[3,1,2]))/sqrt(size(mean_coherences,3)));
        mean_across_subjects = nanmean(mean_coherences,3);
        
    case 'with coherence levels'
        
        for coh = 1:3
            sem_across_subjects(:,:,coh) = squeeze(std(permute(squeeze(mean_coherences(:,:,:,coh)),[3,1,2]))/sqrt(size(mean_coherences,3)));
            mean_across_subjects(:,:,coh) = nanmean(squeeze(mean_coherences(:,:,:,coh)),3);
        end
end
%% plot conditions across subjects for each coh level
thres = round(.05 * 100); 
repetitions = 100; 
nS = 28; 
anova_samples = 28; 
plt = 1; 
color = [1,0,0; 0.5 1 1; 0.5 1 0.5];
for coh = 1:3

    [sig_tps,p,tbl] = shuffled_permutation_test_Fscore(squeeze(mean_coherences(:,:,:,coh)), lags, thres, repetitions, anova_samples, plt,nS);
        figure (10)
    subplot(1,3,coh)
    for i = 1:4
        %subplot(2,2,i)
        hold on
        plot(mean_across_subjects(:,i,coh), 'LineWidth', 3,'Color',cl(i,:));
        hold on
        h = shadedErrorBar(1:lags+1,mean_across_subjects(:,i,coh),(sem_across_subjects(:,i,coh)), 'lineprops', '-k');
        h.patch.FaceColor = cl(i,:);
        h.mainLine.Color = cl(i,:);
        
        title(['integration kernel at ',num2str(coherence(coh))])
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
    
    y_level = [0.6 0.55 0.5];
   
   figure (10)
   subplot(1,3,coh)
   hold on
    for id = 1:3
    if any(sig_tps{id})
plot(sig_tps{id},ones(length(sig_tps{id}),1).* y_level(id),'.','MarkerSize',10)

    end
    end
ylim([-0.2 0.7])
    hold off
    keyboard;
  %  legend(condition)
    
end


%%


repetitions = 100; 
thres = round(.05 * repetitions); 
nS = 28; 
anova_samples = 28; 
plt = 1; 
color = [1,0,0; 0.5 1 1; 0.5 1 0.5];
[sig_tps,p,tbl] = shuffled_permutation_test_Fscore(mean_coherences, lags, thres, repetitions, anova_samples, plt,nS);
for i = 1:4
    %subplot(2,2,i)
    figure (11)
    hold on
   b(i) = plot(mean_across_subjects(:,i), 'LineWidth', 3,'Color',cl(i,:));
    hold on
    h = shadedErrorBar(1:lags+1,mean_across_subjects(:,i),(sem_across_subjects(:,i)), 'lineprops', '-k');
    h.patch.FaceColor = cl(i,:);
    h.mainLine.Color = cl(i,:);
    
    
    title('mean coherence leading to a button press participants')
    xlim([0 501])
    
    %tidyfig
    xticks([0:100:500])
    xticklabels([-5 -4 -3 -2 -1 0])
    if i == 1
        xlabel('time to button press [s]')
        ylabel('mean coherence')
    end
end
plot(1:501,zeros(501,1),'k--','LineWidth',3)
hold on 
y_level = [0.6 0.55 0.5];
for id = 1%:3
    if any(sig_tps{id})
plot(sig_tps{id},ones(length(sig_tps{id}),1).* y_level(id),'.','MarkerSize',2)

    end
end
hold off
ylim([-0.1 0.65])
%plot(sig_tps{1},ones(length(sig_tps{1}),1).* 0.5,'k*','MarkerSize',2)
legend(b,condition)
%%
% plot each subject separately for conditions

for bl = 1:4
    figure(bl+1)
    for i = 1:size(mean_coherences,3)
        subplot(4,7,i)
        plot(mean_coherences(:,bl,i), 'LineWidth', 1,'Color',cl(bl,:));
        hold on
        h = shadedErrorBar(1:lags+1,mean_coherences(:,bl,i),(sem_coherence(:,bl,i)), 'lineprops', '-k');
        h.patch.FaceColor = cl(bl,:);
        h.mainLine.Color = cl(bl,:);
        %plot(1:501,zeros(501,1),'k--','LineWidth',3)
        xlim([0 501])
        %ylim([0 0.8])
        title(['subject',' ',num2str(i)])
        tidyfig
        xticks([0:100:500])
        xticklabels([-5 -4 -3 -2 -1 0])
        if i == 1
            xlabel('time to FA [s]')
            ylabel('mean coherence')
            title(condition{bl})
        end
    end
end
%% %% fit exponential (with method from Murray 2014 Nature Neuroscience: A hierarchy of intrinsic timescales across primate cortex)
for con = 1:4
    data = mean_across_subjects(:,con); % data to fit exp model to
    
    % find the peak of the data as starting point
    [val,idx_max] = max(data);
    data_new{con} = data(1:idx_max);
    % time steps
    dt = 0.01;
    num_steps = length(data_new{con});
    
    t{con} = dt : dt : num_steps * dt;
    
    
    % initial param guesses for exp model
    pstart(1) = 1; % Amplitude
    pstart(2) = 1; % 1/tau
    
    
    fun = @(p)exp_residual(p,data_new{con},t{con}); % this is the correct cost function that works
    pnew(con,:) = fminsearch(fun,pstart);
    
end
%% plot the results
for con = 1:4
    model_mean{con} = eval_exp(pnew(con,:),t{con});
    subplot(2,2,con)
    title(condition{con})
    hold on
    plot(data_new{con},'Color',cl(con,:),'LineWidth',3);
    plot(model_mean{con},'k', 'LineWidth', 3);
    xticks([0:100:500])
    xticklabels([-5 -4 -3 -2 -1 0])
    ylabel('mean coherence')
    xlabel('t to FA [s]')
    tidyfig
    
    hold off
end
%% plot all modeled curves in one graph
figure;
hold on
for con = 1:4
    plot([nan(500-length(model_mean{con}),1); model_mean{con}'],'-','Color',cl(con,:),'LineWidth',3)
    
end
plot([1:500],zeros(500,1),'k--','LineWidth',3)
xticks([0:100:500])
xticklabels([-5 -4 -3 -2 -1 0])
legend(condition)
title('model fit averaged across subjects for each condition')
ylabel('mean coherence')
xlabel('t to FA [s]')
tidyfig;
%% fit curves do invidual subject data
options = optimset('MaxFunEvals',1000000,'MaxIter',100000);
for con = 1 : 4
    for sj = 1 : size(mean_coherences,3)
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

% find subject with less than 20 responses to exclude from further analysis
[r,c] = find(num_false_alarms<=20);
r = unique(r);

switch which_responses
    
    case 'trials'
        save('pnew_sj_TR', 'pnew_sj','r')
        
    case 'false alarms'
        
        save('pnew_sj_FA', 'pnew_sj','r')
        
    case '3sec_rts'
        save('pnew_sj_3sec_rts', 'pnew_sj', 'r')      
end

%% plot single subject fits
for con = 1:4
    figure;
    for sj = 1:size(mean_coherences,3)
        model{con,sj} = eval_exp(pnew_sj(con,:,sj),t_sj{con,sj});
        subplot(7,4,sj)
        title(['subject', ' ', num2str(sj)])
        if sj == 1
            title(condition{con})
            ylabel('mean coherence')
            xlabel('t to FA [s]')
        end
        hold on
        plot(data_new{con,sj},'Color',cl(con,:),'LineWidth',3);
        plot(model{con,sj},'k', 'LineWidth', 3);
        
        txt = ['FA: ', num2str(num_false_alarms(sj,con))];
        text(20,0.5,txt,'FontSize',14)
        xticks([0:100:500])
        xticklabels([ 5000 4000 3000 2000 1000 0])
        set(gca, 'XDir','reverse')
        tidyfig
        
        hold off
    end
end
%% figure for laurence 

figure;
i = 1; 
for sj = [9,15,24]
for con = 1:4
    
    subplot(3,4,i)
    i = i+1;
        model{con,sj} = eval_exp(pnew_sj(con,:,sj),t_sj{con,sj});
     
%         title(['subject', ' ', num2str(sj)])
%        
%             title(condition{con})
%             ylabel('mean coherence')
%             xlabel('t to FA [s]')
        
        hold on
        plot(data_new{con,sj},'Color',cl(con,:),'LineWidth',3);
        plot(model{con,sj},'k', 'LineWidth', 3);
%         
%         txt = ['FA: ', num2str(num_false_alarms(sj,con))];
%         text(20,0.5,txt,'FontSize',14)
        xticks([0:100:500])
        xticklabels([ 5 4 3 2 1 0])
       % set(gca, 'XDir','reverse')
        tidyfig
        
        hold off
    end
end



%% plot fits for each subject in one graph
figure
for sj = 1:size(mean_coherences,3)
    subplot(5,3,sj)
    for con = 1:4
        hold on
        plot([nan(500-length(model{con,sj}),1); model{con,sj}'],'Color',cl(con,:),'LineWidth',3)
    end
    hold off
    xticks([0:100:500])
    xticklabels([-5 -4 -3 -2 -1 0])
    title(['subject', ' ', num2str(sj)])
    if sj == 1
        ylabel('mean coherence')
        xlabel('t to FA [s]')
        legend(condition)
    end
    tidyfig;
end

%% plot curves of one condition on top of each other
figure
for con = 1:4
    subplot(2,2,con)
    title(condition{con})
    hold on
    for sj = 1:size(mean_coherences,3)
        prepared_data = [nan(500-length(model{con,sj}),1); model{con,sj}'];
        plot(prepared_data,'Color',cl(con,:),'LineWidth',3)
    end
    hold off
    
end




%% swarm plot of paramter values for each condition/ parameter - here only for tau invers
tau = pnew_sj(:,2,:);
tau = squeeze(permute(tau,[2,3,1]));
tau(tau>=3.9e05)=nan;

% generate x values
a = 0.6;
b = 1.4;
figure;
for i = 1:4
    x(:,i) = (b-a).*rand(15,1) + a;
    b = b+1;
    a = a+1;
    
%     hold on
%     plot(x(:,i),tau(:,i),'x','Color',cl(i,:))
    
end

xticks([0:1:5])
xticklabels({'','ITIS INTS','ITIS INTL','ITIL INTS','ITIL INTL'})
tidyfig;
ylabel('inverse tau')
%%


%% plot tau across conditions for each subject

% only include taus for estimates with more than 20 responses 
[r,c] = find(num_false_alarms<=20);
r = unique(r);
pnew_sj(:,:,r) = [];
%%
tau = pnew_sj(:,2,:);
tau = squeeze(permute(tau,[2,3,1]));

mean_tau = mean(tau(:));
sd_tau = 2 * std(tau(:));
[sj_id,con_id] = find(tau >= (mean_tau + sd_tau));

tau(sj_id,:) = [];


%tau = tau./tau(:,3);
cl1 = cbrewer('qual','Set1', 28);

figure
for sj = 1:size(tau(:,1),1)
    
    hold on
    
    plot([1:4],tau(sj,:),'x-','Color',cl1(sj,:),'LineWidth',3)
    
end
xlim([0 5])
%ylim([-0.5 2.4])
xticks([0:1:5])
xticklabels({'','ITIS INTS','ITIS INTL','ITIL INTS','ITIL INTL',''})
tidyfig;
ylabel('tau')
legend(num2str([1:size(tau(:,1),1)]'));
title('amp across conditions for each subject trials 3sec_rt')
%%

taus = squeeze(pnew_sj(:,2,:));
% normal distributed? check with histogram - for trials this could be
% trueish

histogram(taus(:),10);

% calculate mean and stds for each conditoin and plot
mean_taus = mean(taus,2);
std_taus = std(taus');
x = 1:length(mean_taus);
figure
errorbar(x,mean_taus,std_taus,'kx','LineWidth',3,'MarkerSize',10)
ylabel('1/tau')
xlim([0,5])
xticks([1:4])
xticklabels(condition)
xlabel('condition')
title ('mean tau per condition')
tidyfig
%% anovan to test for difference (not sure this is correct)

taus = squeeze(pnew_sj(:,2,:));

taus_reorderd = [[taus(1,:)';taus(2,:)'],[taus(3,:)';taus(4,:)']];

[p,tbl] = anova2(taus_reorderd,16);

%% investigate sig difference by calculating a GLM 

% regressors to regress out mean across conditions for a subject 
sj_mean_reg = zeros(length(tau)*4,length(tau));

for sj = 1:length(tau)
    if sj == 1 
       start_idx = 1; 
    else 
         start_idx = (sj-1)*4 + 1;
    end 
       
    end_idx = sj * 4; 
    
    sj_mean_reg(start_idx:end_idx,sj) = ones(4,1); 
    
    
end 

% regressor frequency cond - frequent = 1, rare = -1
freq_reg = zeros(length(tau)*4,1); 
freq_reg(1:4:end) = 1; 
freq_reg(2:4:end) = 1;
freq_reg(3:4:end) = -1;
freq_reg(4:4:end) = -1;


% regressor long vs short - short = 1; long = -1; 
len_reg = zeros(length(tau)*4,1); 
len_reg(1:4:end) = 1; 
len_reg(2:4:end) = -1; 
len_reg(3:4:end) = 1; 
len_reg(4:4:end) = -1; 

% interaction long vs short 

interact_reg = len_reg .* freq_reg; 


% regressors combined 
reg_comb = [sj_mean_reg freq_reg len_reg interact_reg];


taus = tau(:); 

[b,~,stats] = glmfit(reg_comb,taus,'normal','constant','off');


figure; 

bar(stats.t(length(tau)+1:end))
xticklabels({'frequent trials','short trials','interaction'})
ylabel('t stats')
xlabel('regressor')
tidyfig;


%% correlate FA taus with TR taus
FA_exp_fits = load('pnew_sj_FA.mat');
TR_exp_fits = load('pnew_sj_3sec_rts.mat');

FA_taus = squeeze(FA_exp_fits.pnew_sj(:,1,:))';
TR_taus = squeeze(TR_exp_fits.pnew_sj(:,1,:))';

num_tr = unique([FA_exp_fits.r;TR_exp_fits.r]);

FA_taus(num_tr,:) = [];
TR_taus(num_tr,:) = [];

mean_FA_taus = mean(FA_taus(:));
sd_FA_taus = 2 * std(FA_taus(:));
[sj_ls_FA,con_id] = find(FA_taus >= (mean_FA_taus + sd_FA_taus));

mean_TR_taus = mean(TR_taus(:));
sd_TR_taus = 2 * std(TR_taus(:));
[sj_ls_TR,con_id] = find(TR_taus >= (mean_TR_taus + sd_TR_taus));

sj_id = unique([sj_ls_FA;sj_ls_TR]); 
 
FA_taus(sj_id,:) = []; 
TR_taus(sj_id,:) = [];

figure

for con = 1:4
    
    [R,P] = corrcoef(FA_taus(:,con),TR_taus(:,con));

    subplot(2,2,con)
    plot(FA_taus(:,con),TR_taus(:,con),'*','Color',cl(con,:),'LineWidth',1,'MarkerSize',5)
    title([condition{con},' ','R= ',num2str(round(R(1,2),2)),' ','P= ', num2str(round(P(1,2),2))])
    xlabel('FA amplitude')
    ylabel('TR amplitude')
    tidyfig;
    
    
end

%% correlate taus and amplitudes across conditions for false alarms to go with single subject integration kernels
FA_exp_fits = load('pnew_sj_FA.mat');

sj_fa_out = FA_exp_fits.r; 

tau_short = squeeze(FA_exp_fits.pnew_sj([1,3],2,:)); 
tau_short_mean = mean(tau_short);
tau_short_mean(sj_fa_out) = [];

tau_long = squeeze(FA_exp_fits.pnew_sj([2,4],2,:));
tau_long_mean = mean(tau_long);
tau_long_mean(sj_fa_out) = [];

tau_freq = squeeze(FA_exp_fits.pnew_sj([1,2],2,:));
tau_freq_mean = mean(tau_freq);
tau_freq_mean(sj_fa_out) = [];

tau_rare = squeeze(FA_exp_fits.pnew_sj([3,4],2,:)); 
tau_rare_mean = mean(tau_rare);
tau_rare_mean(sj_fa_out) = [];


amp_short = squeeze(FA_exp_fits.pnew_sj([1,3],1,:));
amp_short_mean = mean(amp_short);
amp_short_mean(sj_fa_out) = [];
 

amp_long = squeeze(FA_exp_fits.pnew_sj([2,4],1,:));
amp_long_mean = mean(amp_long); 
amp_long_mean(sj_fa_out) = [];

amp_freq = squeeze(FA_exp_fits.pnew_sj([1,2],1,:));
amp_freq_mean = mean(amp_freq); 
amp_freq_mean(sj_fa_out) = [];

amp_rare = squeeze(FA_exp_fits.pnew_sj([3,4],1,:)); 
amp_rare_mean = mean(amp_rare); 
amp_rare_mean(sj_fa_out) = [];


[R_tau_len, P_tau_len] = corrcoef(tau_short_mean, tau_long_mean); 
[R_amp_len, P_amp_len] = corrcoef(amp_short_mean, amp_long_mean);

[R_tau_freq, P_tau_freq] = corrcoef(tau_freq_mean, tau_rare_mean); 
[R_amp_freq, P_amp_freq] = corrcoef(amp_freq_mean, amp_rare_mean);


figure

    subplot(2,2,1)
    plot(tau_short_mean,tau_long_mean,'kd','LineWidth',1,'MarkerSize',5)
    title([condition{con},' ','R= ',num2str(round(R_tau_len(1,2),2)),' ','P= ', num2str(P_tau_len(1,2),2)])
    xlabel('short')
    ylabel('long')
    tidyfig;
    
    subplot(2,2,2)
    plot(amp_short_mean,amp_long_mean,'kd','LineWidth',1,'MarkerSize',5)
    title([condition{con},' ','R= ',num2str(round(R_amp_len(1,2),2)),' ','P= ', num2str(P_amp_len(1,2),2)])
    xlabel('short')
    ylabel('long')
    tidyfig;
    
    subplot(2,2,3)
    plot(tau_freq_mean,tau_rare_mean,'kd','LineWidth',1,'MarkerSize',5)
    title([condition{con},' ','R= ',num2str(round(R_tau_freq(1,2),2)),' ','P= ', num2str(P_tau_freq(1,2),2)])
    xlabel('frequent')
    ylabel('rare')
    tidyfig;

    subplot(2,2,4)
    plot(amp_freq_mean,amp_rare_mean,'kd','LineWidth',1,'MarkerSize',5)
    title([condition{con},' ','R= ',num2str(round(R_amp_freq(1,2),2)),' ','P= ', num2str(P_amp_freq(1,2),2)])
    xlabel('frequent')
    ylabel('rare')
    tidyfig;

%% helper function for exp fit

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
