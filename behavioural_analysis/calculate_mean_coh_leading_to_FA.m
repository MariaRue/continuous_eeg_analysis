% this script takes in streams and triggers from the matrixes genarated
% with 'read_in_behav_data'

%% add paths to plotting packages/data

addpath(genpath('/Users/maria/Documents/Matlab/cbrewer'))
addpath('/Users/maria/Documents/MATLAB/raacampbell-shadedErrorBar-9b19a7b')
EEGpreproc = '/Volumes/LaCie 1/data_preproc';  % path to behav data all subjs

cl = cbrewer('div','RdBu', 12);
cl = cl([2 4 12 10],:);
condition = {'ITIS INTS', 'ITIS INTL','ITIL INTS', 'ITIL INTL'};
%% load data
load_name = fullfile(EEGpreproc,'behav_data_all_subjs');
load(load_name)

lags = 500; % frames back in time that are leading up to FA
%%
% loop through subjects and calculate mean streams of coherences leading up
% to button press for each condition
for sj =  1:length(stim_streams(1,1,1,:))
    all_streams = [];
    all_streams = permute(squeeze(stim_streams(:,:,:,sj)),[1 3 2]); % permute order of entries so that we can cat across sessions for one subject
    all_streams = reshape(all_streams, [], size(stim_streams,2),1); % reshape is actually doing the concatenating
    
    
    % repeat the same for the triggers
    all_triggers = [];
    all_triggers = permute(squeeze(trigger_streams(:,:,:,sj)),[1 3 2]);
    all_triggers = reshape(all_triggers, [], size(trigger_streams,2),1);
    
    
    for bl  = 1:4 % loop through the conditions
        
        
        % find all triggers for a left and right button press FA
        % triggers for FA = 202 right button press, 206 = left button press
        triggers_right =  find(all_triggers(:,bl) == 202);
        triggers_left =  find(all_triggers(:,bl) == 206);
        
        triggers_right(triggers_right <= lags) = [];
        triggers_left(triggers_left <= lags) = [];
        
        matrix_coherences_right = [];
        for i = 1:length(triggers_right)
            matrix_cohrences_right(i,:) = all_streams(triggers_right(i) - lags : triggers_right(i),bl);
        end
        
        matrix_coherences_left = [];
        for ii = 1:length(triggers_left)
            matrix_coherences_left(ii,:) = all_streams(triggers_left(ii) - lags : triggers_left(ii),bl) .* -1;     % multiply by -1 to cmopare to button press to the riht 
        end
        
        num_false_alarms(sj,bl) = size([matrix_coherences_right; matrix_coherences_left],1);
        
        % calculate mean and sem
        mean_coherence(:,bl,sj) = mean([matrix_coherences_right; matrix_coherences_left]);
        sem_coherence(:,bl,sj) = std([matrix_coherences_right; matrix_coherences_left])/sqrt(size(mean_coherence,1));
        
        
    end
    
end

% calculate mean across subjects for each block

sem_across_subjects = squeeze(std(permute(mean_coherence,[3,1,2]))/sqrt(size(mean_coherence,3)));
mean_across_subjects = mean(mean_coherence,3);

%% plot

% prepare colours
cl = cbrewer('div','RdBu', 12);
cl = cl([2 4 12 10],:);

% plot average across subjects
figure(1)

for i = 1:4
    subplot(2,2,i)
    plot(mean_across_subjects(:,i), 'LineWidth', 3,'Color',cl(i,:));
    hold on
    h = shadedErrorBar(1:lags+1,mean_across_subjects(:,i),(sem_across_subjects(:,i)), 'lineprops', '-k');
    h.patch.FaceColor = cl(i,:);
    h.mainLine.Color = cl(i,:);
    %plot(1:501,zeros(501,1),'k--','LineWidth',3)
    title(condition{i})
    %xlim([0 501])
    ylim([0 0.8])
    tidyfig
    %xticks([0:100:500])
    %xticklabels([-5 -4 -3 -2 -1 0])
    if i == 1
        xlabel('time to FA [s]')
        ylabel('mean coherence')
    end
end
hold off

%%
% plot each subject separately for conditions
for bl = 1:4
    figure(bl+1)
    for i = 1:size(mean_coherence,3)
        subplot(3,5,i)
        plot(mean_coherence(:,bl,i), 'LineWidth', 2,'Color',cl(bl,:));
        hold on
        h = shadedErrorBar(1:lags+1,mean_coherence(:,bl,i),(sem_coherence(:,bl,i)), 'lineprops', '-k');
        h.patch.FaceColor = cl(bl,:);
        h.mainLine.Color = cl(bl,:);
        %plot(1:501,zeros(501,1),'k--','LineWidth',3)
        %xlim([0 501])
        %ylim([0 0.8])
        title(['subject',' ',num2str(i)])
        tidyfig
        %xticks([0:100:500])
        %xticklabels([-5 -4 -3 -2 -1 0])
        if i == 1
            xlabel('time to FA [s]')
            ylabel('mean coherence')
            title(condition{bl})
        end
    end
end

%% plot difference curve for each subject and condition and then average - difference curves are ITI long - ITI short and for both ITI conditions long trials - short trials 
% difference 1:ITI long INT long(4) - ITI short INT long(2)
%            2:ITI long INT short(3) - ITI short INT short(1)
%            3:ITI long INT long(4) - ITI long INT short(3)
%            4:ITI short INT long(2) - ITI short INT short(1)

    difference_curves(:,:,1) = squeeze(mean_coherence(:,4,:))' - squeeze(mean_coherence(:,2,:))';
    difference_curves(:,:,2) = squeeze(mean_coherence(:,3,:))' - squeeze(mean_coherence(:,1,:))';
    difference_curves(:,:,3) = squeeze(mean_coherence(:,4,:))' - squeeze(mean_coherence(:,3,:))';
    difference_curves(:,:,4) = squeeze(mean_coherence(:,2,:))' - squeeze(mean_coherence(:,1,:))';

    mean_difference_curves = squeeze(mean(difference_curves));
    sem_difference_curves = squeeze(std(difference_curves)/sqrt(size(mean_difference_curves,2)));

    figure 
    diff_con = {'ITIL INTL - ITIS INTL','ITIL INTS - ITIS INTS', 'ITIL INTL - ITIL INTS','ITIS INTL - ITIS INTS'};
    for con =  1:4
    subplot(2,2,con)
    plot(mean_difference_curves(:,con), 'LineWidth', 2,'Color',cl(con,:));
      hold on
        h = shadedErrorBar(1:lags+1,mean_difference_curves(:,con),(sem_difference_curves(:,con)), 'lineprops', '-k');
        h.patch.FaceColor = cl(con,:);
        h.mainLine.Color = cl(con,:);  

      
        xticks([0:100:500])
        xticklabels([-5 -4 -3 -2 -1 0])    
        
                if con == 1
            xlabel('time to FA [s]')
            ylabel('mean difference coherence')
            
        end
        title(diff_con{con})
          tidyfig
    end 

%% fit exponential (with method from Murray 2014 Nature Neuroscience: A hierarchy of intrinsic timescales across primate cortex)
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
    %pstart(3) = 1; % offset
    
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
    for sj = 1 : size(mean_coherence,3)
        data = mean_coherence(:,con,sj); % data to fit exp model to
        
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
        %pstart(3) = 1; % offset
        
        fun = @(p)exp_residual(p,data_new{con,sj},t_sj{con,sj}); % this is the correct cost function that works
        [pnew_sj(con,:,sj),~,exitflag(con,sj)] = fminsearch(fun,pstart,options);
        
    end
end
%% plot single subject fits
for con = 1:4
    figure;
    for sj = 1:size(mean_coherence,3)
        model{con,sj} = eval_exp(pnew_sj(con,:,sj),t_sj{con,sj});
        subplot(5,3,sj)
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
        xticklabels([-5 -4 -3 -2 -1 0])
  
        tidyfig
        
        hold off
    end
end

%% plot fits for each subject in one graph
figure
for sj = 1:size(mean_coherence,3)
    %subplot(5,3,sj)
    figure
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
    for sj = 1:size(mean_coherence,3)
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
    
    hold on
    plot(x(:,i),tau(:,i),'x','Color',cl(i,:))
    
end

xticks([0:1:5])
xticklabels({'','ITIS INTS','ITIS INTL','ITIL INTS','ITIL INTL'})
tidyfig;
ylabel('inverse tau')
%% 


%% plot tau across conditions for each subject

tau = pnew_sj(:,1,:);
tau = squeeze(permute(tau,[2,3,1]));
%tau = 1./tau;
%tau = tau ./ tau(:,3);
cl1 = cbrewer('qual','Set1', 15);

figure
for sj = 1:size(mean_coherence,3)
    
    hold on 
    
    plot([1:4],tau(sj,:),'x-','Color',cl1(sj,:),'LineWidth',3)

end
xlim([0 5])
%ylim([-0.5 2.4])
xticks([0:1:5])
xticklabels({'','ITIS INTS','ITIS INTL','ITIL INTS','ITIL INTL',''})
tidyfig;
ylabel('amplitude')
legend(num2str([1:size(mean_coherence,3)]'));
title('amplitude across conditions for each subject')


%% helper function for exp fit

function [residual_reg] = exp_residual(params, data, t)
% this function computes the residuals for the exponential fit which is the
% sum of errers data - exp_model
% params is a vector with
%  amplitude, tau, offset


Amp = params(1);
tau = params(2);
%offset = params(3);

%model = Amp .* (exp(t/tau)+ offset); % compute the model
model = Amp .* exp(t/tau);


model(isinf(model)) = 0;
residual = sum((data - model').^2); % compute the error
residual_reg = residual + sum(params.^2) .* 0.01; 

end

%% calculate model for estimated paramters

function model = eval_exp(params,t)
Amp = params(1);
tau = params(2);
%offset = params(3);

%model = Amp .* (exp(t/tau)+ offset); % compute the model
model = Amp .* exp(t/tau); % compute the model
end