%% add paths to plotting packages/data
clear all
addpath(genpath('/Users/maria/Documents/Matlab/cbrewer'))
addpath('/Users/maria/Documents/MATLAB/raacampbell-shadedErrorBar-9b19a7b')
EEGpreproc = '/Volumes/crdkData/preprocessedData/behaviour';   % path to behav data all subjs

cl = cbrewer('div','RdBu', 12);
cl = cl([2 4 12 10],:);
condition = {'ITIS INTS', 'ITIS INTL','ITIL INTS', 'ITIL INTL'};
condition_bar = {'ITIS INTS', '','ITIS INTL','','ITIL INTS','','ITIL INTL','single sj'};
%% load data
load_name = fullfile(EEGpreproc,'behav_data_all_subjs_all3');
load(load_name)

nS = max(all_responses(:,11)); % num subjects 

coherence_list = [0.3, 0.4, 0.5];

%% 

for sj = 1:nS 
    
    % plot rts for different conditions and coherences for each subject
    
    figure 
    
    idx_fig = 0; 
    for coh = 1:3 
        idx_fig = idx_fig + 1; 
        subplot(3,2,idx_fig)
        
        idx_rts = abs(all_responses(:,4)) == coherence_list(coh) & all_responses(:,7) == 1 & all_responses(:,11) == sj; % idx to all trials with a correct response, a given coherence and a specific subject 
         
        rts = all_responses(idx_rts,[2,9]); % select rt and condition for idx_rts 
        
        hold on 
        histogram(rts(rts(:,2) == 1,1), 10, 'FaceColor', cl(1,:))
        histogram(rts(rts(:,2) == 2,1), 10, 'FaceColor', cl(2,:))

        hold off 
        legend(condition{1:2})
        title(['subject: ', num2str(sj),' coherence: ', num2str(coherence_list(coh))])
        xlabel('rt sec')
        tidyfig;
        
        idx_fig = idx_fig + 1;
        subplot(3,2,idx_fig)
        hold on
        histogram(rts(rts(:,2) == 3,1), 7, 'FaceColor', cl(3,:))
        histogram(rts(rts(:,2) == 4,1), 7, 'FaceColor', cl(4,:))
        hold off 
        xlabel('rt sec')
        legend(condition{3:4})
        tidyfig;
    end 
 
end 

%% only include rts up to 3.5 seconds 



for sj = 1:nS 
    
    % plot rts for different conditions and coherences for each subject
    
    figure 
    
    idx_fig = 0; 
    for coh = 1:3 
        idx_fig = idx_fig + 1; 
        subplot(3,2,idx_fig)
        
        idx_rts = abs(all_responses(:,4)) == coherence_list(coh) & all_responses(:,7) == 1 & all_responses(:,11) == sj & all_responses(:,2) <= 3.5; % idx to all trials with a correct response, a given coherence and a specific subject 
         
        rts = all_responses(idx_rts,[2,9]); % select rt and condition for idx_rts 
        
        hold on 
        histogram(rts(rts(:,2) == 1,1), 10, 'FaceColor', cl(1,:))
        histogram(rts(rts(:,2) == 2,1), 10, 'FaceColor', cl(2,:))

        hold off 
        legend(condition{1:2})
        title(['subject: ', num2str(sj),' coherence: ', num2str(coherence_list(coh))])
         xlabel('rt sec')
         tidyfig;
        idx_fig = idx_fig + 1;
        subplot(3,2,idx_fig)
        hold on
        histogram(rts(rts(:,2) == 3,1), 7, 'FaceColor', cl(3,:))
        histogram(rts(rts(:,2) == 4,1), 7, 'FaceColor', cl(4,:))
        hold off 
        legend(condition{3:4})
        xlabel('rt sec')
        tidyfig;
    end 
 
end 
%% plot rts ITI short vs ITI long

for sj = 1:nS 
    
    % plot rts for different conditions and coherences for each subject
    
    figure 
    
    idx_fig = 0; 
    for coh = 1:3 
        idx_fig = idx_fig + 1; 
        subplot(3,2,idx_fig)
        
        idx_rts = abs(all_responses(:,4)) == coherence_list(coh) & all_responses(:,7) == 1 & all_responses(:,11) == sj ; % idx to all trials with a correct response, a given coherence and a specific subject 
         
        rts = all_responses(idx_rts,[2,9]); % select rt and condition for idx_rts 
        
        hold on 
        histogram(rts(rts(:,2) == 1,1), 10,'Normalization','probability', 'FaceColor', cl(1,:))
        histogram(rts(rts(:,2) == 3,1), 10, 'Normalization','probability', 'FaceColor', cl(2,:))

        hold off 
        legend(condition{[1,3]})
        title(['subject: ', num2str(sj),' coherence: ', num2str(coherence_list(coh))])
         xlabel('rt sec')
         ylabel('count per numel(X)')
        idx_fig = idx_fig + 1;
        subplot(3,2,idx_fig)
        hold on
        histogram(rts(rts(:,2) == 2,1), 7, 'Normalization','probability', 'FaceColor', cl(3,:))
        histogram(rts(rts(:,2) == 4,1), 7, 'Normalization','probability','FaceColor', cl(4,:))
        hold off 
        legend(condition{[2,4]})
        xlabel('rt sec')
        
    end 
 
end 

%% short ITI vs long ITI across trial lengths



for sj = 1:nS 
    
    % plot rts for different conditions and coherences for each subject
    
    figure 
    
    idx_fig = 0; 
    for coh = 1:3 
        idx_fig = idx_fig + 1; 
        subplot(3,1,idx_fig)
        
        idx_rts = abs(all_responses(:,4)) == coherence_list(coh) & all_responses(:,7) == 1 & all_responses(:,11) == sj ; % idx to all trials with a correct response, a given coherence and a specific subject 
         
        rts = all_responses(idx_rts,[2,9]); % select rt and condition for idx_rts 
        
        hold on 
        histogram(rts(rts(:,2) == 1 | rts(:,2) == 2,1), 10,'Normalization','probability', 'FaceColor', cl(1,:))
        histogram(rts(rts(:,2) == 3 | rts(:,2) == 4,1), 10, 'Normalization','probability', 'FaceColor', cl(3,:))

        hold off 
        legend({'ITI short', 'ITI long'})
        title(['subject: ', num2str(sj),' coherence: ', num2str(coherence_list(coh))])
         xlabel('rt sec')
         ylabel('count per numel(X)')
tidyfig;
        
    end 
 
end 

%%  median rt analysis 
clear rts
% median rt for each coherence per subject and condition
for sj = 1:nS
    for con = 1:4
        for coh = 1:3
            
    idx_rts = abs(all_responses(:,4)) == coherence_list(coh) & all_responses(:,9) == con & all_responses(:,7) == 1 & all_responses(:,11) == sj & all_responses(:,2) <= 3.5; % idx to all trials with a  a given coherence, condition , correct response,and a specific subject 
    
   
    rts(coh,con,sj)  = nanmedian(all_responses(idx_rts,2));
        end 
    end 

end 

rts_sum = nanmedian(rts,3);

figure 
hold on 
x = 1; 
for coh = 1:3
    
    for con = 1:4 
    
       b(con) = bar(x,rts_sum(coh,con),'FaceColor',cl(con,:),'EdgeColor',cl(con,:)); 
        ln = plot([ones(nS,1).*x],squeeze(rts(coh,con,:)),'k.','MarkerSize',7,'LineWidth',3);
        
        
    x = x+1;     
    end 
    
    x = x + 1.5; 
end 
hold off

 xticks([2.5, 8, 13.5])
 xticklabels([0.3 0.4 0.5])
 xlabel('coherence')
 ylabel('median rt')
 legend(b(1:4),condition)
 tidyfig;
%% regression rt analysis 

clear X
idx = 0; 


sj_idx = 0; % for regressors for mean for each sj 
for sj = 1:nS
    
    for coh = 1:3
        
        for con = 1:4
            
          
            idx = idx+1;
            Y(idx) = rts(coh,con,sj);
            X(idx,1) = coherence_list(coh);
            if con == 1 || con == 2
            X(idx,2) = 1;
            else 
                X(idx,2) =-1; 
            end
            
            if con == 1 || con == 3
                X(idx,3) = 1; 
            else 
                X(idx,3) = -1; 
            end 
            
            
            X(idx,4) = X(idx,3) * X(idx,2); 
            
          
            X(idx,5 + sj_idx) = 1;
    
            
        end 
        
    end 
    
     sj_idx = sj_idx + 1;
    
end 
%zscore

    for reg = 1:3
        
        
    X(:,reg) = zscore(X(:,reg));
    end
    

    idx = isnan(Y);
    Y(idx) = [];
    X(idx,:) = [];
    
[b,~,stats] = glmfit(X,Y,'normal','constant','off');

figure; 

bar(stats.t(1:4))
xticklabels({'coherence','frequent trials','short trials','interaction','sj1','sj2','sj3','sj4','sj5','sj6','sj7','sj8','sj9','sj10','sj11','sj12','sj13','sj14'})
ylabel('t stats')
xlabel('regressor')
tidyfig;
ylim([-20 5])
% why is constant so high? Shouldn't that be in sj? 
% why X over parameterised when I split up coherences, because technically 
% they are also conditions and all the variance that isn't explained by 0.3 or 0.4 should be 0.5? 
% also is trial length vs frequence ok? Looks like they are not correlated
%% rt analysis regression analysis for each subject and test betas in t-test 


clear rts
% rts for each subject 
for sj = 1:nS
    
    idx_rt = all_responses(:,11) == sj & all_responses(:,2)<= 3.5; 
    rts  = all_responses(idx_rt,2);

    idx_rt_tr_long = all_responses(idx_rt,9) == 2 |  all_responses(idx_rt,9) == 4;
    idx_rt_tr_short = all_responses(idx_rt,9) == 1 |  all_responses(idx_rt,9) == 3;
    
    idx_rt_tr_freq = all_responses(idx_rt,9) == 1 |  all_responses(idx_rt,9) == 2;
    idx_rt_tr_rare = all_responses(idx_rt,9) == 3 |  all_responses(idx_rt,9) == 4;
    
    tr_len_reg = zeros(length(rts),1); 
    tr_len_reg(idx_rt_tr_long) = -1; 
    tr_len_reg(idx_rt_tr_short) = 1;
    
    tr_freg_reg = zeros(length(rts),1); 
    tr_freg_reg(idx_rt_tr_rare) = -1; 
    tr_freg_reg(idx_rt_tr_freq) = 1;
    
    tr_interaction_reg = tr_freg_reg .* tr_len_reg; 
    
    tr_coh_reg = abs(all_responses(idx_rt,4)); 
    
    rt_nan = isnan(rts); 
    
    
    rts(rt_nan) = []; 
    tr_coh_reg(rt_nan) = []; 
    tr_len_reg(rt_nan) = []; 
    tr_freg_reg(rt_nan) = []; 
    tr_interaction_reg(rt_nan) = []; 
    

  X = [tr_coh_reg tr_len_reg  tr_freg_reg tr_interaction_reg]; 
  
  X(:,1) = tr_coh_reg - mean(tr_coh_reg);
  %X = [tr_coh_reg tr_len_reg  tr_freg_reg]; 
     %for reg = 1:size(X,2)
        
        
    %X(:,reg) = zscore(X(:,reg));
     %end
    
     [b(:,sj),~,stats] = glmfit(X,rts,'normal','constant','off');
     
end 


% ttest
for reg = 1:size(X,2)
[H,P(reg),~,stat(reg)] = ttest(b(reg,:)); 
end 

% 3.5 rts 
% tstat - p-Value 
% coherence -11.62 1e-5 
% length  -4.77 0.0001 
% frequency 13.31 1e-5 


%% correlate median rts with taus from trial periods - mean coherences leading up to button press - step_2 script

% calculate rts across sjs and conditions but without respecting coherences
clear rts 
for sj = 1:nS
    for con = 1:4
       
            
    idx_rts = all_responses(:,9) == con & all_responses(:,7) == 1 & all_responses(:,11) == sj ; % idx to all trials with a  a given coherence, condition , correct response,and a specific subject 
    
   
    rts(con,sj)  = nanmedian(all_responses(idx_rts,2));
       
    end 

end 

% load taus 
load('pnew_sj_TR');
taus = squeeze(pnew_sj(:,2,:))';
mean_taus = mean(taus(:));
sd_taus = 2 * std(taus(:));
[sj_ls,con_id] = find(taus >= (mean_taus + sd_taus));



taus(sj_ls,:) = [];

rts(:,sj_ls) = [];

for con = 1:4
    
    [R,P] = corrcoef(taus(:,con),rts(con,:));
    subplot(2,2,con)
    plot(rts(con,:),taus(:,con),'x','Color',cl(con,:),'LineWidth',3,'MarkerSize',10)
    title([condition{1},' ','R= ',num2str(round(R(1,2),2)),' ','P= ', num2str(round(P(1,2),2))])
    xlabel('median rt (sec)')
    ylabel('tau')
    tidyfig; 
    
   
end 


%% correlate median rts with taus from FAs 

% calculate rts across sjs and conditions but without respecting coherences
%
for sj = 1:nS
    for con = 1:4
       
            
    idx_rts = all_responses(:,9) == con & all_responses(:,7) == 1 & all_responses(:,11) == sj & all_responses(:,2)<= 3.5; % idx to all trials with a  a given coherence, condition , correct response,and a specific subject 
    
   
    rts(con,sj)  = median(all_responses(idx_rts,2));
       
    end 

end 

% load taus 
FA_exp_fits = load('pnew_sj_FA.mat');
FA_taus = squeeze(FA_exp_fits.pnew_sj(:,1,:))';

num_tr = unique(FA_exp_fits.r);
FA_taus(num_tr,:) = [];
rts(:,num_tr) = [];

mean_taus = mean(FA_taus(:));
sd_taus = 2 * std(FA_taus(:));
[sj_ls,con_id] = find(FA_taus >= (mean_taus + sd_taus));

FA_taus(sj_ls,:) = [];

rts(:,sj_ls) = [];

figure;
for con = 1:4
    
    [R,P] = corrcoef(FA_taus(:,con),rts(con,:));
    subplot(2,2,con)
    plot(rts(con,:),FA_taus(:,con),'d','Color',cl(con,:),'LineWidth',2,'MarkerSize',8)
    title([condition{1},' ','R= ',num2str(round(R(1,2),2)),' ','P= ', num2str(round(P(1,2),2))])
    xlabel('median rt (sec)')
    ylabel('tau')
    tidyfig; 
    
   
end 

for con = 1:4
[R{con},P{con},RLO{con},RUP{con}] = corrcoef(FA_taus(:,con),rts(con,:));

end 
