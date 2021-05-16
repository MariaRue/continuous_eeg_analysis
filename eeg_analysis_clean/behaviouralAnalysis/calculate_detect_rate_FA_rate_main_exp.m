%% add paths
addpath(genpath('/Users/maria/Documents/Matlab/cbrewer'))
addpath('/Users/maria/Documents/MATLAB/raacampbell-shadedErrorBar-9b19a7b')
EEGpreproc = '/Volumes/crdkData/preprocessedData/behaviour'; % path to behav data all subjs

condition = {'Tr frequent TR short', 'Tr frequent Tr long','Tr rare Tr short','Tr rare Tr long'};
% plot each subject separately for conditions
cl = cbrewer('div','RdBu', 12);
cl = cl([2 4 12 10],:);

%% load data
load_name = fullfile(EEGpreproc,'behav_data_all_subjs_all3');
load(load_name)

lags = 500; % frames back in time that are leading up to FA

nS = max(all_responses(:,11)); % number of subjects

coherence = [0.3 0.4 0.5];

%% calculate FA rate
for sj = 1:nS
    se_total = 6; 
    if sj == 26 
        se_total = 5; 
    end 
    for se = 1:se_total
        for c = 1:4
            
            
            time_ITI_sec =   sum(mean_stim_streams{sj,se}(:,c) == 0)/100;
            
            
            FA_num = sum(all_responses(:,7) == 2 & all_responses(:,9) == c & all_responses(:,10) == se & all_responses(:,11) == sj);
            
            FA_rate(se, sj, c) = FA_num/time_ITI_sec;
            
            
        end
    end
    
    for c = 1:4
    % mean per subject
FA_rate_sj(c,sj) = mean(squeeze(FA_rate(:,sj,c))); 
    end
end
mean_FA_rate_total = mean(FA_rate_sj(:)); 
sd_FA_rate_total = 2 * std(FA_rate_sj(:));

[~,sj_id] = find(FA_rate_sj>= mean_FA_rate_total + sd_FA_rate_total | FA_rate_sj<= mean_FA_rate_total - sd_FA_rate_total );
FA_rate_sj(:,sj_id) = [];

% mean across subjects
for c = 1:4
    
FA_rate_all(c) = mean(FA_rate_sj(c,:)); 
se_FA_rate_all(c) = std(FA_rate_sj(c,:))/sqrt(nS);

end 
% figure
% bar(FA_rate_all)
% hold on
% errorbar(FA_rate_all,se_FA_rate_all,'.');
% xticklabels(condition)
% tidyfig;
%%% 2 way anova 
y = [[FA_rate_sj(1,:)';FA_rate_sj(2,:)'],[FA_rate_sj(3,:)';FA_rate_sj(4,:)']];

[p,tbl,stats] = anova2(y,length(FA_rate_sj(1,:)));

figure 
hold on 
for condition = 1:4

 bar(condition, FA_rate_all(condition), 'FaceColor',cl(condition,:),'EdgeColor',cl(condition,:))
 ln = plot([ones(27,1).*condition],squeeze(FA_rate_sj(condition,:)),'k.','MarkerSize',7,'LineWidth',3);
 
end 

%% detect rate
subj_list_mistake = [];
for sj = 1:nS
    
    for c = 1:4
        se_total = 6;
        
        if sj == 26
            se_total = 5; 
        end 
        
        
        for se = 1:se_total

            all_trials_idx = all_responses(:,9) == c & all_responses(:,10) == se & all_responses(:,11) == sj;
            all_trials = all_responses(all_trials_idx,[4,3,6,7]);
            
            
         
            for coh = 1:3
         
                num_of_trials = sum(abs(all_trials(:,1)) == coherence(coh));
                
                num_of_trials_1 = sum(abs(mean_stim_streams{sj,se}(2:end,c)) == coherence(coh) & abs(mean_stim_streams{sj,se}(1:end-1,c)) ==0 );

%                
                

                hit_trials =  sum(abs(all_trials(:,1)) == coherence(coh) & all_trials(:,4) ==1);
                 
                
                detect_rate_se(coh,se) = hit_trials/num_of_trials;
                
                ntr(sj,c,se,coh) = num_of_trials; %from the response matrix
                ntr1(sj,c,se,coh) = num_of_trials_1; %from the stimulus stream
                
                if num_of_trials~=num_of_trials_1
                    keyboard;
                end
                
            end
            
            
            
            
        end
        
        detect_rate_sj(1,c,sj) = nanmean(detect_rate_se(1,:));
        detect_rate_sj(2,c,sj) = nanmean(detect_rate_se(2,:));
        detect_rate_sj(3,c,sj) = nanmean(detect_rate_se(3,:));
        
    end
    
    
    
end

for c = 1:4
    for coh = 1:3
        
        
        detect_rate(coh,c) = nanmean(squeeze(detect_rate_sj(coh,c,:)));
        se_detect_rate(coh,c) = std(squeeze(detect_rate_sj(coh,c,:)))/sqrt(nS);
        
    end
end

figure (1)
hold on

errorbar(coherence,detect_rate(:,1), se_detect_rate(:,1),'x-','Color',cl(1,:),'LineWidth',3)
errorbar(coherence,detect_rate(:,2), se_detect_rate(:,2),'x-','Color',cl(2,:),'LineWidth',3)
errorbar(coherence,detect_rate(:,3), se_detect_rate(:,3),'x-','Color',cl(3,:),'LineWidth',3)
errorbar(coherence,detect_rate(:,4), se_detect_rate(:,4),'x-','Color',cl(4,:),'LineWidth',3)
xlim([0.2 0.6])

legend(condition)
tidyfig;



%% GLM for detect rate 

Y = detect_rate_sj(:);

coh_reg = zeros(length(Y),1);
coh_reg(1:3:end) = 0.3; 
coh_reg(2:3:end) = 0.4; 
coh_reg(3:3:end) = 0.5; 

tr_len_reg = zeros(size(coh_reg));
tr_len_reg(1:6:end) = 1; 
tr_len_reg(2:6:end) = 1; 
tr_len_reg(3:6:end) = 1; 
tr_len_reg(4:6:end) = -1; 
tr_len_reg(5:6:end) = -1; 
tr_len_reg(6:6:end) = -1; 

freq_reg = zeros(size(coh_reg));
freq_reg(1:12:end) = 1; 
freq_reg(2:12:end) = 1; 
freq_reg(3:12:end) = 1; 
freq_reg(4:12:end) = 1; 
freq_reg(5:12:end) = 1; 
freq_reg(6:12:end) = 1; 

freq_reg(7:12:end) = -1; 
freq_reg(8:12:end) = -1; 
freq_reg(9:12:end) = -1; 
freq_reg(10:12:end) = -1; 
freq_reg(11:12:end) = -1; 
freq_reg(12:12:end) = -1; 



sj_mean = zeros(length(coh_reg),nS);
r_start = 1;
for sj = 1:nS
    sj_mean(r_start : r_start + 11,sj) = 1; 
    r_start = r_start + 12; 
    
end 


X = [coh_reg, freq_reg, tr_len_reg, sj_mean];

[b,~,stats] = glmfit(X,Y,'normal','constant','off');
%% random effects analysis detect rate 

for sj = 1:nS
    
    dt_rate = detect_rate_sj(:,:,sj);
Y = dt_rate(:);

coh_reg = zeros(length(Y),1);
coh_reg(1:3:end) = 0.3; 
coh_reg(2:3:end) = 0.4; 
coh_reg(3:3:end) = 0.5; 

tr_len_reg = zeros(size(coh_reg));
tr_len_reg(1:6:end) = 1; 
tr_len_reg(2:6:end) = 1; 
tr_len_reg(3:6:end) = 1; 
tr_len_reg(4:6:end) = -1; 
tr_len_reg(5:6:end) = -1; 
tr_len_reg(6:6:end) = -1; 

freq_reg = zeros(size(coh_reg));
freq_reg(1:12:end) = 1; 
freq_reg(2:12:end) = 1; 
freq_reg(3:12:end) = 1; 
freq_reg(4:12:end) = 1; 
freq_reg(5:12:end) = 1; 
freq_reg(6:12:end) = 1; 

freq_reg(7:12:end) = -1; 
freq_reg(8:12:end) = -1; 
freq_reg(9:12:end) = -1; 
freq_reg(10:12:end) = -1; 
freq_reg(11:12:end) = -1; 
freq_reg(12:12:end) = -1; 



tr_interaction_reg = freq_reg .* tr_len_reg; 


X = [coh_reg, freq_reg, tr_len_reg, tr_interaction_reg];

[b(:,sj),~,stats(sj)] = glmfit(X,Y,'normal');
end 
%%

for reg = 1:5
    [Hs,Ps(reg),~,statss(reg)] = ttest(b(reg,:));
end 



%% plot detect rate per subject 
for sj = 1:nS 
    
    subplot(7,4,sj)
    for c = 1:4
    
     hold on 
     
    plot(coherence, squeeze(detect_rate_sj(:,c,sj)),'-x', 'Color', cl(c,:), 'LineWidth',3)
    xlim([0.2 0.6])

    hold off 
    title(sj)
    end 
    
    if sj == 1
        legend(condition)
        ylabel('detect rate')
        xlabel('coherence')
    end 
    
    
end 

%%  dectrate exclude trials with rts < 3.5 sec 


for sj = 1:nS
    
    for c = 1:4
        
        for se = 1:6
            
            if c == 2 || c == 4
            all_trials_idx = all_responses(:,2) <= 3.5 & all_responses(:,9) == c & all_responses(:,10) == se & all_responses(:,11) == sj;
            all_trials = all_responses(all_trials_idx,[4,7]);
            else 
            all_trials_idx = all_responses(:,9) == c & all_responses(:,10) == se & all_responses(:,11) == sj;
            all_trials = all_responses(all_trials_idx,[4,7]);                
            end 
            

            for coh = 1:3
                
                num_of_trials = sum(abs(all_trials(:,1)) == coherence(coh));
                hit_trials =  sum(abs(all_trials(:,1)) == coherence(coh) & all_trials(:,2) ==1);
                
                
                detect_rate_se(coh,se) = hit_trials/num_of_trials;
                
                
                
            end
            
            
            
            
        end
        
        detect_rate_sj(1,c,sj) = nanmean(detect_rate_se(1,:));
        detect_rate_sj(2,c,sj) = nanmean(detect_rate_se(2,:));
        detect_rate_sj(3,c,sj) = nanmean(detect_rate_se(3,:));
        
    end
    
    
    
end

for c = 1:4
    for coh = 1:3
        
        
        detect_rate(coh,c) = nanmean(squeeze(detect_rate_sj(coh,c,:)));
         se_detect_rate(coh,c) = std(squeeze(detect_rate_sj(coh,c,:)))/sqrt(nS);
        
    end
end

figure

hold on

errorbar(coherence,detect_rate(:,1), se_detect_rate(:,1),'x-','Color',cl(1,:),'LineWidth',3)
errorbar(coherence,detect_rate(:,2), se_detect_rate(:,2),'x-','Color',cl(2,:),'LineWidth',3)
errorbar(coherence,detect_rate(:,3), se_detect_rate(:,3),'x-','Color',cl(3,:),'LineWidth',3)
errorbar(coherence,detect_rate(:,4), se_detect_rate(:,4),'x-','Color',cl(4,:),'LineWidth',3)
xlim([0.25 0.55])

title('detect rate truncated')
