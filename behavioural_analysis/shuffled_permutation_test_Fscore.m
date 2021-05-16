
function [sig_tps,p,tbl] = shuffled_permutation_test_Fscore(mean_coherences, lags, thres, repetitions, anova_samples, plt,nS)
% loop through time points and calculate F-scores for original collect
% labelling

Fscore_org = zeros(lags+1,3);

labels = {'Fscore resampled frequency','Fscore resampled trial length', 'Fscore resampled interaction'};

for t = 1:lags + 1

    data =  squeeze(mean_coherences(t,:,:))'; % get data at one time point
    
    
    y = [[data(:,1);data(:,2)],[data(:,3);data(:,4)]]; % arrange for 2-way anova - columns = frequent, fare/ first nS rows = short trials, second nS rows = long trials
    
    %%% 2 way anova
    
    [p,tbl] = anova2(y,anova_samples,'off');
    
    Fscore_org(t,1) = tbl{2,5}; % interaction columns here between frequent and rare
    Fscore_org(t,2) = tbl{3,5}; % interaction short vs long
    Fscore_org(t,3) = tbl{4,5}; % interaction of trial length with frequency
    
end



%% now switch labels



Fscore_resampled = zeros(lags+1,3,repetitions);

for rep = 1:repetitions % reshuffle labels for each participant 100 times per each timepoint and repeat anova and save F-score
    
    if rep == 500 
        disp('half of the repetitions done')
    end 
    
    for t = 1:lags+1
        

        perm_vec = zeros(size(data));
        
        data =  squeeze(mean_coherences(t,:,:))'; % get data at one time point
        
        for sj = 1 : nS
            perm_vec(sj,:) = data(sj,randperm(4,4));
            
        end
        
        y =  [[perm_vec(:,1);perm_vec(:,2)],[perm_vec(:,3);perm_vec(:,4)]];
        [p,tbl] = anova2(y,anova_samples,'off');
        
        Fscore_resampled(t,1,rep) = tbl{2,5};
        Fscore_resampled(t,2,rep) = tbl{3,5};
        Fscore_resampled(t,3,rep) = tbl{4,5};
        
        
    end
end

disp('repetitions done')

%% select highest F for each time point and make new distribution for each F

for t = 1:lags+1
    highest_Fscore(t,1) = max(max(Fscore_resampled(t,1,:)));
    highest_Fscore(t,2) = max(max(Fscore_resampled(t,2,:)));
    highest_Fscore(t,3) = max(max(Fscore_resampled(t,3,:)));
end 




%% transform Fscore_resampled into distributions and test whether original Fscore is 


for id = 1:3

 



if plt
    figure (1)
    subplot(3,1,id)
    h = histogram(highest_Fscore(:,id),12);
    ylabel('count')
    xlabel(labels{id})

end 



[pd,edges] = histcounts(highest_Fscore(:,id),12,'Normalization','pdf');
binCenters = edges + (h.BinWidth/2);
F_sort = sort(squeeze(highest_Fscore(:,id))); 

try 
thres_F = F_sort(lags+1-thres);
catch 
    
    keyboard; 
end


if plt
figure (2)
subplot(3,1,id)
hold on 
plot(binCenters(1:end-1), pd, 'r-')
hold off 
end 

sig_timepoints = [];
for t = 1:lags+1

if plt
    figure (2)
subplot(3,1,id)
hold on
plot([ones(2,1)*thres_F],[0 0.25],'k-')
plot(Fscore_org(t,id),0.2,'xg','MarkerSize',10)

ylabel('pdf','FontSize',14)
xlabel('F score','FontSize',14)
title(labels{id})
hold off
end 


if Fscore_org(t,id) > thres_F
    
    sig_timepoints = [sig_timepoints,t];
    
end 


end

subplot(3,1,id)


sig_tps{id} = sig_timepoints; 
end



end 