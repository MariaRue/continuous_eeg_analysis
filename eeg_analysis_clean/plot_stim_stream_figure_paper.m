load('/Volumes/crdkData/rawData/experiment/sub016/stim/sub016_sess001_stim.mat')

mean_coherence = S.mean_coherence_org{4};

coherence_noise = S.coherence_frame{4};
coherence_noise(coherence_noise>1) = 1;
coherence_noise(coherence_noise<-1) = -1;
%% 
figure
plot(coherence_noise,'-','LineWidth',1,'Color',[0.5 0.5 0.5])
hold on 
plot(mean_coherence,'k-','LineWidth',3)


minute = length(mean_coherence)/5; 
half_minute = minute/2; 

xticks([0:half_minute:length(mean_coherence)])
xticklabels({'0','0.5','1','1.5','2','2.5','3','3.5','4','4.5','5'})
xlabel('minutes')