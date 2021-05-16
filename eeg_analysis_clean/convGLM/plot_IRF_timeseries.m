
[dat1_vs_dat2_avg] = calculate_difference_waveform(cond_tr_freq_cp, cond_tr_rare_cp);
[dat1_vs_dat2_avgLen] = calculate_difference_waveform(cond_tr_short_cp, cond_tr_long_cp);

figure;
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [200 100]);
set(gcf, 'Position',  [500, 500, 800, 1290])


cl = cbrewer('qual','Set1',3);


subplot(2,1,1)
hold on
% channel = {'CPz'};
channel = [39,40, 41]; % 22 - FCZ, 40 = CPZ
% [t1,t2] = match_str( stat.label, channel);
hold on
h = shadedErrorBar(cond_tr_freq_avg_cp.time,cond_tr_freq_avg_cp.avg,dat1_vs_dat2_avg.se, 'lineprops', '-k');
h.patch.FaceColor = cl(1,:);
h.mainLine.Color = cl(1,:);

h = shadedErrorBar(cond_tr_rare_avg_cp.time,cond_tr_rare_avg_cp.avg, dat1_vs_dat2_avg.se, 'lineprops', '-k');
h.patch.FaceColor = cl(2,:);
h.mainLine.Color = cl(2,:);
title('difference between freq and rare conditions correct trial CPZ')
legend({'frequent','rare'})
%time points in which either or all of the channes cp1 cp2 cpz are
%significant
xlim([0 0.8])
%ylim([-0.5 1])

%idx = pos>0;
% %
%plot(stat.time(idx),ones(sum(idx),1).*1,'kx')
xlabel('time(s) 0 start of event')
% hold off
tidyfig;


subplot(2,1,2)

hold on
% channel = {'CPz'};
channel = [21,22, 23]; % 22 - FCZ, 40 = CPZ
% [t1,t2] = match_str( stat.label, channel);
hold on
h = shadedErrorBar(cond_tr_short_avg_cp.time,cond_tr_short_avg_cp.avg,dat1_vs_dat2_avgLen.se, 'lineprops', '-k');
h.patch.FaceColor = cl(1,:);
h.mainLine.Color = cl(1,:);

h = shadedErrorBar(cond_tr_long_avg_cp.time,cond_tr_long_avg_cp.avg, dat1_vs_dat2_avgLen.se, 'lineprops', '-k');
h.patch.FaceColor = cl(2,:);
h.mainLine.Color = cl(2,:);
title('difference between freq and rare conditions correct trial FCZ')
legend({'short','long'})
%xlim([-4 0])
xlim([0 0.8])
%ylim([-0.5 1])
%time points in which either or all of the channes cp1 cp2 cpz are
%significant

% idx = pos>0;
% % %
% plot(stat_len.time(idx),ones(sum(idx),1).*0.93,'kx')
% xlabel('time(s) 0 start of event')
% hold off
tidyfig;


%% 


figure 
% 
% figure;
% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperSize', [200 100]);
% set(gcf, 'Position',  [500, 500, 800, 1290])


cl = cbrewer('qual','Set1',3);


hold on
% channel = {'CPz'};
channel = [39,40, 41]; % 22 - FCZ, 40 = CPZ
% [t1,t2] = match_str( stat.label, channel);
hold on
h = shadedErrorBar(cond_tr_all_avg_cp1.time,cond_tr_all_avg_cp1.avg,cond_tr_all_avg_cp1.se, 'lineprops', '-k');
h.patch.FaceColor = cl(1,:);
h.mainLine.Color = cl(1,:);

h = shadedErrorBar(cond_tr_all_avg_cp2.time,cond_tr_all_avg_cp2.avg, cond_tr_all_avg_cp2.se, 'lineprops', '-k');
h.patch.FaceColor = cl(2,:);
h.mainLine.Color = cl(2,:);
title('difference between freq and rare conditions correct trial CPZ')
legend({'decision relevant','decision irrelevant'})
%time points in which either or all of the channes cp1 cp2 cpz are
%significant
xlim([-0.1 0.8])
%ylim([-0.5 1])

%idx = pos>0;
% %
%plot(stat.time(idx),ones(sum(idx),1).*1,'kx')
xlabel('time(s) 0 start of event')
% hold off
tidyfig;
