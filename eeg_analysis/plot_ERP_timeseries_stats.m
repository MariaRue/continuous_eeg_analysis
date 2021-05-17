function plot_ERP_timeseries_stats(pos,stat,data_cp1,data_cp2,data_cp_all1,data_cp_all2) %%% plot significance for single channels

figure


cl = cbrewer('qual','Set1',3);

[dat1_vs_dat2_avg] = calculate_difference_waveform(data_cp_all1, data_cp_all2);


subplot(2,1,1)
hold on 
%channel = { 'CP1','CPz','CP2'};
%channel = { 'FP1','FPz','FP2'};
%[t1,t2] = match_str( stat.label, channel);
hold on
h = shadedErrorBar(data_cp1.time,data_cp1.avg,dat1_vs_dat2_avg.se, 'lineprops', '-k');
h.patch.FaceColor = cl(1,:);
h.mainLine.Color = cl(1,:);
h.mainLine.LineWidth = 0.5;
h.patch.FaceAlpha = 0.3;

h = shadedErrorBar(data_cp2.time,data_cp2.avg, dat1_vs_dat2_avg.se, 'lineprops', '-k');
h.patch.FaceColor = cl(2,:);
h.mainLine.Color = cl(2,:);
h.mainLine.LineWidth = 0.5;
h.patch.FaceAlpha = 0.3;
%title('difference between freq and rare conditions CPz')
title('difference between short and long conditions CPz')
legend({'short','long'})
%legend({'frequent','rare'})
% time points in which either or all of the channes cp1 cp2 cpz are
% significant
% ylim([-6 6])

%idx = sum(s)>0;
idx = pos>0;
plot(stat.time(idx),ones(sum(idx),1).*1,'k*','MarkerSize',2)
xlabel('time(s) 0 start of event')
tidyfig;
hold off



subplot(2,1,2)
channel = {'FC1', 'FCz', 'FC2'};
%channel = {'AF1', 'AF2', 'AFz'};
[t1,t2] = match_str( stat.label, channel);
hold on
h = shadedErrorBar(dat1_vs_dat2_avg.time,dat1_vs_dat2_avg.avg,dat1_vs_dat2_avg.se, 'lineprops', '-k');
h.patch.FaceColor = cl(1,:);
h.mainLine.Color = cl(1,:);
h.mainLine.LineWidth = 0.5;
h.patch.FaceAlpha = 0.3;


% time points in which either or all of the channes cp1 cp2 cpz are
% significant

idx = pos>0;
plot(stat.time(idx),ones(sum(idx),1).* 2,'k*','MarkerSize',2)
%title('difference freq vs rare FCZ')
title('difference wave cpz')
xlabel('time(s) 0 start of event')
%legend({'frequent','rare'})

tidyfig; 
hold off
end