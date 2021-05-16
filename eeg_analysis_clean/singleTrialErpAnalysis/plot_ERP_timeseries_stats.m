function plot_ERP_timeseries_stats(pos,stat,data_cp1,data_cp2,data_cp_all1,data_cp_all2) %%% plot significance for single channels


figure

    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [200 100]);    
    set(gcf, 'Position',  [500, 500, 800, 1290])
 
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
xlim([-0.5 8])
plot(zeros(10,1),[linspace(-3,3,10)],'k-')
title('difference between freq and rare conditions')
%title('difference between short and long conditions CPz')
legend({'short','long'})
%legend({'frequent','rare'})
% time points in which either or all of the channes cp1 cp2 cpz are
% significant
ylim([-4 4])
xlim([-4 0.1])
%ylim([-3e-4 3e-4])
%idx = sum(s)>0;
idx = pos>0;
plot(stat.time(idx),ones(sum(idx),1).*4e-4,'kx','MarkerSize',6)
xlabel('time(s) 0 start of event')
tidyfig;
hold off
set(gca,'FontName','Arial')


subplot(2,1,2)
channel = {'FC1', 'FCz', 'FC2'};
%channel = {'AF1', 'AF2', 'AFz'};
%[t1,t2] = match_str( stat.label, channel);
hold on
h = shadedErrorBar(dat1_vs_dat2_avg.time,dat1_vs_dat2_avg.avg,dat1_vs_dat2_avg.se, 'lineprops', '-k');
% h = shadedErrorBar(cond_tr_all_avg.time,cond_tr_all_avg.avg,cond_tr_all_avg.se, 'lineprops', '-k');
% h = shadedErrorBar(cond_avg_ALL.time,cond_avg_ALL.avg,cond_avg_ALL.se, 'lineprops', '-k');
h.patch.FaceColor = cl(1,:);
h.mainLine.Color = cl(1,:);
h.mainLine.LineWidth = 0.5;
h.patch.FaceAlpha = 0.3;
xlim([-0.5 8])
%ylim([-2.1 2])
 ylim([-2 2])
%ylim([-6e-4 3e-4])
% time points in which either or all of the channes cp1 cp2 cpz are
% significant
xlim([-4 0.1])
idx = pos>0;
plot(stat.time(idx),ones(sum(idx),1).* 2e-4,'kx','MarkerSize',6)
plot(zeros(10,1),[linspace(-2,2,10)],'k-')
%title('difference freq vs rare FCZ')
title('difference wave')
xlabel('time(s) 0 start of event')
%legend({'frequent','rare'})
set(gca,'FontName','Arial')
tidyfig; 
hold off
end