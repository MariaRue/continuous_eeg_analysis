function make_figure4c(plotVariables, options)



%behavioural figure for trial periods
subjectListEEG = [16 18:21, 24, 26, 27, 28, 29, 31,  33, 34, 35, 42, 43, 47, 50, 51, 52, 54, 55, 57, 58]; %32 taken out

lags = 500;
EEGpreproc = options.path.preproc.behaviour;  % path to behav data all subjs
load_name = fullfile(EEGpreproc,'behav_data_all_subjs_all3');
load(load_name)
completeSubjectListBehaviour = unique(all_responses(:,12));
%get behavioural subjects who are in subjectList from EEG
[~,~,SubjectListBehaviourEEG] = intersect(subjectListEEG,completeSubjectListBehaviour');
nS = length(SubjectListBehaviourEEG);



[ExpParameters,Correlation] = calculate_collapsed_integration_kernels(all_responses,SubjectListBehaviourEEG,nS, mean_stim_streams, stim_streams, trigger_streams,lags);



fh = figure;
subplot(2,2,1)
line([0 4],[0 4], 'color', 'k', 'linewidth', 1);
hold on;
plot(ExpParameters.sorted.short.tau,ExpParameters.sorted.long.tau, ...
    'o', 'LineWidth', 1, 'MarkerSize', 10, 'MarkerFaceColor', [0.8 0.8 0.8], 'Color', [1 150 200]/255);
xlabel('\tau short')
ylabel('\tau long')
xlim([0 4]);ylim([0 4]);
xticks(0:4); yticks(0:4);
tidyfig;
hAx = gca; hAx.LineWidth = 1; hAx.TickLength = [0.025 0.025];
box off

subplot(2,2,2)
line([0 1],[0 1], 'color', 'k', 'linewidth', 1);
hold on;
plot(ExpParameters.sorted.short.amplitude,ExpParameters.sorted.long.amplitude, ...
    'o', 'LineWidth', 1, 'MarkerSize', 10, 'MarkerFaceColor', [0.8 0.8 0.8], 'Color', [1 150 200]/255);

xlabel('A short')
ylabel('A long')
xlim([0 1]);ylim([0 1]);
xticks(0:0.5:1); yticks(0:0.5:1);
tidyfig;
hAx = gca; hAx.LineWidth = 1; hAx.TickLength = [0.025 0.025];
box off

subplot(2,2,3)
line([0 4],[0 4], 'color', 'k', 'linewidth', 1);
hold on;
plot(ExpParameters.sorted.frequent.tau,ExpParameters.sorted.rare.tau, ...
    'o', 'LineWidth', 1, 'MarkerSize', 10, 'MarkerFaceColor', [0.8 0.8 0.8], 'Color', [1 150 200]/255);
xlabel('\tau frequent')
ylabel('\tau rare')
xlim([0 4]);ylim([0 4]);
xticks(0:4); yticks(0:4);
tidyfig;
hAx = gca; hAx.LineWidth = 1; hAx.TickLength = [0.025 0.025];
box off;

subplot(2,2,4)
line([0 1],[0 1], 'color', 'k', 'linewidth', 1);
hold on;
plot(ExpParameters.sorted.frequent.amplitude,ExpParameters.sorted.rare.amplitude, ...
    'o', 'LineWidth', 1, 'MarkerSize', 10, 'MarkerFaceColor', [0.8 0.8 0.8], 'Color', [1 150 200]/255);
xlabel('A frequent')
ylabel('A rare')
xlim([0 1]);ylim([0 1]);
xticks(0:0.5:1); yticks(0:0.5:1);
tidyfig;
hAx = gca; hAx.LineWidth = 1; hAx.TickLength = [0.025 0.025];
box off

fh.Position = [1831 243 884 725];

%% STATS TESTING
% frequent vs. rare tau
[h,p,ci,stats] = ttest(ExpParameters.sorted.frequent.tau-ExpParameters.sorted.rare.tau);

% frequent vs. rare A
[h,p,ci,stats] = ttest(ExpParameters.sorted.frequent.amplitude-ExpParameters.sorted.rare.amplitude);

% long vs. short tau
[h,p,ci,stats] = ttest(ExpParameters.sorted.long.tau-ExpParameters.sorted.short.tau);

% long vs. short A
[h,p,ci,stats] = ttest(ExpParameters.sorted.long.amplitude-ExpParameters.sorted.short.amplitude);

keyboard;

end