% load plot variables
addpath(genpath(pwd));
options = continuous_RDK_set_options('iMac');
plotStructure; %returns plotVariables to workspace

%% figure 1: plots prediction error effect sorted by long vs. short, rare vs. frequent
make_figure1(plotVariables,options);

%% figure 2: plots ramp to response, sorted by long vs. short, rare vs. frequent
make_figure2(plotVariables,options);

%% figure 3: split by horizontal and vertical motion (all different regressors; only n=6 subjects)
make_figure3_LOCAL_5567(plotVariables,options);

%% figure 4: cross-subject correlation with Tau?
make_figure4(plotVariables,options);

%% figure 5: think this is just the main effect across all conditions, to show at start of paper
make_figure5(plotVariables, options);

%% figure 6: behavioural figure detection rate, FAs, RTs, Integration kernels
make_figure6(plotVariables, options)

%% figure 7: signal period behav figure
make_figure7(plotVariables, options)

%% figure 8: baseline period behav figure
make_figure8(plotVariables, options)