% load plot variables
addpath(genpath(pwd));
options = continuous_RDK_set_options('LTHiMac');
plotStructure; %returns plotVariables to workspace

%% figure 1: plots prediction error effect sorted by long vs. short, rare vs. frequent
make_figure1(plotVariables,options);

%% figure 2: plots ramp to response, sorted by long vs. short, rare vs. frequent
make_figure2(plotVariables,options);

%% figure 3: split by horizontal and vertical motion (all different regressors; only n=6 subjects)
make_figure3(plotVariables,options);

%% figure 4: cross-subject correlation with Tau?
make_figure4(plotVariables,options);

%% figure 5: think this is just the main effect across all conditions, to show at start of paper
make_figure5(plotVariables, options);

