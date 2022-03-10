% load plot variables
addpath(genpath(pwd));
options = continuous_RDK_set_options('LTHiMac');
plotStructure; %returns plotVariables to workspace

%% figure 1: plots prediction error effect sorted by long vs. short, rare vs. frequent
make_figure1(plotVariables,options);

%%
make_figure2(plotVariables,options);

%%
make_figure3(plotVariables,options);

%%
make_figure4(plotVariables,options);

%%
make_figure5(plotVariables, options);

%%
