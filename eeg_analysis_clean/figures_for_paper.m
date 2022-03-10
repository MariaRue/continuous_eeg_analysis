% load plot variables 
addpath(genpath(pwd));
options = continuous_RDK_set_options('LTHiMac');
plotStructure; %returns plotVariables to workspace

%% figure 1: plots prediction error effect sorted by long vs. short, rare vs. frequent
make_figure1(plotVariables,options); 

%% figure 2: plots different response-locked ERPs in run-up to making a final choice

make_figure2(plotVariables,options); 

