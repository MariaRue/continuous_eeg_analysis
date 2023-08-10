% load plot variables
addpath(genpath(pwd));
options = continuous_RDK_set_options('LTHiMac');
plotStructure; %returns plotVariables to workspace

%% figure 2: just the main effect across all conditions, to show at start of paper
make_figure2(plotVariables, options); %n.b. used to be make_figure5.m

%% figure 3a: figure to show example of 2*2 interaction design

make_figure3a(plotVariables, options);

%% figure 3 bcd: behavioural figure during response preiods: detection rate, RTs, Integration kernels
make_figure3bcd(plotVariables, options);

%% figure 4 ab:  behavioural figure during false alarms: detection rate, RTs, Integration kernels
make_figure4ab(plotVariables, options);

%% figure 4c: tau and amplitude fits for each condition
make_figure4c(plotVariables, options);

%% figure 5: plots prediction error effect sorted by long vs. short, rare vs. frequent
make_figure5(plotVariables,options);

%% figure 6: plots ramp to response, sorted by long vs. short, rare vs. frequent
make_figure6(plotVariables,options); 

%for response to reviewers, see also:
%make_figure6_tf(plotVariables,options); 
%make_figure6_tf_lateralised(plotVariables,options); 

%% figure 7: split by horizontal and vertical motion (all different regressors; only n=6 subjects)
make_figure7(plotVariables,options);

%% figure 8: cross-subject correlation with Tau?
make_figure8(plotVariables,options);

%% OLD figure 7: signal period behav figure
%make_figure7(plotVariables, options)

%% OLD figure 8: baseline period behav figure collapsed false alarm rate only 
%make_figure8(plotVariables, options)

%% OLD figure 9: baseline period behav figure fa rate and integration kernels not collapsed 
%make_figure9(plotVariables, options)

%% OLD figure 10: single subject integration kernels 
%make_figure10(plotVariables, options)

%% OLD figure 10b: single subject Integration kernels correlation of model fit parameters 
%make_figure10b(plotVariables, options)

