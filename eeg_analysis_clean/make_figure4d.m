function make_figure4d(plotVariables, options)

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


% single subject kernels with exponential fits
[GroupIntegrationKernels, SubjectIntegrationKernels, SignificantTimePoints, ExParameters] = calculate_integration_kernels(all_responses,SubjectListBehaviourEEG,nS, mean_stim_streams, stim_streams, trigger_streams,lags);

[Model] = calculate_integration_kernel_based_on_model(ExParameters, nS);


% correlation of exponential fits

plotVariables.conditions;
figure;set(gcf,'Position',[1000         527        1025         810]);
subplotCounter = 1;
for subject = [19 13 22]
    for condition = 1:4
        
        subplot(3,4,subplotCounter)
        
        if  subplotCounter == 1
            
            ylabel('Mean coherence subject 1')
            title(plotVariables.originalConditions{1})
        elseif subplotCounter == 2
            title(plotVariables.originalConditions{2})
            
        elseif subplotCounter == 3
            title(plotVariables.originalConditions{3})
            
        elseif subplotCounter == 4
            title(plotVariables.originalConditions{4})
            
        elseif subplotCounter == 5
            ylabel('Mean coherence subject 2')
            
        elseif subplotCounter == 9
            ylabel('Mean coherence subject 3')
            
        end
        if subplotCounter>=9
            xlabel('time to False Alarm (s)')
        end
        
        
        hold on
        plot(SubjectIntegrationKernels.DataForModelFit{condition,subject},'Color',plotVariables.figure6.colours(condition,:),'LineWidth',3);
        plot(Model{subject,condition},'k', 'LineWidth', 3);
        Asj   = ExParameters.parameters(condition,1,subject);
        tausj = ExParameters.parameters(condition,2,subject);
        text(50,0.6,sprintf('%c = %0.2fs\nA = %0.2f', 964,tausj,Asj))
        %
        %         txt = ['FA: ', num2str(num_false_alarms(sj,con))];
        %         text(20,0.5,txt,'FontSize',14)
        xticks([0:100:500])
        xticklabels([ 5 4 3 2 1 0])
        % set(gca, 'XDir','reverse')
        tidyfig
        xlim([0 500]);
        ylim([-0.2 0.7]);
        
        hold off
        subplotCounter = subplotCounter+1;
    end
end

end