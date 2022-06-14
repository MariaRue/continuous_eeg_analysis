function make_figure3a(plotVariables,options)

% figure to illustrate the 2*2 design from an example subject.

subID = 16;

reference = 'LMRM';
csdFlag = 0;

[details,paths] =  conrdk_subjects(subID,options,reference,csdFlag); 
file_ID_behaviour = paths.behaviour(1).sessionList;

response_mat= load(file_ID_behaviour,'respMat','B'); % responses

figure;
set(gcf,'Position',[138         300        1387         770]);
for i = 1:4; 
    subplot(2,2,i); 
    ss = response_mat.B.mean_coherence{i};
    dd = [diff(ss); 0];
    
    %now, some code to correct for the fact that the stimulus in ss is
    %truncated whenever the subject presses the button
    stim_onsets = find(abs(dd)>0&ss==0);
    stim_offsets = find(abs(dd)>0&abs(ss)>0);
    stim_dur = max(stim_offsets-stim_onsets); disp(stim_dur);
    for so = 1:length(stim_onsets)
        ss(stim_onsets(so):(stim_onsets(so)+stim_dur)) = ss(stim_onsets(so)+1);
    end
    
    time = [1:length(ss)]/100/60;
    ll = plot(time,ss); 
    set(ll,'LineWidth',2,'Color','k'); box off;
    set(gca,'XTick',0:1:5);
    xlim([0 5]); ylim([-1 1]); 
    if i==4
        xlabel('Time (minutes)');
    elseif i==3
        xlabel('Time (minutes)');
        ylabel('Response periods SHORT');
    elseif i==1
        title('Repsonse periods FREQUENT');
        ylabel('Response periods LONG');
    elseif i==2
        title('Response periods RARE');
    end
    
    hold on;
    tidyfig;
    %plot(stim_onsets/100,zeros(size(stim_onsets)),'r.');
end