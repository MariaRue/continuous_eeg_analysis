function ReactionTimes = calculate_reaction_times(responses,SubjectList,nS,coherence) 

for subject = 1:nS
    subjectID = SubjectList(subject);
    for con = 1:4
        for coh = 1:3
            
            idx_rts = abs(responses(:,4)) == coherence(coh) & responses(:,9) == con & responses(:,7) == 1 & responses(:,11) == subjectID ; % idx to all trials with a  a given coherence, condition , correct response,and a specific subject and only reactiones <= 3.5s
            
            ReactionTimes.subjectLevel(coh,con,subject)  = nanmedian(responses(idx_rts,2));
        end
    end
    
end

ReactionTimes.groupLevel = nanmean(ReactionTimes.subjectLevel,3); %may be nans if subject failed to make any responses in a given coherence level
ReactionTimes.groupLevelSe = nanstd(ReactionTimes.subjectLevel,[],3)./sqrt(nS); %may be nans if subject failed to make any responses in a given coherence level


%% now perform 3-way ANOVA for factors of FREQUENCY, LENGTH, COHERENCE

slDataforANOVA = permute(ReactionTimes.subjectLevel,[3 2 1]);

varNames = {'Y1' 'Y2' 'Y3' 'Y4' 'Y5' 'Y6' 'Y7' 'Y8' 'Y9' 'Y10' 'Y11' 'Y12'};
t = array2table(slDataforANOVA(:,:),'VariableNames',varNames);

factorNames = {'Freq','Length','Coherence'};
within = table({'F';'F';'R';'R';'F';'F';'R';'R';'F';'F';'R';'R'},...
               {'S';'L';'S';'L';'S';'L';'S';'L';'S';'L';'S';'L'},...
               {'30';'30';'30';'30';'40';'40';'40';'40';'50';'50';'50';'50'},'VariableNames',factorNames); %F = frequent, R = Rare, S = short, L = long

% fit the repeated measures model
rm = fitrm(t,'Y1-Y12~1','WithinDesign',within);
ReactionTimes.ranovatbl = ranova(rm, 'WithinModel','Freq*Length*Coherence');

end