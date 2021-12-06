function [betas] = run_subject_level_convolutional_glm(laggedDesignMatrix, EEGdata, badsamples, VEOG_indx,chanlabels)
%RUN_SUBJECT_LEVEL_CONVOLUTIONAL_GLM
% runs covnolutional GLM usind the pseudoinverse to compute betas for the
% regressors in the lagged design matrix for one block of one subject 
%
%INPUT: 
% laggedDesignMatrix = [nRegressos X EEG samples]
% EEGdata = [chan x time] that is matched to stimulus 
% badsamples = [can x time] indicating badsamples across channels - used
%              for eyeblink and artifact removal in GLM 
% VEOG_indx = channel index into badsamples to remove eyeblink artifacts 
%
% OUTPUT:
% betas = [can x num Regressors]


numChannels = length(chanlabels);

% pre-allocate space
betas = zeros(numChannels,length(laggedDesignMatrix(:,1)));

pDM = geninv(laggedDesignMatrix'); %pseudoinvesrse of design matrix

% calculate betas for each regressor and channel 
for ch = 1:numChannels
   
    % indices of badsamples for this channel and eyeblinks
    removed_samples = badsamples(ch,:) | badsamples(VEOG_indx,:);
    
    betas(ch,:) = pDM(:,~removed_samples)*EEGdata(ch,~removed_samples)';
end




end