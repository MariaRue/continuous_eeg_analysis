function [betas] = run_subject_level_convolutional_glm(laggedDesignMatrix, EEGdata, badsamples, VEOG_indx,chanlabels,doReg,lambda,breakPoints)
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
% chanlabels = list of channels
% doReg = 1 if regularisation applied, 0 if not,
% lambda = regularisation hyperparameter
%
% OUTPUT:
% betas = [can x num Regressors]

if nargin<6
    doReg = false;
end
if nargin<7
    warning('need to specify lambda if regularisation is switched on! Regularisation is turned OFF')
    lambda = 0; %equivalent to turning off regularisation
end
if nargin<8
    warning('should specify breakpoints in design matrix when using regularisation!! none specified');
    breakPoints = [];
end


numChannels = length(chanlabels);

% pre-allocate space
betas = zeros(numChannels,length(laggedDesignMatrix(:,1)));

if doReg
    DM = laggedDesignMatrix';
    DM = (DM - mean(DM))./std(DM); %normalise regressors so that a single regularisation parameter can be applied to all
    pDM = pinv_reg(DM,lambda,'onediff',breakPoints);
else
    pDM = geninv(laggedDesignMatrix'); %pseudoinvesrse of design matrix
end

% calculate betas for each regressor and channel
for ch = 1:numChannels
    
    % indices of badsamples for this channel and eyeblinks
    removed_samples = badsamples(ch,:) | badsamples(VEOG_indx,:);
    
    betas(ch,:) = pDM(:,~removed_samples)*EEGdata(ch,~removed_samples)';
end




end