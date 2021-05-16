function [lagged_design_matrix] = create_subject_level_design_matrix_for_convolutional_glm(all_regressors,options,glmFlag)
%CREATE_SUBJECT_LEVEL_DESIGN_MATRIX_FOR_CONVOLUTIONAL_GLM 
% creates a design matrix of lagged regressors for the convolutional GLM 
% 
% INPUTS:
% all_regressors is a structure with fields for each regressor containing
% the value [1 * nEEGSamples vector of regressor's value] that is lagged 
%
% options with information: 
%   .nLagsBack [how many lags backward does this regressor go?]
%   .nLagsForward [how many lags forward does this regressor go?]
%   .name [what is this regressor's name?]
%
% glmFlag is a string that specifies the GLM model we want to use and the
% regressors we need 
%
%
% OUTPUTS:
% lagged_design_matrix
% .regressors is a design matrix for the EEG data, with rows (nLagsBack+nLagsForward+1) for each regressor and columns nEEGSamples
% .dm_row_idx contains row indices of each regressor type - this is needed
% to later save the betas from the convolutional glm according to the
% regressor type 



% select correct GLM model 

fields = fieldnames(options.subjectLevelGLM);
glmIdx = strcmp(fields,glmFlag);

% loop through regressors of the GLM 
lengthRegressorList = length(options.subjectLevelGLM.(fields{glmIdx}).regressors);

% loop through regressors of the GLM and make sure that data is in the
% correct format
for regressor = 1:lengthRegressorList 
    regressorName = options.subjectLevelGLM.(fields{glmIdx}).regressors(regressor).name;
    nSamples = length(all_regressors.(regressorName));
    
     if any(size(all_regressors.(regressorName))~=[nSamples 1])
        error('allRegressors must be nEEGsamples * 1 for all regressors');
    elseif round(options.subjectLevelGLM.(fields{glmIdx}).regressors(regressor).nLagsBack)~=options.subjectLevelGLM.(fields{glmIdx}).regressors(regressor).nLagsBack || options.subjectLevelGLM.(fields{glmIdx}).regressors(regressor).nLagsBack<0
        error('nLagsBack must be positive integer or 0');
    elseif round(options.subjectLevelGLM.(fields{glmIdx}).regressors(regressor).nLagsForward)~=options.subjectLevelGLM.(fields{glmIdx}).regressors(regressor).nLagsForward|| options.subjectLevelGLM.(fields{glmIdx}).regressors(regressor).nLagsForward<0
        error('nLagsForward must be positive integer or 0');
    end

end 


%now build design matrix 
row_start = 1; %indicator for where we are in the design matrix rows

for regressor = 1:lengthRegressorList 
    regressorName = options.subjectLevelGLM.(fields{glmIdx}).regressors(regressor).name;
    nSamples = length(all_regressors.(regressorName));
    
    % do we want to keep this???? 
    lagged_design_matrix.regressor_indices(regressor).dm_row_idx = row_start:row_start+options.subjectLevelGLM.(fields{glmIdx}).regressors(regressor).nLagsBack + options.subjectLevelGLM.(fields{glmIdx}).regressors(regressor).nLagsForward;
    
    row_start =  lagged_design_matrix.regressor_indices(regressor).dm_row_idx(end) + 1;
    
    nRows = length( lagged_design_matrix.regressor_indices(regressor).dm_row_idx); %number of rows for this regressor
    
    reg_zp = [zeros(options.subjectLevelGLM.(fields{glmIdx}).regressors(regressor).nLagsForward,1); all_regressors.(regressorName); ...
        zeros(options.subjectLevelGLM.(fields{glmIdx}).regressors(regressor).nLagsBack,1)]; %zero pad regressor
   
    % fill designmatrix for one regressor
    lagged_regressor = nan(nRows,nSamples);
    for j = 1:nRows
        rowidx = nRows - j + 1;
        lagged_regressor(rowidx,:) = reg_zp(j:(j-1+nSamples));
    end
    
    lagged_design_matrix.regressors_matrix( lagged_design_matrix.regressor_indices(regressor).dm_row_idx,:) = lagged_regressor;
    
end 
end 