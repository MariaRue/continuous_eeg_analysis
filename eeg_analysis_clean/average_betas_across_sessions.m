function betas_all_subjects_sessavg = average_betas_across_sessions (betas_all_subjects, glmFlag, options)

for r = 1:length(options.subjectLevelGLM.(glmFlag).regressors) % loop over regressors
    betas_all_subjects_sessavg{r} = squeeze(mean(betas_all_subjects{r},4));
end


end 