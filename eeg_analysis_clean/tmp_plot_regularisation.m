%mini script to plot effects of switching on regularisation

dat = load(paths.(reference).subjectLevelGLM.(glmFlag).saveName);
dd = dat.betas_subject{4}(40,:,:,:);
ddm = mean(mean(dd,4),3);

[pp,nn,ee] = fileparts(paths.(reference).subjectLevelGLM.(glmFlag).saveName);
nn = [nn '_regularised'];
datr = load((fullfile(pp,[nn ee])));
dd_reg = datr.betas_subject{4}(40,:,:,:);
ddm_reg = mean(mean(dd_reg,4),3);

figure;plot(ddm_reg,'r');hold on;plot(ddm,'k');