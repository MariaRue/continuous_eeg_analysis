%-- eye blink detection ----------------------------------------------%
% subject-specific EB detection thresholds
options.preproc.eyeblinkthreshold = details.eyeblinkthreshold;

[Dm, trialStats.numEyeblinks] = tnueeg_eyeblink_detection_spm(D, options);
savefig(details.ebdetectfig);
fprintf('\nEye blink detection done.\n\n');

%-- eye blink epoching and confounds ----------------------------------------------------------%
Db = tnueeg_epoch_blinks(Dm, options);
fprintf('\nEpoching to eye blinks done.\n\n');

% preclean the EB epochs to get better confound estimates
Dbc = dprst_preclean_eyeblink_epochs(Db, details, options);

Dbc = tnueeg_get_spatial_confounds(Dbc, options);
savefig(details.ebspatialfig);
fprintf('\nSpatial confounds done.\n\n');

%-- eye blink correction & headmodel ---------------------------------%
if options.preproc.preclean.doBadChannels
    badPrecleanChannels = badchannels(Dbc);
    if ~isempty(badPrecleanChannels)
        fprintf('\nMarking %s channels as bad due to precleaning.', ...
            num2str(numel(badPrecleanChannels)));
        D = badchannels(D, badPrecleanChannels, ones(1, numel(badPrecleanChannels)));
    end
end
D = tnueeg_add_spatial_confounds(D, Dbc, options);
fprintf('\nAdding spatial confounds done.\n\n');

Dc = tnueeg_eyeblink_correction_simple(D, options);
fprintf('\nSimple eye blink correction done.\n\n');
