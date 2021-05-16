function create_set_file(pathSubjectRawEEG,pathSetFile)

	EEG = loadcurry(pathSubjectRawEEG, 'CurryLocations', 'False');

pop_saveset(EEG,'filename',pathSetFile);
end 