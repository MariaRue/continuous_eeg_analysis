function [ D ] = tnueeg_sensor_locations( D )
%TNUEEG_SENSOR_LOCATIONS Loads a file indicating the locations of channels
%in 3D space
%   IN:     D           - EEG data set
%           elecFile    - filename of the sensor locations file incl. path
%   OUT:    D           - EEG data set with 3D sensor location information

S = [];

S.D = D;
S.task = 'defaulteegsens';
%S.source = 'mat';
%S.sensfile = options.preproc.path.layout;
S.save = 1;

D = spm_eeg_prep(S);


end

