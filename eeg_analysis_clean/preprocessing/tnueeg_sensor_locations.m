function [ D ] = tnueeg_sensor_locations( D )
%TNUEEG_SENSOR_LOCATIONS Assigns the default locations of channels
%in 3D space
%   IN:     D           - EEG data set
%   OUT:    D           - EEG data set with 3D sensor location information

S = [];

S.D = D;
S.task = 'defaulteegsens';
S.save = 1;

D = spm_eeg_prep(S);


end

