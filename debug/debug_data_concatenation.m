clearvars; clc;

subject = 'MA25VE';

pattern     = '.mi.';
modality    = '';

experiment  = 'cybathlon';
datapath    = [pwd '/analysis/'];

%% Get datafiles
[Files, NumFiles] = cnbiutil_getfile(datapath, '.mat', [subject '*' modality '*' pattern]);

%% Concatenate data
[F, events, labels, settings] = cnbiutil_concatenate_data(Files);