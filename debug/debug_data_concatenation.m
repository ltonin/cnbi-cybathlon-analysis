clearvars; clc;

subject = 'AN14VE';

pattern     = '.mi.';
modality    = 'offline';

experiment  = 'cybathlon';
datapath    = [pwd '/analysis/'];

%% Get datafiles
[Files, NumFiles] = cnbiutil_getfile(datapath, '.mat', [subject '*' modality '*' pattern]);

%% Concatenate data
[F, events, labels, settings] = cnbiutil_concatenate_data(Files(1:10));
NumSamples = size(F, 1);

%% Extract events
CueEventTypes = [769 770 771 773 783];
[DbgCueLabels, CueEvents] = cnbiproc_get_event(CueEventTypes, NumSamples, events.POS, events.TYP, events.DUR);

CFbType = 781;
[DbgCFbLabels, CFbEvents] = cnbiproc_get_event(CFbType, NumSamples, events.POS, events.TYP, events.DUR);

[DataLabels, DataEvents] = cnbiproc_get_event(CueEventTypes, NumSamples, CueEvents.POS, CueEvents.TYP, CFbEvents.DUR + CueEvents.DUR - 1, CueEvents.DUR);





%% Plotting
fig = figure;

NumRows = 1;
NumCols = 1;

subplot(NumRows, NumCols, 1);
hold on;
plot(DbgCueLabels); 
plot(DbgCFbLabels); 
plot(DataLabels); 
hold off
title('Debug - Comparison labels');
legend('Cue', 'CFeedback', 'Trial');