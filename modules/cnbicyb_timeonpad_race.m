clearvars; clc; 

subject = 'AN14VE';
%subject = 'MA25VE';

pattern     = '.mi.';
modality    = 'race';

experiment  = 'cybathlon';
datapath    = [pwd '/analysis/'];
figuredir  = './figures/';

PadTypeId = [768 769 770 771 773 783];
PadTypeLb = {'Slide', 'Slide', 'Speed', 'Jump', 'Speed', 'Rest'};
PadTypeInd = [3 3 1 2 1 4];

CommTypeId = PadTypeId + hex2dec('6000');
CommTypeLb = {'Slide', 'Slide', 'Speed', 'Jump', 'Speed', 'Rest'};

%% Get datafiles
[Files, NumFiles] = cnbiutil_getfile(datapath, '.mat', [subject '*' modality '*' pattern]);

%% Concatenate data
cnbiutil_bdisp(['[io] - Import psd datafiles (' modality ')']);
[U, events, labels, settings] = cnbiutil_concatenate_data(Files);

DataLength  = size(U, 1);

Dk = labels.Dk;
Days    = unique(Dk);
NumDays = length(Days);
Dl = labels.Dl;

Mk = labels.Mk;

%% Compute the overall time-on-pad distributions
[TOPPad TOPTask] = cnbiproc_timeonpad(events.extra.pad, PadTypeId, PadTypeLb);

%% Compute time-on-pad per day

for dId = 1:NumDays
    KeepInd = [];
    for tr=1:length(events.extra.pad.TYP)
        if(Dk(events.extra.pad.POS(tr))==dId)
            KeepInd = [KeepInd; tr];
        end
    end
    thisSessionTOP.POS = events.extra.pad.POS(KeepInd);
    thisSessionTOP.TYP = events.extra.pad.TYP(KeepInd);
    thisSessionTOP.DUR = events.extra.pad.DUR(KeepInd);
    [TOPPadSes{dId} TOPTaskSes{dId}] = cnbiproc_timeonpad(thisSessionTOP, PadTypeId, PadTypeLb);
    
end


%% Plotting
fig1 = figure;
cnbifig_set_position(fig1, 'All');
plotpad = [];
plotlbl = [];
for p=1:5
    plotpad = [plotpad TOPPad{p}'];
    plotlbl = [plotlbl p*ones(1,length(TOPPad{p}))];
end
boxplot(plotpad,plotlbl,'labels',{'Speed','Jump','Slide','Rest','Start/End'});
xlabel('Pad Type','FontSize',20,'LineWidth',3);
ylabel('Time On Pad (sec)','FontSize',20,'LineWidth',3);
title(subject);
axis([0 6 0 21]);
set(gca,'FontSize',20,'LineWidth',3);
cnbifig_export(fig1, [figuredir '/' subject '.commacc.' modality '.png'], '-png');