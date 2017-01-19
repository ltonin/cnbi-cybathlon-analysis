clearvars; clc;

%subject = 'AN14VE';
subject = 'MA25VE';

pattern     = '.mi.';
modality    = 'online';

experiment  = 'cybathlon';
datapath    = [pwd '/analysis/'];
figuredir  = './figures/';

CueTypeId = [769 770 771 773 774 775 783];
CueTypeLb = {'LeftHand', 'RightHand', 'BothFeet', 'BothHands', 'Boh1', 'Boh2', 'Rest'};

CFbTypeId = 781;
CFbTypeLb = {'Continous Feedback'};

SelectedClassId = [773 771];
SelectedClassLb = {'BothHands', 'BothFeet'};
NumClasses = length(SelectedClassId);

SelFreqs = 4:2:32;

KFolds = 10;
NFeats = 5;

%% Get datafiles
[Files, NumFiles] = cnbiutil_getfile(datapath, '.mat', [subject '*' modality '*' pattern]);

%% Concatenate data
cnbiutil_bdisp(['[io] - Import psd datafiles (' modality ')']);
[U, events, labels, settings] = cnbiutil_concatenate_data(Files);
[FreqGrid, SelFreqIds] = intersect(settings.spectrogram.freqgrid, SelFreqs);

F = log(U(:, SelFreqIds, :));
DataLength  = size(F, 1);
NumFreqs = size(F, 2);
NumChans = size(F, 3);


%% Extract events
cnbiutil_bdisp('[proc] - Extract events');
[~, CueEvents] = cnbiproc_get_event(CueTypeId, DataLength, events.POS, events.TYP, events.DUR);
[~, CFbEvents] = cnbiproc_get_event(CFbTypeId, DataLength, events.POS, events.TYP, events.DUR);

NumTrials = length(CFbEvents.POS);

%% Labeling data

TrialCls = zeros(DataLength, 1);
TrialInd = zeros(DataLength, 1);
TrialRun = zeros(DataLength, 1);
TrialDay = zeros(DataLength, 1);
for trId = 1:NumTrials
   cstart = CFbEvents.POS(trId);
   cstop  = cstart + CFbEvents.DUR(trId) - 1;
   cclass = CueEvents.TYP(trId);
   
   TrialCls(cstart:cstop) = cclass;
   TrialInd(cstart:cstop) =trId;
end
Days    = unique(labels.Dk);
NumDays = length(Days);

%% Selecting classes
Ck = zeros(length(TrialCls), 1);
for cId = 1:NumClasses
    Ck(TrialCls == SelectedClassId(cId)) = SelectedClassId(cId);
end

%% Compute the overall simulated singl
Lbl = Ck(Ck > 0);
Lbl(Lbl==SelectedClassId(1))=1;
Lbl(Lbl==SelectedClassId(2))=2;
SimAcc =  cvk(F(Ck > 0, :, :), Lbl, TrialInd(Ck > 0), KFolds, NFeats);

%% Compute discriminancy per day
SimAccSes = [];
SimAccSesLb = {};
for dId = 1:NumDays
    cindex = labels.Dk == Days(dId) & Ck > 0;
    if length(unique(Ck(cindex))) == 2
        Lbl = Ck(cindex);
        Lbl(Lbl==SelectedClassId(1))=1;
        Lbl(Lbl==SelectedClassId(2))=2;
        SimAccSes = cat(1, SimAccSes, cvk(F(cindex, :, :), Lbl, TrialInd(cindex), KFolds, NFeats));
        SimAccSesLb = cat(1, SimAccSesLb, labels.Dl(dId, :));
    end
end

%% Plotting simulated single-sample accuracy
fig1 = figure;
cnbifig_set_position(fig1, 'All');
plot(1:length(SimAccSes),SimAccSes,'LineWidth',3);
xlabel('Feedback Session','FontSize',20,'LineWidth',3);
ylabel('Simulated single-sample Accuracy (%)','FontSize',20,'LineWidth',3);
title(subject);
axis([0 length(SimAccSes)+1 45 100]);
set(gca,'FontSize',20,'LineWidth',3);
set(gca,'XTick',[1:length(SimAccSes)]);
set(gca,'XTickLabel',SimAccSesLb);
xticklabel_rotate([],45,[])

cnbifig_export(fig1, [figuredir '/' subject '.simacc.' modality '.png'], '-png');