% clearvars; clc; 
% 
subject = 'MA25VE';
% %subject = 'AN14VE';

pattern     = '.mi.';
modality    = 'race';

experiment  = 'cybathlon';
datapath    = [pwd '/analysis/'];
figuredir   = './figures/';
savedir     = [pwd '/analysis/'];
AnalysisType = 'all'; % 'correct'       % Select if all trials are taken into account or only the correct ones (no big changes)


TrialTypeId = [768 769 770 771 773 783];
TrialTypeLb = {'Slide', 'Slide', 'Speed', 'Jump', 'Speed', 'Rest'};

SelectedClassId = [773 771];
SelectedClassLb = {'BothFeet', 'BothHands'};
NumClasses = length(SelectedClassId);

CommandId = SelectedClassId + hex2dec('6000');

SelFreqs = 4:2:32;

KFolds = 3;
NFeats = 5;

%% Get datafiles
[Files, NumFiles] = cnbiutil_getfile(datapath, '.mat', [subject '*' modality '*' pattern]);

%% Concatenate data
cnbiutil_bdisp(['[io] - Import psd datafiles (' modality ')']);
[U, events, labels, settings] = cnbiutil_concatenate_data(Files);

[FreqGrid, SelFreqIds] = intersect(settings.spectrogram.freqgrid, SelFreqs);
F = log(U(:, SelFreqIds, :));
% F = U(:, SelFreqIds, :);

DataLength  = size(F, 1);
NumFreqs = size(F, 2);
NumChans = size(F, 3);

Dk = labels.Dk;
Days    = unique(Dk);
NumDays = length(Days);

Mk = labels.Mk;

%% Extract events
cnbiutil_bdisp('[proc] - Extract events');
[TrialLb, TrialEvents] = cnbiproc_get_event(TrialTypeId, DataLength, events.POS, events.TYP, events.DUR);
[~, CmdEvents]         = cnbiproc_get_event(CommandId, DataLength, events.POS, events.TYP, events.DUR);
[EyeLb, EyeEvents]     = cnbiproc_get_event(267, DataLength, events.POS, events.TYP, events.DUR);
ArtifactFree = EyeLb == 0;

%% Create correct command label (only for jump and speed pad)
CorrectLb = false(DataLength, 1);
TrialInd = zeros(DataLength, 1);
for trId = 1:length(TrialEvents.TYP)
   cstart = TrialEvents.POS(trId);
   cstop  = cstart + TrialEvents.DUR(trId) + 1;
   ctype  = TrialEvents.TYP(trId);
   
   index = ctype == (CmdEvents.TYP - hex2dec('6000')) & CmdEvents.POS >= cstart & CmdEvents.POS <= cstop;
   
   if sum(index) > 0
       CorrectLb(cstart:cstop) = true;
   end
   
   TrialInd(cstart:cstop)=trId;
end

%% Selecting classes
Ck = zeros(DataLength, 1);
for cId = 1:NumClasses
    Ck(TrialLb == SelectedClassId(cId)) = SelectedClassId(cId);
end

%% Generic condition
GenericCondition = Ck > 0 & ArtifactFree;

if strcmpi(AnalysisType, 'correct')
    GenericCondition = GenericCondition & CorrectLb;
end

%% Compute the overall simulated accuracy
Lbl = Ck(GenericCondition);
Lbl(Lbl==SelectedClassId(1))=1;
Lbl(Lbl==SelectedClassId(2))=2;
[SimAcc, SimCM] =  cvk(F(GenericCondition, :, :), Lbl, TrialInd(GenericCondition), KFolds, NFeats);
SimAccClass = [SimCM(1,1) SimCM(2,2)];

%% Compute simulated accuracy per day
SimAccSes = [];
SimAccClassSes = [];
SimCMSes = [];
SimAccSesLb = {};
for dId = 1:NumDays
    cindex = labels.Dk == Days(dId) & Ck > 0;
    if length(unique(Ck(cindex))) == 2
        Lbl = Ck(cindex);
        Lbl(Lbl==SelectedClassId(1))=1;
        Lbl(Lbl==SelectedClassId(2))=2;
        [tmpSimAcc, tmpCM] = cvk(F(cindex, :, :), Lbl, TrialInd(cindex), KFolds, NFeats);
        SimAccSes = cat(1, SimAccSes, tmpSimAcc);
        SimCMSes = cat(3, SimCMSes, tmpCM);
        SimAccSesLb = cat(1, SimAccSesLb, labels.Dl(dId, :));
        SimAccClassSes = cat(3,SimAccClassSes,[tmpCM(1,1) tmpCM(2,2)]);
    end
end
SimAccClassSes = squeeze(SimAccClassSes)';
%% Compute simulated accuracy per race for competition day

CybRaceId = unique(labels.Rk(labels.Mk == 3));
NumCybRaces = length(CybRaceId);
SimAccCyb = zeros(2,1);
SimCMCyb = zeros(2,2,2);
for rId = 1:NumCybRaces
    craceid = CybRaceId(rId);
    cindex = labels.Rk == craceid & GenericCondition;
    Lbl = Ck(cindex);
    Lbl(Lbl==SelectedClassId(1))=1;
    Lbl(Lbl==SelectedClassId(2))=2;
    [SimAccCyb(rId) SimCMCyb(rId,:,:)] = cvk(F(cindex, :, :), Lbl, TrialInd(cindex), KFolds, NFeats);
end

%% Plotting simulated single-sample accuracy
fig1 = figure;
cnbifig_set_position(fig1, 'All');
plot(1:length(SimAccSes),SimAccSes,'LineWidth',3);
xlabel('Race Session','FontSize',20,'LineWidth',3);
ylabel('Simulated single-sample Accuracy (%)','FontSize',20,'LineWidth',3);
title(subject);
axis([0 length(SimAccSes)+1 45 100]);
set(gca,'FontSize',20,'LineWidth',3);
set(gca,'XTick',[1:length(SimAccSes)]);
set(gca,'XTickLabel',SimAccSesLb);
xticklabel_rotate([],45,[])

cnbifig_export(fig1, [figuredir '/' subject '.simaccrace.' modality '.png'], '-png');

%% Saving metadata

% Grouping results
simulated.accuracy = SimAccSes;
simulated.label    = cell2mat(SimAccSesLb);

savefile = [savedir '/' subject '.metadata.mat'];
if exist(savefile, 'file')
    cnbiutil_bdisp(['Loading metadata from: ' savefile]);
    load(savefile);
end

metadata.race.simulated = simulated;

cnbiutil_bdisp(['Saving simulated accuracy (race) results in: ' savefile]);
save(savefile, 'metadata');
