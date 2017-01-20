clearvars; clc; 

subject = 'AN14VE';

pattern     = '.mi.';
modality    = 'race';

experiment  = 'cybathlon';
datapath    = [pwd '/analysis/'];
figuredir  = './figures/';

PadTypeId = [768 769 770 771 773 783];
PadTypeLb = {'Slide', 'Slide', 'Speed', 'Jump', 'Speed', 'Rest'};

SelectedClassId = [771 773];
SelectedClassLb = {'BothFeet', 'BothHands'};
NumClasses = length(SelectedClassId);

SelFreqs = 4:2:32;


%% Get racetime
racedata = load([datapath '/' subject '.race.time.mat']); 

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
[TrialLb, TrialEvents] = cnbiproc_get_event(PadTypeId, DataLength, events.POS, events.TYP, events.DUR);
[EyeLb, EyeEvents] = cnbiproc_get_event(267, DataLength, events.POS, events.TYP, events.DUR);
ArtifactFree = EyeLb == 0;

%% Selecting classes
Ck = zeros(DataLength, 1);
for cId = 1:NumClasses
    Ck(TrialLb == SelectedClassId(cId)) = SelectedClassId(cId);
end

%% Generic condition
GenericCondition = Ck > 0 & ArtifactFree;

%% Compute discriminancy per day
NSigma = [];
discrday = [];
discrdlb = [];
for dId = 1:NumDays
    cindex = labels.Dk == Days(dId) & GenericCondition;
    
    if length(unique(Ck(cindex))) == 2
        discrday = cat(2, discrday, cnbiproc_fisher(F(cindex, :, :), Ck(cindex), NSigma));
        discrdlb = cat(1, discrdlb, labels.Dl(dId, :));
    end
end

%% Extract race time according to the daylabels of DP maps
NumDPDays = size(discrday, 2);
dtime = zeros(2, NumDPDays);
for dId = 1:NumDPDays
    cdaylabel = discrdlb(dId, :);
    cdayid = find(ismember(racedata.race.Dl, cdaylabel, 'rows'));
    cindex = racedata.race.Dk == cdayid;
    dtime(:, dId) = [mean(racedata.race.time(cindex)) std(racedata.race.time(cindex))];
    
end

%% Get correlations

AlphaBand = 8:12;
BetaBand  = 14:30;

[~, AlphaBandId] = intersect(FreqGrid, AlphaBand);
[~, BetaBandId]  = intersect(FreqGrid, BetaBand);

rdiscr = reshape(discrday, [NumFreqs NumChans NumDPDays]);

[calpha, palpha] = corr(squeeze(sum(sum(rdiscr(AlphaBandId, 9, :), 1), 2)), dtime(1, :)');
[cbeta, pbeta] = corr(squeeze(sum(sum(rdiscr(BetaBandId, 9, :), 1), 2)), dtime(1, :)');




