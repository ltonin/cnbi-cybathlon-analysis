clearvars; clc; 

% subject = 'MA25VE';
subject = 'AN14VE';

pattern     = '.mi.';
modality    = 'race';

experiment  = 'cybathlon';
datapath    = [pwd '/analysis/'];
figuredir   = './figures/';
savedir     = [pwd '/analysis/'];

TaskTypeId = [769 770 771 773 774 775 783];
TaskTypeLb = {'LeftHand', 'RightHand', 'BothFeet', 'BothHands', 'Boh1', 'Boh2', 'Rest'};

CFbTypeId = 781;
CFbTypeLb = {'Continous Feedback'};

FixTypeId = 786;
FixTypeLb = {'Fixation'};

SelectedClassId = [773 771];
SelectedClassLb = {'BothFeet', 'BothHands'};
NumClasses = length(SelectedClassId);

CommandId = SelectedClassId + hex2dec('6000');

SelFreqs = 4:2:32;

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

Dk = labels.Dk;
Days    = unique(Dk);
NumDays = length(Days);

Mk = labels.Mk;

%% Extract events
cnbiutil_bdisp(['[proc] - Extract events (' modality ')']);
if strcmpi(modality, 'online') || strcmpi(modality, 'offline') 
    disp('         Based on cue and continous feedback');
    [~, CueEvents] = cnbiproc_get_event(TaskTypeId, DataLength, events.POS, events.TYP, events.DUR);
    [~, CFbEvents] = cnbiproc_get_event(CFbTypeId,  DataLength, events.POS, events.TYP, events.DUR);
    NumTrials = length(CFbEvents.POS);
    ArtifactFree = true(DataLength, 1);
elseif strcmpi(modality, 'race')
    disp('         Based on pad and eye artifacts');
    [TrialLb, TrialEvents] = cnbiproc_get_event(TaskTypeId, DataLength, events.POS, events.TYP, events.DUR);
    [EyeLb, EyeEvents]     = cnbiproc_get_event(267, DataLength, events.POS, events.TYP, events.DUR);
    ArtifactFree = EyeLb == 0;
end

%% Labeling data
cnbiutil_bdisp(['[proc] - Labeling data (' modality ')']);
if strcmpi(modality, 'online') || strcmpi(modality, 'offline') 
    TrialLb = zeros(DataLength, 1);
    TrialRun = zeros(DataLength, 1);
    TrialDay = zeros(DataLength, 1);
    for trId = 1:NumTrials
       cstart = CFbEvents.POS(trId);
       cstop  = cstart + CFbEvents.DUR(trId) - 1;
       cclass = CueEvents.TYP(trId);
       TrialLb(cstart:cstop) = cclass;
    end
elseif strcmpi(modality, 'race')
    for trId = 1:length(TrialEvents.TYP)
       cstart = TrialEvents.POS(trId);
       cstop  = cstart + TrialEvents.DUR(trId) + 1;
       ctype  = TrialEvents.TYP(trId);
    end
end

%% Selecting classes
cnbiutil_bdisp(['[proc] - Selecting required classes: ' num2str(SelectedClassId)]);
Ck = zeros(length(TrialLb), 1);
for cId = 1:NumClasses
    Ck(TrialLb == SelectedClassId(cId)) = SelectedClassId(cId);
end

%% Generic condition
cnbiutil_bdisp('[proc] - Set generic condition mask');
GenericCondition = Ck > 0 & ArtifactFree;

%% Compute discriminancy per day
cnbiutil_bdisp('[proc] - Compute discriminancy per day');
NSigma   = [];
discrlb  = cell(NumDays, 1);
discrday = zeros(NumFreqs*NumChans, NumDays);
for dId = 1:NumDays
    cindex = labels.Dk == Days(dId) & Ck > 0;
    if length(unique(Ck(cindex))) == 2
        discrday(:, dId) = cnbiproc_fisher(F(cindex, :, :), Ck(cindex), NSigma);
    end
    discrlb{dId}     = char(labels.Dl(dId, :));
end

%% Plotting discriminancy map
cnbiutil_bdisp('[proc] - Plotting');
load('chanlocs64.mat');
fig1 = figure;
cnbifig_set_position(fig1, 'All');

AlphaBand = 8:12;
BetaBand  = 14:32;

[~, AlphaBandId] = intersect(FreqGrid, AlphaBand);
[~, BetaBandId]  = intersect(FreqGrid, BetaBand);

NumRows = 3;
NumCols = size(discrday, 2);

% Discriminancy maps
for dId = 1:size(discrday, 2)
    subplot(NumRows, NumCols, dId);
    cdata = reshape(discrday(:, dId), [NumFreqs NumChans]);
    imagesc(FreqGrid, 1:NumChans, cdata');
    
    if dId == 1
        ylabel('Channel');
    end
    xlabel('[Hz]');
    title(discrlb{dId});
end

% Topoplots
for dId = 1:size(discrday, 2)
    subplot(NumRows, NumCols, dId + NumCols);
    cdata = reshape(discrday(:, dId), [NumFreqs NumChans]);
    
    tdata = convChans(mean(cdata(AlphaBandId, :), 1));
    topoplot(tdata, chanlocs, 'headrad', 'rim');
    axis image;
    
    if dId == 1
         h = axes('Position', get(gca, 'Position'), 'Visible', 'off');
        set(h.YLabel, 'Visible', 'on');
        ylabel('Alpha band');
    end
    
    if dId == size(discrday, 2)
        colorbar(gca, 'EastOutside')
    end
end

% Topoplots
for dId = 1:size(discrday, 2)
    subplot(NumRows, NumCols, dId + 2*NumCols);
    cdata = reshape(discrday(:, dId), [NumFreqs NumChans]);
    
    tdata = convChans(mean(cdata(BetaBandId, :), 1));
    topoplot(tdata, chanlocs, 'headrad', 'rim');
    axis image;
    title('');
    if dId == 1
         h = axes('Position', get(gca, 'Position'), 'Visible', 'off');
        set(h.YLabel, 'Visible', 'on');
        ylabel('Beta band');
    end
    if dId == size(discrday, 2)
        colorbar(gca, 'EastOutside')
    end
end

suptitle([subject ' - DP - ' modality ' - Classes: ' num2str(SelectedClassId)]);

cnbifig_export(fig1, [figuredir '/' subject '.discriminancy.' modality '.png'], '-png');

%% Saving metadata

% Grouping results
discriminancy.fisherscore = discrday;
discriminancy.label       = discrlb;

savefile = [savedir '/' subject '.discriminancy.' modality '.mat'];

cnbiutil_bdisp(['Saving discriminancy (' modality ') results in: ' savefile]);
save(savefile, 'discriminancy');
