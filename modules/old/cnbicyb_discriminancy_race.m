% clearvars; clc; 
% 
% subject = 'MA25VE';
% %subject = 'AN14VE';

pattern     = '.mi.';
modality    = 'race';

experiment  = 'cybathlon';
datapath    = [pwd '/analysis/'];
figuredir   = './figures/';
savedir     = [pwd '/analysis/'];

TrialTypeId = [768 769 770 771 773 783];
TrialTypeLb = {'Slide', 'Slide', 'Speed', 'Jump', 'Speed', 'Rest'};

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
cnbiutil_bdisp('[proc] - Extract events');
[TrialLb, TrialEvents] = cnbiproc_get_event(TrialTypeId, DataLength, events.POS, events.TYP, events.DUR);
[~, CmdEvents]         = cnbiproc_get_event(CommandId, DataLength, events.POS, events.TYP, events.DUR);
[EyeLb, EyeEvents]     = cnbiproc_get_event(267, DataLength, events.POS, events.TYP, events.DUR);
ArtifactFree = EyeLb == 0;

%% Create correct command label (only for jump and speed pad)
CorrectLb = false(DataLength, 1);
for trId = 1:length(TrialEvents.TYP)
   cstart = TrialEvents.POS(trId);
   cstop  = cstart + TrialEvents.DUR(trId) + 1;
   ctype  = TrialEvents.TYP(trId);
   
   index = ctype == (CmdEvents.TYP - hex2dec('6000')) & CmdEvents.POS >= cstart & CmdEvents.POS <= cstop;
   
   if sum(index) > 0
       CorrectLb(cstart:cstop) = true;
   end
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

%% Compute the overall discriminancy
NSigma = [];
discrovl = cnbiproc_fisher(F(GenericCondition, :, :), Ck(GenericCondition), NSigma);

%% Compute discriminancy per day

discrday = [];
discrdlb = [];
for dId = 1:NumDays
    cindex = labels.Dk == Days(dId) & GenericCondition;
    
    if length(unique(Ck(cindex))) == 2
        discrday = cat(2, discrday, cnbiproc_fisher(F(cindex, :, :), Ck(cindex), NSigma));
        discrdlb = cat(1, discrdlb, labels.Dl(dId, :));
    end
end

%% Compute discriminancy per race for competition day

CybRaceId = unique(labels.Rk(labels.Mk == 3));
NumCybRaces = length(CybRaceId);
discrcyb = zeros(NumFreqs*NumChans, NumCybRaces);
for rId = 1:NumCybRaces
    craceid = CybRaceId(rId);
    cindex = labels.Rk == craceid & GenericCondition;
    discrcyb(:, rId) = cnbiproc_fisher(F(cindex, :, :), Ck(cindex), NSigma);
end

%% Plotting discriminancy map per races
load('chanlocs64.mat');
fig1 = figure;
cnbifig_set_position(fig1, 'All');

AlphaBand = 8:12;
BetaBand  = 14:30;

[~, AlphaBandId] = intersect(FreqGrid, AlphaBand);
[~, BetaBandId]  = intersect(FreqGrid, BetaBand);

NumRows = 3;
NumCols = size(discrday, 2);

% Discriminancy maps
for dId = 1:size(discrday, 2)
    subplot(NumRows, NumCols, dId);
    cdata = reshape(discrday(:, dId), [NumFreqs NumChans]);
    imagesc(FreqGrid, 1:NumChans, cdata',[0 0.3]);
    
    if dId == 1
        ylabel('Channel');
    end
    xlabel('[Hz]');
    title(discrdlb(dId, :));

end

% Topoplots
for dId = 1:size(discrday, 2)
    subplot(NumRows, NumCols, dId + NumCols);
    cdata = reshape(discrday(:, dId), [NumFreqs NumChans]);
    
    tdata = convChans(mean(cdata(AlphaBandId, :), 1));
    topoplot(tdata, chanlocs, 'headrad', 'rim', 'maplimits', [0 0.3]);
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
    topoplot(tdata, chanlocs, 'headrad', 'rim', 'maplimits', [0 0.3]);
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

suptitle([subject ' - DP - ' modality ' (' AnalysisType ')']);
cnbifig_export(fig1, [figuredir '/' subject '.discriminancy.' modality '.png'], '-png');


%% Plotting DP maps per race for competition day
fig2 = figure;
cnbifig_set_position(fig2, 'All');
NumRows = 3;
NumCols = NumCybRaces;
discrcyblb = {'Qualifier', 'Final'};

% Discriminancy maps
for rId = 1:NumCybRaces
    subplot(NumRows, NumCols, rId);
    cdata = reshape(discrcyb(:, rId), [NumFreqs NumChans]);
    imagesc(FreqGrid, 1:NumChans, cdata',[0 1]);
    
    if rId == 1
        ylabel('Channel');
    end
    xlabel('[Hz]');
    title(discrcyblb{rId});
end

% Topoplots
for rId = 1:NumCybRaces
    subplot(NumRows, NumCols, rId + NumCols);
    cdata = reshape(discrcyb(:, rId), [NumFreqs NumChans]);
    
    tdata = convChans(mean(cdata(AlphaBandId, :), 1));
    topoplot(tdata, chanlocs, 'headrad', 'rim', 'maplimits', [0 0.6]);
    axis image;
    if dId == 1
        h = axes('Position', get(gca, 'Position'), 'Visible', 'off');
        set(h.YLabel, 'Visible', 'on');
        ylabel('Alpha band');
    end
    if rId == NumCybRaces
        colorbar(gca, 'EastOutside')
    end
end

% Topoplots
for rId = 1:NumCybRaces
    subplot(NumRows, NumCols, rId + 2*NumCols);
    cdata = reshape(discrcyb(:, rId), [NumFreqs NumChans]);
    
    tdata = convChans(mean(cdata(BetaBandId, :), 1));
    topoplot(tdata, chanlocs, 'headrad', 'rim', 'maplimits', [0 0.6]);
    axis image;
    title('');
    if dId == 1
        h = axes('Position', get(gca, 'Position'), 'Visible', 'off');
        set(h.YLabel, 'Visible', 'on');
        ylabel('Beta band');
    end
    if rId == NumCybRaces
        colorbar(gca, 'EastOutside')
    end
end

suptitle([subject ' - DP - Cybathlon (' AnalysisType ')']);

%% Saving metadata

% Grouping results
discriminancy.fisherscore = discrday;
discriminancy.label       = discrdlb;

savefile = [savedir '/' subject '.metadata.mat'];
if exist(savefile, 'file')
    cnbiutil_bdisp(['Loading metadata from: ' savefile]);
    load(savefile);
end

metadata.race.discriminancy = discriminancy;

cnbiutil_bdisp(['Saving discriminancy (race) results in: ' savefile]);
save(savefile, 'metadata');