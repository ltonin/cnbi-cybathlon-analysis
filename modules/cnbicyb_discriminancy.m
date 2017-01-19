% clearvars; clc;
% 
% subject = 'MA25VE';

pattern     = '.mi.';
modality    = 'online';

experiment  = 'cybathlon';
datapath    = [pwd '/analysis/'];
figuredir   = './figures/';
savedir     = [pwd '/analysis/'];

CueTypeId = [769 770 771 773 774 775 783];
CueTypeLb = {'LeftHand', 'RightHand', 'BothFeet', 'BothHands', 'Boh1', 'Boh2', 'Rest'};

CFbTypeId = 781;
CFbTypeLb = {'Continous Feedback'};

FixTypeId = 786;
FixTypeLb = {'Fixation'};

SelectedClassId = [771 773];
SelectedClassLb = {'BothFeet', 'BothHands'};
NumClasses = length(SelectedClassId);

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


%% Extract events
cnbiutil_bdisp('[proc] - Extract events');
[~, CueEvents] = cnbiproc_get_event(CueTypeId, DataLength, events.POS, events.TYP, events.DUR);
[~, CFbEvents] = cnbiproc_get_event(CFbTypeId, DataLength, events.POS, events.TYP, events.DUR);

NumTrials = length(CFbEvents.POS);

%% Labeling data

TrialCls = zeros(DataLength, 1);
TrialRun = zeros(DataLength, 1);
TrialDay = zeros(DataLength, 1);
for trId = 1:NumTrials
   cstart = CFbEvents.POS(trId);
   cstop  = cstart + CFbEvents.DUR(trId) - 1;
   cclass = CueEvents.TYP(trId);
   
   TrialCls(cstart:cstop) = cclass;
end
Days    = unique(labels.Dk);
NumDays = length(Days);

%% Selecting classes
Ck = zeros(length(TrialCls), 1);
for cId = 1:NumClasses
    Ck(TrialCls == SelectedClassId(cId)) = SelectedClassId(cId);
end


%% Compute the overall discriminancy
NSigma = [];
discrovl = cnbiproc_fisher(F(Ck > 0, :, :), Ck(Ck > 0), NSigma);

%% Compute discriminancy per day

discrday = [];
discrdlb = [];
for dId = 1:NumDays
    cindex = labels.Dk == Days(dId) & Ck > 0;
    if length(unique(Ck(cindex))) == 2
        discrday = cat(2, discrday, cnbiproc_fisher(F(cindex, :, :), Ck(cindex), NSigma));
        discrdlb = cat(1, discrdlb, labels.Dl(dId, :));
    end
end


%% Plotting discriminancy map
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
    imagesc(FreqGrid, 1:NumChans, cdata', [0 0.5]);
    
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
    topoplot(tdata, chanlocs, 'headrad', 'rim', 'maplimits', [0 0.2]);
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
    topoplot(tdata, chanlocs, 'headrad', 'rim', 'maplimits', [0 0.2]);
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

suptitle([subject ' - DP - ' modality]);

cnbifig_export(fig1, [figuredir '/' subject '.discriminancy.' modality '.png'], '-png');

%% Saving metadata

% Grouping results
discriminancy.fisherscore = discrday;
discriminancy.label       = discrdlb;

savefile = [savedir '/' subject '.metadata.mat'];
if exist(savefile, 'file')
    cnbiutil_bdisp(['Loading metadata from: ' savefile]);
    load(savefile);
end

metadata.online.discriminancy = discriminancy;

cnbiutil_bdisp(['Saving discriminancy (online) results in: ' savefile]);
save(savefile, 'metadata');
