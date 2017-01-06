clearvars; clc; 

subject = 'AN14VE';

pattern     = '.mi.';
modality    = 'race';

experiment  = 'cybathlon';
datapath    = [pwd '/analysis/'];
figuredir  = '/figures/';

PadTypeId = [768 771 773 783];
PadTypeLb = {'Slide', 'Jump', 'Speed', 'Rest'};

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

%% Compute the overall discriminancy
NSigma = 5;
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


%% Plotting discriminancy map
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
    imagesc(FreqGrid, 1:NumChans, cdata');
    
    if dId == 1
        ylabel('Channel');
    end
    xlabel('Frequency [Hz]');
    title(discrdlb(dId, :));

end

% Topoplots
for dId = 1:size(discrday, 2)
    subplot(NumRows, NumCols, dId + NumCols);
    cdata = reshape(discrday(:, dId), [NumFreqs NumChans]);
    
    tdata = convChans(mean(cdata(AlphaBandId, :), 1));
    topoplot(tdata, chanlocs, 'headrad', 'rim', 'maplimits', [0 0.3]);
    axis image;
    
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
    if dId == size(discrday, 2)
        colorbar(gca, 'EastOutside')
    end
end

suptitle([subject ' - DP - ' modality]);
