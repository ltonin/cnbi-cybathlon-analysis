clearvars; clc; 

% subject = 'MA25VE';
subject = 'AN14VE';

pattern     = '.mi.';
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

%% Get datafiles
[Files, NumFiles] = cnbiutil_getfile(datapath, '.mat', [subject '*' pattern]);

%% Concatenate data
cnbiutil_bdisp('[io] - Import psd all datafiles');
[U, events, labels, settings] = cnbiutil_concatenate_data(Files);
F = log(U);

NumSamples = size(F, 1);
NumFreqs   = size(F, 2);
NumChans   = size(F, 3);

Dk      = labels.Dk;
Days    = unique(Dk);
NumDays = length(Days);

Rk      = labels.Rk;
Runs    = unique(Rk);
NumRuns = length(Runs);

Mk         = labels.Mk;
Modalities = unique(Mk);
NumMods    = length(Modalities);

%% Get general events
[TrialLb, TrialEvents] = cnbiproc_get_event([771 773 783 781 786], NumSamples, events.POS, events.TYP, events.DUR);

%% Iterate across each run

TrialLb = zeros(NumSamples, 1);
for rId = 1:NumRuns
    cindex    = Rk == Runs(rId);
    cmodality = unique(Mk(cindex));
    crunstart = find(cindex, 1, 'first');
    crunstop  = find(cindex, 1, 'last');
    
    
    switch(cmodality)
        case {0, 1}          % Offline/Online
            [~, CueEvents] = cnbiproc_get_event(TaskTypeId, NumSamples, events.POS, events.TYP, events.DUR);
            [~, CFbEvents] = cnbiproc_get_event(781,  NumSamples, events.POS, events.TYP, events.DUR);
            
            evtcueId   = CueEvents.POS >= crunstart & CueEvents.POS <= crunstop;
            evtcue.POS = CueEvents.POS(evtcueId);
            evtcue.TYP = CueEvents.TYP(evtcueId);
            evtcue.DUR = CueEvents.DUR(evtcueId);
            
            evtcfbId   = CFbEvents.POS >= crunstart & CFbEvents.POS <= crunstop;
            evtcfb.POS = CFbEvents.POS(evtcfbId);
            evtcfb.TYP = CFbEvents.TYP(evtcfbId);
            evtcfb.DUR = CFbEvents.DUR(evtcfbId);
            cnumtrials = length(evtcfb.POS);
            
            for trId = 1:cnumtrials
                cstart = evtcfb.POS(trId);
                cstop  = cstart + evtcfb.DUR(trId) - 1;
                cclass = evtcue.TYP(trId);
                TrialLb(cstart:cstop) = cclass;
            end
        case {2, 3}         % Race/Competition
            clabels = cnbiproc_get_event(TaskTypeId, NumSamples, events.POS, events.TYP, events.DUR);
            TrialLb(crunstart:crunstop) = clabels(crunstart:crunstop);
    end
    
    cclasses = unique(TrialLb(crunstart:crunstop));
    
    if ismember(783, cclasses) == 0
        fixlabels = cnbiproc_get_event(786, NumSamples, events.POS, events.TYP, events.DUR);
        fixlabels(fixlabels > 0) = 783;
        TrialLb(crunstart:crunstop) = TrialLb(crunstart:crunstop) + fixlabels(crunstart:crunstop);
    end
end

[EyeLb, EyeEvents] = cnbiproc_get_event(267, NumSamples, events.POS, events.TYP, events.DUR);
ArtifactFree = EyeLb == 0;


%% Generic condition
cnbiutil_bdisp('[proc] - Set generic condition mask');
GenericCondition = ArtifactFree;

%% Generate selected pairs

SelectedTasks = [771 773 783];

combinations = nchoosek(SelectedTasks, 2);
NumCombinations = size(combinations, 1);

%% Compute discriminancy per run
cnbiutil_bdisp('[proc] - Compute discriminancy per run');

FisherScores = nan(NumFreqs*NumChans, NumRuns, NumCombinations);
rDk = zeros(NumRuns, 1);
rMk = zeros(NumRuns, 1);
for rId = 1:NumRuns
    for cId = 1:NumCombinations
        ctasks = combinations(cId, :);
        cindex = Rk == Runs(rId) & GenericCondition &(TrialLb == ctasks(1) | TrialLb == ctasks(2));
       
        if length(unique(TrialLb(cindex))) == 2
            FisherScores(:, rId, cId) = cnbiproc_fisher(F(cindex, :, :), TrialLb(cindex));
        end 
    end
    rDk(rId) = unique(Dk(Rk == Runs(rId)));
    rMk(rId) = unique(Mk(Rk == Runs(rId)));
end

%% Saving data

% Grouping results
discriminancy.run.fisherscore  = FisherScores;
discriminancy.run.label.Dk     = rDk;
discriminancy.run.label.Dl     = labels.Dl;
discriminancy.run.label.Mk     = rMk;
discriminancy.run.combinations = combinations;
discriminancy.freqs            = settings.spectrogram.freqgrid;

savefile = [savedir '/' subject '.discriminancy.maps.mat'];

cnbiutil_bdisp(['Saving discriminancy results in: ' savefile]);
save(savefile, 'discriminancy');



%% Plotting

fig1 = figure;
fig_set_position(fig1, 'All');

SelFreqs = 4:2:32;
[FreqGrid, SelFreqIds] = intersect(settings.spectrogram.freqgrid, SelFreqs);

Modlb = {'Offline', 'Online', 'Race', 'Competition'};
NumRows = NumMods;
NumCols = NumCombinations;

for mId = 1:NumMods
    for cId = 1:NumCombinations
        subplot(NumRows, NumCols, cId + NumCols*(mId - 1));
        cdata = reshape(nanmean(FisherScores(:, rMk == Modalities(mId), cId), 2), [NumFreqs NumChans]);
        
        cdata = cdata(SelFreqIds, :);
        imagesc(FreqGrid, 1:NumChans, cdata');
        
        [~, comblb] = intersect(TaskTypeId, combinations(cId, :));
        
        if cId == 1
            ylabel('Channel');
        end
        
        if mId == NumMods
            xlabel('[Hz]');
        end
        
        title([Modlb{mId} ' - ' cell2mat(TaskTypeLb(comblb))]);
    end
    
        
end

suptitle(subject);
cnbifig_export(fig1, [figuredir '/' subject '.discriminancy.maps.overall.png'], '-png');


fig2 = figure;
fig_set_position(fig2, 'All');
load('chanlocs64.mat');
BetaFreqs = 22:2:32;
[~, SelBetaFreqIds] = intersect(settings.spectrogram.freqgrid, BetaFreqs);

tModlb = {'Offline', 'Online', 'Race'};
tMk = rMk;
% tMk(tMk == 0) = 1;
tMk(tMk == 3) = 2;
tMods = unique(tMk);
tNumMods = length(tMods);
NumRows = tNumMods;
NumCols = NumCombinations;

for mId = 1:tNumMods
    for cId = 1:NumCombinations
        subplot(NumRows, NumCols, mId + NumCols*(cId -1));
        
        cdata = reshape(nanmean(FisherScores(:, tMk == tMods(mId), cId), 2), [NumFreqs NumChans]);
        tdata = convChans(mean(cdata(SelBetaFreqIds, :), 1));
        if strcmp(subject, 'AN14VE')
            maplimits = [0 0.4];
        elseif strcmp(subject, 'MA25VE')
            maplimits = [0 0.5];
        end
        topoplot(tdata, chanlocs, 'headrad', 'rim', 'maplimits', [maplimits]);
        
        if mId == tNumMods
            colorbar('location', 'EastOutside');
        end
        axis image;
        [~, comblb] = intersect(TaskTypeId, combinations(cId, :));
        title([tModlb{mId} ' - ' cell2mat(TaskTypeLb(comblb))]);
    end
end

suptitle([subject ' - Beta Band']);
cnbifig_export(fig2, [figuredir '/' subject '.discriminancy.topoplots.beta.overall.png'], '-png');