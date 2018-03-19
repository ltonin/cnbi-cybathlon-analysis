clearvars; clc; 

% subject = 'MA25VE';
subject = 'AN14VE';

pattern     = '.mi.';
experiment  = 'cybathlon';
datapath    = [pwd '/analysisforce1_strict/'];
figuredir   = './analysisforce1_strict/';
savedir     = [pwd '/analysisforce1_strict/'];

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

Pk          = labels.Pk;
Paradigms   = unique(Pk);
NumParadigms = length(Paradigms);

% Month labels
MonthStr  = labels.Dl(:, 5:6);
MonthsLb  = str2double(MonthStr);
Months    = unique(MonthsLb);
NumMonths = length(unique(Months));

Mnk = Dk;
for dId = 1:NumDays
    cday = Days(dId);
    cmonth = str2double(labels.Dl(dId, 5:6));
    
    Mnk(Dk == cday) = cmonth;
end

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

SelectedTasks = [770 771 773 783];

combinations = nchoosek(SelectedTasks, 2);
NumCombinations = size(combinations, 1);

%% Compute discriminancy per run
cnbiutil_bdisp('[proc] - Compute discriminancy per run');

FisherScores = nan(NumFreqs*NumChans, NumRuns, NumCombinations);
rDk = zeros(NumRuns, 1);
rMk = zeros(NumRuns, 1);
rPk = zeros(NumRuns, 1);
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
    rPk(rId) = unique(Pk(Rk == Runs(rId)));
end

%% Compute discriminancy per month
freqs = settings.spectrogram.freqgrid;
cnbiutil_bdisp('[proc] - Compute discriminancy per month');
FisherScoresMonth = nan(NumFreqs*NumChans, NumMonths, NumCombinations);
excludedCompetition = true(size(Mnk));
excludedCompetition(Mk == 3) = false;
ntrials = zeros(NumMonths, NumCombinations);
nruns   = zeros(NumMonths, NumCombinations);
for mId = 1:NumMonths
    for cId = 1:NumCombinations
        ctasks = combinations(cId, :);
        cindex = Mnk == Months(mId) & GenericCondition & (TrialLb == ctasks(1) | TrialLb == ctasks(2));
%         cindex = Mnk == Months(mId) & GenericCondition & excludedCompetition & (TrialLb == ctasks(1) | TrialLb == ctasks(2));
 
        if length(unique(TrialLb(cindex))) == 2
            FisherScoresMonth(:, mId, cId) = cnbiproc_fisher(F(cindex, :, :), TrialLb(cindex));
        end 
        ntrials(mId, cId) = sum(cindex);
        nruns(mId, cId) = length(unique(Rk(cindex)));
    end
end

%% Saving data

% Grouping results
discriminancy.run.fisherscore  = FisherScores;
discriminancy.run.label.Dk     = rDk;
discriminancy.run.label.Dl     = labels.Dl;
discriminancy.run.label.Mk     = rMk;
discriminancy.run.label.Pk     = rPk;
discriminancy.run.label.Pl     = labels.Pl;
discriminancy.run.combinations = combinations;
discriminancy.month.fisherscore = FisherScoresMonth;
discriminancy.month.label.Mnk  = unique(Mnk);
discriminancy.freqs            = settings.spectrogram.freqgrid;

savefile = [savedir '/' subject '.discriminancy.maps.mat'];

cnbiutil_bdisp(['Saving discriminancy results in: ' savefile]);
save(savefile, 'discriminancy');



%% Plotting




fig1 = figure;
%fig_set_position(fig1, 'All');

SelFreqs = 4:2:96;

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
%fig_set_position(fig2, 'All');
load('chanlocs64.mat');
BetaFreqs = 22:2:96;
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
        subplot(NumRows, NumCols, cId + NumCols*(mId -1));
        
        cdata = reshape(nanmean(FisherScores(:, tMk == tMods(mId), cId), 2), [NumFreqs NumChans]);
        tdata = convChans(mean(cdata(SelBetaFreqIds, :), 1));
        if strcmp(subject, 'AN14VE')
            maplimits = [0 0.4];
        elseif strcmp(subject, 'MA25VE')
            maplimits = [0 0.5];
        end
        topoplot(tdata, chanlocs, 'headrad', 'rim', 'maplimits', [maplimits]);
        
        

        
        axis image;
        [~, comblb] = intersect(TaskTypeId, combinations(cId, :));
        title([cell2mat(TaskTypeLb(comblb))]);
        if cId == 1
            h = axes('Position', get(gca, 'Position'), 'Visible', 'off');
            set(h.YLabel, 'Visible', 'on');
            ylabel(tModlb{mId});
        end
        if cId == NumCombinations
            colorbar('location', 'EastOutside');
        end
    end
end

suptitle([subject ' - Beta Band']);
cnbifig_export(fig2, [figuredir '/' subject '.discriminancy.topoplots.beta.overall.png'], '-png');


%% Figure dp per month (dp computed per month)
fig3 = figure;
fig_set_position(fig3, 'Top');
SelectedClassId = [771 773];
SelectedCombination = find(ismember(combinations, SelectedClassId, 'rows'));
MonthsLabels = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
for m = 1:NumMonths
    subplot(1, 4, m); 
    maplimits = [0 0.6];
    imagesc(freqs, 1:16, reshape(FisherScoresMonth(:, m, SelectedCombination), [NumFreqs NumChans])', [0 0.6]); 
    title([MonthsLabels{Months(m)} ' - nruns=' num2str(nruns(m, SelectedCombination))]); 
end

suptitle([subject ' - dp per month']);

%% Figure dp per run per month
SelectedClassId = [771 773];
SelectedCombination = find(ismember(combinations, SelectedClassId, 'rows'));

SelectedMonth = [8 9 10];

for mId = 1:length(SelectedMonth)
    figure;
    fig_set_position(gcf, 'All');
   
    runmonth = unique(Rk(Mnk == SelectedMonth(mId)));
    nrunmonth = length(runmonth);
    maplimits = [0 1.5];
    for rId = 1:nrunmonth
        subplot(6, ceil(nrunmonth/6), rId);
      
        imagesc(freqs, 1:16, reshape(FisherScores(:, Runs == runmonth(rId), SelectedCombination), [NumFreqs NumChans])', maplimits); 
        title(['Run= ' num2str(runmonth(rId))]);
    end
    
    suptitle([subject ' - Month ' MonthsLabels{SelectedMonth(mId)}]);
end


%% Figure dp per month (average)
fig4 = figure;
fig_set_position(fig4, 'Top');

ExcludedRunIds = [43 45 51 53 59 61 62 63 66 73];


for mId = 1:NumMonths
    subplot(1, NumMonths, mId);
    cmonthId = Months(mId);
    runmonth = unique(Rk(Mnk == cmonthId));

    runindex = setdiff(runmonth, ExcludedRunIds);
    cdata = nanmean(FisherScores(:, runindex, SelectedCombination), 2);
    
    maplimits = [0 0.5];
    imagesc(freqs, 1:16, reshape(cdata, [NumFreqs NumChans])', maplimits); 
    title([MonthsLabels{Months(m)} ' - nruns=' num2str(length(runindex))]);
    
end


