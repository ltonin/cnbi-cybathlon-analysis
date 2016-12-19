clearvars; clc; close all;

subject = 'AN14VE';

pattern     = '.mi.';
modality    = 'offline';

experiment  = 'cybathlon';
datapath    = [pwd '/analysis/'];
figuredir  = '/figures/';


CueTypeId = [769 770 771 773 774 775 783];
CueTypeLb = {'LeftHand', 'RightHand', 'BothFeet', 'BothHands', 'Boh1', 'Boh2', 'Rest'};

CFbTypeId = 781;
CFbTypeLb = {'Continous Feedback'};

FixTypeId = 786;
FixTypeLb = {'Fixation'};

TrialPeriod = [-2 4];       % With respect to continous feedback (781)

%% Get datafiles
[Files, NumFiles] = cnbiutil_getfile(datapath, '.mat', [subject '*' modality '*' pattern]);

%% Concatenate data
cnbiutil_bdisp('[io] - Import psd datafiles');
[F, events, labels, settings] = cnbiutil_concatenate_data(Files);
DataLength = size(F, 1);

%% Extract events

[~, CueEvents] = cnbiproc_get_event(CueTypeId, DataLength, events.POS, events.TYP, events.DUR);
[~, CFbEvents] = cnbiproc_get_event(CFbTypeId, DataLength, events.POS, events.TYP, events.DUR);

% [DataLabels, DataEvents] = cnbiproc_get_event(CueTypeId, DataLength, CueEvents.POS, CueEvents.TYP, CFbEvents.DUR + CueEvents.DUR - 1, CueEvents.DUR);
[FixLabels, FixEvents]   = cnbiproc_get_event(FixTypeId, DataLength, events.POS, events.TYP, events.DUR);

%% Compute ERD/ERS for each trial

TrialSize = TrialPeriod/settings.spectrogram.wshift;            % <- (TrialPeriod*settings.data.samplerate)/(settings.spectrogram.wshift*settings.data.samplerate)
FreqGrid  = settings.spectrogram.freqgrid;


U = log(F);

NumSamples = length(TrialSize(1):TrialSize(2)) - 1;
NumFreqs   = size(F, 2);
NumChans   = size(F, 3);
NumTrials  = length(CueEvents.POS);

ERSP = zeros(NumSamples, NumFreqs, NumChans, NumTrials);
Ck   = zeros(NumTrials, 1);
Rk   = zeros(NumTrials, 1);
Mk   = zeros(NumTrials, 1);
Dk   = zeros(NumTrials, 1);

cnbiutil_bdisp('[proc] - Compute ERSP for each trial');

for trId = 1:NumTrials
    cnbiutil_disp_progress(trId, NumTrials, '    ');
    
    % Reference period [Fixation]
    Rstart = FixEvents.POS(trId);
    Rstop  = Rstart + FixEvents.DUR(trId) - 1;

    
    % Trial period [with respect to continous feedback event]
    Tstart = CFbEvents.POS(trId) + TrialSize(1);
    Tstop  = CFbEvents.POS(trId) + TrialSize(2) - 1;

    creference = U(Rstart:Rstop, :, :);
    ctrial     = U(Tstart:Tstop, :, :);

   
    for fId = 1:NumFreqs
        for chId = 1:NumChans
            ERSP(:, fId, chId, trId) = cnbiproc_ersp(creference(:, fId, chId), ctrial(:, fId, chId), 3);
        end
    end

    
    Ck(trId) = CueEvents.TYP(trId);
    Rk(trId) = labels.Rk(Tstart);
    Mk(trId) = labels.Mk(Tstart);
    Dk(trId) = labels.Dk(Tstart);
end

Days    = unique(Dk);
NumDays = length(Days);

%% Plotting


load('extra/chanlocs64.mat');
ChannelsLb  = {'Fz', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4', 'C3', 'C1', 'Cz', 'C2', 'C4', 'CP3', 'CP1', 'CPz', 'CP2', 'CP4'};
Layout      = [0 0 1 0 0; 2 3 4 5 6; 7 8 9 10 11; 12 13 14 15 16];
LayoutPos   = find(Layout' > 0);
SelClassId  = [771 773];
SelClassLb   = {'BothFeet', 'BothHands'}; 

CFPosition  = abs(TrialSize(1));
CuePosition = CFPosition - min(unique(CueEvents.DUR));

t = TrialPeriod(1):settings.spectrogram.wshift:TrialPeriod(2) - settings.spectrogram.wshift;

NumRows = 4;
NumCols = 5;
fig = zeros(length(SelClassId), 1);

for cId = 1:length(SelClassId)
    fig(cId) = figure;
    cnbifig_set_position(fig(cId), 'All');
    for chId = 1:NumChans
        
        subplot(NumRows, NumCols, LayoutPos(chId));
        cdata = mean(ERSP(:, :, chId, Ck == SelClassId(cId)), 4)';
        imagesc(1:size(ERSP, 1), FreqGrid, cdata,[-1 0.5])
        set(gca,'YDir','normal')
        plot_vline(CuePosition, 'k', 'cue');
        plot_vline(CFPosition, 'k--', 'feedback');
        set(gca, 'XTick', [CuePosition CFPosition size(ERSP, 1)]);
        set(gca, 'XTickLabel', ceil(t([CuePosition CFPosition size(ERSP, 1)])));
        title(ChannelsLb{chId});
        
        [ccol, crow] = ind2sub([NumCols NumRows], LayoutPos(chId));
        
        if (crow == NumRows)
            xlabel('Time [s]');
        else
            set(gca, 'XTick', []);
        end
            
        
        if (ccol == 1 || (ccol == 3 && crow == 1))
            ylabel('Frequency [Hz]');
        else
            set(gca, 'YTick', []);
        end
        
        if (ccol == 3 && crow == 1)
           colorbar;
        end
        
    end
    suptitle([subject ' - ' modality ' - ' SelClassLb{cId}]);
end

%% Saving figures
[~, figurepath] = cnbiutil_mkdir(pwd, figuredir);
for cId = 1:length(SelClassId)
    cnbifig_figure2pdf(fig(cId), [figurepath '/' subject '_' modality '_' SelClassLb{cId} '_ersp.pdf']);
end
