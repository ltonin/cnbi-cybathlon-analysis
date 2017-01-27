clearvars; clc; close all;

SubList = {'AN14VE', 'MA25VE'};
NumSubjects = length(SubList);

datapath  = [pwd '/analysis/'];
figuredir = './figures/';

SelectedClassId = [771 773];
SelectedClassLb = {'BothFeet', 'BothHands'};

AltSelectedClassId = [770 771];
AltSelectedClassLb = {'RightHand', 'BothFeet'};

PatternLocationsId = {[4 9 14], [2 7 12 6 11 16]};
PatternLocationsLb = {'FCz-Cz-CPz', 'FC3-C3-CP3-CP4-C4-FC4'}; 
PatternLabels = {'Middle locs', 'Lateral locs'};
NumPatterns = length(PatternLocationsId);
PatternCorr = zeros(NumPatterns, NumSubjects);
PatternPVal = zeros(NumPatterns, NumSubjects);



fisherscore = [];
altfisherscore = [];
Mk = [];
Dk = [];
Sk = [];
Dm = [];
Dml = [];
Dl = cell(NumSubjects, 1);
for sId = 1:NumSubjects
    csubject = SubList{sId};
    cfilepath = [datapath '/' csubject '.discriminancy.maps.mat'];
    
    cdata = load(cfilepath);     
    
    cfisherscore = cdata.discriminancy.run.fisherscore;
    FreqGrid = cdata.discriminancy.freqs;
    
    NumFreqs = length(FreqGrid);
    NumChans = size(cfisherscore, 1)./NumFreqs;
    
    combinations = cdata.discriminancy.run.combinations;
    NumCombinations = size(combinations, 1);
    
    cMk = cdata.discriminancy.run.label.Mk;
    cMk(cMk == 3) = 2;
    Mk = cat(1, Mk, cMk);
    
    cDk = cdata.discriminancy.run.label.Dk;
    Dk = cat(1, Dk, cDk);
    Dl{sId} = cdata.discriminancy.run.label.Dl;
    
    cDmk = zeros(length(cDk), 1);
    cDml = str2double(Dl{sId}(cDk, 5:6));
    
    cmonthId = 0;
    cmonthLb = [];
    for dId = 1:length(cDml)
        if isequal(cmonthLb, cDml(dId)) == false
            cmonthLb = cDml(dId);
            cmonthId = cmonthId + 1;
        end
        cDmk(dId) = cmonthId;
    end
    
    Dm  = cat(1, Dm, cDmk);
    Dml = cat(1, Dml, cDml);
    
    Sk = cat(1, Sk, sId*ones(length(cdata.discriminancy.run.label.Mk), 1));

    SelectedCombination = find(ismember(combinations, SelectedClassId, 'rows'));
    
    BetaFreqs = 22:2:32;
    [~, SelBetaFreqIds] = intersect(FreqGrid, BetaFreqs);
    
    

    fisherscore = cat(3, fisherscore, reshape(cfisherscore(:, :, SelectedCombination), [NumFreqs NumChans size(cfisherscore, 2)]));
    cindex = Sk == sId;
    for pId = 1:length(PatternLocationsId)
        cvalues = squeeze(nanmean(nanmean(fisherscore(SelBetaFreqIds, PatternLocationsId{pId}, cindex), 1), 2))';
        [PatternCorr(pId, sId), PatternPVal(pId, sId)] = corr((1:length(cvalues))', cvalues', 'rows', 'pairwise', 'type', 'spearman');
        
    end
    
    AltSelectedCombination = find(ismember(combinations, AltSelectedClassId, 'rows'));
    altfisherscore = cat(3, altfisherscore, reshape(cfisherscore(:, :, AltSelectedCombination), [NumFreqs NumChans size(cfisherscore, 2)]));
end

Modalities   = unique(Mk);
NumModalities = length(Modalities);
ModalitiesLb = {'Offline', 'Online', 'Race'};

 %% Plotting

fig1 = figure;
fig_set_position(fig1, 'Top');
load('chanlocs64.mat');

NumRows = NumSubjects;
NumCols = NumModalities;

for sId = 1:NumSubjects
    csubject = SubList{sId};
    for mId = 1:NumModalities
        subplot(NumRows, NumCols, mId + NumCols*(sId -1));

        cindex = Sk == sId & Mk == Modalities(mId);
        cdata = squeeze(nanmean(nanmean(fisherscore(SelBetaFreqIds, :, cindex), 1), 3));
        tdata = convChans(cdata);

        if strcmp(csubject, 'AN14VE')
            maplimits = [0 0.4];
        elseif strcmp(csubject, 'MA25VE')
            maplimits = [0 0.5];
        end
        topoplot(tdata, chanlocs, 'headrad', 'rim', 'maplimits', maplimits);
        axis image;

        title(ModalitiesLb{mId});
        if mId == 1
            h = axes('Position', get(gca, 'Position'), 'Visible', 'off');
            set(h.YLabel, 'Visible', 'on');
            ylabel(csubject);
        end
    end
end
suptitle(['Discriminancy - Pattern stability - Beta Band (' SelectedClassLb{1} '-' SelectedClassLb{2} ')']);
cnbifig_export(fig1, [figuredir '/cybathlon.journal.discriminancy.stability.topoplot.png'], '-png');

fig2 = figure;
fig_set_position(fig2, 'Top');

NumRows = 1;
NumCols = NumSubjects;

Marker = {'v', 's', 'o'};

ax = gca;
Colors = ax.ColorOrder;

for sId = 1:NumSubjects
    csubject = SubList{sId};
    subplot(NumRows, NumCols, sId);
    cindex = Sk == sId;
    hold on;
    XThicks = 1:sum(cindex);
    cdata = [];
    for pId = 1:length(PatternLocationsId)
        cpatterns = PatternLocationsId{pId};
        cdata = cat(2, cdata, squeeze(nanmean(nanmean(fisherscore(SelBetaFreqIds, cpatterns, cindex), 1), 2)));
    end
    ax = plot(XThicks,  cdata, '.');
    hl = lsline;
    set(hl, 'LineWidth', 3);
    cDk = Dk(cindex);
    cticks = XThicks(isnan(cdata(:, 1)) == 0);
    cdates = cDk(cticks);
    changeId = find(diff(cdates)) + 1;
    set(gca, 'XTick', cticks(changeId));
    set(gca, 'XTickLabel', Dl{sId}(cdates(changeId), :));
    xticklabel_rotate([],90,[])
    set(ax, 'Visible', 'off');
    
        
    for pId = 1:length(PatternLocationsId)
        cpatterns = PatternLocationsId{pId};
        for mId = 1:NumModalities
            cdata = squeeze(nanmean(nanmean(fisherscore(SelBetaFreqIds, cpatterns, cindex & Mk == Modalities(mId)), 1), 2));
            plot(XThicks(Mk(cindex) == Modalities(mId)), cdata, 'LineStyle', 'none', 'Marker', Marker{mId}, 'MarkerEdgeColor', Colors(pId, :), 'MarkerFaceColor', Colors(pId, :), 'MarkerSize', 4);
        end
    end
    
    hold off;
    ylim([0 0.8]);

    h = gca;
    legend(h.Children([end-3 end-2 6 5 4]), [PatternLabels ModalitiesLb])

    
    grid on;
    ylabel('Discriminancy');
    xlabel('Runs');
    title(csubject);

    cpos = get(gca,'Position');
    for pId = 1:length(PatternLocationsId)
        apos = cpos;
        apos(2) = (cpos(2) - 0.05) - (pId-1)*0.035;
        textcolor = 'k';
        annotation('textbox', apos, 'String', ['rho=' num2str(PatternCorr(pId, sId), '%3.2f') ', p=' num2str(PatternPVal(pId, sId), '%3.3f')], 'LineStyle', 'none', 'Color', Colors(pId, :), 'FontWeight', 'bold')
    end
    
    
end

suptitle(['Discriminancy - Emerging patterns correlation - Beta Band (' SelectedClassLb{1} '-' SelectedClassLb{2} ')']);
cnbifig_export(fig2, [figuredir '/cybathlon.journal.discriminancy.emerging.correlation.png'], '-png');

%% Plot 3
fig3 = figure;
fig_set_position(fig3, 'Top');


NumRows = NumSubjects;
NumCols = max(Dml) - min(Dml) + 1;
PlotLoc = min(Dml):max(Dml);

Months = unique(Dml);
NumMonths = length(Months);

for sId = 1:NumSubjects
   csubject = SubList{sId};
   cmonths = unique(Dml(Sk == sId));
   for dmId = 1:length(cmonths)
        crowloc = find(PlotLoc == cmonths(dmId));
        subplot(NumRows, NumCols, crowloc + NumCols*(sId -1));
        cindex = Sk == sId & Dml == cmonths(dmId);
        cdata = squeeze(nanmean(nanmean(fisherscore(SelBetaFreqIds, :, cindex), 1), 3));
        if sum(isnan(cdata) > 0)
            cdata = squeeze(nanmean(nanmean(altfisherscore(SelBetaFreqIds, :, cindex), 1), 3));
        end
        
        tdata = convChans(cdata);

        if strcmp(csubject, 'AN14VE')
            maplimits = [0 0.5];
        elseif strcmp(csubject, 'MA25VE')
            maplimits = [0 0.5];
        end
        topoplot(tdata, chanlocs, 'headrad', 'rim', 'maplimits', maplimits);
        
        [~, cmonthname] = month(num2str(unique(Dml(Dml == cmonths(dmId)))), 'mm');
        title(cmonthname);       
        
   end
           
    u = subplot(NumRows, NumCols, 1 + NumCols*(sId -1));
    h = axes('Position', get(gca, 'Position'), 'Visible', 'off');
    set(u, 'Visible', 'off');
    set(h.YLabel, 'Visible', 'on');
    ylabel(csubject);
   
end

suptitle(['Discriminancy - Emerging patterns topoplot - Beta Band (' SelectedClassLb{1} '-' SelectedClassLb{2} ')'])
cnbifig_export(fig3, [figuredir '/cybathlon.journal.discriminancy.emerging.topoplot.png'], '-png');