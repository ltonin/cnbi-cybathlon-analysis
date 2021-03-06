clearvars; clc; close all;

SubList = {'AN14VE', 'MA25VE'};
NumSubjects = length(SubList);

datapath  = [pwd '/analysisforce1_strict/'];
figuredir = './figures/';

SelectedClassId = [771 773];
SelectedClassLb = {'Both feet', 'Both hands'};

AltSelectedClassId = [770 771];
AltSelectedClassLb = {'RightHand', 'BothFeet'};

PatternLocationsId = {[4 9 14], [2 7 12 6 11 16]};
PatternLocationsLb = {'FCz-Cz-CPz', 'FC3-C3-CP3-CP4-C4-FC4'}; 
PatternLabels = {'Medial locs', 'Lateral locs'};
NumPatterns = length(PatternLocationsId);
PatternCorr = zeros(NumPatterns, NumSubjects);
PatternPVal = zeros(NumPatterns, NumSubjects);

[~, plotdatapath] = cnbiutil_mkdir(pwd, '/plotdata');
data = struct([]);

CompetitionDay = '20161008';

fisherscore = [];
altfisherscore = [];
Mk = [];
Dk = [];
Sk = [];
Dm = [];
Dml = [];
Dl = cell(NumSubjects, 1);
Cyk = [];
Pk = [];
Pl = [];
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
    
    % Competition day
    cCyId = find(ismember(cdata.discriminancy.run.label.Dl, CompetitionDay, 'rows'));
    cCyk = cDk == cCyId;
    Cyk = cat(1, Cyk, cCyk); 
    
    cDmk = zeros(length(cDk), 1);
    %cDml = str2double(Dl{sId}(cDk, 5:6));
    cDml = str2num(Dl{sId}(cDk, 5:6));
    
    cPk = cdata.discriminancy.run.label.Pk;
    Pk = cat(1, Pk, cPk);
    Pl{sId} = cdata.discriminancy.run.label.Pl;
    
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

    SelectedCombination    = find(ismember(combinations, SelectedClassId, 'rows'));
    AltSelectedCombination = find(ismember(combinations, AltSelectedClassId, 'rows'));
    
    BetaFreqs = 22:2:32;
    [~, SelBetaFreqIds] = intersect(FreqGrid, BetaFreqs);
    
    % Get fisherscore for the main combination (both hand vs. both feet)
    sfisherscore = reshape(cfisherscore(:, :, SelectedCombination), [NumFreqs NumChans size(cfisherscore, 2)]);
    
    % Get fisherscore for the alternative combination (right hand vs. both feet)
    afisherscore = reshape(cfisherscore(:, :, AltSelectedCombination), [NumFreqs NumChans size(cfisherscore, 2)]);
    
    % If fisherscore for the main combination is nan, then use the
    % alternative combination
    ffisherscore = sfisherscore;
    for vId = 1:size(ffisherscore, 3)
        if sum(sum(isnan(ffisherscore(:, :, vId)))) > 0
            ffisherscore(:, :, vId) = afisherscore(:, :, vId);
        end
    end
            
    
    fisherscore = cat(3, fisherscore, ffisherscore);
    cindex = Sk == sId;
    for pId = 1:length(PatternLocationsId)
        cvalues = squeeze(nanmean(nanmean(fisherscore(SelBetaFreqIds, PatternLocationsId{pId}, cindex), 1), 2))';
        [PatternCorr(pId, sId), PatternPVal(pId, sId)] = corr((1:length(cvalues))', cvalues', 'rows', 'pairwise', 'type', 'pearson');
    end
    
%     AltSelectedCombination = find(ismember(combinations, AltSelectedClassId, 'rows'));
%     altfisherscore = cat(3, altfisherscore, );
end

for sId = 1:NumSubjects
    csubject = SubList{sId};
    cnbiutil_bdisp(['Correlation (discriminancy vs. run index) for ' csubject ':']);
    for pId = 1:length(PatternLabels)
        disp([PatternLabels{pId} ': r= ' num2str(PatternCorr(pId, sId)) ', pval= ' num2str(PatternPVal(pId, sId), '%3.3e')]);
    end
end

Modalities   = unique(Mk);
NumModalities = length(Modalities);
ModalitiesLb = {'Offline', 'Online', 'Race'};

 %% Plotting

fig1 = figure;
cnbifig_set_position(fig1, 'All');
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
            maplimits = [0 0.3];
        elseif strcmp(csubject, 'MA25VE')
            maplimits = [0 0.3];
        end
        topoplot(tdata, chanlocs, 'headrad', 'rim', 'maplimits', maplimits, 'electrodes', 'off');
        axis image;

        title(ModalitiesLb{mId});
        if mId == 1
            h = axes('Position', get(gca, 'Position'), 'Visible', 'off');
            set(h.YLabel, 'Visible', 'on');
            ylabel(csubject);
        end
    end
    
      
    % Storing plot data
    data(sId).fisherscore          = fisherscore;
    data(sId).beta_frequency_id    = SelBetaFreqIds;
    data(sId).chanlocs             = chanlocs;
    data(sId).labels.subject_id    = Sk;
    data(sId).labels.modality_id   = Mk;
    data(sId).labels.modality_name = ModalitiesLb;
end
suptitle(['Discriminancy - Modality - Beta Band - ' SelectedClassLb{1} '/' SelectedClassLb{2}]);

%% Saving plot data
cnbiutil_bdisp(['Saving plot data FigS2 in ' plotdatapath]);
save([plotdatapath '/FigS2.mat'], 'data');
data = struct([]);
tmpdata = struct([]);
%% Saving plots

cnbifig_export(fig1, [figuredir '/cybathlon.journal.discriminancy.modality.topoplot.png'], '-png');
cnbifig_export(fig1, [figuredir '/cybathlon.journal.discriminancy.modality.topoplot.pdf'], '-pdf');
%% Figure 2
fig2 = figure;
cnbifig_set_position(fig2, 'Top');

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
    xlim([XThicks(1) XThicks(end) + 1]);
    cDk = Dk(cindex);
    cticks = XThicks(isnan(cdata(:, 1)) == 0);
    cdates = cDk(cticks);
    changeId = find(diff(cdates)) + 1;
    set(gca, 'XTick', cticks(changeId));
    xlabelfontsize = 10;
    set(gca, 'XTickLabel', [repmat(['\fontsize{' num2str(xlabelfontsize) '} '], size(Dl{sId}(cdates(changeId), :), 1), 1) Dl{sId}(cdates(changeId), :)]);
    
    xticklabel_rotate([],90,[])
    set(ax, 'Visible', 'off');
    
    MeanVec{sId} = [];
    StdVec{sId} = [];
    NVec{sId} = [];
    pVal{sId} = [];
    for pId = 1:length(PatternLocationsId)
        cpatterns = PatternLocationsId{pId};
        odata = [];
        oDk = [];
        for mId = 1:NumModalities
            cdata = squeeze(nanmean(nanmean(fisherscore(SelBetaFreqIds, cpatterns, cindex & Mk == Modalities(mId)), 1), 2));
            plot(XThicks(Mk(cindex) == Modalities(mId)), cdata, 'LineStyle', 'none', 'Marker', Marker{mId}, 'MarkerEdgeColor', Colors(pId, :), 'MarkerFaceColor', Colors(pId, :), 'MarkerSize', 4);
            odata = [odata; cdata];
            oDk = [oDk; Dk(cindex & Mk == Modalities(mId))];
        end
        % Put them in chronological order
        [oDk, sortIndoDk] = sort(oDk,'ascend');
        odata = odata(sortIndoDk);
        oDays = length(unique(oDk));
        
        % Save pre-post (first and last four sessions for later use)
        MeanVec{sId} = [MeanVec{sId} [nanmean(odata(find(ismember(oDk,[1:4])))) ; nanmean(odata(find(ismember(oDk,[max(oDays)-3:max(oDays)]))))]];
        StdVec{sId} = [StdVec{sId} [nanstd(odata(find(ismember(oDk,[1:4])))) ; nanstd(odata(find(ismember(oDk,[max(oDays)-3:max(oDays)]))))]];
        NVec{sId} = [NVec{sId} [size(odata(find(ismember(oDk,[1:4]))),1) ; size(odata(find(ismember(oDk,[max(oDays)-3:max(oDays)]))),1)]];
        pVal{sId} = [pVal{sId} ranksum(odata(find(ismember(oDk,[1:4]))), odata(find(ismember(oDk,[max(oDays)-3:max(oDays)]))))];

        % Storing for after
        tmpdata(sId).fisherscore(:, pId) = odata;
        tmpdata(sId).labels.day_id = oDk;
    end
    
    hold off;
    
    ylim([0 0.8]);

    h = gca;
    [~, hobjl] = legend(h.Children([end-3 end-2 6 5 4]), [PatternLabels ModalitiesLb]);
    set(hobjl(end), 'MarkerFaceColor', 'none');   set(hobjl(end),   'MarkerEdgeColor', 'k')
    set(hobjl(end-2), 'MarkerFaceColor', 'none'); set(hobjl(end-2), 'MarkerEdgeColor', 'k')
    set(hobjl(end-4), 'MarkerFaceColor', 'none'); set(hobjl(end-4), 'MarkerEdgeColor', 'k')
    
    
    grid on;
    ylabel('Discriminancy');
    hxlabel = xlabel('Session');
    set(hxlabel, 'Position', get(hxlabel, 'Position') - [0 0.02 0])
    title(csubject);

    cpos = get(gca,'Position');
    for pId = 1:length(PatternLocationsId)
        apos = cpos;
        apos(2) = (cpos(2) - 0.05) - (pId-1)*0.035;
        textcolor = 'k';
        annotation('textbox', apos, 'String', ['r=' num2str(PatternCorr(pId, sId), '%3.2f') ', p=' num2str(PatternPVal(pId, sId), '%3.3f')], 'LineStyle', 'none', 'Color', Colors(pId, :), 'FontWeight', 'bold')
    end
    
    % Storing plot data
    data(sId).fisherscore          = fisherscore;
    data(sId).beta_frequency_id    = SelBetaFreqIds;
    data(sId).chanlocs             = chanlocs;
    data(sId).labels.subject_id    = Sk;
    data(sId).labels.modality_id   = Mk;
    data(sId).labels.modality_name = ModalitiesLb;
    data(sId).labels.pattern_id    = PatternLocationsId;
    data(sId).labels.pattern_name  = PatternLabels;
    data(sId).labels.day_id        = cDk;
    data(sId).labels.day_date      = Dl{sId};
    
    % Storing for after
    tmpdata(sId).labels.pattern_name  = PatternLabels;
   
end
suptitle(['Discriminancy - Emerging patterns - correlation - Beta Band - ' SelectedClassLb{1} '/' SelectedClassLb{2}]);

%% Saving plot data
cnbiutil_bdisp(['Saving plot data Fig4B in ' plotdatapath]);
save([plotdatapath '/Fig4B.mat'], 'data');
data = tmpdata;
cnbiutil_bdisp(['Saving plot data Fig4C in ' plotdatapath]);
save([plotdatapath '/Fig4C.mat'], 'data');
data = struct([]);
tmpdata = struct([]);

%% Saving plots
cnbifig_export(fig2, [figuredir '/cybathlon.journal.discriminancy.emerging.correlation.png'], '-png');
cnbifig_export(fig2, [figuredir '/cybathlon.journal.discriminancy.emerging.correlation.pdf'], '-pdf');

%% Plot 2bis
fig2b = figure;
cnbifig_set_position(fig2b, 'Top');


ax = gca;
Colors = ax.ColorOrder;
RaceModality = 2;

csubject = 'AN14VE';
csubjectId = cell2mat(strfind(SubList, csubject));

cColors = {'b', 'r'};
SelParadigmId = [3 5 6];

if isempty(csubjectId)
    warning([csubject ' not found']);
else
    
    cindex = Sk == csubjectId & Mk == RaceModality;
    
    hold on;
    XThicks = 1:sum(cindex);
    cdata = [];
    MeanF = zeros(1, length(SelParadigmId), length(PatternLocationsId));
    StdF = zeros(2, length(SelParadigmId), length(PatternLocationsId));
    for pId = 1:length(PatternLocationsId)
        
        subplot(1, 3, pId);
        cpatterns = PatternLocationsId{pId};
        cdata = squeeze(nanmean(nanmean(fisherscore(SelBetaFreqIds, cpatterns, cindex), 1), 2));
        
        
        
        hold on;
        for i = 1:length(SelParadigmId)
            ccthick = Pk(cindex) == SelParadigmId(i);
            ax = plot(XThicks(ccthick),  cdata(ccthick), '.', 'Color', Colors(SelParadigmId(i), :));
            [tmpcorr, tmppval] = corr(XThicks(ccthick)', cdata(ccthick) , 'rows', 'pairwise', 'type', 'pearson');
            disp([PatternLabels{pId} '- Paradigm' num2str(i) ': r= ' num2str(tmpcorr) ', pval= ' num2str(tmppval)]);
            
            ccindex = Pk(cindex) == SelParadigmId(i);
            cvalues = cdata(ccindex);
            MeanF(1, i, pId) = nanmean(cvalues(1:4));
            StdF(1, i, pId) = nanstd(cvalues(1:4));
            MeanF(2, i, pId) = nanmean(cvalues(end-3:end));
            StdF(2, i, pId) = nanstd(cvalues(end-3:end));
        end
        hold off;
        hl = lsline;
        set(hl, 'LineWidth', 3);
        xlim([XThicks(1) XThicks(end) + 1]);
    
        grid on;
        
        cDk = Dk(cindex);
        cticks = XThicks(isnan(cdata(:, 1)) == 0);
        cdates = cDk(cticks);
        changeId = find(diff(cdates)) + 1;
        set(gca, 'XTick', cticks(changeId));
        xlabelfontsize = 10;
        set(gca, 'XTickLabel', [repmat(['\fontsize{' num2str(xlabelfontsize) '} '], size(Dl{csubjectId}(cdates(changeId), :), 1), 1) Dl{csubjectId}(cdates(changeId), :)]);

        xticklabel_rotate([],90,[])
        ylim([0 0.8]);
        title(PatternLabels{pId});
    end
    
    subplot(1, 3, 3)
    bar(cat(2, MeanF(:, :, 1), MeanF(:, :, 2))');
    %plot_vline(3.5, 'k');
    ylim([0 0.4]);
    title('Average discriminancy');
    legend('First runs', 'Last runs');
    xlabel('Paradigm/Locations');
    grid on;
    set(gca, 'XTickLabel', [1 2 3 1 2 3]);
    
    
    
    
%     cDk = Dk(cindex);
%     cticks = XThicks(isnan(cdata(:, 1)) == 0);
%     cdates = cDk(cticks);
%     changeId = find(diff(cdates)) + 1;
%     set(gca, 'XTick', cticks(changeId));
%     xlabelfontsize = 10;
%     set(gca, 'XTickLabel', [repmat(['\fontsize{' num2str(xlabelfontsize) '} '], size(Dl{csubjectId}(cdates(changeId), :), 1), 1) Dl{csubjectId}(cdates(changeId), :)]);
%     
%     xticklabel_rotate([],90,[])
%     set(ax, 'Visible', 'off');
    
%     MeanVec{sId} = [];
%     StdVec{sId} = [];
%     NVec{sId} = [];
%     pVal{sId} = [];
%     for pId = 1:length(PatternLocationsId)
%         cpatterns = PatternLocationsId{pId};
%         odata = [];
%         oDk = [];
%         
%         cdata = squeeze(nanmean(nanmean(fisherscore(SelBetaFreqIds, cpatterns, cindex & Mk == RaceModality), 1), 2));
%         plot(XThicks(Mk(cindex) == RaceModality), cdata, 'LineStyle', 'none', 'Marker', Marker{mId}, 'MarkerEdgeColor', Colors(pId, :), 'MarkerFaceColor', Colors(pId, :), 'MarkerSize', 4);
%         odata = [odata; cdata];
%         oDk = [oDk; Dk(cindex & Mk == RaceModality)];
%         
%         % Put them in chronological order
%         [oDk, sortIndoDk] = sort(oDk,'ascend');
%         odata = odata(sortIndoDk);
%         oDays = length(unique(oDk));
%         
%         % Save pre-post (first and last four sessions for later use)
%         MeanVec{sId} = [MeanVec{sId} [nanmean(odata(find(ismember(oDk,[1:4])))) ; nanmean(odata(find(ismember(oDk,[max(oDays)-3:max(oDays)]))))]];
%         StdVec{sId} = [StdVec{sId} [nanstd(odata(find(ismember(oDk,[1:4])))) ; nanstd(odata(find(ismember(oDk,[max(oDays)-3:max(oDays)]))))]];
%         NVec{sId} = [NVec{sId} [size(odata(find(ismember(oDk,[1:4]))),1) ; size(odata(find(ismember(oDk,[max(oDays)-3:max(oDays)]))),1)]];
%         pVal{sId} = [pVal{sId} ranksum(odata(find(ismember(oDk,[1:4]))), odata(find(ismember(oDk,[max(oDays)-3:max(oDays)]))))];
% 
%     end
    

    
    

%     h = gca;
%     [~, hobjl] = legend(h.Children([end-3 end-2 6 5 4]), [PatternLabels 'Race']);
%     set(hobjl(end), 'MarkerFaceColor', 'none');   set(hobjl(end),   'MarkerEdgeColor', 'k')
%     set(hobjl(end-2), 'MarkerFaceColor', 'none'); set(hobjl(end-2), 'MarkerEdgeColor', 'k')
%     set(hobjl(end-4), 'MarkerFaceColor', 'none'); set(hobjl(end-4), 'MarkerEdgeColor', 'k')
%     
%     

%     ylabel('Discriminancy');
%     hxlabel = xlabel('Session');
%     set(hxlabel, 'Position', get(hxlabel, 'Position') - [0 0.02 0])
%     title(csubject);
% 
%     cpos = get(gca,'Position');
%     for pId = 1:length(PatternLocationsId)
%         apos = cpos;
%         apos(2) = (cpos(2) - 0.05) - (pId-1)*0.035;
%         textcolor = 'k';
%         annotation('textbox', apos, 'String', ['r=' num2str(PatternCorr(pId, sId), '%3.2f') ', p=' num2str(PatternPVal(pId, sId), '%3.3f')], 'LineStyle', 'none', 'Color', Colors(pId, :), 'FontWeight', 'bold')
%     end
    
    
end

suptitle('P1 - Discriminancy per control paradigm');
cnbifig_export(fig2b, [figuredir '/cybathlon.journal.discriminancy.paradigm.png'], '-png');
cnbifig_export(fig2b, [figuredir '/cybathlon.journal.discriminancy.paradigm.pdf'], '-pdf');

%% Plot 3
fig3 = figure;
cnbifig_set_position(fig3, 'Top');

for sId = 1:NumSubjects
    csubject = SubList{sId};
    subplot(NumRows, NumCols, sId);
    hold on;
    offset = 0.2;
    cellsigstart = {};
    for b=1:2 % Location condition (medial, lateral)
        hBar(b,:) = bar([b-offset b+offset],MeanVec{sId}(:,b),'LineWidth',2,'FaceColor',Colors(b,:));
    end
    % The separate for loops is for the legend to show correctly
    for b=1:2
        hEBar(b,:) = errorbar([b-offset b+offset],MeanVec{sId}(:,b),zeros(1,2),StdVec{sId}(:,b),'.','LineWidth',2,'Color',Colors(b,:));
        cellsigstart{b} = [b-offset,b+offset];
    end
    H = sigstar(cellsigstart,pVal{sId}); 
    hold off;
    
    grid on;    
    ylim([0 0.6]);
    xlim([0.5 2.5]);
   
    set(gca, 'XTick', cell2mat(cellsigstart)-0.01);
    set(gca, 'XTickLabel', repmat({'First 4 sessions','Last 4 sessions'},1,2));
    xticklabel_rotate([],90,[]);
    %hxlabel = xlabel('Session');
    %set(hxlabel, 'Position', get(hxlabel, 'Position') - [0 0.02 0])
    ylabel('Discriminancy');
    title(csubject);
    
    legend(PatternLabels, 'location', 'SouthEast');    
   
    
end

cnbifig_export(fig3, [figuredir '/cybathlon.journal.discriminancy.emerging.firstlast.png'], '-png');
cnbifig_export(fig3, [figuredir '/cybathlon.journal.discriminancy.emerging.firstlast.pdf'], '-pdf');





%% Plot 4
fig4 = figure;
cnbifig_set_position(fig4, 'All');


NumRows = NumSubjects;
NumCols = max(Dml) - min(Dml) + 1 + 1;
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
%         if sum(isnan(cdata) > 0)
%             cdata = squeeze(nanmean(nanmean(altfisherscore(SelBetaFreqIds, :, cindex), 1), 3));
%         end
        
        tdata = convChans(cdata);

        if strcmp(csubject, 'AN14VE')
            maplimits = [0 0.4];
        elseif strcmp(csubject, 'MA25VE')
            maplimits = [0 0.4];
        end
        topoplot(tdata, chanlocs, 'headrad', 'rim', 'maplimits', maplimits, 'electrodes', 'off');
        
        [~, cmonthname] = month(num2str(unique(Dml(Dml == cmonths(dmId)))), 'mm');
        title(cmonthname);       
        
   end
    
    % Competition day
    subplot(NumRows, NumCols, NumCols + NumCols*(sId -1));
    cindex = Sk == sId & Cyk == true;
    cdata = squeeze(nanmean(nanmean(fisherscore(SelBetaFreqIds, :, cindex), 1), 3));
    tdata = convChans(cdata);

    if strcmp(csubject, 'AN14VE')
        maplimits = [0 0.4];
    elseif strcmp(csubject, 'MA25VE')
        maplimits = [0 0.4];
    end
    topoplot(tdata, chanlocs, 'headrad', 'rim', 'maplimits', maplimits, 'electrodes', 'off');
    title('Competition');
   
    u = subplot(NumRows, NumCols, 1 + NumCols*(sId -1));
    h = axes('Position', get(gca, 'Position'), 'Visible', 'off');
    set(u, 'Visible', 'off');
    set(h.YLabel, 'Visible', 'on');
    ylabel(csubject);
   
          
    % Storing plot data
    data(sId).fisherscore          = fisherscore;
    data(sId).beta_frequency_id    = SelBetaFreqIds;
    data(sId).chanlocs             = chanlocs;
    data(sId).labels.subject_id    = Sk;
    data(sId).labels.month_id      = Dml;
end
suptitle(['Discriminancy - Emerging patterns - topoplot - Beta Band - ' SelectedClassLb{1} '/' SelectedClassLb{2}])

%% Saving plot data
cnbiutil_bdisp(['Saving plot data Fig4A in ' plotdatapath]);
save([plotdatapath '/Fig4A.mat'], 'data');
data = struct([]);

%% Saving plots
cnbifig_export(fig4, [figuredir '/cybathlon.journal.discriminancy.emerging.topoplot.png'], '-png');
cnbifig_export(fig4, [figuredir '/cybathlon.journal.discriminancy.emerging.topoplot.pdf'], '-pdf');

%% Competition
fig5 = figure;
cnbifig_set_position(fig5, 'All');
maplimits = [0 0.5];
for sId = 1:NumSubjects
    csubject = SubList{sId};
    competRunId = max(Dk(Sk == sId));
    cindex = Sk == sId & Dk == competRunId;

    cdata = squeeze(nanmean(fisherscore(SelBetaFreqIds, :, cindex), 1));
    
   
    u = subplot(2, 2, 1 + 2*(sId -1));
    tdata = convChans(cdata(:, 1));
    topoplot(tdata, chanlocs, 'headrad', 'rim', 'maplimits', maplimits, 'electrodes', 'off');
    title('Qualifier');
    
    h = axes('Position', get(gca, 'Position'), 'Visible', 'off');
    set(u, 'Visible', 'off');
    set(h.YLabel, 'Visible', 'on');
    ylabel(csubject);
    
    subplot(2, 2, 2 + 2*(sId -1));
    tdata = convChans(cdata(:, 2));
    topoplot(tdata, chanlocs, 'headrad', 'rim', 'maplimits', maplimits, 'electrodes', 'off');
    title('Final');
    
        % Storing plot data
    data(sId).fisherscore          = fisherscore;
    data(sId).beta_frequency_id    = SelBetaFreqIds;
    data(sId).chanlocs             = chanlocs;
    data(sId).labels.subject_id    = Sk;
    data(sId).labels.day_id        = Dk;
end
suptitle('Discriminancy - Competition');

%% Saving plot data
cnbiutil_bdisp(['Saving plot data Fig5 in ' plotdatapath]);
save([plotdatapath '/Fig5.mat'], 'data');
data = struct([]);

%% Saving plots
cnbifig_export(fig5, [figuredir '/cybathlon.journal.discriminancy.competition.topoplot.png'], '-png');
cnbifig_export(fig5, [figuredir '/cybathlon.journal.discriminancy.competition.topoplot.pdf'], '-pdf');
