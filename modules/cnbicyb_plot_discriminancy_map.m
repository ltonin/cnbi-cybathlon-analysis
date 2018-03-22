clearvars; clc; close all;

% SubList = {'AN14VE', 'MA25VE'};
SubList = {'MA25VE'};
NumSubjects = length(SubList);

%datapath  = [pwd '/analysisforce1_strict/'];
datapath  = [pwd '/analysis/'];
%figuredir = './analysisforce1_strict/';
figuredir = './analysis/';

SelectedClassId = [771 773];
SelectedClassLb = {'Both feet', 'Both hands'};

AltSelectedClassId = [770 771];
AltSelectedClassLb = {'RightHand', 'BothFeet'};

[~, plotdatapath] = cnbiutil_mkdir(pwd, '/plotdata');
data = struct([]);

CompetitionDay = '20161008';
freqs = 4:2:30;
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
    
    %BetaFreqs = 22:2:32;
    BetaFreqs = 22:2:48;
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

end

Modalities   = unique(Mk);
NumModalities = length(Modalities);
ModalitiesLb = {'Offline', 'Online', 'Race'};




%% Plot 4
fig4 = figure;
cnbifig_set_position(fig4, 'All');

SelFreqs = 4:2:48;
[~, SelFreqIds] = intersect(FreqGrid, SelFreqs);

NumRows = NumSubjects;
NumCols = max(Dml) - min(Dml) + 1 + 1;
PlotLoc = min(Dml):max(Dml);

Months = unique(Dml);
NumMonths = length(Months);
% colormap(flipud(bone));
for sId = 1:NumSubjects
   csubject = SubList{sId};
   cmonths = unique(Dml(Sk == sId));
   for dmId = 1:length(cmonths)
        crowloc = find(PlotLoc == cmonths(dmId));
        subplot(NumRows, NumCols, crowloc + NumCols*(sId -1));
        cindex = Sk == sId & Dml == cmonths(dmId) ;%& Cyk == false;
        cnruns = sum(cindex);
   
        cdata = nanmean(fisherscore(:, :, cindex), 3);
        
        maplimits = [0 0.6];
        if(cnruns < 5)
            maplimits = [0 1];
        end
        
        
        imagesc(SelFreqs(SelFreqIds), 1:16, cdata(SelFreqIds, :)', maplimits);
        axis image;

        [~, cmonthname] = month(num2str(unique(Dml(Dml == cmonths(dmId)))), 'mm');
        title([cmonthname ' (' num2str(cnruns) ')']);       
        
        
        
        disp([csubject ' - ' cmonthname ' (N=' num2str(cnruns) ')' ]);
   end
    
%     % Competition day
%     subplot(NumRows, NumCols, NumCols + NumCols*(sId -1));
%     cindex = Sk == sId & Cyk == true;
%     cnruns = sum(cindex);
% 
%     cdata = nanmean(fisherscore(:, :, cindex), 3);
% 
% 
%     if(cnruns < 5)
%         maplimits = [0 1.5];
%     else
%         maplimits = [0 0.8];
%     end
%     
%     imagesc(freqs(1:15), 1:16, cdata(1:15, :)', maplimits);
%     axis image;
%     title('Competition');

    u = subplot(NumRows, NumCols, 1 + NumCols*(sId -1));
    h = axes('Position', get(gca, 'Position'), 'Visible', 'off');
    set(u, 'Visible', 'off');
    set(h.YLabel, 'Visible', 'on');
    ylabel(csubject);
   
    
    % Storing plot data
    data(sId).fisherscore          =  fisherscore;
    data(sId).frequency_id_selected = SelFreqIds;
    data(sId).frequency_grid       = SelFreqs;
    data(sId).labels.subject_id    = Sk;
    data(sId).labels.month_id      = Dml;
    
end

% suptitle(['Discriminancy - Emerging patterns - Maps - ' SelectedClassLb{1} '/' SelectedClassLb{2}])

%% Saving plot data
cnbiutil_bdisp(['Saving plot data FigS1 in ' plotdatapath]);
save([plotdatapath '/FigS1.mat'], 'data');

%% Saving plots
cnbifig_export(fig4, [figuredir '/cybathlon.journal.discriminancy.emerging.map.png'], '-png');
cnbifig_export(fig4, [figuredir '/cybathlon.journal.discriminancy.emerging.map.pdf'], '-pdf');



% Simis quick run
%r=27;
% r=41;
%r=70;
r=80;
r=r-1;
close all;
%imagesc(squeeze(fisherscore(1:23,:,r))',[0 0.9]);
imagesc(squeeze(fisherscore(1:23,:,r))');
set(gca,'XTick',[1:23]);
set(gca,'XTickLabel',[4:2:48]);
set(gca,'YTick',[1:16]);
set(gca,'YTickLabel',{'Fz','FC3','FC1','FCz','FC2','FC4','C3','C1','Cz','C2','C4','CP3','CP1','CPz','CP2','CP4'});
title(['P2, ' Dl{1}(Dk(r),:)],'FontSize',30);
set(gca,'FontSize',30);
xlabel('Frequency [Hz]','FontSize',30);
ylabel('Channel','FontSize',30);


% Simis quick session
% s=1;
% s=41;
% s=14;
s=s-1;
close all;
%imagesc(squeeze(fisherscore(1:23,:,r))',[0 0.9]);
imagesc(squeeze(mean(fisherscore(1:23,:,Dk==s),3))');
set(gca,'XTick',[1:23]);
set(gca,'XTickLabel',[4:2:48]);
set(gca,'YTick',[1:16]);
set(gca,'YTickLabel',{'Fz','FC3','FC1','FCz','FC2','FC4','C3','C1','Cz','C2','C4','CP3','CP1','CPz','CP2','CP4'});
title(['P2, ' Dl{1}(s,:)],'FontSize',30);
set(gca,'FontSize',30);
xlabel('Frequency [Hz]','FontSize',30);
ylabel('Channel','FontSize',30);