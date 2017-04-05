clearvars; clc; close all;

SubList = {'AN14VE', 'MA25VE'};
NumSubjects = length(SubList);

datapath  = [pwd '/analysis/'];
figuredir = './figures/';

PadNames = {'Spin', 'Jump', 'Slide'};

ax = gca;
Colors = ax.ColorOrder;
close all;
% Speed, Jump, Slide, Rest
PadColors = {Colors(6, :), Colors(4, :), Colors(3, :), 'k', 'k', 'k'};

NumRows = 1;
NumCols = NumSubjects;

fig1 = figure;
cnbifig_set_position(fig1, 'Top');

for sId = 1:NumSubjects
    
    csubject = SubList{sId};
    cfilepath = [datapath '/' csubject '.padaccuracy.race.mat'];
    cdata = load(cfilepath);
    
    cRk = cdata.pad.label.run.Rk;
    cruns = unique(cRk);
    cnruns = length(cruns);
    
    cDk = cdata.pad.label.run.Dk;
    cDays = unique(cDk);
    cndays = length(cDays);
    
    caccuracyday = zeros(cndays, size(cdata.pad.accuracy.run.values, 1) - 1);
    for dId = 1:cndays
        cvaluesday = cdata.pad.accuracy.run.values(1:3, cDk == cDays(dId));
        cvaluesday(isnan(cvaluesday)) = 0;
        caccuracyday(dId, :) = mean(cvaluesday, 2);
    end
    
    caccuracyrun = zeros(cnruns, 4);
    for rId = 1:cnruns
       cvaluesrun =  cdata.pad.accuracy.run.values(1:3, rId);
       cvaluesrun(isnan(cvaluesrun)) = 0;
       caccuracyrun(rId,:) = [cvaluesrun ; mean(cvaluesrun)]';
    end
    
    
    subplot(NumRows, NumCols, sId);
    

    hold on;
    
    % Plotting bars of first and last
    MeanVec = [mean(caccuracyrun(find(ismember(cDk,[1:4])),:)) ; mean(caccuracyrun(find(ismember(cDk,[max(cDays)-3:max(cDays)])),:))];
    StdVec = [std(caccuracyrun(find(ismember(cDk,[1:4])),:)) ; std(caccuracyrun(find(ismember(cDk,[max(cDays)-3:max(cDays)])),:))];
    NVec = [size(caccuracyrun(find(ismember(cDk,[1:4]))),1) ; size(caccuracyrun(find(ismember(cDk,[max(cDays)-3:max(cDays)]))),1)];
    pVal = ranksum(caccuracyrun(find(ismember(cDk,[1:4]))), caccuracyrun(find(ismember(cDk,[max(cDays)-3:max(cDays)]))));
    
    offset = 0.2;
    cellsigstart = {};
    for b=1:4
        hBar(b,:) = bar([b-offset b+offset],MeanVec(:,b),'LineWidth',2,'FaceColor',PadColors{b});
    end
    % The separate for loops is for the legend to show correctly
    for b=1:4
        hEBar(b,:) = errorbar([b-offset b+offset],MeanVec(:,b),zeros(1,2),StdVec(:,b),'.','LineWidth',2,'Color',PadColors{b});
        p(b) = ranksum(caccuracyrun(find(ismember(cDk,[1:4])),b), caccuracyrun(find(ismember(cDk,[max(cDays)-3:max(cDays)])),b));
        cellsigstart{b} = [b-offset,b+offset];
    end
    H = sigstar(cellsigstart,p); 
    hold off;
    
    grid on;    
    ylim([0 120]);
    xlim([0 5]);
    ytick = get(gca, 'YTick');
    set(gca, 'YTick', intersect(ytick, 1:100));
    
    set(gca, 'XTick', cell2mat(cellsigstart)-0.1);
    set(gca, 'XTickLabel', repmat({'First 4 sessions','Last 4 sessions'},1,4));
    xticklabel_rotate([],90,[]);
    %hxlabel = xlabel('Session');
    %set(hxlabel, 'Position', get(hxlabel, 'Position') - [0 0.02 0])
    ylabel('Accuracy [%]');
    title(csubject);
    
    legend([PadNames, 'Average'], 'location', 'SouthEast');
    
end

suptitle('Command accuracy');
cnbifig_export(fig1, [figuredir '/cybathlon.journal.accuracy.pad.png'], '-png');