clearvars; clc; close all;

SubList = {'AN14VE', 'MA25VE'};
NumSubjects = length(SubList);

datapath  = [pwd '/analysis/'];
figuredir = './figures/';

PadNames = {'Spin', 'Jump', 'Slide'};

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
    cdays = unique(cDk);
    cndays = length(cdays);
    
    caccuracyday = zeros(cndays, size(cdata.pad.accuracy.run.values, 1) - 1);
    for dId = 1:cndays
        cvaluesday = cdata.pad.accuracy.run.values(1:3, cDk == cdays(dId));
        cvaluesday(isnan(cvaluesday)) = 0;
        caccuracyday(dId, :) = mean(cvaluesday, 2);
    end
    
    caccuracyrun = zeros(cnruns, 1);
    for rId = 1:cnruns
       cvaluesrun =  cdata.pad.accuracy.run.values(1:3, rId);
       cvaluesrun(isnan(cvaluesrun)) = 0;
       caccuracyrun(rId) = mean(cvaluesrun);
    end
    
    [ccorr, cpval] = corr(cruns, caccuracyrun, 'rows', 'pairwise');
    cnbiutil_bdisp(['Correlation (accuracy vs. run index) for ' csubject ':']);
    disp(['r= ' num2str(ccorr) ', pval= ' num2str(cpval, '%3.3e') ' - N=' num2str(max(cruns))]);
    
    
    
    subplot(NumRows, NumCols, sId);
    

    hold on;
    plot(caccuracyday, '.', 'MarkerSize', 20);
    plot(nanmean(caccuracyday, 2), ':g', 'LineWidth', 2);
    hold off;
    
    ylim([0 120]);
    xlim([1 cndays]);
    ytick = get(gca, 'YTick');
    set(gca, 'YTick', intersect(ytick, 1:100));
    xtick = [1; cDk(find(diff(cDk))+1)] ;
    set(gca, 'XTick', xtick);
    set(gca, 'XTickLabel', cdata.pad.label.session.Dl(xtick, :));
    xticklabel_rotate([],90,[])
    grid on;
    
    hxlabel = xlabel('Session');
    set(hxlabel, 'Position', get(hxlabel, 'Position') - [0 0.02 0])
    ylabel('Accuracy [%]');
    title(csubject);
    
    cpos = get(gca,'Position'); 
    cpos(2) = 0.15; 
    textcolor = 'g';
    annotation('textbox', cpos, 'String', ['r=' num2str(ccorr, '%3.2f') ', p=' num2str(cpval, '%3.3f')], 'LineStyle', 'none', 'FontWeight', 'bold')
    legend([PadNames, 'Average'], 'location', 'SouthEast');
    
end

suptitle('Command accuracy');
cnbifig_export(fig1, [figuredir '/cybathlon.journal.accuracy.pad.png'], '-png');
cnbifig_export(fig1, [figuredir '/cybathlon.journal.accuracy.pad.pdf'], '-pdf');