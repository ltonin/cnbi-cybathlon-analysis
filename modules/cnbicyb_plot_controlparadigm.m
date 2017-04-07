clearvars; clc; close all;

SubList = {'AN14VE'};
NumSubjects = length(SubList);

datapath  = [pwd '/analysis/'];
figuredir = './figures/';

ParadigmNames = {'Control paradigm 3', 'Control paradigm 4'};
PadNames = {'Spin', 'Jump', 'Slide'};

ax = gca;
Colors = ax.ColorOrder;
close all;
% Speed, Jump, Slide, Rest
PadColors = {Colors(6, :), Colors(4, :), Colors(3, :), 'k', 'k', 'k'};
ParadigmColors = {Colors(5, :), Colors(6, :)};

fig1 = figure;
cnbifig_set_position(fig1, 'Top');

for sId = 1:NumSubjects
    
    csubject = SubList{sId};
    cfilepath = [datapath '/' csubject '.padaccuracy.race.mat'];
    cdatapadacc = load(cfilepath);
    
    cfilepath = [datapath '/' csubject '.time.race.mat'];
    cdatatimerace = load(cfilepath);
    
    
    % Race time, paradigm 3 vs 4
    MeanVecRT = [mean(cdatatimerace.time.run.values(cdatatimerace.time.run.label.Pk==5)) ; mean(cdatatimerace.time.run.values(cdatatimerace.time.run.label.Pk==6))];
    StdVecRT = [std(cdatatimerace.time.run.values(cdatatimerace.time.run.label.Pk==5)) ; std(cdatatimerace.time.run.values(cdatatimerace.time.run.label.Pk==6))];
    NVecRT = [size(cdatatimerace.time.run.values(cdatatimerace.time.run.label.Pk==5),2) ; size(cdatatimerace.time.run.values(cdatatimerace.time.run.label.Pk==6),2)];
    pValRT = ranksum(cdatatimerace.time.run.values(cdatatimerace.time.run.label.Pk==5), cdatatimerace.time.run.values(cdatatimerace.time.run.label.Pk==6));    
    
    % Pad accuracy of slide, paradigm 3 vs 4
    MeanVecPA = [mean(cdatapadacc.pad.accuracy.run.values(1:3,cdatapadacc.pad.label.run.Pk==5),2)' ; mean(cdatapadacc.pad.accuracy.run.values(1:3,cdatapadacc.pad.label.run.Pk==6),2)'];
    StdVecPA = [std(cdatapadacc.pad.accuracy.run.values(1:3,cdatapadacc.pad.label.run.Pk==5)')' std(cdatapadacc.pad.accuracy.run.values(1:3,cdatapadacc.pad.label.run.Pk==6)')']';
    NVecPA = [size(cdatapadacc.pad.accuracy.run.values(1:3,cdatapadacc.pad.label.run.Pk==5),2) ; size(cdatapadacc.pad.accuracy.run.values(1:3,cdatapadacc.pad.label.run.Pk==6),2)];
    for t=1:3
        pValPA(t) = ranksum(cdatapadacc.pad.accuracy.run.values(t,cdatapadacc.pad.label.run.Pk==5), cdatapadacc.pad.accuracy.run.values(t,cdatapadacc.pad.label.run.Pk==6));
    end
    
    offset = 0.2;
    subplot(1,2,1);hold on;
    cellsigstart = {};
    
    hBar{1} = bar([1:3]-offset,MeanVecPA(1,:),0.25,'LineWidth',2,'FaceColor',ParadigmColors{1});
    hBar{2} = bar([1:3]+offset,MeanVecPA(2,:),0.25,'LineWidth',2,'FaceColor',ParadigmColors{2});
    hEBar{1} = errorbar([1:3]-offset,MeanVecPA(1,:),zeros(1,3),StdVecPA(1,:),'.','LineWidth',2,'Color',ParadigmColors{1});
    hEBar{2} = errorbar([1:3]+offset,MeanVecPA(2,:),zeros(1,3),StdVecPA(2,:),'.','LineWidth',2,'Color',ParadigmColors{2});
    
    for b=1:3    
        cellsigstart{b} = [b-offset,b+offset];
    end
    H = sigstar(cellsigstart,pValPA); 
    hold off;
    
    ylim([0 120]);
    xlim([0.5 3.5]);
    ytick = get(gca, 'YTick');
    set(gca, 'YTick', intersect(ytick, 1:100));
    xtick = [1:3] ;
    set(gca, 'XTick', xtick);
    set(gca, 'XTickLabel', PadNames);
    %xticklabel_rotate([],90,[])
    grid on;
    
    legend(ParadigmNames, 'location', 'SouthEast');
    hxlabel = xlabel('Command Type');
    set(hxlabel, 'Position', get(hxlabel, 'Position') - [0 0.02 0])
    ylabel('Accuracy [%]');
    title('Command Accuracy');

    
    subplot(1,2,2);hold on;
    cellsigstart = {};
    
    hBar{1} = bar(1-offset,MeanVecRT(1),0.25,'LineWidth',2,'FaceColor',ParadigmColors{1});
    hBar{2} = bar(1+offset,MeanVecRT(2),0.25,'LineWidth',2,'FaceColor',ParadigmColors{2});
    hEBar{1} = errorbar(1-offset,MeanVecRT(1),0,StdVecRT(1),'.','LineWidth',2,'Color',ParadigmColors{1});
    hEBar{2} = errorbar(1+offset,MeanVecRT(2),0,StdVecRT(2),'.','LineWidth',2,'Color',ParadigmColors{2});
    
    cellsigstart{1} = [1-offset,1+offset];
    
    H = sigstar(cellsigstart,pValRT); 
    hold off;
    
    ylim([80 155]);
    xlim([0.5 1.5]);
    ytick = get(gca, 'YTick');
    set(gca, 'YTick', intersect(ytick, 1:240));
    %xtick = 1 ;
    %set(gca, 'XTick', xtick);
    set(gca, 'XTickLabel', {});
    %xticklabel_rotate([],90,[])
    grid on;
    
    legend(ParadigmNames, 'location', 'SouthEast');
    %hxlabel = xlabel('Command Type');
    %set(hxlabel, 'Position', get(hxlabel, 'Position') - [0 0.02 0])
    ylabel('Time [s]');
    title('Race Time');    

end

suptitle('P1 control paradigm effects');
cnbifig_export(fig1, [figuredir '/cybathlon.journal.controlparadigm.png'], '-png');
cnbifig_export(fig1, [figuredir '/cybathlon.journal.controlparadigm.pdf'], '-pdf');