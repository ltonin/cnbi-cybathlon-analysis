clearvars; clc; close all;

SubList = {'AN14VE'};
NumSubjects = length(SubList);

datapath  = [pwd '/analysis/'];
figuredir = './figures/';

ParadigmNames = {'Control paradigm 1', 'Control paradigm 3', 'Control paradigm 4'};
PadNames = {'Spin', 'Jump', 'Slide'};

ax = gca;
Colors = ax.ColorOrder;
close all;
% Speed, Jump, Slide, Rest
PadColors = {Colors(6, :), Colors(4, :), Colors(3, :), 'k', 'k', 'k'};
ParadigmColors = {Colors(3, :) Colors(5, :), Colors(6, :)};

fig1 = figure;
cnbifig_set_position(fig1, 'All');

[~, plotdatapath] = cnbiutil_mkdir(pwd, '/plotdata');

for sId = 1:NumSubjects
    
    csubject = SubList{sId};
    cfilepath = [datapath '/' csubject '.padaccuracy.race.mat'];
    cdatapadacc = load(cfilepath);
    
    cfilepath = [datapath '/' csubject '.time.race.mat'];
    cdatatimerace = load(cfilepath);
    
    cfilepath = [datapath '/' csubject '.timeonpad.race.mat'];
    cdatatop = load(cfilepath);    
    
    % Race time, paradigm 1 vs 3 vs 4
    MeanVecRT = [mean(cdatatimerace.time.run.values(cdatatimerace.time.run.label.Pk==3)) ; mean(cdatatimerace.time.run.values(cdatatimerace.time.run.label.Pk==5)) ; mean(cdatatimerace.time.run.values(cdatatimerace.time.run.label.Pk==6))];
    StdVecRT = [std(cdatatimerace.time.run.values(cdatatimerace.time.run.label.Pk==3)) ; std(cdatatimerace.time.run.values(cdatatimerace.time.run.label.Pk==5)) ; std(cdatatimerace.time.run.values(cdatatimerace.time.run.label.Pk==6))];
    NVecRT = [size(cdatatimerace.time.run.values(cdatatimerace.time.run.label.Pk==3),1) ; size(cdatatimerace.time.run.values(cdatatimerace.time.run.label.Pk==5),1) ; size(cdatatimerace.time.run.values(cdatatimerace.time.run.label.Pk==6),1)];
    pValRT3_4 = ranksum(cdatatimerace.time.run.values(cdatatimerace.time.run.label.Pk==5), cdatatimerace.time.run.values(cdatatimerace.time.run.label.Pk==6));    
    pValRT1_3 = ranksum(cdatatimerace.time.run.values(cdatatimerace.time.run.label.Pk==3), cdatatimerace.time.run.values(cdatatimerace.time.run.label.Pk==5));
    pValRT1_4 = ranksum(cdatatimerace.time.run.values(cdatatimerace.time.run.label.Pk==3), cdatatimerace.time.run.values(cdatatimerace.time.run.label.Pk==6));
    

    % Pad accuracy of slide, paradigm 1 vs 3 vs 4
    MeanVecPA = [nanmean(cdatapadacc.pad.accuracy.run.values(1:3,cdatapadacc.pad.label.run.Pk==3),2)' ; mean(cdatapadacc.pad.accuracy.run.values(1:3,cdatapadacc.pad.label.run.Pk==5),2)' ; mean(cdatapadacc.pad.accuracy.run.values(1:3,cdatapadacc.pad.label.run.Pk==6),2)'];
    StdVecPA = [nanstd(cdatapadacc.pad.accuracy.run.values(1:3,cdatapadacc.pad.label.run.Pk==3)')' std(cdatapadacc.pad.accuracy.run.values(1:3,cdatapadacc.pad.label.run.Pk==5)')' std(cdatapadacc.pad.accuracy.run.values(1:3,cdatapadacc.pad.label.run.Pk==6)')']';
    NVecPA = [size(cdatapadacc.pad.accuracy.run.values(1:3,cdatapadacc.pad.label.run.Pk==3),2) ; size(cdatapadacc.pad.accuracy.run.values(1:3,cdatapadacc.pad.label.run.Pk==5),2) ; size(cdatapadacc.pad.accuracy.run.values(1:3,cdatapadacc.pad.label.run.Pk==6),2)];
    for t=1:3
        pValPA3_4(t) = ranksum(cdatapadacc.pad.accuracy.run.values(t,cdatapadacc.pad.label.run.Pk==5), cdatapadacc.pad.accuracy.run.values(t,cdatapadacc.pad.label.run.Pk==6));
    end
    
    for t=1:3
        pValPA1_3(t) = ranksum(cdatapadacc.pad.accuracy.run.values(t,cdatapadacc.pad.label.run.Pk==3), cdatapadacc.pad.accuracy.run.values(t,cdatapadacc.pad.label.run.Pk==5));
    end
    
    for t=1:3
        pValPA1_4(t) = ranksum(cdatapadacc.pad.accuracy.run.values(t,cdatapadacc.pad.label.run.Pk==3), cdatapadacc.pad.accuracy.run.values(t,cdatapadacc.pad.label.run.Pk==6));
    end
    
    % Time on pad, paradigm 1 vs 3 vs 4
    IndProt1 = [(find(cdatapadacc.pad.label.run.Pk==3)-1)*18+1 find(cdatapadacc.pad.label.run.Pk==3)*18];
    IndFProt1 = [];
    for i=1:size(IndProt1,1)
        IndFProt1 = [IndFProt1 IndProt1(i,1):IndProt1(i,2)];
    end
    
    IndProt3 = [(find(cdatapadacc.pad.label.run.Pk==5)-1)*18+1 find(cdatapadacc.pad.label.run.Pk==5)*18];
    IndFProt3 = [];
    for i=1:size(IndProt3,1)
        IndFProt3 = [IndFProt3 IndProt3(i,1):IndProt3(i,2)];
    end
    IndProt4 = [(find(cdatapadacc.pad.label.run.Pk==6)-1)*18+1 find(cdatapadacc.pad.label.run.Pk==6)*18];
    IndFProt4 = [];
    for i=1:size(IndProt4,1)
        IndFProt4 = [IndFProt4 IndProt4(i,1):IndProt4(i,2)];
    end    
    
    
    bpdata = [];
    bplbl = [];
    xpos = [];
    offset = 0.2;
    for t=1:3
        %MeanVecTOP(1,t) = median(cdatatop.timeonpad.topall.values(intersect(find(cdatatop.timeonpad.topall.labels==t),IndFProt3)));
        %MeanVecTOP(1,t) = median(cdatatop.timeonpad.topall.values(intersect(find(cdatatop.timeonpad.topall.labels==t),IndFProt3)));
        bpdata = [ bpdata ; cdatatop.timeonpad.topall.values(intersect(find(cdatatop.timeonpad.topall.labels==t),IndFProt1))];
        bplbl = [ bplbl ; (3*(t-1)+1)*ones(length(cdatatop.timeonpad.topall.values(intersect(find(cdatatop.timeonpad.topall.labels==t),IndFProt1))),1)];
        xpos = [xpos t-offset];
        bpdata = [ bpdata ; cdatatop.timeonpad.topall.values(intersect(find(cdatatop.timeonpad.topall.labels==t),IndFProt3))];
        bplbl = [ bplbl ; (3*(t-1)+2)*ones(length(cdatatop.timeonpad.topall.values(intersect(find(cdatatop.timeonpad.topall.labels==t),IndFProt3))),1)];
        xpos = [xpos t];
        bpdata = [ bpdata ; cdatatop.timeonpad.topall.values(intersect(find(cdatatop.timeonpad.topall.labels==t),IndFProt4))];
        bplbl = [ bplbl ; (3*(t-1)+3)*ones(length(cdatatop.timeonpad.topall.values(intersect(find(cdatatop.timeonpad.topall.labels==t),IndFProt4))),1)];        
        xpos = [xpos t+offset];
        pValTOP3_4(t) = ranksum(cdatatop.timeonpad.topall.values(intersect(find(cdatatop.timeonpad.topall.labels==t),IndFProt3)),cdatatop.timeonpad.topall.values(intersect(find(cdatatop.timeonpad.topall.labels==t),IndFProt4)));
        pValTOP1_3(t) = ranksum(cdatatop.timeonpad.topall.values(intersect(find(cdatatop.timeonpad.topall.labels==t),IndFProt1)),cdatatop.timeonpad.topall.values(intersect(find(cdatatop.timeonpad.topall.labels==t),IndFProt3)));
        pValTOP1_4(t) = ranksum(cdatatop.timeonpad.topall.values(intersect(find(cdatatop.timeonpad.topall.labels==t),IndFProt1)),cdatatop.timeonpad.topall.values(intersect(find(cdatatop.timeonpad.topall.labels==t),IndFProt4)));
    end

    
    % Time on pad per paradigm
    subplot(2,2,[1 2]);hold on;
    hbp = boxplot(bpdata,bplbl,'Positions',xpos,'Colors',[ParadigmColors{1} ; ParadigmColors{2} ; ParadigmColors{3}]);
    set(hbp,{'linew'},{2});
%     tmppos = cell2mat(get(hbp(1,:),'XData')); % positions of whiskers
%     tmppos = tmppos(:,1)';
%     cellsigstart = {tmppos([1,2]),tmppos([3,4]),tmppos([5,6])};
%     H = sigstar(cellsigstart,pValTOP); 
    
    hold off;
    
    ylim([0 22]);
%     xlim([tmppos(1)-0.3 tmppos(end)+0.3]);
    ytick = get(gca, 'YTick');
    set(gca, 'YTick', intersect(ytick, 1:22));
    %xticklabel_rotate([],90,[])
    grid on;
    
    hl = findall(gca,'Tag','Box');
    hl = flipud(hl);
    hLegend = legend(hl, {'Control paradigm 1', 'Control paradigm 3', 'Control paradigm 4'});
    xtick = [1:3] ;
    set(gca, 'XTick', xtick);
    set(gca, 'XTickLabel', PadNames);    
    %hxlabel = xlabel('Command Type');
    %set(hxlabel, 'Position', get(hxlabel, 'Position') - [0 0.02 0])
    ylabel('Time [s]');
    title('Time On Pad');    

    % Storing plot data
    data.time_on_pad          = bpdata;
    data.labels.group_id      = bplbl;
    data.labels.paradigm_name = {'Control paradigm 1', 'Control paradigm 3', 'Control paradigm 4'};
    data.labels.pad_name      = PadNames;
    
    
    % BCI command accuracy per paradigm
    offset = 0.3;
    subplot(2,2,3);hold on;
%     cellsigstart = {};
    
    hBar{1} = bar([1:3]-offset,MeanVecPA(1,:),0.25,'LineWidth',2,'FaceColor',ParadigmColors{1});
    hBar{2} = bar([1:3],MeanVecPA(2,:),0.25,'LineWidth',2,'FaceColor',ParadigmColors{2});
    hBar{3} = bar([1:3]+offset,MeanVecPA(3,:),0.25,'LineWidth',2,'FaceColor',ParadigmColors{3});
    hEBar{1} = errorbar([1:3]-offset,MeanVecPA(1,:),zeros(1,3),StdVecPA(1,:),'.','LineWidth',2,'Color',ParadigmColors{1});
    hEBar{2} = errorbar([1:3],MeanVecPA(2,:),zeros(1,3),StdVecPA(2,:),'.','LineWidth',2,'Color',ParadigmColors{2});
    hEBar{3} = errorbar([1:3]+offset,MeanVecPA(3,:),zeros(1,3),StdVecPA(3,:),'.','LineWidth',2,'Color',ParadigmColors{3});
    
%     for b=1:3    
%         cellsigstart{b} = [b-offset,b+offset];
%     end
%     H = sigstar(cellsigstart,pValPA); 
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
    %hxlabel = xlabel('Command Type');
    %set(hxlabel, 'Position', get(hxlabel, 'Position') - [0 0.02 0])
    ylabel('Accuracy [%]');
    title('Command Accuracy');

    % Storing plot data
    data2.accuracy_on_pad      = cdatapadacc.pad.accuracy.run.values(1:3,:);
    data2.labels.paradigm_id   = cdatapadacc.pad.label.run.Pk;
    data2.labels.paradigm_sel  = [3 5 6];
    data2.labels.paradigm_name = {'Control paradigm 1', 'Control paradigm 3', 'Control paradigm 4'};
    data2.labels.pad_name      = PadNames;
    
    % Race time per paradigm
    subplot(2,2,4);hold on;
%     cellsigstart = {};
    
    hBar{1} = bar(1-offset,MeanVecRT(1),0.25,'LineWidth',2,'FaceColor',ParadigmColors{1});
    hBar{2} = bar(1,MeanVecRT(2),0.25,'LineWidth',2,'FaceColor',ParadigmColors{2});
    hBar{3} = bar(1+offset,MeanVecRT(3),0.25,'LineWidth',2,'FaceColor',ParadigmColors{3});
    hEBar{1} = errorbar(1-offset,MeanVecRT(1),0,StdVecRT(1),'.','LineWidth',2,'Color', ParadigmColors{1});
    hEBar{2} = errorbar(1,MeanVecRT(2),0,StdVecRT(2),'.','LineWidth',2,'Color', ParadigmColors{2});
    hEBar{3} = errorbar(1+offset,MeanVecRT(3),0,StdVecRT(3),'.','LineWidth',2,'Color', ParadigmColors{3});
    
%     cellsigstart{1} = [1-offset,1+offset];
    
%     H = sigstar(cellsigstart,pValRT);
    hold off;
    
    ylim([110 180]);
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
    
    % Storing plot data
    data3.accuracy_on_pad      = cdatatimerace.time.run.values;
    data3.labels.paradigm_id   = cdatatimerace.time.run.label.Pk;
    data3.labels.paradigm_sel  = [3 5 6];
    data3.labels.paradigm_name = {'Control paradigm 1', 'Control paradigm 3', 'Control paradigm 4'};
    data3.labels.pad_name      = PadNames;
    
    % Statistical tests
    disp('Time on pad - Spin:');
    disp(['           - Paradigm 1 vs 3: p=' num2str(pValTOP1_3(1))]);
    disp(['           - Paradigm 1 vs 4: p=' num2str(pValTOP1_4(1))]);
    disp(['           - Paradigm 3 vs 4: p=' num2str(pValTOP3_4(1))]);
    
    disp('Time on pad - Jump:');
    disp(['           - Paradigm 1 vs 3: p=' num2str(pValTOP1_3(2))]);
    disp(['           - Paradigm 1 vs 4: p=' num2str(pValTOP1_4(2))]);
    disp(['           - Paradigm 3 vs 4: p=' num2str(pValTOP3_4(2))]);
    
    disp('Time on pad - Slide:');
    disp(['           - Paradigm 1 vs 3: p=' num2str(pValTOP1_3(3))]);
    disp(['           - Paradigm 1 vs 4: p=' num2str(pValTOP1_4(3))]);
    disp(['           - Paradigm 3 vs 4: p=' num2str(pValTOP3_4(3))]);
    
    disp('Accuracy on pad - Spin:');
    disp(['           - Paradigm 1 vs 3: p=' num2str(pValPA1_3(1))]);
    disp(['           - Paradigm 1 vs 4: p=' num2str(pValPA1_4(1))]);
    disp(['           - Paradigm 3 vs 4: p=' num2str(pValPA3_4(1))]);
    
    disp('Accuracy on pad - Jump:');
    disp(['           - Paradigm 1 vs 3: p=' num2str(pValPA1_3(2))]);
    disp(['           - Paradigm 1 vs 4: p=' num2str(pValPA1_4(2))]);
    disp(['           - Paradigm 3 vs 4: p=' num2str(pValPA3_4(2))]);
    
    disp('Accuracy on pad - Slide:');
    disp(['           - Paradigm 1 vs 3: p=' num2str(pValPA1_3(3))]);
    disp(['           - Paradigm 1 vs 4: p=' num2str(pValPA1_4(3))]);
    disp(['           - Paradigm 3 vs 4: p=' num2str(pValPA3_4(3))]);
    
    disp('Race Time:');
    disp(['      - Paradigm 1 vs 3: p=' num2str(pValRT1_3)]);
    disp(['      - Paradigm 1 vs 4: p=' num2str(pValRT1_4)]);    
    disp(['      - Paradigm 3 vs 4: p=' num2str(pValRT3_4)]);

end

suptitle('P1 control paradigm effects');

%% Saving plot data
cnbiutil_bdisp(['Saving plot data Fig6A in ' plotdatapath]);
save([plotdatapath '/Fig6A.mat'], 'data');
cnbiutil_bdisp(['Saving plot data Fig6B in ' plotdatapath]);
save([plotdatapath '/Fig6B.mat'], 'data2');
cnbiutil_bdisp(['Saving plot data Fig6C in ' plotdatapath]);
save([plotdatapath '/Fig6C.mat'], 'data3');

cnbifig_export(fig1, [figuredir '/cybathlon.journal.controlparadigm3.png'], '-png');
cnbifig_export(fig1, [figuredir '/cybathlon.journal.controlparadigm3.pdf'], '-pdf');