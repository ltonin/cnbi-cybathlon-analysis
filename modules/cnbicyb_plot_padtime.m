clearvars; clc; close all;

SubList = {'AN14VE', 'MA25VE'};
NumSubjects = length(SubList);

datapath  = [pwd '/analysis/'];
figuredir = './figures/';

MinOnPad = 1.9;
MaxOnPad = 19;
NoCmdOnPad = 11;
MinOnRest = 5.5;


ax = gca;
Colors = ax.ColorOrder;
close all;
% Speed, Jump, Slide, Rest
PadColors = {Colors(6, :), Colors(4, :), Colors(3, :), 'k', 'k', 'k'};

TP = [];
TPk    = [];
Sk      = [];
for sId = 1:NumSubjects
    csubject = SubList{sId};
    cfilepath = [datapath '/' csubject '.timeonpad.race.mat'];
    
    cdata = load(cfilepath);
    
    TP  = cat(1, TP, cdata.timeonpad.values');
    TPk = cat(1, TPk, cdata.timeonpad.label.pad');
    Sk  = cat(1, Sk, sId*ones(length(cdata.timeonpad.values), 1));
end


LSk = zeros(length(Sk), length(SubList{1}));
for sId = 1:NumSubjects
    LSk(Sk == sId, :) = repmat(SubList{sId}, sum(Sk == sId), 1);
end

fig1 = figure;
cnbifig_set_position(fig1, 'Top');

SelectedIndex = TPk <= 4;   % Speed, Jump, Slide, Rest

boxplot(TP(SelectedIndex), {TPk(SelectedIndex) Sk(SelectedIndex)}, 'FactorGap', 'auto', 'FactorSeparator', 1, 'labels', char(LSk(SelectedIndex, :)));
h = findobj(gca, 'Tag', 'Box');
h = flipud(h);

for pId = 1:length(unique(TPk(SelectedIndex)))
    set(h(2*pId-1:2*pId), 'Color', PadColors{pId}, 'LineWidth', 2);
end

colorgrey = [0.3 0.3 0.3];
dotstyle = ':';
line([0 11],[MaxOnPad MaxOnPad], 'LineStyle', dotstyle, 'Color', colorgrey);
text(0.6,MaxOnPad + 0.5,'MaxTime')
text(6.9,MaxOnPad + 0.5,'MaxTime')

line([0 6.85],[NoCmdOnPad NoCmdOnPad], 'LineStyle', dotstyle, 'Color', colorgrey);
text(0.6,NoCmdOnPad + 0.5,'NoDelivery')

line([0 6.85],[MinOnPad MinOnPad], 'LineStyle', dotstyle, 'Color', colorgrey);
text(0.6,MinOnPad + 0.5,'MinTime')

line([6.85 11],[MinOnRest MinOnRest], 'LineStyle', dotstyle, 'Color', colorgrey);
text(6.9,MinOnRest + 0.5,'MinTime')

legend(h(1:2:end), 'Speed', 'Jump', 'Slide', 'Rest')
ylim([1 20]);
ylabel('Time [s]');
xlabel('Pad/Pilot');
title('Action pads');


suptitle('Time distribution on pads');

cnbifig_export(fig1, [figuredir '/cybathlon.journal.padtime.png'], '-png');