clearvars; clc; 

subject = 'AN14VE';
%subject = 'MA25VE';

pattern     = '.mi.';
modality    = 'race';

experiment  = 'cybathlon';
datapath    = [pwd '/analysis/'];
figuredir   = './figures/';
savedir     = [pwd '/analysis/'];

PadTypeId = [768 769 770 771 773 783];
PadTypeLb = {'Slide', 'Slide', 'Speed', 'Jump', 'Speed', 'Rest'};
PadTypeInd = [3 3 1 2 1 4];

CommTypeId = PadTypeId + hex2dec('6000');
CommTypeLb = {'Slide', 'Slide', 'Speed', 'Jump', 'Speed', 'Rest'};

%% Get datafiles
[Files, NumFiles] = cnbiutil_getfile(datapath, '.mat', [subject '*' modality '*' pattern]);

%% Concatenate data
cnbiutil_bdisp(['[io] - Import psd datafiles (' modality ')']);
[U, events, labels, settings] = cnbiutil_concatenate_data(Files);

DataLength  = size(U, 1);

Dk = labels.Dk;
Days    = unique(Dk);
NumDays = length(Days);
Dl = labels.Dl;

Mk = labels.Mk;

%% Compute the overall time-on-pad distributions
[TOPPad TOPTask] = cnbiproc_timeonpad(events.extra.pad, PadTypeId, PadTypeLb);

%% Compute time-on-pad per day

for dId = 1:NumDays
    KeepInd = [];
    for tr=1:length(events.extra.pad.TYP)
        if(Dk(events.extra.pad.POS(tr))==dId)
            KeepInd = [KeepInd; tr];
        end
    end
    thisSessionTOP.POS = events.extra.pad.POS(KeepInd);
    thisSessionTOP.TYP = events.extra.pad.TYP(KeepInd);
    thisSessionTOP.DUR = events.extra.pad.DUR(KeepInd);
    [TOPPadSes{dId} TOPTaskSes{dId}] = cnbiproc_timeonpad(thisSessionTOP, PadTypeId, PadTypeLb);
    
    for t=1:6
        MTOPPad(dId,t) = mean(TOPPadSes{dId}{t});
        STOPPad(dId,t) = std2(TOPPadSes{dId}{t});
        MedTOPPad(dId,t) = median(TOPPadSes{dId}{t});
        PLowTOPPad(dId,t) = prctile(TOPPadSes{dId}{t},25);
        PUpTOPPad(dId,t) = prctile(TOPPadSes{dId}{t},75);
    end
end


%% Plotting
fig1 = figure;
cnbifig_set_position(fig1, 'All');
plotpad = [];
plotlbl = [];
for p=1:6
    plotpad = [plotpad TOPPad{p}'];
    plotlbl = [plotlbl p*ones(1,length(TOPPad{p}))];
end
boxplot(plotpad,plotlbl,'labels',{'Speed','Jump','Slide','Rest','Start','End'});
xlabel('Pad Type','FontSize',20,'LineWidth',3);
ylabel('Time On Pad (sec)','FontSize',20,'LineWidth',3);
title(subject);
axis([0 7 0 21]);
set(gca,'FontSize',20,'LineWidth',3);
cnbifig_export(fig1, [figuredir '/' subject '.timeonpad.' modality '.png'], '-png');

% Per session results, mean std
fig2 = figure;
cnbifig_set_position(fig2, 'All');
plot(1:NumDays,MTOPPad(:,1),'c',1:NumDays,MTOPPad(:,2),'m',1:NumDays,MTOPPad(:,3),'y',1:NumDays,MTOPPad(:,4),'k',1:NumDays,MTOPPad(:,5),'r',1:NumDays,MTOPPad(:,6),'b','LineWidth',3);
legend({'Speed','Jump','Slide','Rest','Start','End'});
% hold on;
% shadedErrorBar(1:NumDays,squeeze(MTOPPad(:,1))',squeeze(STOPPad(:,1))','c',1);
% shadedErrorBar(1:NumDays,squeeze(MTOPPad(:,2))',squeeze(MTOPPad(:,2))','m',1);
% shadedErrorBar(1:NumDays,squeeze(MTOPPad(:,3))',squeeze(MTOPPad(:,3))','y',1);
% shadedErrorBar(1:NumDays,squeeze(MTOPPad(:,4))',squeeze(MTOPPad(:,4))','k',1);
% shadedErrorBar(1:NumDays,squeeze(MTOPPad(:,5))',squeeze(MTOPPad(:,5))','r',1);
% shadedErrorBar(1:NumDays,squeeze(MTOPPad(:,6))',squeeze(MTOPPad(:,6))','b',1);
% hold off;
xlabel('Pad Type','FontSize',20,'LineWidth',3);
ylabel('Average Time On Pad (sec)','FontSize',20,'LineWidth',3);
title(subject);
%axis([1 max(Dk) min(min(MTOPPad - STOPPad)) max(max(MTOPPad + STOPPad))]);
axis([1 max(Dk) 0 21]);
set(gca,'FontSize',20,'LineWidth',3);
set(gca,'XTick',unique(Dk));
set(gca,'XTickLabel',Dl);
xticklabel_rotate([],45,[])
cnbifig_export(fig2, [figuredir '/' subject '.racetimesession.' modality '.png'], '-png');

% Per session results, median and percentiles
fig3 = figure;
cnbifig_set_position(fig3, 'All');
plot(1:NumDays,MedTOPPad(:,1),'c',1:NumDays,MedTOPPad(:,2),'m',1:NumDays,MedTOPPad(:,3),'y',1:NumDays,MedTOPPad(:,4),'k',1:NumDays,MedTOPPad(:,5),'r',1:NumDays,MedTOPPad(:,6),'b','LineWidth',3);
legend({'Speed','Jump','Slide','Rest','Start','End'});
% hold on;
% shadedErrorBar(1:NumDays,squeeze(MedTOPPad(:,1))',[PUpTOPPad(:,1)' ; PLowTOPPad(:,1)'],'c',1);
% shadedErrorBar(1:NumDays,squeeze(MedTOPPad(:,2))',[PUpTOPPad(:,2)' ; PLowTOPPad(:,2)'],'m',1);
% shadedErrorBar(1:NumDays,squeeze(MedTOPPad(:,3))',[PUpTOPPad(:,3)' ; PLowTOPPad(:,3)'],'y',1);
% shadedErrorBar(1:NumDays,squeeze(MedTOPPad(:,4))',[PUpTOPPad(:,4)' ; PLowTOPPad(:,4)'],'k',1);
% shadedErrorBar(1:NumDays,squeeze(MedTOPPad(:,5))',[PUpTOPPad(:,5)' ; PLowTOPPad(:,5)'],'r',1);
% shadedErrorBar(1:NumDays,squeeze(MedTOPPad(:,6))',[PUpTOPPad(:,6)' ; PLowTOPPad(:,6)'],'b',1);
% hold off;
xlabel('Pad Type','FontSize',20,'LineWidth',3);
ylabel('Average Time On Pad (sec)','FontSize',20,'LineWidth',3);
title(subject);
%axis([1 max(Dk) min(min(MTOPPad - STOPPad)) max(max(MTOPPad + STOPPad))]);
axis([1 max(Dk) 0 21]);
set(gca,'FontSize',20,'LineWidth',3);
set(gca,'XTick',unique(Dk));
set(gca,'XTickLabel',Dl);
xticklabel_rotate([],45,[])
cnbifig_export(fig3, [figuredir '/' subject '.racetimesessionmedian.' modality '.png'], '-png');

%% Saving metadata

% Grouping results
timeonpad.values    = plotpad;
timeonpad.label.pad = plotlbl;
timeonpad.label.padlb = {'Speed','Jump','Slide','Rest','Start','End'};

savefile = [savedir '/' subject '.timeonpad.' modality '.mat'];

cnbiutil_bdisp(['Saving timings on pad (' modality ') results in: ' savefile]);
save(savefile, 'timeonpad');
