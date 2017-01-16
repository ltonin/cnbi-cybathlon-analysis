clearvars; clc; 

%subject = 'AN14VE';
subject = 'MA25VE';

pattern     = '.mi.';
modality    = 'race';

experiment  = 'cybathlon';
datapath    = [pwd '/analysis/'];
figuredir  = './figures/';

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

%% Extract pad events
cnbiutil_bdisp('[proc] - Extract events');
[TrialLb, TrialEvents] = cnbiproc_get_event(PadTypeId, DataLength, events.POS, events.TYP, events.DUR);

%% Extract command events
cnbiutil_bdisp('[proc] - Extract commands');
[CommLb, CommEvents] = cnbiproc_get_event(CommTypeId, DataLength, events.POS, events.TYP, events.DUR);

%% Compute the overall accuracies (one decision per pad)
[TPFPPad, TPFPTask, SpeedPad, SpeedTask] = cnbiproc_commacc(CommLb, TrialEvents, PadTypeId, PadTypeLb, PadTypeInd);

%% Compute the overall accuracies (one decision per command/trial)
[TPFPPad, TPFPTask, SpeedPad, SpeedTask] = cnbiproc_commacc(CommLb, TrialEvents, PadTypeId, PadTypeLb, PadTypeInd);

%% Compute accuracies per day
SessionTrialEvents = [];
MSpeedPadSes = [];
for dId = 1:NumDays
    KeepInd = [];
    for tr=1:length(TrialEvents.DUR)
        if(Dk(TrialEvents.POS(tr))==dId)
            KeepInd = [KeepInd; tr];
        end
    end
    thisSessionTrialEvents.POS = TrialEvents.POS(KeepInd);
    thisSessionTrialEvents.TYP = TrialEvents.TYP(KeepInd);
    thisSessionTrialEvents.DUR = TrialEvents.DUR(KeepInd);
    [TPFPPadSes{dId}, TPFPTaskSes{dId}, SpeedPadSes, SpeedTaskSes] = cnbiproc_commacc(CommLb, thisSessionTrialEvents, PadTypeId, PadTypeLb, PadTypeInd);
    
    for t=1:size(TPFPPad,1)
        AccPad(t,dId) = TPFPPadSes{dId}(t,1);
    end
    for t=1:size(TPFPTask,1)
        AccTask(t,dId) = TPFPTaskSes{dId}(t,1);
    end
    
    MSpeedPadSes(dId,1) = mean(SpeedPadSes{1});
    SSpeedPadSes(dId,1) = std2(SpeedPadSes{1});
    MSpeedPadSes(dId,2) = mean(SpeedPadSes{2});
    SSpeedPadSes(dId,2) = std2(SpeedPadSes{2});
    MSpeedPadSes(dId,3) = mean(SpeedPadSes{3});
    SSpeedPadSes(dId,3) = std2(SpeedPadSes{3});
end


%% Plotting
fig1 = figure;
cnbifig_set_position(fig1, 'All');
plot(1:NumDays,AccPad(1,:),'c',1:NumDays,AccPad(2,:),'m',1:NumDays,AccPad(3,:),'y',1:NumDays,AccPad(4,:),'k',1:NumDays,nanmean(AccPad(1:3,:)),'--g','LineWidth',3);
legend({'Speed','Jump','Slide','Rest','Average (active)'});
xlabel('Race Session','FontSize',20,'LineWidth',3);
ylabel('Command Accuracy (sec)','FontSize',20,'LineWidth',3);
title(subject);
axis([0 NumDays+1 0 109]);
set(gca,'FontSize',20,'LineWidth',3);
set(gca,'XTick',unique(Dk));
set(gca,'XTickLabel',Dl);
xticklabel_rotate([],45,[])

cnbifig_export(fig1, [figuredir '/' subject '.commacc.' modality '.png'], '-png');


fig2 = figure;
cnbifig_set_position(fig2, 'All');
plot(1:NumDays,MSpeedPadSes(:,1),'c',1:NumDays,MSpeedPadSes(:,2),'m',1:NumDays,MSpeedPadSes(:,3),'y','LineWidth',3, 'MarkerSize',10,'LineWidth',3);
legend({'Speed','Jump','Slide'});
hold on;
shadedErrorBar(1:NumDays,MSpeedPadSes(:,1),SSpeedPadSes(:,1),'*c',1);
shadedErrorBar(1:NumDays,MSpeedPadSes(:,2),SSpeedPadSes(:,2),'*m',1);
shadedErrorBar(1:NumDays,MSpeedPadSes(:,3),SSpeedPadSes(:,3),'*y',1);
hold off;
xlabel('Race Session','FontSize',20,'LineWidth',3);
ylabel('Average Command Delivery Time (sec)','FontSize',20,'LineWidth',3);
title(subject);
axis([1 max(Dk) 0 11]);
set(gca,'FontSize',20,'LineWidth',3);
set(gca,'XTick',unique(Dk));
set(gca,'XTickLabel',Dl);
xticklabel_rotate([],45,[])
cnbifig_export(fig2, [figuredir '/' subject '.commtime.' modality '.png'], '-png');