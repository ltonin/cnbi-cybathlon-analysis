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

%% Extract pad events
cnbiutil_bdisp('[proc] - Extract events');
[TrialLb, TrialEvents] = cnbiproc_get_event(PadTypeId, DataLength, events.POS, events.TYP, events.DUR);

%% Extract command events
cnbiutil_bdisp('[proc] - Extract commands');
[CommLb, CommEvents] = cnbiproc_get_event(CommTypeId, DataLength, events.POS, events.TYP, events.DUR);

%% Compute the overall accuracies (one decision per pad)
[TPFPPad, TPFPTask, SpeedPad, SpeedTask, FPPerSec] = cnbiproc_padacc(CommLb, TrialEvents, PadTypeId, PadTypeLb, PadTypeInd);

%% Compute the overall accuracies (one decision per command/trial)
[TPFPPad, TPFPTask, SpeedPad, SpeedTask] = cnbiproc_padacc(CommLb, TrialEvents, PadTypeId, PadTypeLb, PadTypeInd);

%% Compute accuracies per day
SessionTrialEvents = [];
MSpeedPadSes = [];
FPPerSecSes = [];
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
    [TPFPPadSes{dId}, TPFPTaskSes{dId}, SpeedPadSes, SpeedTaskSes, FPPerSecSes(dId)] = cnbiproc_padacc(CommLb, thisSessionTrialEvents, PadTypeId, PadTypeLb, PadTypeInd);
    
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


%% Compute accuracies per race
RaceTrialEvents = [];
MSpeedPadRace = [];
FPPerSecRun = [];
RkSess = zeros(length(TrialEvents.DUR)/18, 1);
DkSess = zeros(length(TrialEvents.DUR)/18, 1);
for rId = 1:(length(TrialEvents.DUR)/18)
    KeepInd = [];
    for tr=1:length(TrialEvents.DUR)
        if(labels.Rk(TrialEvents.POS(tr))==rId)
            KeepInd = [KeepInd; tr];
        end
    end
    thisRaceTrialEvents.POS = TrialEvents.POS(KeepInd);
    thisRaceTrialEvents.TYP = TrialEvents.TYP(KeepInd);
    thisRaceTrialEvents.DUR = TrialEvents.DUR(KeepInd);
    [TPFPPadRace{rId}, TPFPTaskRace{rId}, SpeedPadRace, SpeedTaskRace, FPPerSecRace(rId)] = cnbiproc_padacc(CommLb, thisRaceTrialEvents, PadTypeId, PadTypeLb, PadTypeInd);
    
    for t=1:size(TPFPPad,1)
        AccPadRace(t,rId) = TPFPPadRace{rId}(t,1);
    end
    for t=1:size(TPFPTask,1)
        AccTaskRace(t,dId) = TPFPTaskRace{rId}(t,1);
    end
    
    MSpeedPadRace(rId,1) = mean(SpeedPadRace{1});
    SSpeedPadRace(rId,1) = std2(SpeedPadRace{1});
    MSpeedPadRace(rId,2) = mean(SpeedPadRace{2});
    SSpeedPadRace(rId,2) = std2(SpeedPadRace{2});
    MSpeedPadRace(rId,3) = mean(SpeedPadRace{3});
    SSpeedPadRace(rId,3) = std2(SpeedPadRace{3});
    
    RkSess(rId) = rId;
    DkSess(rId) = unique(labels.Dk(labels.Rk == rId));
    
end

[r pvalPearson] = corr([1:(length(TrialEvents.DUR)/18)]',nanmean(AccPadRace(1:3,:))','type','Pearson');
disp(['Per run, Pearson Correlation r = ' num2str(r) ' with pval = ' num2str(pvalPearson)]);
[rho pvalSperaman] = corr([1:(length(TrialEvents.DUR)/18)]',nanmean(AccPadRace(1:3,:))','type','Spearman');
disp(['Per run, Spearman Correlation r = ' num2str(rho) ' with pval = ' num2str(pvalSperaman)]);

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

cnbifig_export(fig1, [figuredir '/' subject '.padacc.' modality '.png'], '-png');


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
cnbifig_export(fig2, [figuredir '/' subject '.padtime.' modality '.png'], '-png');

fig3 = figure;
cnbifig_set_position(fig3, 'All');
plot(1:NumDays,FPPerSecSes,'k','MarkerSize',10,'LineWidth',3);
xlabel('Race Session','FontSize',20,'LineWidth',3);
ylabel('False Positives per second','FontSize',20,'LineWidth',3);
title(subject);
axis([1 max(Dk) 0 max(FPPerSecSes)+0.1*max(FPPerSecSes)]);
set(gca,'FontSize',20,'LineWidth',3);
set(gca,'XTick',unique(Dk));
set(gca,'XTickLabel',Dl);
xticklabel_rotate([],45,[])
cnbifig_export(fig3, [figuredir '/' subject '.fppersec.' modality '.png'], '-png');


fig4 = figure;
cnbifig_set_position(fig4, 'All');
plot(1:length(labels.Rl),AccPadRace(1,:),'c',1:length(labels.Rl),AccPadRace(2,:),'m',1:length(labels.Rl),AccPadRace(3,:),'y',1:length(labels.Rl),nanmean(AccPadRace(1:3,:)),'--g','MarkerSize',10,'LineWidth',3);
legend({'Speed','Jump','Slide','Average'});
xlabel('Race Index (chronological)','FontSize',20,'LineWidth',3);
ylabel('Command Accuracy (%)','FontSize',20,'LineWidth',3);
title(subject);
axis([1 max(labels.Rk) 0 105]);
set(gca,'FontSize',20,'LineWidth',3);
set(gca,'XTick',unique(labels.Rk));
%set(gca,'XTickLabel',Dl);
%xticklabel_rotate([],45,[])
cnbifig_export(fig4, [figuredir '/' subject '.padaccrace.' modality '.png'], '-png');

%% Saving metadata

% Grouping results
pad.accuracy.run.values     = AccPadRace;
pad.accuracy.run.names      = {'Speed','Jump','Slide','Rest'};
pad.accuracy.session.values = AccPad;
pad.speed.session.values    = MSpeedPadSes;
pad.label.session.Dl        = Dl;
pad.label.run.Rk            = RkSess;
pad.label.run.Dk            = DkSess;

savefile = [savedir '/' subject '.padaccuracy.' modality '.mat'];

cnbiutil_bdisp(['Saving pad accuracy/speed (' modality ') results in: ' savefile]);
save(savefile, 'pad');
