% clearvars; clc;
% 
% subject = 'AN14VE';
% % subject = 'MA25VE';

pattern     = '.mi.';
modality    = 'online';

experiment  = 'cybathlon';
datapath    = [pwd '/analysis/'];
figuredir   = './figures/';
savedir     = [pwd '/analysis/'];

CueTypeId = [769 770 771 773 783];
CueTypeLb = {'LeftHand', 'RightHand', 'BothFeet', 'BothHands', 'Rest'};

HitMissTypeId = [897 898 899];
HitMissTypeLb = {'Hit','Miss', 'Timeout'};

%% Get datafiles
[Files, NumFiles] = cnbiutil_getfile(datapath, '.mat', [subject '*' modality '*' pattern]);

%% Concatenate data
cnbiutil_bdisp(['[io] - Import psd datafiles (' modality ')']);
[U, events, labels, settings] = cnbiutil_concatenate_data(Files);
DataLength  = size(U, 1);

Dk = labels.Dk;
Dl = labels.Dl;
Days    = unique(Dk);
NumDays = length(Days);

%% Extract events
cnbiutil_bdisp('[proc] - Extract events');
[CueLb, CueEvents] = cnbiproc_get_event(CueTypeId, DataLength, events.POS, events.TYP, events.DUR);
[HitMissLb, HitMissEvents] = cnbiproc_get_event(HitMissTypeId, DataLength, events.POS, events.TYP, events.DUR);

%% Compute overall command accuracy
for c=1:length(CueTypeId)
    CommAcc(c) = 100*sum(CueEvents.TYP==CueTypeId(c) & HitMissEvents.TYP==897)/sum(CueEvents.TYP==CueTypeId(c));
end

%% Compute discriminancy per day

for dId = 1:NumDays
    
    SesTrialInd = find(Dk(CueEvents.POS)==dId);
    for c=1:length(CueTypeId)
        CommAccSes(c,dId) = 100*sum(CueEvents.TYP(SesTrialInd)==CueTypeId(c) & HitMissEvents.TYP(SesTrialInd)==897)/sum(CueEvents.TYP(SesTrialInd)==CueTypeId(c));
    end
    
end

%% Plotting discriminancy map
fig1 = figure;
cnbifig_set_position(fig1, 'All');
plot(1:NumDays,CommAccSes(1,:),'y',1:NumDays,CommAccSes(2,:),'b',1:NumDays,CommAccSes(3,:),'m',1:NumDays,CommAccSes(4,:),'c',1:NumDays,CommAccSes(5,:),'k',1:NumDays,nanmean(CommAccSes,1),'--g','LineWidth',3);
legend({'Left Hand','Right Hand','Both Feet','Both Hands','Rest','Average'});
xlabel('Feedback Session','FontSize',20,'LineWidth',3);
ylabel('Command Accuracy (%)','FontSize',20,'LineWidth',3);
title(subject);
axis([0 NumDays+1 0 109]);
set(gca,'FontSize',20,'LineWidth',3);
set(gca,'XTick',unique(Dk));
set(gca,'XTickLabel',Dl);
xticklabel_rotate([],45,[])

cnbifig_export(fig1, [figuredir '/' subject '.commacc.' modality '.png'], '-png');

%% Saving metadata

% Grouping results
command.accuracy  = CommAccSes;
command.label     = Dl;

savefile = [savedir '/' subject '.metadata.mat'];
if exist(savefile, 'file')
    cnbiutil_bdisp(['Loading metadata from: ' savefile]);
    load(savefile);
end

metadata.online.command = command;

cnbiutil_bdisp(['Saving command accuracy (online) results in: ' savefile]);
save(savefile, 'metadata');
