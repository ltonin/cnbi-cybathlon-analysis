clearvars; clc; 

subject = 'AN14VE';
% subject = 'MA25VE';

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
[CommRL, CommRLEvents] = cnbiproc_get_event(CommTypeId, DataLength, events.extra.cmd.POS, events.extra.cmd.TYP, events.extra.cmd.DUR);
[CommAll, CommAllEvents] = cnbiproc_get_event(CommTypeId, DataLength, events.extra.cmdg.POS, events.extra.cmdg.TYP, events.extra.cmdg.DUR);

%% Compute the overall accuracies
%[CommAcc] = cnbiproc_commacc(CommRL, CommAll, TrialEvents, events, PadTypeId, PadTypeLb, PadTypeInd);

% Find all commands
IndComm = find(CommAll~=0);
tmpval = max(TrialLb(IndComm),TrialLb(IndComm-1));
CommGT = zeros(length(TrialLb),1);
CommGT(IndComm) = tmpval;

remIndSample = setdiff(find(CommAll~=0),find(CommGT~=0));
remIndTrial = find(ismember(CommRLEvents.POS,remIndSample));

CommRL(remIndSample)=0;
CommAll(remIndSample)=0;
CommGT(remIndSample)=0;

CommRLEvents.TYP(remIndTrial)=[];
CommRLEvents.POS(remIndTrial)=[];
CommRLEvents.DUR(remIndTrial)=[];
CommAllEvents.TYP(remIndTrial)=[];
CommAllEvents.POS(remIndTrial)=[];
CommAllEvents.DUR(remIndTrial)=[];

CommProt = nan(length(CommAll),1);
for p=1:length(events.extra.protocol.TYP)
    CommProt(events.extra.protocol.POS(p):events.extra.protocol.POS(p)+events.extra.protocol.DUR(p)-1) = events.extra.protocol.TYP(p);
end

% The problem is only with the 768 labels, that can have a number of
% interpretations depending on the protocol and order. Reassign them
RealGT = CommGT;
RealComm = CommRL;

% Case of protocol 5 (mi_cybathlon2), where slide is essentially rest
IndPad786Prot5 = intersect(find(CommProt==5), find(CommGT==768));
RealGT(IndPad786Prot5) = 783;
RealComm(IndPad786Prot5(find(RealComm(IndPad786Prot5)==25344))) = 783+hex2dec('6000');


% Case of protocol 6 (mi_cybathlon3/controller), where reverse + timeout
% sends a command
IndRealComm = find(RealComm~=0);
Prot = [];
GT = [];
PossibleClasses = [771 773];
for tr=1:length(events.extra.pad.TYP)
    Prot = events.extra.protocol.TYP(floor((tr-1)/18)+1);
    BiasedClass = mode(RealComm(IndRealComm(intersect(find(IndRealComm >= events.extra.race.POS(floor((tr-1)/18)+1)),...
            find(IndRealComm < events.extra.race.POS(floor((tr-1)/18)+1) + ...
            events.extra.race.DUR(floor((tr-1)/18)+1))))))-hex2dec('6000');
    BiasClassInd = find(PossibleClasses==BiasedClass);
    GT = events.extra.pad.TYP(tr);
    BeginPad = events.extra.pad.POS(tr);
    if( (Prot==6) && (GT==768) )
        % Find the indices of commands inside there
        CommIndInTrial = intersect(find(IndRealComm >= events.extra.pad.POS(tr)),...
            find(IndRealComm < events.extra.pad.POS(tr) + events.extra.pad.DUR(tr)));
        if(isempty(CommIndInTrial))
            continue;
        end
        
        % Find time since immediately previous command
        TimeFromPrevious = (BeginPad-IndRealComm(CommIndInTrial(1)-1))*0.0625;
        TypePrevious = RealComm(IndRealComm(CommIndInTrial(1)-1));
        PreviousInd = find(PossibleClasses==BiasedClass);
        
        if(TimeFromPrevious <= 3.0)
            % If small amount of time has elapsed since the last command
            % when the pad starts, we can almost safely assume that he was
            % trying to reach the inverse, so that he can slide.
            PredictedInd = 3-PreviousInd;
            PredictedType = PossibleClasses(PredictedInd);
        else
            % If too much time has passed, itis most likely they wanted to
            % deliver their "easy"/biased command first
            PredictedInd = BiasClassInd;
            PredictedType = PossibleClasses(BiasClassInd);
        end
        
        for c=1:length(CommIndInTrial)
            if(CommAll(IndRealComm(CommIndInTrial(c))) == 25344)
                % Never mind the prediction, if we are lucky and there is a
                % slide command, we are sure the real command is correct.
                % So, change the 768 label to the one corresponding to this
                % command
                RealGT(IndRealComm(CommIndInTrial(c))) = RealComm(IndRealComm(CommIndInTrial(c)))-hex2dec('6000');
            else
                % Change to the predicted one
                RealGT(IndRealComm(CommIndInTrial(c))) = PredictedType;
                % Update the predicted to the opposite of this one, now we
                % can be pretty sure that he simply tried for the inverse
                IndCurrent = find(PossibleClasses == (RealComm(IndRealComm(CommIndInTrial(c)))-hex2dec('6000')));
                PredictedInd = 3-IndCurrent;
                PredictedType = PossibleClasses(PredictedInd);
            end
        end        
    end
end

% Allrighty, now it is easy to calculate a confusion matrix
ConfMat = cnbiproc_commacc(RealGT, RealComm, PadTypeId, ones(size(Dk)));


%% Compute accuracies per day
CommAcc = [];
for dId = 1:NumDays
    UseInd = (Dk==dId);
    ConfMatSes{dId} = cnbiproc_commacc(RealGT, RealComm, PadTypeId, UseInd);
    % We are mainly interested in 771 and 773
    CommAcc(1,dId) = ConfMatSes{dId}(4,4); % 771, feet
    CommAcc(2,dId) = ConfMatSes{dId}(5,5); % 773, hands
    CommAcc(3,dId) = ConfMatSes{dId}(3,3); % 770, right hand
    CommAcc(4,dId) = ConfMatSes{dId}(2,2); % 769, left hand
    %CommAcc(3,dId) = sum([ ConfMatSes{dId}(3,3) ConfMatSes{dId}(5,5) ]); % 770 or 773
end


%% Plotting
fig1 = figure;
cnbifig_set_position(fig1, 'All');
plot(1:NumDays,CommAcc(1,:),'m',1:NumDays,CommAcc(2,:),'c',1:NumDays,CommAcc(3,:),'b',1:NumDays,CommAcc(4,:),'y',1:NumDays,nanmean(CommAcc,1),'k','LineWidth',3);
legend({'Both Feet','Both Hands','Right Hand','Left Hand','Average'});
xlabel('Race Session','FontSize',20,'LineWidth',3);
ylabel('Command Accuracy (sec)','FontSize',20,'LineWidth',3);
title(subject);
axis([0 NumDays+1 0 109]);
set(gca,'FontSize',20,'LineWidth',3);
set(gca,'XTick',unique(Dk));
set(gca,'XTickLabel',Dl);
xticklabel_rotate([],45,[])

cnbifig_export(fig1, [figuredir '/' subject '.commaccrace.' modality '.png'], '-png');

%% Saving metadata

% Grouping results
command.accuracy  = CommAcc;
command.label     = Dl;

savefile = [savedir '/' subject '.cmdaccuracy.' modality '.mat'];

cnbiutil_bdisp(['Saving command accuracy (' modality ') results in: ' savefile]);
save(savefile, 'command');
