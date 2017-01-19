% clearvars; clc; 
% 
% subject = 'MA25VE';

pattern     = '.mi.';
modality    = 'online';

experiment  = 'cybathlon';
datapath    = [pwd '/analysis/'];
figuredir   = './figures/';
savedir     = [pwd '/analysis/'];

CueTypeId = [769 770 771 773 774 775 783];
CueTypeLb = {'LeftHand', 'RightHand', 'BothFeet', 'BothHands', 'Boh1', 'Boh2', 'Rest'};

CFbTypeId = 781;
CFbTypeLb = {'Continous Feedback'};

FixTypeId = 786;
FixTypeLb = {'Fixation'};

SelectedClassId = [773 771];
SelectedClassLb = {'BothFeet', 'BothHands'};
NumClasses = length(SelectedClassId);



%% Get datafiles
[Files, NumFiles] = cnbiutil_getfile(datapath, '.mat', [subject '*' modality '*' pattern]);

%% Concatenate data
cnbiutil_bdisp(['[io] - Import psd datafiles (' modality ')']);
[U, events, labels, settings, classifiers] = cnbiutil_concatenate_data(Files);
F = log(U);

NumSamples = size(F, 1);
NumFreqs   = size(F, 2);
NumChans   = size(F, 3);

Rk      = labels.Rk;
Runs    = unique(Rk);
NumRuns = length(Runs);

Dk      = labels.Dk;
Days    = unique(Dk);
NumDays = length(Days);

Mk            = labels.Mk;
Modalities    = unique(Mk);
NumModalities = length(Modalities);

Xk             = labels.Xk;
ClassifierIds  = unique(Xk);
NumClassifiers = length(ClassifierIds);

%% Extract events
cnbiutil_bdisp('[proc] - Extract events');
[~, CueEvents] = cnbiproc_get_event(CueTypeId, NumSamples, events.POS, events.TYP, events.DUR);
[~, CFbEvents] = cnbiproc_get_event(CFbTypeId, NumSamples, events.POS, events.TYP, events.DUR);

NumTrials = length(CFbEvents.POS);
  
%% Labeling data
TrialCls = zeros(NumSamples, 1);
TrialRun = zeros(NumSamples, 1);
TrialDay = zeros(NumSamples, 1);
for trId = 1:NumTrials
   cstart = CFbEvents.POS(trId);
   cstop  = cstart + CFbEvents.DUR(trId) - 1;
   cclass = CueEvents.TYP(trId);
   
   TrialCls(cstart:cstop) = cclass;
end
Days    = unique(labels.Dk);
NumDays = length(Days);

%% Selecting classes
Ck = zeros(length(TrialCls), 1);
for cId = 1:NumClasses
    Ck(TrialCls == SelectedClassId(cId)) = SelectedClassId(cId);
end

%% Generic condition
GenericCondition = Ck > 0;

%% Reshaping features
rF = cnbiproc_reshape_ts_bc(F);

%% Computing single sample accuracy for each online

ConfMatr = zeros(2, 3, NumRuns);
for rId = 1:NumRuns
    
    cindex = Rk == Runs(rId) & Ck > 0;
    
    if length(unique(Ck(cindex))) ~= 2
        continue;
    end
    
    % Classifier identification
    cclassifierId = unique(Xk(cindex));
    
    if cclassifierId == 0
        continue;
    end
    
    if(length(cclassifierId) ~= 1)
        error('chk:cls', [num2str(length(cclassifierId)) ' classifiers found in the same run (' num2str(Runs(rId)) ')']);
    end
    cclassifier = classifiers(cclassifierId);
    
    % Feature extraction
    findices  = cnbiproc_features2indices(cclassifier.features, settings.spectrogram.freqgrid);
    cfeatures = rF(cindex, findices);
    reallb    = Ck(cindex);
    
    cprobs = zeros(size(cfeatures, 1), 2);
    for sId = 1:size(cfeatures, 1)
        [~, cprobs(sId, :)] = gauClassifier(cclassifier.gau.M, cclassifier.gau.C, cfeatures(sId, :));
    end
    
    ConfMatr(:, :, rId) = cnbiproc_confusionmat(reallb, cprobs, 'classes', cclassifier.classes, 'rejection', cclassifier.rejection);
end

%% Computing overall accuracy/rejection for each run
RunAccuracy  = nan(NumRuns, 1);
RunRejection = nan(NumRuns, 1);
for rId = 1:NumRuns
    RunAccuracy(rId) = trace(ConfMatr(:, 1:NumClasses, rId))/NumClasses;
    RunRejection(rId) = mean(ConfMatr(:, 3, rId));
end

%% Computing overall accuracy/rejection for each day

DayAccuracy  = nan(NumDays, 2);
DayRejection = nan(NumDays, 2);
for dId = 1:NumDays
    cruns = unique(Rk(Dk == Days(dId)));
    
    DayAccuracy(dId, 1) = mean(RunAccuracy(cruns));
    DayAccuracy(dId, 2) = std(RunAccuracy(cruns));
    
    DayRejection(dId, 1) = mean(RunRejection(cruns));
    DayRejection(dId, 2) = std(RunRejection(cruns));
end
    
%% Plotting

fig1 = figure;
cnbifig_set_position(fig1, 'Top');

NumRows = 1;
NumCols = 2;

% Accuracy/Rejection per run
subplot(NumRows, NumCols, 1);
cindex = RunAccuracy > 0;
hold on;
plot(Runs(cindex), RunAccuracy(cindex), 'k', 'LineWidth', 2);
plot(Runs(cindex), RunRejection(cindex), 'r');
hold off;
xlim([min(Runs(cindex)) max(Runs(cindex))]);
ylim([0 100]);
cnbiplot_vline(Rk(diff(Xk) > 0) + 0.5, '--k');
grid on;
xlabel('Runs');
ylabel('[%]');
legend('Accuracy', 'Rejection', 'Location', 'Best');
title('Single sample accuracy/rejection (per runs)');

% Accuracy/Rejection per day
subplot(NumRows, NumCols, 2);
cindex = DayAccuracy(:, 1) > 0;
hold on;
plot(Days(cindex), DayAccuracy(cindex, 1), 'k', 'LineWidth', 2);
hea = errorbar(Days(cindex), DayAccuracy(cindex, 1), DayAccuracy(cindex, 2), 'k', 'LineWidth', 2);
plot(Days(cindex), DayRejection(cindex, 1), 'r');
her = errorbar(Days(cindex), DayRejection(cindex, 1), DayRejection(cindex, 2), 'r');
hold off;
xlim([min(Days(cindex)) max(Days(cindex))]);
ylim([0 100]);
cnbiplot_vline(Dk(diff(Xk) > 0) + 0.5, '--k');
grid on;
set(gca, 'XTickLabel', labels.Dl(cindex, :));
xlabel('Days');
ylabel('[%]');
title('Single sample accuracy/rejection (per day)');

suptitle([subject ' - Accuracy/Rejection (' modality ')']);

%% Saving metadata

% Grouping results
singlesample.accuracy  = DayAccuracy;
singlesample.rejection = DayRejection;
singlesample.label     = labels.Dl;

savefile = [savedir '/' subject '.metadata.mat'];
if exist(savefile, 'file')
    cnbiutil_bdisp(['Loading metadata from: ' savefile]);
    load(savefile);
end

metadata.online.singlesample = singlesample;

cnbiutil_bdisp(['Saving singlesample (online) results in: ' savefile]);
save(savefile, 'metadata');
