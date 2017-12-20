clearvars; clc; close all;

subject = 'AN14VE';

datapath  = [pwd '/analysis/'];
figuredir = './figures/';

SelectedClassId = [771 773];
SelectedClassLb = {'Both feet', 'Both hands'};

AltSelectedClassId = [770 771];
AltSelectedClassLb = {'RightHand', 'BothFeet'};

PatternLocationsId = {[4 9 14], [2 7 12 6 11 16]};
PatternLocationsLb = {'FCz-Cz-CPz', 'FC3-C3-CP3-CP4-C4-FC4'}; 
PatternLabels = {'Medial locs', 'Lateral locs'};
NumPatterns = length(PatternLocationsId);

ParadigmId = [3 5 6];
NumParadigm = length(ParadigmId);

BetaFreqs = 22:2:32;

%% Loading pre-computed discriminancy
filepath = [datapath '/' subject '.discriminancy.maps.mat'];
data = load(filepath);     

%% Extracting data
fs = data.discriminancy.run.fisherscore;

FreqGrid = data.discriminancy.freqs;
[~, SelBetaFreqIds] = intersect(FreqGrid, BetaFreqs);

NumFreqs = length(FreqGrid);
NumChans = size(fs, 1)./NumFreqs;
    
combinations = data.discriminancy.run.combinations;
NumCombinations = size(combinations, 1);
    
Mk = data.discriminancy.run.label.Mk;
Modalities   = unique(Mk);
NumModalities = length(Modalities);
Ml = {'Offline', 'Online', 'Race', 'Competition'};
    
Dk = data.discriminancy.run.label.Dk;
Dl = data.discriminancy.run.label.Dl;
    
Pk = data.discriminancy.run.label.Pk;
Pl = data.discriminancy.run.label.Pl;
    
SelectedCombination    = find(ismember(combinations, SelectedClassId, 'rows'));
AltSelectedCombination = find(ismember(combinations, AltSelectedClassId, 'rows'));
    

%% Get fisherscore for selected and alternative task combinations
% Get fisherscore for the main combination (both hand vs. both feet)
sfs = reshape(fs(:, :, SelectedCombination), [NumFreqs NumChans size(fs, 2)]);

% Get fisherscore for the alternative combination (right hand vs. both feet)
afs = reshape(fs(:, :, AltSelectedCombination), [NumFreqs NumChans size(fs, 2)]);

% If fisherscore for the main combination is nan, then use the
% alternative combination
ffs = sfs;
for vId = 1:size(ffs, 3)
    if sum(sum(isnan(ffs(:, :, vId)))) > 0
        ffs(:, :, vId) = afs(:, :, vId);
    end
end

%% Splitting fisherscore per paradigm
RaceModality = 2; % Competition modality is equal to 3
SelNumRuns = 10;
l_avgfs = zeros(2, NumParadigm, NumPatterns);
l_stdfs = zeros(2, NumParadigm, NumPatterns);
o_avgfs = zeros(2, NumParadigm);
o_stdfs = zeros(2, NumParadigm);
o_pval  = zeros(NumParadigm, 1);
l_pval  = zeros(NumParadigm, NumPatterns);

for pId = 1:NumParadigm
    cindex = Mk >= RaceModality & Pk == ParadigmId(pId);
    
    o_cfs = squeeze(nanmean(nanmean(ffs(SelBetaFreqIds, cell2mat(PatternLocationsId), cindex), 1), 2));
    
    rstartid = 1:SelNumRuns;
    rstopid  = length(o_cfs)-(SelNumRuns-1):length(o_cfs);
    
    for lId = 1:NumPatterns
        l_cfs = squeeze(nanmean(nanmean(ffs(SelBetaFreqIds, PatternLocationsId{lId}, cindex), 1), 2));
        l_avgfs(1, pId, lId) = mean(l_cfs(rstartid)); 
        l_avgfs(2, pId, lId) = mean(l_cfs(rstopid)); 
        l_stdfs(1, pId, lId) = std(l_cfs(rstartid)); 
        l_stdfs(2, pId, lId) = std(l_cfs(rstopid)); 
        l_pval(pId, lId) = ranksum(l_cfs(rstartid), l_cfs(rstopid));
    end
    
    o_avgfs(1, pId) = mean(o_cfs(rstartid));
    o_avgfs(2, pId) = mean(o_cfs(rstopid));
    o_stdfs(1, pId) = std(o_cfs(rstartid));
    o_stdfs(2, pId) = std(o_cfs(rstopid));
    
    o_pval(pId) = ranksum(o_cfs(rstartid), o_cfs(rstopid));
end

%% Plotting
fig1 = figure;
cnbifig_set_position(fig1, 'Top');

for lId = 1:NumPatterns
    subplot(1, NumPatterns+1, lId);
    %bar(l_avgfs(:, :, lId)', l_stdfs(:, :, lId)');
    barwitherr(cat(3, zeros(size(l_stdfs(:, :, lId)')), l_stdfs(:, :, lId)'), l_avgfs(:, :, lId)');
    title(PatternLabels{lId});
    ylim([0 0.5]);
    xlabel('Control paradigm');
    ylabel('Discriminancy');
    grid on;
end

subplot(1, NumPatterns+1, NumPatterns+1);
barwitherr(cat(3, zeros(size(o_stdfs')), o_stdfs'), o_avgfs');
title('Average over locations');
ylim([0 0.5]);
xlabel('Control paradigm');
ylabel('Discriminancy');
grid on;

legend(['First ' num2str(SelNumRuns) ' runs'], ['Last ' num2str(SelNumRuns) ' runs']);

 