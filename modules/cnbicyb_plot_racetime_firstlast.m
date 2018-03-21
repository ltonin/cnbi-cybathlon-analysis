clearvars; clc; close all;

SubList = {'AN14VE', 'MA25VE'};
NumSubjects = length(SubList);

datapath  = [pwd '/analysis/'];
figuredir = './figures/';
datasubmission = [pwd '/data/'];

fig1 = figure;
cnbifig_set_position(fig1, 'All');
NumRows = 1;
NumCols = NumSubjects;

ax = gca;
Colors = ax.ColorOrder;

[~, plotdatapath] = cnbiutil_mkdir(pwd, '/plotdata');
data = struct([]);

for sId = 1:NumSubjects
    csubject = SubList{sId};
    cfilepath = [datapath '/' csubject '.time.race.mat'];
    
    cdata = load(cfilepath);
    
    crt = cdata.time.run.values;
    cRk = cdata.time.run.label.Rk;
    cRuns = unique(cRk);
    cDk = cdata.time.run.label.Dk;
    cDl = cdata.time.run.label.Dl;
    cDays = unique(cDk);
    cPk = cdata.time.run.label.Pk;
    cProtocols = unique(cPk);
    cNumProtocols = length(cProtocols);
    
    subplot(NumRows, NumCols, sId);
    
    hold on;
    
    % Plotting bars of first and last
    MeanVec = [mean(crt(find(ismember(cDk,[1:4])))) ; mean(crt(find(ismember(cDk,[max(cDays)-3:max(cDays)]))))];
    StdVec = [std2(crt(find(ismember(cDk,[1:4])))) ; std2(crt(find(ismember(cDk,[max(cDays)-3:max(cDays)]))))];
    NVec = [length(crt(find(ismember(cDk,[1:4])))) ; length(crt(find(ismember(cDk,[max(cDays)-3:max(cDays)]))))];
    pVal = ranksum(crt(find(ismember(cDk,[1:4]))), crt(find(ismember(cDk,[max(cDays)-3:max(cDays)]))));
    
    hBar = bar(MeanVec,'LineWidth',2);
    hEBar = errorbar([1:2],MeanVec,zeros(1,2),StdVec,'.','LineWidth',2);
    H = sigstar({[1,2]},pVal);
    hold off;
    xlim([0 3]);
    ylim([50 200]);
    % Adding x-ticks
    set(gca, 'XTick', [1:2]);
    set(gca, 'XTickLabel', {['First 4 sessions, N=' num2str(NVec(1))],['Last 4 sessions, N=' num2str(NVec(2))]});
    ylabel('Time [sec]');
    title(csubject);
    
    % Storing plot data
    data(sId).racetime      = crt;
    data(sId).labels.day_id = cDk;
end

%% Saving plot data
cnbiutil_bdisp(['Saving plot data in ' plotdatapath]);
save([plotdatapath '/Fig1B.mat'], 'data');

%% Plotting
suptitle('Race Time')
cnbifig_export(fig1, [figuredir '/cybathlon.journal.racetimefirstlast.png'], '-png');
cnbifig_export(fig1, [figuredir '/cybathlon.journal.racetimefirstlast.pdf'], '-pdf');
