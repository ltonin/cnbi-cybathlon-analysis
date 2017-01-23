clearvars; clc; close all;

SubList = {'AN14VE', 'MA25VE'};
NumSubjects = length(SubList);

datapath  = [pwd '/analysis/'];
figuredir = './figures/';

MinRaceTime = 54;
MaxRaceTime = 327;

fig1 = figure;
cnbifig_set_position(fig1, 'All');
NumRows = 1;
NumCols = NumSubjects;

ax = gca;
Colors = ax.ColorOrder;

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
    % Plotting the least-squares fit
    htot = plot(cRuns, crt, 'o');
    csq = lsline(gca);
    set(htot, 'Visible', 'off');
    set(csq, 'LineWidth', 2);
    set(csq, 'LineStyle', '-');
    set(csq, 'Color', 'k');
    
    % Plotting race time for each protocol
    chandles = zeros(cNumProtocols, 1);
    cPl = cell(cNumProtocols, 1);
    for pId = 1:cNumProtocols
        cindex = cPk == cProtocols(pId);
        chandles(pId) = plot(cRuns(cindex), crt(cindex), '.', 'Color', Colors(cProtocols(pId), :), 'LineWidth', 1, 'MarkerSize', 20);
        cnames = unique(cdata.time.run.label.Pl(cindex));
        cPl{pId} = ['protocol' num2str(cProtocols(pId) - 2)];
    end
    
    % Highlight record
    [crecord, crunrecord] = min(crt);
    hrecord = plot(crunrecord, crecord, 'sr', 'LineWidth', 2, 'MarkerSize', 15);
    legend([chandles; hrecord], [cPl; 'Record: ' num2str(crecord, '%2.1f') ' s']);
    
    hold off;
    grid on;
    xlim([1 max(cRuns)]);
    
    
    
    % Adding x-ticks when date changes
    changedateId = find(diff(cDk)) + 1;
    set(gca, 'XTick', changedateId);
    set(gca, 'XTickLabel', cDl(cDk(changedateId), :));
    xticklabel_rotate([],90,[])
    xlabel('Race');
    ylabel('Time [s]');
    title(csubject);
    
    % Compute correlation and add annotation
    [ccorr, cpval] = corr(cRuns, crt);

    cpos = get(gca,'Position');
    %cpos(1) = cpos(1) + cpos(3) - 0.1;
    cpos(2) = cpos(2) - 0.1;
    textcolor = 'k';
    annotation('textbox', cpos, 'String', ['rho=' num2str(ccorr, '%3.2f') ', p=' num2str(cpval, '%3.3f')], 'LineStyle', 'none', 'Color', textcolor, 'FontWeight', 'bold')
    
    % Add max and min race time
    cnbiplot_hline(MinRaceTime, '--k', ['Min race time (' num2str(MinRaceTime) ' s)']);
    cnbiplot_hline(MaxRaceTime, '--k', ['Max race time (' num2str(MaxRaceTime) ' s)']);
end

suptitle('Race Time')

cnbifig_export(fig1, [figuredir '/cybathlon.journal.racetime.png'], '-png');
