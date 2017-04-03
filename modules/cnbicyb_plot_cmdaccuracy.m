clearvars; clc; close all;

SubList = {'AN14VE', 'MA25VE'};
NumSubjects = length(SubList);

datapath  = [pwd '/analysis/'];
figuredir = './figures/';

fig1 = figure;
%fig_set_position(fig1, 'Top');

for sId = 1:NumSubjects
    csubject = SubList{sId};
    cfileonline = [datapath '/' csubject '.cmdaccuracy.online.mat'];
    cfilerace   = [datapath '/' csubject '.cmdaccuracy.race.mat'];
    
    cdata.online = load(cfileonline);
    cdata.race   = load(cfilerace);
    
    oDl = str2double(cdata.online.command.session.label);
    rDl = str2double(cdata.race.command.session.label);
    cDl = [oDl; rDl];
    
    
    
    oTsk = cdata.online.command.session.tasks;
    rTsk = cdata.race.command.session.tasks;
    tTsk = union(oTsk, rTsk);
    
    notsk = length(oTsk);
    nrtsk = length(rTsk);
    nttsk = length(tTsk);
    
    noday = length(oDl);
    nrday = length(rDl);
    c_oMk = ones(noday, 1);
    c_rMk = 2*ones(nrday, 1);
    cMk = [c_oMk; c_rMk];
    
    c_oacc = nan(nttsk, noday);
    c_racc = nan(nttsk, nrday);
    
    [~, oindex] = ismember(tTsk, oTsk);
    [~, rindex] = ismember(tTsk, rTsk);
    
    c_oacc(oindex>0, :) = cdata.online.command.session.accuracy(oindex(oindex > 0),:);
    c_racc(rindex>0, :) = cdata.race.command.session.accuracy(rindex(rindex > 0),:);
    
    tAcc = cat(2, c_oacc, c_racc);
    
    [sDl, iDl] = sort(cDl);
    stAcc = tAcc(:, iDl);
    sMk = cMk(iDl);
    
    DThick = 1:size(stAcc, 2);
    
    subplot(1, 2, sId);
    hold on;
    plot(DThick, stAcc');
    ax = gca;
    ax.ColorOrderIndex = 1;
    plot(DThick(sMk == 1), stAcc(:, sMk == 1)', 'o');
    ax = gca;
    ax.ColorOrderIndex = 1;
    plot(DThick(sMk == 2), stAcc(:, sMk == 2)', 's');
    plot(DThick, nanmean(stAcc), 'k', 'LineWidth', 2);
    hold off;
    legend(tTsk, 'location', 'best');
    grid on;
%     
%     [~, ~, oindex] = intersect(oTsk, totTasks);
%     [~, ~, rindex] = intersect(rTsk, totTasks);
%     
%     oacc = zeros(max(oindex), length(oDl));
%     oacc(oindex, :) = cdata.online.command.session.accuracy;
%     racc = zeros(max(rindex), length(rDl));
%     racc(rindex, :) = cdata.race.command.session.accuracy;
%     
%     nobservations = size(oacc, 2) + size(racc, 2);
%     
%     tAcc = zeros(length(totTasks), nobservations);
%     
%     tAcc(oindex, 1:size(oacc, 2)) = oacc;
%     tAcc(rindex, size(oacc, 2)+1:end) = racc;    
    
end