clearvars; clc; 

subject = 'AN14VE';
%subject = 'MA25VE';

pattern     = '.mi.';
modality    = 'race';

experiment  = 'cybathlon';
%datapath    = [pwd '/analysis/'];
datapath    = '/home/sperdikis/Desktop/tst/AN14VE/AN14VE_RaceMat/';
%datapath    = '/home/sperdikis/Desktop/tst/MA25VE/MA25VE_RaceMat/';
figuredir  = '/figures/';

rejectlim = 240; % Reject races above this limit, 

%% Get datafiles
[Files, NumFiles] = cnbiutil_getfile(datapath, '.mat', [subject '*' modality '*' pattern]);

numfiles = length(Files);
Rk = []; 
Rl = cell(numfiles, 1);
Dk = [];            
Dl = [];            
Mk = [];
RT = [];

lday = [];
nday = 0;
fprintf('[io] - Get race times:\n');
for fId = 1:numfiles
    cnbiutil_disp_progress(fId, numfiles, '        ');
    cfilepath = Files{fId};
    cdata = load(cfilepath);

    % Extract information from filename
    cinfo = cnbiutil_getfile_info(cfilepath);

    % Get modality from filename
    switch lower(cinfo.modality)
        case 'offline'
            modality = 0;
        case 'online'
            modality = 1;
        case 'race'
            modality = 2;
            if strcmp(cinfo.extra{1}, 'competition')
                modality = 3;
            end
        otherwise
            error('chk:mod', ['[' mfilename '] Unknown modality']);
    end
    Mk = cat(1, Mk, modality*ones(1, 1));
    
    % Get day from filename
    if strcmpi(cinfo.date, lday) == false
        nday = nday + 1;
        Dl = cat(1, Dl, cinfo.date);
        lday = cinfo.date;
    end
    Dk = cat(1, Dk, nday*ones(1, 1));

    % Create run vector
    Rk = cat(1, Rk, fId*ones(1, 1));
    Rl{fId} = cinfo.extra;

    % Extract race time
    if modality >= 2    % race runs
        RT = [RT ; cdata.Race.RaceTime];
    else
        error('There is no race time in modality ~= 2');
    end
end

% Reject races
if(rejectlim > 0)
    IndRej = find(RT>rejectlim);
    RT(IndRej) = [];
    Rk(IndRej) = [];
    Rl(IndRej) = [];
    Dk(IndRej) = [];
    Mk(IndRej) = [];
end

%% Plotting per run results
fig1 = figure;
cnbifig_set_position(fig1, 'All');
xx=Rk;
p1 = polyfit(xx,RT,1);
p2 = polyfit(xx,RT,2);
lf1 = @(x)(p1(1)*x +p1(2));
lf2 = @(x)(p2(1)*x.^2 +p2(2)*x+p2(3));
plot(xx,RT,'*r',xx,lf1(xx),'k',xx,lf2(xx),'b','LineWidth',3, 'MarkerSize',10,'LineWidth',3);
xlabel('Race Index (chronologically)','FontSize',20,'LineWidth',3);
ylabel('Race Completion Time (sec)','FontSize',20,'LineWidth',3);
title(subject);
axis([0 max(Rk)+1 min(RT)-5 max(RT)+5]);
legend({'Race Time','Linear Fit','2nd order polynomial fit'});
set(gca,'FontSize',20,'LineWidth',3);
[r pvalPearson] = corr(xx,RT,'type','Pearson');
disp(['Per race, Pearson Correlation r = ' num2str(r) ' with pval = ' num2str(pvalPearson)]);
[rho pvalSperaman] = corr(xx,RT,'type','Spearman');
disp(['Per race, Spearman Correlation r = ' num2str(rho) ' with pval = ' num2str(pvalSperaman)]);


%% Plotting per session results

SMRT = [];
SSRT = [];
for s=1:length(Dl)
    SMRT(s) = mean(RT(Dk==s));
    SSRT(s) = std2(RT(Dk==s));
end
fig2 = figure;
cnbifig_set_position(fig2, 'All');
xx=unique(Dk);
p1 = polyfit(xx',SMRT,1);
p2 = polyfit(xx',SMRT,2);
lf1 = @(x)(p1(1)*x +p1(2));
lf2 = @(x)(p2(1)*x.^2 +p2(2)*x+p2(3));
plot(xx,SMRT,'r',xx,lf1(xx),'k',xx,lf2(xx),'b','LineWidth',3, 'MarkerSize',10,'LineWidth',3);
legend({'Race Time','Linear Fit','2nd order polynomial fit'});
hold on;
shadedErrorBar(xx,SMRT,SSRT,'*r',1);
hold off;
xlabel('Race Session','FontSize',20,'LineWidth',3);
ylabel('Average Race Completion Time (sec)','FontSize',20,'LineWidth',3);
title(subject);
axis([1 max(Dk) min(SMRT-SSRT) max(SMRT+SSRT)]);
set(gca,'FontSize',20,'LineWidth',3);
set(gca,'XTick',unique(Dk));
set(gca,'XTickLabel',Dl);
xticklabel_rotate([],45,[])

[r pvalPearson] = corr(xx,SMRT','type','Pearson');
disp(['Per session, Pearson Correlation r = ' num2str(r) ' with pval = ' num2str(pvalPearson)]);
[rho pvalSperaman] = corr(xx,SMRT','type','Spearman');
disp(['Per session, Spearman Correlation r = ' num2str(rho) ' with pval = ' num2str(pvalSperaman)]);
