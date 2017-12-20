clearvars; clc;

subject = 'AN14VE';
% subject = 'MA25VE';


experiment  = 'cybathlon';
datapath    = ['/mnt/data/Research/' experiment '/' subject '/'];
% datapath    = ['/home/sperdikis/Data/Raw/Cybathlon/' subject '/'];
targetpath  = [datapath '/' subject '_racemat/'];

%% Sync race and gdf
logpath = [datapath '/BrainRunnersLogs/v1326/'];
cnbiutil_bdisp(['[sync] - Syncronizing races between gdf (' datapath ') and logs (' logpath ')']);
syncGDFlog(subject, datapath, logpath, targetpath)

%% Sync race time from log (manually for each subject)
cnbiutil_bdisp('[sync] - Getting only race time from log');
logpath = [datapath '/BrainRunnersLogs/v1326/'];
if strcmpi(subject, 'MA25VE')
    disp('         Date: 20160825');
    RaceTimeFromLog(subject, '20160825', 1, logpath, targetpath);  % First face off (no gdf data)
    disp('         Date: 20160922');
    RaceTimeFromLog(subject, '20160922', 1, logpath, targetpath);  % Second face off (no gdf data)
elseif strcmpi(subject, 'AN14VE')
    disp('         Date: 20160825');
    RaceTimeFromLog(subject, '20160825', 2, logpath, targetpath);  % First face off (no gdf data)
    disp('         Date: 20161006');
    RaceTimeFromLog(subject, '20161006', 1, logpath, targetpath);  
    disp('         RehearsalStyle');
    findRehearsalStyleRaces(subject, datapath, targetpath)
end

%% Sync competition data
logpath_competition = [datapath '/BrainRunnersLogs/vCompetition/'];
cnbiutil_bdisp(['[sync] - Syncronizing competition races between gdf (' datapath ') and logs (' logpath_competition ')']);
syncGDFlog(subject, datapath, logpath_competition, targetpath)

