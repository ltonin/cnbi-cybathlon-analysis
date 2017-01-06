clearvars; clc;

subject = 'AN14VE';


experiment  = 'cybathlon';
datapath    = ['/mnt/data/Research/' experiment '/' subject '/'];
% logpath     = [datapath '/BrainRunnersLogs/v1326/'];
logpath     = [datapath '/BrainRunnersLogs/vCompetition/'];
targetpath  = [datapath '/' subject '_racemat/'];

cnbiutil_bdisp(['[sync] - Syncronizing races between gdf (' datapath ') and logs (' logpath ')']);
syncGDFlog(subject, datapath, logpath, targetpath)