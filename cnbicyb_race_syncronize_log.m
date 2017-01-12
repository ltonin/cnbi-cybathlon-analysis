clearvars; clc;

subject = 'AN14VE';


experiment  = 'cybathlon';
datapath    = ['/mnt/data/Research/' experiment '/' subject '/'];
targetpath  = [datapath '/' subject '_racemat/'];

logpath     = [datapath '/BrainRunnersLogs/v1326/'];
cnbiutil_bdisp(['[sync] - Syncronizing races between gdf (' datapath ') and logs (' logpath ')']);
syncGDFlog(subject, datapath, logpath, targetpath)

logpath     = [datapath '/BrainRunnersLogs/vCompetition/'];
cnbiutil_bdisp(['[sync] - Syncronizing competition races between gdf (' datapath ') and logs (' logpath ')']);
syncGDFlog(subject, datapath, logpath, targetpath)