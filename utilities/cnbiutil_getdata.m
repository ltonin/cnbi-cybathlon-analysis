function [files, numfiles] = cnbiutil_getdata(datapath, dirpattern, filepattern, extension)
% [files, numfiles] = cnbiutil_getdata(datapath, dirpattern, filepattern, extension)
%
% The function get cvsa datafiles from structured directory. It looks for
% directories in DATAPATH that match DIRPATTERN. Then, for each directory
% founded, it looks for all files with EXTENSION that match FILEPATTERN.
% It return a cell array with all files found and the number of files
% found.
%
% Example:
%
% Data folders is structured as follows:
%
% + data/
% |+ AN14VE_20160404/
%  |- AN14VE.20160404.144113.offline.mi.mi_rlbf.gdf
%  |- AN14VE.20160404.145454.offline.mi.mi_rlbf.gdf
%  |- AN14VE.20160404.150545.offline.mi.mi_rlbf.gdf
% |+ AN14VE_20160408/
%  |- AN14VE.20160408.143756.offline.mi.mi_rlsf.gdf
%  |- AN14VE.20160408.150617.offline.mi.mi_rlsf.gdf
%
% datafile = cnbiutil_getdata(datapath, 'AN14VE', '.mi.', '.gdf');
% datafile = 
%
%   '/mnt/data/Research/cybathlon/AN14VE//AN14VE_20160404//AN14VE.20160404.144113.offline.mi.mi_rlbf.gdf'
%   '/mnt/data/Research/cybathlon/AN14VE//AN14VE_20160404//AN14VE.20160404.145454.offline.mi.mi_rlbf.gdf'
%   '/mnt/data/Research/cybathlon/AN14VE//AN14VE_20160404//AN14VE.20160404.150545.offline.mi.mi_rlbf.gdf'
%   '/mnt/data/Research/cybathlon/AN14VE//AN14VE_20160408//AN14VE.20160408.143756.offline.mi.mi_rlsf.gdf'
%   '/mnt/data/Research/cybathlon/AN14VE//AN14VE_20160408//AN14VE.20160408.150617.offline.mi.mi_rlsf.gdf'
%
% Five files are found that match the input arguments.
%
% SEE ALSO: cnbiutil_getdir, cnbiutil_getfile


    folders = cnbiutil_getdir(datapath, dirpattern);
    nfolders = length(folders);
    
    files = {};
  
    for dId = 1:nfolders
        cfiles = cnbiutil_getfile(folders{dId}, extension, filepattern);
        files = cat(1, files, cfiles);
    end
    
    numfiles = length(files);

end