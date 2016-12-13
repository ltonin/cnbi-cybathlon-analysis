function [filenames, nfiles] = cnbiutil_getfile(path, extension, pattern, index)
% [filenames, nfiles] = cnbiutil_getfile(path, extension, [pattern, index])
%
% Look for files with a given extension in a given directory.
% It returns cell array with the absolute path of the files.
% Optional arguments are:
%   pattern     Default value: ''
%   index       Default value: all files
%
% SEE ALSO: cnbiutil_getdir

    if nargin < 2
        error('chk:arg', ['[' mfilename '] Provide at least a path and an extension']);
    end    

    if ischar(path)
        path = cellstr(path);
    end
    
    if iscell(path) == 0
        error('chk:path', ['[' mfilename '] Provide path argument as a string or as a cell array'])
    end
    
    
    NumDirs = length(path);
    
    if nargin == 2
        index = repmat({[]}, NumDirs, 1);
        pattern = '';
    end
    
    if nargin == 3
        index = repmat({[]}, NumDirs, 1);
    end
    
    if iscell(index) == false
        error('chk:index', ['[' mfilename '] Provide file indexes for each directory in cell format'])
    end
    
    if (NumDirs ~= length(index))
       error('chk:index', ['[' mfilename '] Provide file indexes for each directory'])
    end
    
    
    wildcard = ['*' pattern '*' extension];
    
    filenames = [];
    for dId = 1:NumDirs
        cpath  = path{dId};
        cindex = index{dId}; 
        
        cstruct = dir([cpath wildcard]);
        csize = length(cstruct);
        cfilename = cell(csize, 1);
        
        for eId = 1:csize
            cfilename{eId} = [cpath '/' cstruct(eId).name];
        end
           
        
        if isempty(cindex) == false
            cfilename = cfilename(cindex);
        end

        
        filenames = cat(1, filenames, cfilename);
        
    end
    
    nfiles = length(filenames);
    
end
