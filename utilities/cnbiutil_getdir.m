function dirpath = cnbiutil_getdir(path, pattern)
% dirpath = cnbiutil_getdir([path, pattern])
%
% This function returns the list of directories founded in path and that
% match pattern string.
% If no directories are found, error is thrown.
%
% If no agurments are provided, the function will use the default values:
%
%   path    = '.'       (current directory)
%   pattern = ''        (all directories are listed)
%
% SEE ALSO: cnbiutil_getfile


    if nargin < 1
        path = '.';
        pattern = '';
    end
    
    if nargin < 2
        pattern = '';
    end
    
    wildcard = ['*' pattern '*'];
    pathnames = [path '/' wildcard];
    entries = dir(pathnames);    
    NumEntries = length(entries);
    
    if(NumEntries == 0)
        error('chk:dir', ['[' mfilename '] No folders with ' pattern ' found in ' path])
    end
    
    dirpath = cell(NumEntries, 1);
    isdir   = false(NumEntries, 1);
    
    for dId = 1:NumEntries
        cname = entries(dId).name;
        dirpath{dId} = [path '/' cname '/'];
       
        isdir(dId)   = entries(dId).isdir;
        if strcmp(cname, '.') || strcmp(cname, '..')
            isdir(dId) = false;
        end
    end
    
    dirpath = dirpath(isdir);
end
    
    


