function [success, newpath] = cnbiutil_mkdir(parentdir, newdir)
% success = cnbiutil_mkdir(rootpath, dirpath)
%
% The function creates a new directory tree NEWDIR starting from PARENTDIR 
% (if PARENTDIR exists, otherwise it throws an error). 
%
% The function returns:
%       ->  1       In case of success in making the new directory
%       ->  2       If the directory already exists
%       ->  0       If an error occurred in creating the new directory
%
% SEE ALSO: mkdir
    
    if(exist(parentdir, 'dir') ~= 7)
        error('chk:root', ['[' mfilename '] - Parent directory does not exist']);
    end
    
    if (exist([parentdir newdir], 'dir') ~= 7)
        cnbiutil_bdisp(['[' mfilename '] - Creating new directory at: ' parentdir '/' newdir]);
        success = mkdir([parentdir newdir]);
    else
        success = 2;
    end
    
    newpath = [parentdir newdir];

end