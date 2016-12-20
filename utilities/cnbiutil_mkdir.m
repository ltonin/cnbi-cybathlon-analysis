function [success, newpath] = cnbiutil_mkdir(parentdir, newdir)
% [success, NEWDIR] = cnbiutil_mkdir(PARENTDIR, NEWDIR)
%
% The function creates a new directory tree NEWDIR starting from PARENTDIR 
% (if PARENTDIR exists, otherwise it throws an error). 
%
% [success, NEWDIR] = cnbiutil_mkdir(NEWDIR)
%
% The function creates a new directory at NEWDIR 
%
% The function returns:
%       ->  1       In case of success in making the new directory
%       ->  2       If the directory already exists
%       ->  0       If an error occurred in creating the new directory
%
% SEE ALSO: mkdir
    
    if nargin == 1
        NEWDIR = parentdir;
    elseif nargin == 2
        NEWDIR = [parentdir newdir];
    end


    if (exist(NEWDIR, 'dir') ~= 7)
        cnbiutil_bdisp(['[' mfilename '] - Creating new directory at: ' NEWDIR]);
        success = mkdir(NEWDIR);
    else
        success = 2;
    end

    newpath = NEWDIR;
   
    
    

end