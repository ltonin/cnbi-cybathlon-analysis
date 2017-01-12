function cnbifig_export(F, filename, format, orientation, resize)
% cnbifig_export(f, filename, format [, orientation, resize])
%
% The function exports the figure in the given format. 
%
% Input:
%   - F                 Figure handle
%   - filename          Target filename
%   - format            Format to be used to export the figure. Available
%                       formats are:
%                       - '-pdf'
%                       - '-png'
%                       - '-jpg'
%   - orientation       Option to change orientation of page when printing 
%                       figure to PDF ('portrait', 'landscape') 
%                       [Default: 'landscape']
%   - resize            Option to expand figure to fill page when printing 
%                       figure to PDF ('-fillpage', '-bestfit) 
%                       [Default: '-bestfit' for pdf and '-fillpage' for other formats]
%
% SEE ALSO: print
    
    if nargin < 3
        error('chk:arg', 'At least three argument are required: figure handle, filename and format');
    end

    if nargin == 3
        orientation = [];
        resize      = [];
    elseif nargin == 4
        resize      = [];
    end
    
    switch(lower(format))
        case '-pdf'
            dformat = '-dpdf';
        case '-png'
            dformat = '-dpng';
        case '-jpg'
            dformat = '-djpg';
        otherwise
            error('chk:arg', [format ' format not implemented']);
    end
    
    if isempty(orientation)
         if strcmpi(dformat, '-dpdf')
             orientation = 'landscape';
         else
            orientation = 'portrait';
         end
    end
    
    if isempty(resize)
        resize = '-bestfit';
    end
    
    set(F, 'PaperType', 'a4');
    set(F, 'PaperOrientation',   orientation);
    set(F, 'PaperUnits', 'normalized');
    
    if strcmpi(dformat, '-dpdf')
        print(resize, F, dformat, filename);
    else
        print(F, dformat, filename);
    end
end