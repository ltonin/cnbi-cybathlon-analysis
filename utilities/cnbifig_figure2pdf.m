function cnbifig_figure2pdf(F, filename, orientation, resize)
% cnbifig_figure2pdf(f, filename [, orientation, resize])
%
% The function exports the given figure in PDF. 
%
% Input:
%   - F                 Figure handle
%   - filename          PDF target filename
%   - orientation       Option to change orientation of page when printing 
%                       figure to PDF ('portrait', 'landscape') 
%                       [Default: 'landscape']
%   - resize            Option to expand figure to fill page when printing 
%                       figure to PDF ('-fillpage', '-bestfit) 
%                       [Default: '-bestfit']
%
% SEE ALSO: print

    if nargin == 2
        orientation = [];
        resize      = [];
    elseif nargin == 3
        resize      = [];
    end
    
    if isempty(orientation)
        orientation = 'landscape';
    end
    
    if isempty(resize)
        resize = '-bestfit';
    end
    
    set(F, 'PaperType', 'a4');
    set(F, 'PaperOrientation',   orientation);
    set(F, 'PaperUnits', 'normalized');
    print(resize, F, '-dpdf', filename);
end