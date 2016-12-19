function [FIG, POS] = cnbifig_set_position(FIG, TARGETPOS)
% [FIG, POS] = fig_set_position(FIG, TARGETPOS)
%
% Moves and resizes the given figure.
% 
% Input:
%   FIG             Figure handle
%   TARGETPOS       Figure position: 'Top', 'Bottom', 'Left', 'Right', 'All'
%
% Output:
%   FIG             Figure handle
%   POS             New figure position coordinates 
%
% SEE ALSO: set, get

    ScreenSize = get(0, 'ScreenSize');

    if strcmpi(TARGETPOS, 'Top')
        set(FIG, 'Position', [ScreenSize(1) ScreenSize(4)/2 ScreenSize(3) ScreenSize(4)/2]);
    elseif strcmpi(TARGETPOS, 'Bottom')
        set(FIG, 'Position', [ScreenSize(1) ScreenSize(2) ScreenSize(3) ScreenSize(4)/2]);
    elseif strcmpi(TARGETPOS, 'Left')
        set(FIG, 'Position', [ScreenSize(1) ScreenSize(2) ScreenSize(3)/2 ScreenSize(4)]);
    elseif strcmpi(TARGETPOS, 'Right')
        set(FIG, 'Position', [ScreenSize(3)/2 ScreenSize(2) ScreenSize(3)/2 ScreenSize(4)]);
    elseif strcmpi(TARGETPOS, 'All')
        set(FIG, 'Position', ScreenSize);
    else
        warning('chk:arg', 'Unrecognized position required');
    end
    
    POS = get(FIG, 'Position');
end