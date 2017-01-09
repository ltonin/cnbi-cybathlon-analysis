function [evt, extra] = cnbiproc_extract_event_race(events, PadId, ArtId, EndId)
% [evt, extra] = cnbiproc_extract_event_race(events [, PadId, ArtId, EndId])
%
% The function translates cybathlon race events and converts them in the
% standard gdf format (TYP, POS, DUR), according to the following rules:
%
%   Input events format:
%       -   Events equal to PadId correspond to the avatar entering in the 
%           pad (Default values: [768 771 773 783])
%       -   Events equal to ArtId correpond to the onset (ArtId(1)) and 
%           offset (ArtId(2)) of artifact detection
%       -   Event equals to EndId correponds to the end of the race
%           (Default value: 666)
%       -   Command events are assumed to have the format 
%           ['6' dec2hex(PadId)]
%       -   events is a structure in gdf-like format 
%           (events.TYP, events.POS, events.DUR)
%
%   Conversion operations:
%       -   A trial is considered to start from the entering pad event to
%           the first correct command delivered OR the next pad.
%       -   Command events occuring between two artifact events [ArtId(1)
%           ArtId(2)] are discarded.
%       -   Correct duration is added to all events
%   
%   Output events format:
%       -   evt is a gdf-like event structure (evt.TYP, evt.POS, evt.DUR)
%           with inside the following event types: 
%           - 1                     Race event
%           - [768 771 773 783]     Trial events {'Slide', 'Jump', 'Speed', 'Rest'}
%           - ['6' dec2hex(PadId)]  Command events {'Slide', 'Jump', 'Speed'}
%           - 267                   Eye artifact onset
%       -   extra is a gdf-like event structure with all events grouped in
%           subfields. The events in the subfields are: 
%           - race                  Race start
%           - trial                 Trial
%           - pad                   Avatar entering in the pad
%           - bci                   Commands delivered by BCI
%           - commands              Commands implemented in the game
%           - eye                   Eye artifact onset
%
% SEE ALSO: syncGDFlog


    % Default input parameters
    DefaultPadId = [768 769 770 771 773 783];
    DefaultArtId = [267 268];
    DefaultEndId = 666;
    
    if nargin == 1
        PadId = DefaultPadId;
        ArtId = DefaultArtId;
        EndId = DefaultEndId;
    elseif nargin == 2
        ArtId = DefaultArtId;
        EndId = DefaultEndId;
    elseif nargin == 3
        EndId = DefaultEndId;
    end
    
    CmdId = hex2dec(cat(2, repmat('6', length(PadId), 1), dec2hex(PadId)));
    
    % Extracting start/end event position
    startpos = events.POS(1);
    endpos   = events.POS(events.TYP == EndId);
    
    if isempty(endpos)
        error('chk:end', 'No end event found');
    end
    
    if(length(endpos) > 1)
        error('chk:end', 'More than 1 end event');
    end
        
    
    racetyp  = 1;
    racepos  = startpos;
    racedur  = endpos - startpos;
    
    % Extracting pad event types and positions
    index = false(length(events.TYP), 1);
    for pId = 1:length(PadId)
        index = index | events.TYP == PadId(pId);
    end
    
    padtyp = events.TYP(index);
    padpos = events.POS(index);
    
    
    % Extracting command event types and positions
    index = false(length(events.TYP), 1);
    for cId = 1:length(CmdId)
        index = index | events.TYP == CmdId(cId);
    end
    
    cmdtyp_o = events.TYP(index);
    cmdpos_o = events.POS(index);
    
    % Extracting artifact position and duration
    [eyepos, eyedur] = get_artifact_events(events, ArtId(1), ArtId(2));
    eyelabel = false(endpos, 1);
    try
    for eId = 1:length(eyepos)
        ceyepos = eyepos(eId);
        ceyedur = eyedur(eId);
        eyelabel(ceyepos:ceyepos+ceyedur) = true;
    end
    catch
        keyboard
    end
    
    % Removing commands delivered during artifact on
    eyefreeId = eyelabel(cmdpos_o) == false;
    
    if(isempty(eyefreeId) == false)
        cmdtyp = cmdtyp_o(eyefreeId);
        cmdpos = cmdpos_o(eyefreeId);
    else
        cmdtyp = cmdtyp_o;
        cmdpos = cmdpos_o;        
    end
    
    % Computing trial duration
    npads  = length(padtyp);
    paddur = zeros(npads, 1);
    trltyp = zeros(npads, 1);
    trlpos = zeros(npads, 1);
    trldur = zeros(npads, 1);
    for pId = 1:npads
        cpadtyp = padtyp(pId);
        cpadpos = padpos(pId);
        
        % Default duration until the next pad (or finish line for the last pad)
        if(pId < npads)
            cpaddur = padpos(pId + 1) - cpadpos;
        else
            cpaddur = endpos - cpadpos;
        end
        paddur(pId) = cpaddur;
        
        ctrdur = cpaddur;
        
        
        try
        % Find types and positions of commands delivered in the current pad interval
        cindex = find(cmdpos >= cpadpos & cmdpos < (cpadpos + cpaddur)); 
        catch
        disp('a');    
        end
        if (isempty(cindex) == false)
            ccmdtyp = cmdtyp(cindex);
            ccmdpos = cmdpos(cindex);
            
            % find the first command correct (if exists) with respect to
            % the current pad type
            ccmdvalidId = find(ismember(ccmdtyp - hex2dec('6000'), cpadtyp), 1);
            
            if (isempty(ccmdvalidId) == false)
                ctrdur = ccmdpos(ccmdvalidId) - cpadpos;
            end
        end
        trltyp(pId) = cpadtyp;
        trlpos(pId) = cpadpos;
        trldur(pId) = ctrdur;
        
        if ctrdur < 0
            keyboard;
        end
    end
    
    extra.trl.TYP = trltyp;
    extra.trl.POS = trlpos;
    extra.trl.DUR = trldur;
    extra.trl.DUR(extra.trl.DUR == 0) = 1;
    
    extra.pad.TYP = padtyp;
    extra.pad.POS = padpos;
    extra.pad.DUR = paddur;
    extra.pad.DUR(extra.pad.DUR == 0) = 1;
    
    extra.bci.TYP = cmdtyp_o;
    extra.bci.POS = cmdpos_o;
    extra.bci.DUR = ones(length(cmdpos_o), 1);
    extra.bci.DUR(extra.bci.DUR == 0) = 1;
    
    extra.cmd.TYP = cmdtyp;
    extra.cmd.POS = cmdpos;
    extra.cmd.DUR = ones(length(cmdpos), 1);
    extra.cmd.DUR(extra.cmd.DUR == 0) = 1;
    
    extra.eye.TYP = 267*ones(length(eyepos), 1);
    extra.eye.POS = eyepos;
    extra.eye.DUR = eyedur;
    extra.eye.DUR(extra.eye.DUR == 0) = 1;
    
    extra.race.TYP = racetyp;
    extra.race.POS = racepos;
    extra.race.DUR = racedur;
    extra.race.DUR(extra.race.DUR == 0) = 1;
    
    wTYP = [extra.race.TYP' extra.trl.TYP' extra.cmd.TYP' extra.eye.TYP'];
    wPOS = [extra.race.POS' extra.trl.POS' extra.cmd.POS' extra.eye.POS'];
    wDUR = [extra.race.DUR' extra.trl.DUR' extra.cmd.DUR' extra.eye.DUR'];
    
    [~, sindex] = sort(wPOS);
    evt.TYP = wTYP(sindex)';
    evt.POS = wPOS(sindex)';
    evt.DUR = wDUR(sindex)';
end

function [pos, dur] = get_artifact_events(events, onkey, offkey)
try
    % Get all the on/off positions
    onpos  = events.POS(events.TYP == onkey);
    offpos = events.POS(events.TYP == offkey);
    
    if (isempty(onpos) == false && isempty(offpos) == false) 
        % Get start and end positions
        starpos = events.POS(1);
        endpos  = events.POS(end);

        if(offpos(1) < onpos(1))        % Missing 1 onset event, 
            onpos = [starpos; onpos];   % I assume that the race begins after the onset event. 
        end                             % Adding it at the beginning
        
        if(onpos(end) > offpos(end))    % Missing 1 offset event, 
            offpos = [offpos; endpos];  % I assume that the race finish before the offset event.
        end                             % Adding it at the end

        % Create a POS and DUR for onkey and offkey
        pos = [];
        dur = [];
        for eId = 1:length(onpos)
            cstart = onpos(eId);

            pstop = 0;
            if(isempty(pos) == false)
                pstop = pos(end) + dur(end);
            end

            if(cstart >= pstop)
                cstop  = offpos(find(offpos > cstart, 1, 'first'));
                pos = cat(1, pos, cstart);
                dur = cat(1, dur, cstop - cstart);
            end

        end

        if isequal(length(pos), length(dur))  == false
            error('chk:eye', 'Unrecoverable mismatch of onset/offset events for eye artifacts');
        end
    else
        pos = [];
        dur = [];
    end
catch
    keyboard
end
end

