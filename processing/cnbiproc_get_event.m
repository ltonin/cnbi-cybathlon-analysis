function [evtvec, evtstr] = cnbiproc_get_event(evt, nsamples, POS, TYP, DUR, OFFSET)
% [evtvec, evtstr] = cnbiproc_get_event(evt, nsamples, POS, TYP, DUR [, OFFSET])
%
% The function create a vector related to the provided events.
%
% Input:
%   - evt               List of the events to be extracted.
%   - nsamples          Length of the output vector.
%   - POS, TYP, DUR     Vector with positions, types and durations of all
%                       events (e.g., typical fields in gdf EVENT
%                       structure. The three vectors must have the same
%                       length. Optionally dur can be an integer. In this
%                       case, such a fixed duration is applied to all the
%                       extracted events.
%   - OFFSET            Optional argument to set the offset with respect to
%                       the original position of the events. A positive or
%                       negative integer means to start after or before the
%                       original position, respectively. [Default: 0]
%
% Output:
%   - evtvec            Vector of length nsamples. Value of each sample is
%                       the type of the extracted event or zeros otherwise.
%   - evtstr            Structure with the extracted events in the gdf
%                       EVENT structure format (evtstr.POS, evtstr.TYP,
%                       evstr.DUR). This structure refers to the new
%                       position and duration (for instance if an offset is
%                       provided).


    if length(DUR) == 1
        DUR = DUR*ones(size(POS));
    end

    if(isequal(length(POS), length(TYP)) == false)
        error('chk:length', 'POS and TYP must have the same lenght');
    end
    
    if(isequal(length(POS), length(DUR)) == false)
        error('chk:length', 'DUR must have the same length of POS and TYP')
    end
    
    if nargin == 5
        OFFSET = zeros(size(POS));
    end
    
    if length(OFFSET) == 1
        OFFSET = OFFSET*ones(size(POS));
    end
    
    if(isequal(length(POS), length(OFFSET)) == false)
        error('chk:length', 'OFFSET must have the same length of POS, TYP, DUR (or an integer to be applied to all events)')
    end
        
    
    % Compute the selected event index
    evtId = false(length(POS), 1);
    for eId = 1:length(evt)
        evtId = evtId | TYP == evt(eId);
    end
    
    % Extract the selected events
    pos_s    = POS(evtId);
    typ_s    = TYP(evtId);
    dur_s    = DUR(evtId);
    nsevents = length(pos_s);
    
    % Create label vector for the selected events
    evtvec     = zeros(nsamples, 1); 
    evtstr.POS = zeros(nsevents, 1);
    evtstr.TYP = zeros(nsevents, 1);
    evtstr.DUR = zeros(nsevents, 1);
    
    for eId = 1:nsevents
        cpos = pos_s(eId);
        ctyp = typ_s(eId);
        cdur = dur_s(eId);
        coff = OFFSET(eId);
        
        if coff >= cdur
            warning('chk:offset', ['[' mfilename '] Offset (' num2str(coff) ') for ' num2str(eId) 'th event is bigger than the event duration (' num2str(cdur) ')']);
        end
        
        cstart = cpos + coff;
        cstop  = cpos + cdur;
        evtvec(cstart:cstop) = ctyp;
        evtstr.POS(eId) = cstart;
        evtstr.TYP(eId) = ctyp;
        evtstr.DUR(eId) = length(cstart:cstop) - 1;
    end
    
    evtvec = evtvec(1:nsamples);

end