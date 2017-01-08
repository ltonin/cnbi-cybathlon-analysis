function [F, events, labels, settings] = cnbiutil_concatenate_data(filepaths)
% [F, events, labels, settings] = cnbiutil_concatenate_data(filepaths)
%
% The function concatenates preprocessed psd datafiles (done by 
% cnbicyb_processing_spectrogram). It concatenates:
% - psd(s) along 1st dimension (windows)
% - Event types, positions, duration
% It created labels for modality and days. It checks if frequency grid and
% settings of each file are equal (otherwise it raises an error).
%
% Input:
%   - filepaths         Cell array with absolute filepaths
%
% Output:
%   - F                 Feature matrix (windows x frequencies x channels)
%   - events            Structure with TYP, POS and DUR field
%   - events.extra      Structure with extra event in gdf-like format 
%                       (filled only for race runs)
%   - labels            Structure with Mk (modality), Rk (run), Rl(runlabel), Dk (day) and
%                       Dl (daylabel) fields. 
%                       Mk, Rk and Dk size: (windows x 1)
%                       Dl size: (number of days x 1)
%                       Rl size: (number of runs x 1);

    numfiles = length(filepaths);
    
    F     = [];
    TYP  = []; POS  = []; DUR  = [];
    tTYP = []; tPOS = []; tDUR = [];
    pTYP = []; pPOS = []; pDUR = [];
    bTYP = []; bPOS = []; bDUR = [];
    cTYP = []; cPOS = []; cDUR = [];
    eTYP = []; ePOS = []; eDUR = [];
    rTYP = []; rPOS = []; rDUR = [];
    freqs = [];
    settings = [];
    
    Rk = []; 
    Rl = cell(numfiles, 1);
    Dk = [];            
    Dl = [];            
    Mk = [];
    
    lday = [];
    nday = 0;
    fprintf('[io] - Concatenate psd data and events:\n');
    for fId = 1:numfiles
        cnbiutil_disp_progress(fId, numfiles, '        ');
        cfilepath = filepaths{fId};
        cdata = load(cfilepath);
        
        % Extract information from filename
        cinfo = cnbiutil_getfile_info(cfilepath);
        
        % Get modality from filename
        switch lower(cinfo.modality)
            case 'offline'
                modality = 0;
            case 'online'
                modality = 1;
            case 'race'
                modality = 2;
                if strcmp(cinfo.extra{1}, 'competition')
                    modality = 3;
                end
            otherwise
                error('chk:mod', ['[' mfilename '] Unknown modality']);
        end
        Mk = cat(1, Mk, modality*ones(size(cdata.psd, 1), 1));
        
        % Get day from filename
        if strcmpi(cinfo.date, lday) == false
            nday = nday + 1;
            Dl = cat(1, Dl, cinfo.date);
            lday = cinfo.date;
        end
        Dk = cat(1, Dk, nday*ones(size(cdata.psd, 1), 1));
        
        % Create run vector
        Rk = cat(1, Rk, fId*ones(size(cdata.psd, 1), 1));
        Rl{fId} = cinfo.extra;
        
        % Concatenate events
        if modality == 2    % race runs
            [cevents, cevtextra] = cnbiproc_extract_event_race(cdata.events);
            tTYP = cat(1, TYP, cevtextra.trl.TYP);
            tPOS = cat(1, POS, cevtextra.trl.POS + size(F, 1));
            tDUR = cat(1, DUR, cevtextra.trl.DUR);

            pTYP = cat(1, TYP, cevtextra.pad.TYP);
            pPOS = cat(1, POS, cevtextra.pad.POS + size(F, 1));
            pDUR = cat(1, DUR, cevtextra.pad.DUR);

            bTYP = cat(1, TYP, cevtextra.bci.TYP);
            bPOS = cat(1, POS, cevtextra.bci.POS + size(F, 1));
            bDUR = cat(1, DUR, cevtextra.bci.DUR);

            cTYP = cat(1, TYP, cevtextra.cmd.TYP);
            cPOS = cat(1, POS, cevtextra.cmd.POS + size(F, 1));
            cDUR = cat(1, DUR, cevtextra.cmd.DUR);

            eTYP = cat(1, TYP, cevtextra.eye.TYP);
            ePOS = cat(1, POS, cevtextra.eye.POS + size(F, 1));
            eDUR = cat(1, DUR, cevtextra.eye.DUR);

            rTYP = cat(1, TYP, cevtextra.race.TYP);
            rPOS = cat(1, POS, cevtextra.race.POS + size(F, 1));
            rDUR = cat(1, DUR, cevtextra.race.DUR);
        else
            cevents = cdata.events;
        end
        
        TYP = cat(1, TYP, cevents.TYP);
        DUR = cat(1, DUR, cevents.DUR);
        POS = cat(1, POS, cevents.POS + size(F, 1));
      
        % Concatenate features along 1st dimension (samples)
        F = cat(1, F, cdata.psd);
        
        % Save the current frequency grid and check if it is different from
        % the previous one. If so, it raises an error.
        cfreqs = cdata.freqs;
        if(isempty(freqs))
            freqs = cfreqs;
        end
        
        if(isequal(cfreqs, freqs) == false)
            error('chk:frq', ['Different frequency grid for file: ' cfilepath]);
        end
        
        % Save the current settings, remove number of sample and filename,
        % and compare it with the previous settings. If different, it
        % raises an error.
        csettings = cdata.settings;
        csettings.data.nsamples = nan;
        csettings.data.filename = nan;
        
        if(isempty(settings))
            settings = csettings;
        end
        
        if (isequaln(csettings, settings) == false)
             keyboard
            error('chk:stn', ['Different processing settings for file: ' cfilepath]);
           
        end
        
    end
    
    events.TYP = TYP;
    events.POS = POS;
    events.DUR = DUR;
    
    events.extra.trl.TYP = tTYP;
    events.extra.trl.POS = tPOS;
    events.extra.trl.DUR = tDUR;
    
    events.extra.pad.TYP = pTYP;
    events.extra.pad.POS = pPOS;
    events.extra.pad.DUR = pDUR;
    
    events.extra.bci.TYP = bTYP;
    events.extra.bci.POS = bPOS;
    events.extra.bci.DUR = bDUR;
    
    events.extra.cmd.TYP = cTYP;
    events.extra.cmd.POS = cPOS;
    events.extra.cmd.DUR = cDUR;
    
    events.extra.eye.TYP = eTYP;
    events.extra.eye.POS = ePOS;
    events.extra.eye.DUR = eDUR;
    
    events.extra.race.TYP = rTYP;
    events.extra.race.POS = rPOS;
    events.extra.race.DUR = rDUR;
    
    labels.Rk  = Rk;
    labels.Rl  = Rl;
    labels.Mk  = Mk;
    labels.Dl  = Dl;
    labels.Dk  = Dk;

end