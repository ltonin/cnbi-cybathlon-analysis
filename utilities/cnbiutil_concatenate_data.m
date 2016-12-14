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
%   - labels            Structure with Mk (modality), Dk (day) and
%                       Dl (daylabel) fields. 
%                       Mk and Dk size: (windows x 1)
%                       Dl size: (number of days x 1)

    numfiles = length(filepaths);
    
    F     = [];
    TYP   = [];
    POS   = [];
    DUR   = [];
    freqs = [];
    settings = [];
    
    Dk = [];
    Dl = [];
    Mk = [];
    
    lday = [];
    nday = 0;
    fprintf('[io] - Concatenate psd data and events:\n');
    for fId = 1:numfiles
        util_disp_progress(fId, numfiles, '        ');
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
            % TO BE ADDED: race, faceoff, cybathlon
            otherwise
                error('chk:mod', ['[' mfilename '] Unknown modality']);
        end
        Mk = cat(1, Mk, modality*ones(size(cdata.psd, 1), 1));
        
        % Get day from filename
        if strcmpi(cinfo.date, lday) == false
            nday = nday + 1;
            Dk = cat(1, Dk, nday*ones(size(cdata.psd, 1), 1));
            Dl = cat(1, Dl, cinfo.date);
            lday = cinfo.date;
        end
        
        % Concatenate events
        cevents = cdata.events;
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
    
    labels.Mk  = Mk;
    labels.Dl  = Dl;
    labels.Dk  = Dk;

end