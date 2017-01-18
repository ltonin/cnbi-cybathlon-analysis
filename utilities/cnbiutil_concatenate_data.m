function [F, events, labels, settings, classifiers] = cnbiutil_concatenate_data(filepaths)
% [F, events, labels, settings, classifiers] = cnbiutil_concatenate_data(filepaths)
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
%   - classifiers       Structure with size equal to the number of
%                       different classifier
%   - labels            Structure with Mk (modality), Rk (run), Rl(runlabel), Dk (day) and
%                       Dl (daylabel) fields, Xk (classifier id). 
%                       Mk, Rk Dk and Xk size: (windows x 1)
%                       Dl size: (number of days x 1)
%                       Rl size: (number of runs x 1)


    numfiles = length(filepaths);
    
    % Getting size info to allocate memory and speedup the concatenation
    datasize   = get_data_size(filepaths);
    NumSamples = sum(datasize(1, :));
    NumFreqs   = unique(datasize(2, :));
    NumChans   = unique(datasize(3, :));
    
    if length(NumFreqs) > 1
        error('chk:dim', 'Files do not have the same number of frequencies. Cannot concatenate');
    end
    
    if length(NumChans) > 1
        error('chk:dim', 'Files do not have the same number of channels. Cannot concatenate');
    end
    
    % Pre-allocate variables
    F     = zeros(NumSamples, NumFreqs, NumChans);
    Rk    = zeros(NumSamples, 1);
    Dk    = zeros(NumSamples, 1);
    Mk    = zeros(NumSamples, 1);
    Xk    = zeros(NumSamples, 1);
    Rl    = cell(numfiles, 1);
    Dl    = [];
    
    % I don't know the number of events a-priori. Anyway it should not
    % impact the speed too much
    TYP   = [];  POS  = [];  DUR  = [];
    tTYP  = [];  tPOS = [];  tDUR = [];
    pTYP  = [];  pPOS = [];  pDUR = [];
    bTYP  = [];  bPOS = [];  bDUR = [];
    cTYP  = [];  cPOS = [];  cDUR = [];
    cgTYP = []; cgPOS = []; cgDUR = [];
    eTYP  = [];  ePOS = [];  eDUR = [];
    rTYP  = [];  rPOS = [];  rDUR = [];
    prTYP = []; prPOS = []; prDUR = [];
    
    freqs = [];
    settings = [];
    classifiers = [];
    
    lday = [];
    nday = 0;
    fprintf('[io] - Concatenate psd data and events:\n');
    
    fileseek = 1;
    for fId = 1:numfiles
        cnbiutil_disp_progress(fId, numfiles, '        ');
        cfilepath = filepaths{fId};
        cdata = load(cfilepath);
        if(sum(sum(sum(isnan(cdata.psd)))) ~=0)
            keyboard
        end
        % Get current position 
        cstart   = fileseek;
        cstop    = cstart + datasize(1, fId) - 1;
       
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
        Mk(cstart:cstop) = modality;
        
        % Get day from filename
        if strcmpi(cinfo.date, lday) == false
            nday = nday + 1;
            Dl = cat(1, Dl, cinfo.date);
            lday = cinfo.date;
        end
        Dk(cstart:cstop) = nday;
        
        % Create run vector
        Rk(cstart:cstop) = fId;
        Rl{fId} = cinfo.extra;
        
        % Concatenate features along 1st dimension (samples)
        F(cstart:cstop, :, :) = cdata.psd;
        
        % Get events
        cevents = cdata.events;
        
        % Concatenate events
        if modality == 2 || modality == 3    % race runs
            tTYP = cat(1, tTYP, cevents.extra.trl.TYP);
            tPOS = cat(1, tPOS, cevents.extra.trl.POS + fileseek - 1);
            tDUR = cat(1, tDUR, cevents.extra.trl.DUR);

            pTYP = cat(1, pTYP, cevents.extra.pad.TYP);
            pPOS = cat(1, pPOS, cevents.extra.pad.POS + fileseek - 1);
            pDUR = cat(1, pDUR, cevents.extra.pad.DUR);

            prTYP = cat(1, prTYP, cevents.extra.protocol.TYP);
            prPOS = cat(1, prPOS, cevents.extra.protocol.POS + fileseek - 1);
            prDUR = cat(1, prDUR, cevents.extra.protocol.DUR);
            
            bTYP = cat(1, bTYP, cevents.extra.bci.TYP);
            bPOS = cat(1, bPOS, cevents.extra.bci.POS + fileseek - 1);
            bDUR = cat(1, bDUR, cevents.extra.bci.DUR);

            cTYP = cat(1, cTYP, cevents.extra.cmd.TYP);
            cPOS = cat(1, cPOS, cevents.extra.cmd.POS + fileseek - 1);
            cDUR = cat(1, cDUR, cevents.extra.cmd.DUR);

            cgTYP = cat(1, cgTYP, cevents.extra.cmdg.TYP);
            cgPOS = cat(1, cgPOS, cevents.extra.cmdg.POS + fileseek - 1);
            cgDUR = cat(1, cgDUR, cevents.extra.cmdg.DUR);            
            
            eTYP = cat(1, eTYP, cevents.extra.eye.TYP);
            ePOS = cat(1, ePOS, cevents.extra.eye.POS + fileseek - 1);
            eDUR = cat(1, eDUR, cevents.extra.eye.DUR);

            rTYP = cat(1, rTYP, cevents.extra.race.TYP);
            rPOS = cat(1, rPOS, cevents.extra.race.POS + fileseek - 1);
            rDUR = cat(1, rDUR, cevents.extra.race.DUR);
        end
        
        TYP = cat(1, TYP, cevents.TYP);
        DUR = cat(1, DUR, cevents.DUR);
        POS = cat(1, POS, cevents.POS + fileseek - 1);
        if (sum(POS <= 0) > 0)
            error('chk:pos', 'Positions < 0');
        end
        
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
        csettings.protocol = nan;
        
        if(isempty(settings))
            settings = csettings;
        end
        
        if (isequaln(csettings, settings) == false)
            error('chk:stn', ['Different processing settings for file: ' cfilepath]);
        end
        
        % Classifiers concatenation
        if isempty(cdata.classifier) == false
            cclassifier = cdata.classifier;
            if isempty(classifiers) == true
                classifiers = cat(2, classifiers, cclassifier);
            elseif isequal(cclassifier.filename, classifiers(end).filename) == false
                
                classifiers = cat(2, classifiers, cclassifier);
            end
            Xk(cstart:cstop) = length(classifiers);
        else
            Xk(cstart:cstop) = 0;
        end
            
        % Update the fileseek position
        fileseek = cstop + 1;
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
    
    events.extra.cmdg.TYP = cgTYP;
    events.extra.cmdg.POS = cgPOS;
    events.extra.cmdg.DUR = cgDUR;    
    
    events.extra.eye.TYP = eTYP;
    events.extra.eye.POS = ePOS;
    events.extra.eye.DUR = eDUR;
    
    events.extra.race.TYP = rTYP;
    events.extra.race.POS = rPOS;
    events.extra.race.DUR = rDUR;
    
    events.extra.protocol.TYP = prTYP;
    events.extra.protocol.POS = prPOS;
    events.extra.protocol.DUR = prDUR;    
    
    labels.Rk  = Rk;
    labels.Rl  = Rl;
    labels.Mk  = Mk;
    labels.Dl  = Dl;
    labels.Dk  = Dk;
    labels.Xk  = Xk;

end

function dsizes = get_data_size(filepaths)

    nfiles = length(filepaths);
    ndimensions = 3;                            % samples x freqs x chans
    
    dsizes = zeros(ndimensions, nfiles);
    
    for fId = 1:nfiles
        cfilepath = filepaths{fId};
        cinfo = whos('-file', cfilepath, 'psd');
        dsizes(:, fId) = cinfo.size;    
    end

end