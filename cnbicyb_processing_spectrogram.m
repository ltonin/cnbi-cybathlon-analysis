clearvars; clc;

subject = 'AN14VE';

% identifiers = {'.*line.mi.', '.gdf'};
identifiers = {'.race.mi.',  '.mat'};

pattern     = identifiers{1};
extension   = identifiers{2};
experiment  = 'cybathlon';
datapath    = ['/mnt/data/Research/' experiment '/' subject '/'];
savedir     = '/analysis/';

%% Processing parameters
mlength    = 1;
wlength    = 0.5;
pshift     = 0.25;                  
wshift     = 0.0625;                
selfreqs   = 4:2:48;
selchans   = 1:16;                  % <-- Needed for the 2-amplifiers setup
load('extra/laplacian16.mat');              % <-- To be checked if it is the correct one

winconv = 'backward';               % Type of conversion for events from samples to psd windows

%% Get datafiles
[Files, NumFiles] = cnbiutil_getdata(datapath, subject, pattern, extension);

%% Create/Check for savepath
[~, savepath] = cnbiutil_mkdir(pwd, savedir);

%% Processing files

for fId = 1:NumFiles
    cfilename = Files{fId};
    [~, ~, cextension] = fileparts(cfilename);
    cnbiutil_bdisp(['[io] - Loading file ' num2str(fId) '/' num2str(NumFiles)]);
    disp(['       File: ' cfilename]);
    
    % Get information from filename
    cinfo = cnbiutil_getfile_info(cfilename);
    
    % Importing gdf file
    try
        if(strcmp(cextension, '.gdf'))
            [s, h] = sload(cfilename);
        elseif(strcmp(cextension, '.mat'))
            cdata = load(cfilename);
            s = cdata.Race.data;
            if(isnan(s))
                cnbiutil_bdisp(['[io] - Corrupted file (data is NaN), skipping: ' cfilename]);
                continue;
            end
            h.EVENT = cdata.Race.EVENT;
            h.SampleRate = 512;
        else
            error('chk:ext', 'Unknown extension');
        end
            
    catch 
        cnbiutil_bdisp(['[io] - Corrupted file, skipping: ' cfilename]);
        continue;
    end
    
    if isempty(s) 
        cnbiutil_bdisp(['[io] - Corrupted file, skipping: ' cfilename]);
        continue;
    end
    
    s = s(:, selchans);         
    
    % Computed DC removal
    s_dc = s-repmat(mean(s),size(s,1),1);
    
    % Compute Spatial filter
    s_lap = s_dc*lap;
    
    % Compute spectrogram
    [psd, freqgrid] = cnbiproc_spectrogram(s_lap, wlength, wshift, pshift, h.SampleRate, mlength);
    
    % Selecting desired frequencies
    [freqs, idfreqs] = intersect(freqgrid, selfreqs);
    psd = psd(:, idfreqs, :);
    
    % If it is a race, then convert events to standard format (TYP, POS, DUR)
    % and resampling events
    cevents     = h.EVENT;
    cextraevents = [];
    
    if strcmpi(cinfo.modality, 'race')
        disp('       Converting EVENT table from race file');
        [cevents, cextraevents] = cnbiproc_extract_event_race(h.EVENT);
    end
    
    events.TYP = cevents.TYP;
    events.POS = cnbiproc_pos2win(cevents.POS, wshift*h.SampleRate, winconv, mlength*h.SampleRate);
    events.DUR = floor(cevents.DUR/(wshift*h.SampleRate)) + 1;
    events.conversion = winconv;
    
    if isempty(cextraevents) == false
        events.extra.trl.TYP  = cextraevents.trl.TYP;
        events.extra.pad.TYP  = cextraevents.pad.TYP;
        events.extra.bci.TYP  = cextraevents.bci.TYP;
        events.extra.cmd.TYP  = cextraevents.cmd.TYP;
        events.extra.eye.TYP  = cextraevents.eye.TYP;
        events.extra.race.TYP = cextraevents.race.TYP;
        events.extra.trl.POS  = cnbiproc_pos2win(cextraevents.trl.POS,  wshift*h.SampleRate, winconv, mlength*h.SampleRate);
        events.extra.pad.POS  = cnbiproc_pos2win(cextraevents.pad.POS,  wshift*h.SampleRate, winconv, mlength*h.SampleRate);
        events.extra.bci.POS  = cnbiproc_pos2win(cextraevents.bci.POS,  wshift*h.SampleRate, winconv, mlength*h.SampleRate);
        events.extra.cmd.POS  = cnbiproc_pos2win(cextraevents.cmd.POS,  wshift*h.SampleRate, winconv, mlength*h.SampleRate);
        events.extra.eye.POS  = cnbiproc_pos2win(cextraevents.eye.POS,  wshift*h.SampleRate, winconv, mlength*h.SampleRate);
        events.extra.race.POS = cnbiproc_pos2win(cextraevents.race.POS, wshift*h.SampleRate, winconv, mlength*h.SampleRate);
        events.extra.trl.DUR  = floor(cextraevents.trl.DUR/(wshift*h.SampleRate)) + 1;
        events.extra.pad.DUR  = floor(cextraevents.pad.DUR/(wshift*h.SampleRate)) + 1;
        events.extra.bci.DUR  = floor(cextraevents.bci.DUR/(wshift*h.SampleRate)) + 1;
        events.extra.cmd.DUR  = floor(cextraevents.cmd.DUR/(wshift*h.SampleRate)) + 1;
        events.extra.eye.DUR  = floor(cextraevents.eye.DUR/(wshift*h.SampleRate)) + 1;
        events.extra.race.DUR = floor(cextraevents.race.DUR/(wshift*h.SampleRate)) + 1;
    end
    
    % Decided to stop recovering this info (with Simis)
%     % Get classifiers from log file 
%     if strcmpi(cinfo.modality, 'online') || strcmpi(cinfo.modality, 'race')
%         clogfile = [datapath '/' cinfo.subject '_' cinfo.date '/' cinfo.subject '.' cinfo.date '.log'];
%         [~, cfile, cext] = fileparts(cfilename);
%         ctarget = [cfile cext];
%         clogstr = cnbiutil_read_logfile(clogfile, ctarget(1:20));
%        keyboard 
%     end
    
    % Create settings structure
    settings.data.filename          = cfilename;
    settings.data.nsamples          = size(s, 1);
    settings.data.nchannels         = size(s, 2);
    settings.data.samplerate        = h.SampleRate;
    settings.spatial.laplacian      = lap;
    settings.spectrogram.wlength    = wlength;
    settings.spectrogram.wshift     = wshift;
    settings.spectrogram.pshift     = pshift;
    settings.spectrogram.freqgrid   = freqs;
    
    [~, name] = fileparts(cfilename);
    sfilename = [savepath '/' name '.mat'];
    cnbiutil_bdisp(['[out] - Saving psd in: ' sfilename]);
    save(sfilename, 'psd', 'freqs', 'events', 'settings'); 
end
