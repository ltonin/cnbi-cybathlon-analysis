clearvars; clc;


subject = 'AN14VE';
%subject = 'MA25VE_RaceMat';

% identifiers = {'.offline.mi.', '.gdf'};
% identifiers = {'.online.mi.',  '.gdf'};
 identifiers = {'.race.mi.',    '.mat'};

pattern     = identifiers{1};
extension   = identifiers{2};
experiment  = 'cybathlon';
<<<<<<< HEAD
datapath    = ['/mnt/data/Research/' experiment '/' subject '/'];
% datapath    = ['/home/sperdikis/Desktop/tst/AN14VE'];
=======
%datapath    = ['/mnt/data/Research/' experiment '/' subject '/'];
%datapath    = ['/home/sperdikis/Desktop/tst/AN14VE'];
datapath    = ['/home/sperdikis/Desktop/tst/MA25VE'];
>>>>>>> cfc0b6a393340f5e487f5eb7491897aa51cfe392
savedir     = '/analysis/';

%% Processing parameters
wlength    = 0.5;
pshift     = 0.25;                  
wshift     = 0.0625;                
selfreqs   = 4:2:48;
selchans   = 1:16;                  % <-- Needed for the 2-amplifiers setup
load('extra/laplacian16.mat');              % <-- To be checked if it is the correct one

%% Get datafiles
%[Files, NumFiles] = cnbiutil_getdata(datapath, subject, pattern, extension);
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
    [psd, freqgrid] = cnbiproc_spectrogram(s_lap, wlength, wshift, pshift, h.SampleRate);
    
    % Selecting desired frequencies
    [freqs, idfreqs] = intersect(freqgrid, selfreqs);
    psd = psd(:, idfreqs, :);
    
    % Resample events
    events.TYP = h.EVENT.TYP;
    events.POS = floor(h.EVENT.POS/(wshift*h.SampleRate)) + 1;
    events.DUR = floor(h.EVENT.DUR/(wshift*h.SampleRate)) + 1;
    
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
