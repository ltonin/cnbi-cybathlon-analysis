clearvars; clc;

subject = 'MA25VE';

pattern     = '.mi.';

experiment  = 'cybathlon';
datapath    = ['/mnt/data/Research/' experiment '/' subject '/'];
savedir     = '/analysis/';

%% Processing parameters
wlength    = 0.5;
pshift     = 0.25;                  
wshift     = 0.0625;                
selfreqs   = 4:2:48;
selchans   = 1:16;                  % <-- Needed for the 2-amplifiers setup
load('lapmask_16ch.mat');           % <-- To be checked if it is the correct one

%% Get datafiles
[Files, NumFiles] = cnbiutil_getdata(datapath, subject, pattern, '.gdf');

%% Create/Check for savepath
[~, savepath] = util_mkdir(pwd, savedir);

%% Processing files

for fId = 1:NumFiles
    cfilename = Files{fId};
    [~, ~, cextension] = fileparts(cfilename);
    cnbiutil_bdisp(['[io] - Loading file ' num2str(fId) '/' num2str(NumFiles)]);
    disp(['       File: ' cfilename]);
    
    % Importing gdf file
    try
        if(strcmp(cextension, '.gdf'))
            [s, h] = sload(cfilename);
        elseif(strcmp(cextension, '.mat'))
            cdata = load(cfilename);
            s = cdata.data;
            h = cdata.header;
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
    s_lap = s_dc*lapmask;
    
    % Compute spectrogram
    [psd, freqgrid] = cnbiproc_spectrogram(s_lap, wlength, wshift, pshift, h.SampleRate);
    
    % Selecting desired frequencies
    [freqs, idfreqs] = intersect(freqgrid, selfreqs);
    psd = psd(:, idfreqs, :);
    
    % Resample events
    events.TYP = h.EVENT.TYP;
    events.POS = floor(h.EVENT.POS/(wshift*h.SampleRate)) + 1;
    events.DUR = floor(h.EVENT.DUR/(wshift*h.SampleRate)) + 1;
    
    % Create settings structure
    settings.data.filename          = cfilename;
    settings.data.nsamples          = size(s, 1);
    settings.data.nchannels         = size(s, 2);
    settings.data.samplerate        = h.SampleRate;
    settings.spatial.laplacian      = lapmask;
    settings.spectrogram.wlength    = wlength;
    settings.spectrogram.wshift     = wshift;
    settings.spectrogram.pshift     = pshift;
    settings.spectrogram.mavglength = mavglength;
    settings.spectrogram.freqgrid   = freqgrid;
    
    [~, name] = fileparts(cfilename);
    sfilename = [savepath '/' name '.mat'];
    cnbiutil_bdisp(['[out] - Saving psd in: ' sfilename]);
    save(sfilename, 'psd', 'freqs', 'events', 'settings'); 
end
