clearvars; clc;

subject = 'AN14VE';

pattern     = '.mi.';

wlength    = 0.5;
pshift     = 0.25;                  % <-- What is it?
wshift     = 0.0625;                % <-- What is it?
mavglength = 1.0;
selfreqs   = 4:2:48;

load('lapmask_16ch.mat');

experiment  = 'cybathlon';
datapath    = ['/mnt/data/Research/' experiment '/' subject '/'];
savedir     = '/analysis/';

% Get datafiles
[Files, NumFiles] = cnbiutil_getdata(datapath, subject, pattern, '.gdf');

% Create/Check for savepath
[~, savepath] = util_mkdir(pwd, savedir);

%% Processing files

for fId = 1:NumFiles
    cfilename = Files{fId};
    cnbiutil_bdisp(['[io] - Loading file ' num2str(fId) '/' num2str(NumFiles)]);
    disp(['       File: ' cfilename]);
    
    % Importing gdf file
    [s, h] = sload(cfilename);
    s = s(:, 1:end-1);          % Assuming that the last channel is the trigger
    
    % Computed DC removal
    s_dc = s-repmat(mean(s),size(s,1),1);
    
    % Compute Spatial filter
    s_lap = s_dc*lapmask;
    
    % Compute spectrogram
    [psd, freqgrid] = cnbiproc_spectrogram(s_lap, wlength, wshift, pshift, h.SampleRate, mavglength);
    
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
    save(sfilename, 'psd', 'freqs', 'settings'); 
end
