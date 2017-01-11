function [] = analyzeRace(FilePath, lap)

try
    load(FilePath);
    header.EVENT = Race.EVENT;
    data = Race.data;
    
catch
    disp(['Corrupted file, skipping: ' FilePath]);
    return
end

% Remove overall DC
data = data-repmat(mean(data),size(data,1),1);

% Laplacian spatial filtering
data = data(:,1:16);
data = data*lap;

IndPadTrig = find(ismember(header.EVENT.TYP,[768 771 773 783]));
IndComm771 = find(header.EVENT.TYP == 25347);
IndComm773 = find(header.EVENT.TYP == 25349);

usetr = 0;
for tr=1:length(IndPadTrig)
    if(header.EVENT.TYP(IndPadTrig(tr))==771)
        
        if(tr<length(IndPadTrig))
            NextPad = IndPadTrig(tr+1);
        else
            NextPad = find(header.EVENT.TYP==666);
        end
        Next771 = find(IndComm771 > IndPadTrig(tr));
        if(isempty(Next771))
            Next771 = NaN;
        else
            Next771 = IndComm771(Next771(1));
        end
        usetr = usetr + 1;
        trials(usetr,1) = header.EVENT.POS(IndPadTrig(tr));
        trials(usetr,2) = header.EVENT.POS(min(NextPad,Next771));
        labels(usetr) = 771;
        
    elseif(header.EVENT.TYP(IndPadTrig(tr))==773)
        if(tr<length(IndPadTrig))
            NextPad = IndPadTrig(tr+1);
        else
            NextPad = find(header.EVENT.TYP==666);
        end
        Next773 = find(IndComm773 > IndPadTrig(tr));
        if(isempty(Next773))
            Next773 = NaN;
        else
            Next773 = IndComm773(Next773(1));
        end
        
        usetr = usetr + 1;
        trials(usetr,1) = header.EVENT.POS(IndPadTrig(tr));
        trials(usetr,2) = header.EVENT.POS(min(NextPad,Next773));
        labels(usetr) = 773;
    else
        disp('Ignore');
    end
    
end


% Force classes available
ca = [771 773];
desclass = [771 773];


% Useful params for PSD extraction with the fast algorithm
psdshift = 512*0.5*0.5;
winshift = 512*0.0625;

if((mod(psdshift,winshift) ~=0) && (mod(winshift,psdshift) ~=0))
    disp(['[eegc3_smr_simloop_fast] The fast PSD method cannot be applied with the current settings!']);
    disp(['[eegc3_smr_simloop_fast] The internal welch window shift must be a multiple of the overall feature window shift (or vice versa)!']);
    return;
end

% Create arguments for spectrogram
spec_win = 512*0.5;
% Careful here: The overlapping depends on whether the winshift or the
% psdshift is smaller. Some calculated internal windows will be redundant,
% but the speed is much faster anyway

if(psdshift <= winshift)
    spec_ovl = spec_win - psdshift;
else
    spec_ovl = spec_win - winshift;
end

% Calculate all the internal PSD windows
for ch=1:16
    %disp(['[eegc3_smr_simloop_fast] Internal PSDs on electrode ' num2str(ch)]);
    [~,f,t,p(:,:,ch)] = spectrogram(data(:,ch), spec_win, spec_ovl, [], 512);
end

% Keep only desired frequencies
freqs = [8:2:30];
p = p(find(ismember(f,freqs)),:,:);

% Setup moving average filter parameters
FiltA = 1;
if(winshift >= psdshift)
    % Case where internal windows are shifted according to psdshift
    MAsize = (512*1.00)/psdshift - 1;   
    FiltB = (1/MAsize)*ones(1,MAsize);
    MAstep = winshift/psdshift;
else
    % Case where internal windows are shifted according to winshift
    FiltStep = psdshift/winshift;
    MAsize = (512*1.00)/winshift - FiltStep;   
    FiltB = zeros(1,MAsize);
    FiltB(1:FiltStep:end-1) = 1;
    FiltB = FiltB/sum(FiltB);
    MAstep = 1;
end

StartInd = find(FiltB~=0);
StartInd = StartInd(end);

afeats = filter(FiltB,FiltA,p,[],2);
afeats = permute(afeats, [2 1 3]);

% Get rid of initial filter byproducts
afeats = afeats(StartInd:end,:,:);

% In case of psdshift, there will be redundant windows. Remove them
if(MAstep > 1)
   afeats = afeats(1:MAstep:end,:,:);
end

% Take the log as final feature values
afeats = log(afeats);

% Remap trials to PSD space
ftrials(:,1) = floor(trials(:,1)/32)-15;
ftrials(:,2) = floor(trials(:,2)/32)-15;

slabels = zeros(size(afeats,1),1);
for i=1:size(ftrials,1)
    slabels(ftrials(i,1):ftrials(i,2)) = labels(i);
end

nf = length(freqs);

D{1} = afeats(slabels==desclass(1),:,:);
L{1} = slabels(slabels==desclass(1),:,:);

D{2} = afeats(slabels==desclass(2),:,:);
L{2} = slabels(slabels==desclass(2),:,:);

% DP maps of classx vs class y
data = [D{1}; D{2}];
data = reshape(data,size(data,1),nf*16);
lbl = [L{1}; L{2}];
%[COM,PWGR,V,vp,DISC]=cva_tun_opt(data,lbl);
FS = eegc3_fs2(data(lbl==desclass(1),:),data(lbl==desclass(2),:));
DPA = reshape(FS,nf,16)';

imagesc(DPA,[0 1.0]);