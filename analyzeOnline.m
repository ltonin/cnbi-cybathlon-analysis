function [DPA DPS DPE MDPA MDPS MDPE MA CA MS CS ME CE UsedFeat...
    MDA MDS MDE SAccA SAccS SAccE CAccA CAccS CAccE task matfile ...
    SimAccWithINC06 SimAccWithINC07] = ...
    analyzeOnline(FilePath, task, tNA, fs, lap, MATPath)

try
    [data, header] = sload(FilePath);
    
    % Check of there is 772 (both feet) and convert to 770 (right hand)
    if(ismember(772,header.EVENT.TYP))
        header.EVENT.TYP(find(header.EVENT.TYP==772))=770;
    end
    % Check of there is 781, if not (Jerome, c1) add it 1 sec after cue
    if(~ismember(781,header.EVENT.TYP))
        i=1;
        while(i<=length(header.EVENT.TYP))
            if((header.EVENT.TYP(i)==769) || (header.EVENT.TYP(i)==770) || ...
                    (header.EVENT.TYP(i)==771))
                header.EVENT.TYP = [header.EVENT.TYP(1:i); 781; header.EVENT.TYP(i+1:end)];
                header.EVENT.POS = [header.EVENT.POS(1:i); header.EVENT.POS(i)+512 ; header.EVENT.POS(i+1:end)];
            end
            i=i+1;
        end
    end
catch
    disp(['Corrupted file, skipping: ' FilePath]);
    DPA = [];
    DPS = [];
    DPE = [];
    MDPA = [];
    MDPS = [];
    MDPE = [];
    MA = [];
    CA = [];
    MS = [];
    CS = [];
    ME = [];
    CE = [];  
    MDA = [];    
    MDS = [];    
    MDE = [];  
    SAccA = [];
    SAccS = [];
    SAccE = [];
    CAccA = [];
    CAccS = [];
    CAccE = [];    
    UsedFeat = [];
    matfile = [];
    return
end

% Remove overall DC
data = data-repmat(mean(data),size(data,1),1);

% Laplacian spatial filtering
data = data(:,1:16);
data = laplacianSP(data,lap);

% Trial extraction
pos = header.EVENT.POS;
cf = find(header.EVENT.TYP==781);
cue = cf-1;
tend = cf+1;

% Correction for runs with trigger problems
check = tend;
if((sum(header.EVENT.TYP(check(1:end-1))==897) + sum(header.EVENT.TYP(check(1:end-1))==898))==0)
    tend = tend+1;
end
% Extract trials, no-actiity and activity intervals
%trials = [pos(cue)+1-tNA*fs pos(cue) pos(cf)+1 pos(cf)+tA*fs];
try
    trials = [pos(cue)-tNA*fs pos(cue) pos(cf)-fs pos(tend)];
catch
    
    % Last trial corrupted, remove it
    cue = cue(1:end-1);
    cf = cf(1:end-1);
    tend = tend(1:end-1);
    trials = [pos(cue)-tNA*fs pos(cue) pos(cf)-fs pos(tend)];
end

if(strcmp(task,'RHLH'))
    desclass = [770 769];
elseif(strcmp(task,'RHBF'))
    desclass = [770 771];    
elseif(strcmp(task,'BFLH'))
    desclass = [771 769];    
end

% Find classes available
ca = unique(header.EVENT.TYP(cue));

if(sum(ismember(ca,desclass))~=2)
    disp(['One or both classes missing, skipping: ' FilePath]);
    DPA = [];
    DPS = [];
    DPE = [];
    MDPA = [];
    MDPS = [];
    MDPE = [];
    MA = [];
    CA = [];
    MS = [];
    CS = [];
    ME = [];
    CE = [];   
    MDA = [];    
    MDS = [];    
    MDE = [];   
    SAccA = [];
    SAccS = [];
    SAccE = [];  
    CAccA = [];
    CAccS = [];
    CAccE = [];        
    UsedFeat = [];    
    matfile = [];    
    return;    
end

labels = zeros(1,length(cue));
for c=1:length(ca)
    labels(header.EVENT.TYP(cue) == ca(c)) = ca(c);
    N(c) = sum(labels==ca(c));
end

% Exclude runs with too few trials
if((N(1) < 5) || (N(2) < 5))
    disp(['Too few trials, skipping: ' FilePath]);
    DPA = [];
    DPS = [];
    DPE = [];
    MDPA = [];
    MDPS = [];
    MDPE = [];
    MA = [];
    CA = [];
    MS = [];
    CS = [];
    ME = [];
    CE = [];   
    MDA = [];    
    MDS = [];    
    MDE = [];   
    SAccA = [];
    SAccS = [];
    SAccE = [];  
    CAccA = [];
    CAccS = [];
    CAccE = [];        
    UsedFeat = [];    
    matfile = [];    
    return;        
end

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
ftrials(:,1) = floor(trials(:,1)/32)+1;
ftrials(:,3) = floor(trials(:,3)/32)+1;
ftrials(:,2) = floor(trials(:,2)/32)-15;
ftrials(:,4) = floor(trials(:,4)/32)-15;


slabels = zeros(size(afeats,1),1);
for i=1:size(ftrials,1)
    slabels(ftrials(i,1):ftrials(i,2)) = -1;
    slabels(ftrials(i,3):ftrials(i,4)) = labels(i);
end

nf = length(freqs);

BD = afeats(slabels==-1,:,:);
BL = slabels(slabels==-1,:,:);



BDs = BD(1:floor(size(BD,1)/2),:,:);
BDe = BD(floor(size(BD,1)/2)+1:end,:,:);
BDa = [BDs;BDe];

MA{3} = mean(BDa(:,:));
CA{3} = cov(BDa(:,:));
MS{3} = mean(BDs(:,:));
CS{3} = cov(BDs(:,:));    
ME{3} = mean(BDe(:,:));
CE{3} = cov(BDe(:,:));

if(sum(ismember(ca,desclass(1)))==1)
    
    D{1} = afeats(slabels==desclass(1),:,:);
    L{1} = slabels(slabels==desclass(1),:,:);
    
        
    D1s = D{1}(1:floor(size(D{1},1)/2),:,:);
    D1e = D{1}(floor(size(D{1},1)/2)+1:end,:,:);
    D1a = [D1s;D1e];

    MA{1} = mean(D1a(:,:));
    CA{1} = cov(D1a(:,:));
    MS{1} = mean(D1s(:,:));
    CS{1} = cov(D1s(:,:));    
    ME{1} = mean(D1e(:,:));
    CE{1} = cov(D1e(:,:));

else
    MA{1} = [];
    CA{1} = [];
    MS{1} = [];
    CS{1} = [];
    ME{1} = [];
    CE{1} = [];    
end

if(sum(ismember(ca,desclass(2)))==1)
    D{2} = afeats(slabels==desclass(2),:,:);
    L{2} = slabels(slabels==desclass(2),:,:);

    D2s = D{2}(1:floor(size(D{2},1)/2),:,:);
    D2e = D{2}(floor(size(D{2},1)/2)+1:end,:,:);
    D2a = [D2s;D2e];

    MA{2} = mean(D2a(:,:));
    CA{2} = cov(D2a(:,:));
    MS{2} = mean(D2s(:,:));
    CS{2} = cov(D2s(:,:));    
    ME{2} = mean(D2e(:,:));
    CE{2} = cov(D2e(:,:));

else
    MA{2} = [];
    CA{2} = [];
    MS{2} = [];
    CS{2} = [];
    ME{2} = [];
    CE{2} = [];    
end

% DP maps of classx vs class y
if(sum(ismember(ca,desclass))==2)
    
    data = [D{1}; D{2}];
    data = reshape(data,size(data,1),nf*16);
    lbl = [L{1}; L{2}];
    %[COM,PWGR,V,vp,DISC]=cva_tun_opt(data,lbl);
    FS = eegc3_fs2(data(lbl==desclass(1),:),data(lbl==desclass(2),:));
    DPA = reshape(FS,nf,16)';
    
    data = [D{1}(1:floor(size(D{1},1)/2),:,:); D{2}(1:floor(size(D{2},1)/2),:,:)];
    data = reshape(data,size(data,1),nf*16);
    lbl = [L{1}(1:floor(size(L{1},1)/2)); L{2}(1:floor(size(L{2},1)/2))];
    %[COM,PWGR,V,vp,DISC]=cva_tun_opt(data,lbl);
    FS = eegc3_fs2(data(lbl==desclass(1),:),data(lbl==desclass(2),:));    
    DPS = reshape(FS,nf,16)';    
    
    data = [D{1}(floor(size(D{1},1)/2)+1:end,:,:); D{2}(floor(size(D{2},1)/2)+1:end,:,:)];
    data = reshape(data,size(data,1),nf*16);
    lbl = [L{1}(floor(size(L{1},1)/2)+1:end); L{2}(floor(size(L{2},1)/2)+1:end)];
    %[COM,PWGR,V,vp,DISC]=cva_tun_opt(data,lbl);
    FS = eegc3_fs2(data(lbl==desclass(1),:),data(lbl==desclass(2),:));        
    DPE = reshape(FS,nf,16)';    
    
else
    DPA = [];
    DPS = [];
    DPE = [];
end

% DP maps of classx vs baseline
if(sum(ismember(ca,desclass(1)))==1)
    
    data = [D{1}; BD];
    data = reshape(data,size(data,1),nf*16);
    lbl = [L{1}; BL];
    %[COM,PWGR,V,vp,DISC]=cva_tun_opt(data,lbl);
    FS = eegc3_fs2(data(lbl==desclass(1),:),data(lbl==-1,:));                
    MDPA{1} = reshape(FS,nf,16)';
    
    data = [D{1}(1:floor(size(D{1},1)/2),:,:); BD(1:floor(size(BD,1)/2),:,:)];
    data = reshape(data,size(data,1),nf*16);
    lbl = [L{1}(1:floor(size(L{1},1)/2)); BL(1:floor(size(BL,1)/2))];
    %[COM,PWGR,V,vp,DISC]=cva_tun_opt(data,lbl);
    FS = eegc3_fs2(data(lbl==desclass(1),:),data(lbl==-1,:));                    
    MDPS{1} = reshape(FS,nf,16)';    
    
    data = [D{1}(floor(size(D{1},1)/2)+1:end,:,:); BD(floor(size(BD,1)/2)+1:end,:,:)];
    data = reshape(data,size(data,1),nf*16);
    lbl = [L{1}(floor(size(L{1},1)/2)+1:end); BL(floor(size(BL,1)/2)+1:end)];
    %[COM,PWGR,V,vp,DISC]=cva_tun_opt(data,lbl);
    FS = eegc3_fs2(data(lbl==desclass(1),:),data(lbl==-1,:));                    
    MDPE{1} = reshape(FS,nf,16)';    
else
    MDPA{1} = [];
    MDPS{1} = [];
    MDPE{1} = [];
end

if(sum(ismember(ca,desclass(2)))==1)

    data = [D{2}; BD];
    data = reshape(data,size(data,1),nf*16);
    lbl = [L{2}; BL];
    %[COM,PWGR,V,vp,DISC]=cva_tun_opt(data,lbl);
    FS = eegc3_fs2(data(lbl==desclass(2),:),data(lbl==-1,:));                    
    MDPA{2} = reshape(FS,nf,16)';
    
    data = [D{2}(1:floor(size(D{2},1)/2),:,:); BD(1:floor(size(BD,1)/2),:,:)];
    data = reshape(data,size(data,1),nf*16);
    lbl = [L{2}(1:floor(size(L{2},1)/2)); BL(1:floor(size(BL,1)/2))];
    %[COM,PWGR,V,vp,DISC]=cva_tun_opt(data,lbl);
    FS = eegc3_fs2(data(lbl==desclass(2),:),data(lbl==-1,:));                        
    MDPS{2} = reshape(FS,nf,16)';    
    
    data = [D{2}(floor(size(D{2},1)/2)+1:end,:,:); BD(floor(size(BD,1)/2)+1:end,:,:)];
    data = reshape(data,size(data,1),nf*16);
    lbl = [L{2}(floor(size(L{2},1)/2)+1:end); BL(floor(size(BL,1)/2)+1:end)];
    %[COM,PWGR,V,vp,DISC]=cva_tun_opt(data,lbl);
    FS = eegc3_fs2(data(lbl==desclass(2),:),data(lbl==-1,:));                        
    MDPE{2} = reshape(FS,nf,16)';  
else
    MDPA{2} = [];
    MDPS{2} = [];
    MDPE{2} = [];
end

% Classifier dependent metrics

% Load classifier used
analysis = load(MATPath);
analysis = analysis.analysis;

matfile = analysis;

% Find if normalization was used
if(sum(abs(analysis.tools.net.gau.M(:))>1) == 0)
    Norm = 1;
else
    Norm = 0;
end

% Find features used, exclude freqs lower than 8 Hz
UsedFeat = [];
fInd = 0;
remInd = [];
for ch=1:length(analysis.tools.features.channels)
    for fr=1:length(analysis.tools.features.bands{analysis.tools.features.channels(ch)})
        fInd = fInd + 1;
        if((analysis.tools.features.bands{analysis.tools.features.channels(ch)}(fr) >= 8) && ...
                (analysis.tools.features.bands{analysis.tools.features.channels(ch)}(fr) <= 30))
            UsedFeat = [UsedFeat ; analysis.tools.features.channels(ch) ...
                analysis.tools.features.bands{analysis.tools.features.channels(ch)}(fr)];
        else
            remInd = [remInd; fInd];
        end
    end
end

% Remove bad features from used classifier
if(~isempty(remInd))
    analysis.tools.net.gau.M(:,:,remInd) = [];
    analysis.tools.net.gau.C(:,:,remInd) = [];
end

% Crop all data to selected features
sfeats = zeros(size(afeats,1),size(UsedFeat,1));
for s=1:size(afeats,1)
    for f=1:size(UsedFeat,1)
        % Features are in 192d space, not 368d
        sfeats(s,f) = afeats(s,(UsedFeat(f,2)-8)/2+1,UsedFeat(f,1));
    end
end

if(Norm)
    sfeats = eegc3_normalize(sfeats);
end


% Average Mahalanobis distance between used classifier and generated samples
if(sum(ismember(ca,desclass(1)))==1)
    
    D{1} = sfeats(slabels==desclass(1),:);
    L{1} = ones(length(slabels(slabels==desclass(1))),1);
      
    MDA{1} = avgMD(analysis,D{1},1);
    MDS{1} = avgMD(analysis,D{1}(1:floor(size(D{1},1)/2),:),1);
    MDE{1} = avgMD(analysis,D{1}(floor(size(D{1},1)/2)+1:end,:),1);

else
    MDA{1} = [];
    MDS{1} = [];
    MDE{1} = [];
end   

if(sum(ismember(ca,desclass(2)))==1)
    
    D{2} = sfeats(slabels==desclass(2),:);
    L{2} = 2*ones(length(slabels(slabels==desclass(2))),1);
      
    MDA{2} = avgMD(analysis,D{2},2);
    MDS{2} = avgMD(analysis,D{2}(1:floor(size(D{2},1)/2),:),2);
    MDE{2} = avgMD(analysis,D{2}(floor(size(D{2},1)/2)+1:end,:),2);

else
    MDA{2} = [];
    MDS{2} = [];
    MDE{2} = [];
end   

% Single-sample accuracies
cdataA = [D{1} ; D{2}];
clblA = [L{1} ; L{2}];
[conf_matrix, perf] = gauEval(analysis.tools.net.gau.M, analysis.tools.net.gau.C,...
    [cdataA clblA], 0.5);
SAccA = 100*sum(diag(conf_matrix(1:2,1:2)))/sum(sum(conf_matrix(1:2,1:2)));


cdataS = [D{1}(1:floor(size(D{1},1)/2),:); D{2}(1:floor(size(D{2},1)/2),:)];
clblS = [L{1}(1:floor(size(L{1},1)/2)) ; L{2}(1:floor(size(L{2},1)/2))];
[conf_matrix, perf] = gauEval(analysis.tools.net.gau.M, analysis.tools.net.gau.C,...
    [cdataS clblS], 0.5);
SAccS = 100*sum(diag(conf_matrix(1:2,1:2)))/sum(sum(conf_matrix(1:2,1:2)));


cdataE = [D{1}(floor(size(D{1},1)/2)+1:end,:); D{2}(floor(size(D{2},1)/2)+1:end,:)];
clblE = [L{1}(floor(size(L{1},1)/2)+1:end) ; L{2}(floor(size(L{2},1)/2)+1:end)];
[conf_matrix, perf] = gauEval(analysis.tools.net.gau.M, analysis.tools.net.gau.C,...
    [cdataE clblE], 0.5);
SAccE = 100*sum(diag(conf_matrix(1:2,1:2)))/sum(sum(conf_matrix(1:2,1:2)));

slabels(slabels==-1)=3;
BD = sfeats(slabels==3,:);
BL = 3*ones(length(slabels(slabels==3)),1);
cdataA = [D{1} ; D{2} ; BD];
clblA = [L{1} ; L{2} ; BL];
for sample=1:size(cdataA,1)
    [mplamplou prob(sample,:)] = ...
        gauClassifier(analysis.tools.net.gau.M, analysis.tools.net.gau.C,cdataA(sample,:));
end

% Single-sample accuracies with Rejection 0.6
[maxV06 classV06] = max(prob');
classV06(find(maxV06<0.6))=3;
decINC06 = classV06;
[CMClass06 CMAll06 SimAccWithINC06 Error] = ...
    eegc3_confusion_matrix(clblA, decINC06);


% Single-sample accuracies with Rejection 0.7
[maxV07 classV07] = max(prob');
classV07(find(maxV07<0.7))=3;
decINC07 = classV07;
[CMClass07 CMAll07 SimAccWithINC07 Error] = ...
    eegc3_confusion_matrix(clblA, decINC07);
disp('a');

% Command Accuracies (for now, as happened online)
CommCue = header.EVENT.TYP(cue);
CommRes = header.EVENT.TYP(tend);
comm1A = find(CommCue==desclass(1));
comm2A = find(CommCue==desclass(2));
comm1S = comm1A(1:floor(length(comm1A)/2));
comm2S = comm2A(1:floor(length(comm2A)/2));
comm1E = comm1A(floor(length(comm1A)/2)+1:end);
comm2E = comm2A(floor(length(comm2A)/2)+1:end);
CAccA{1} = 100*sum(CommRes(comm1A)==897)/length(comm1A);
CAccA{2} = 100*sum(CommRes(comm2A)==897)/length(comm2A);
CAccS{1} = 100*sum(CommRes(comm1S)==897)/length(comm1S);
CAccS{2} = 100*sum(CommRes(comm2S)==897)/length(comm2S);
CAccE{1} = 100*sum(CommRes(comm1E)==897)/length(comm1E);
CAccE{2} = 100*sum(CommRes(comm2E)==897)/length(comm2E);
CAccA{3} = (length(comm1A)*CAccA{1} + length(comm2A)*CAccA{2})/(length(comm1A)+length(comm2A));
CAccS{3} = (length(comm1S)*CAccS{1} + length(comm2S)*CAccS{2})/(length(comm1S)+length(comm2S));
CAccE{3} = (length(comm1E)*CAccE{1} + length(comm2E)*CAccE{2})/(length(comm1E)+length(comm2E));