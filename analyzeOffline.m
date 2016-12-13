function [DPA DPS DPE MDPA MDPS MDPE MA CA MS CS ME CE SAcc] = ...
    analyzeOffline(FilePath, tNA, tA, fs, lap)

try
    [data, header] = sload(FilePath);
    
    % Check of there is 772 (both feet) and convert to 770 (right hand)
    if(ismember(772,header.EVENT.TYP))
        header.EVENT.TYP(find(header.EVENT.TYP==772))=770;
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
    SAcc = [];
    return
end

% Remove overall DC
data = data-repmat(mean(data),size(data,1),1);

% Laplacian spatial filtering
data = data(:,1:16);
data = laplacianSP(data,lap);

if(sum(ismember(unique(header.EVENT.TYP),781))==0)
    % Treatment of old offline
    % Insert "fake" 781s 1 second after the cues (which I have)
    nheader.EVENT.TYP = zeros(2*length(header.EVENT.TYP),1);
    nheader.EVENT.POS = zeros(2*length(header.EVENT.TYP),1);
    for i=1:length(header.EVENT.TYP)
        nheader.EVENT.TYP(2*i-1) = header.EVENT.TYP(i);
        nheader.EVENT.TYP(2*i) = 781;
        nheader.EVENT.POS(2*i-1) = header.EVENT.POS(i);
        nheader.EVENT.POS(2*i) = header.EVENT.POS(i)+513;
    end
    header.EVENT.TYP = nheader.EVENT.TYP;
    header.EVENT.POS = nheader.EVENT.POS;
end

% Trial extraction
pos = header.EVENT.POS;
cf = find(header.EVENT.TYP==781);
cue = cf-1;

% Extract trials, no-activity and activity intervals
%trials = [pos(cue)+1-tNA*fs pos(cue) pos(cf)+1 pos(cf)+tA*fs];
try
    trials = [pos(cue)-tNA*fs pos(cue) pos(cf) pos(cf)+tA*fs];
catch
    
    % Last trial corrupted, remove it
    cue = cue(1:end-1);
    cf = cf(1:end-1);
    trials = [pos(cue)-tNA*fs pos(cue) pos(cf) pos(cf+1)];
end

% Find classes available
ca = unique(header.EVENT.TYP(cue));

labels = zeros(1,length(cue));
for c=1:length(ca)
    labels(header.EVENT.TYP(cue) == ca(c)) = ca(c);
    N(c) = sum(labels==ca(c));
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
%freqs = [8:2:30];
freqs = [4:2:48];
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
BL = slabels(slabels==-1);


BDs = BD(1:floor(size(BD,1)/2),:,:);
BDe = BD(floor(size(BD,1)/2)+1:end,:,:);
BDa = [BDs;BDe];

MA{4} = mean(BDa(:,:));
CA{4} = cov(BDa(:,:));
MS{4} = mean(BDs(:,:));
CS{4} = cov(BDs(:,:));    
ME{4} = mean(BDe(:,:));
CE{4} = cov(BDe(:,:));

 

if(sum(ismember(ca,770))==1)
    D{1} = afeats(slabels==770,:,:);
    L{1} = slabels(slabels==770,:,:);
    
    
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

if(sum(ismember(ca,769))==1)
    D{2} = afeats(slabels==769,:,:);
    L{2} = slabels(slabels==769,:,:);
    

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



if(sum(ismember(ca,771))==1)
    D{3} = afeats(slabels==771,:,:);
    L{3} = slabels(slabels==771,:,:);

    D3s = D{3}(1:floor(size(D{3},1)/2),:,:);
    D3e = D{3}(floor(size(D{3},1)/2)+1:end,:,:);
    D3a = [D3s;D3e];

    MA{3} = mean(D3a(:,:));
    CA{3} = cov(D3a(:,:));
    MS{3} = mean(D3s(:,:));
    CS{3} = cov(D3s(:,:));    
    ME{3} = mean(D3e(:,:));
    CE{3} = cov(D3e(:,:));

else
    MA{3} = [];
    CA{3} = [];
    MS{3} = [];
    CS{3} = [];
    ME{3} = [];
    CE{3} = [];    
end   

% DP maps of classx vs class y
if(sum(ismember(ca,[770 769]))==2)
    % RHLH
    data = [D{1}; D{2}];
    data = reshape(data,size(data,1),nf*16);
    lbl = [L{1}; L{2}];
    %[COM,PWGR,V,vp,DISC]=cva_tun_opt(data,lbl);
    FS = eegc3_fs2(data(lbl==770,:),data(lbl==769,:));
    DPA{1} = reshape(FS,nf,16)';
    
    data = [D{1}(1:floor(size(D{1},1)/2),:,:); D{2}(1:floor(size(D{2},1)/2),:,:)];
    data = reshape(data,size(data,1),nf*16);
    lbl = [L{1}(1:floor(size(L{1},1)/2)); L{2}(1:floor(size(L{2},1)/2))];
    %[COM,PWGR,V,vp,DISC]=cva_tun_opt(data,lbl);
    FS = eegc3_fs2(data(lbl==770,:),data(lbl==769,:));    
    DPS{1} = reshape(FS,nf,16)';    
    
    data = [D{1}(floor(size(D{1},1)/2)+1:end,:,:); D{2}(floor(size(D{2},1)/2)+1:end,:,:)];
    data = reshape(data,size(data,1),nf*16);
    lbl = [L{1}(floor(size(L{1},1)/2)+1:end); L{2}(floor(size(L{2},1)/2)+1:end)];
    %[COM,PWGR,V,vp,DISC]=cva_tun_opt(data,lbl);
    FS = eegc3_fs2(data(lbl==770,:),data(lbl==769,:));        
    DPE{1} = reshape(FS,nf,16)';    
    
    % Automatic feature selection and classification
    tmpDPA = DPA{1}(:);
    [usl sortedDPA] = sort(tmpDPA,'descend');
    SelFeat = sortedDPA(1:5);
    cdata = data(:,SelFeat);
    [usl Err] = classify(cdata,cdata,lbl,'diaglinear');
    SAcc{1} = 100- 100*Err;
else
    DPA{1} = [];
    DPS{1} = [];
    DPE{1} = [];
    SAcc{1} = [];
end

if(sum(ismember(ca,[770 771]))==2)
    % RHBF
    data = [D{1}; D{3}];
    data = reshape(data,size(data,1),nf*16);
    lbl = [L{1}; L{3}];
    %[COM,PWGR,V,vp,DISC]=cva_tun_opt(data,lbl);
    FS = eegc3_fs2(data(lbl==770,:),data(lbl==771,:));        
    DPA{2} = reshape(FS,nf,16)';
    
    data = [D{1}(1:floor(size(D{1},1)/2),:,:); D{3}(1:floor(size(D{3},1)/2),:,:)];
    data = reshape(data,size(data,1),nf*16);
    lbl = [L{1}(1:floor(size(L{1},1)/2)); L{3}(1:floor(size(L{3},1)/2))];
    %[COM,PWGR,V,vp,DISC]=cva_tun_opt(data,lbl);
    FS = eegc3_fs2(data(lbl==770,:),data(lbl==771,:));            
    DPS{2} = reshape(FS,nf,16)';    
    
    data = [D{1}(floor(size(D{1},1)/2)+1:end,:,:); D{3}(floor(size(D{3},1)/2)+1:end,:,:)];
    data = reshape(data,size(data,1),nf*16);
    lbl = [L{1}(floor(size(L{1},1)/2)+1:end); L{3}(floor(size(L{3},1)/2)+1:end)];
    %[COM,PWGR,V,vp,DISC]=cva_tun_opt(data,lbl);
    FS = eegc3_fs2(data(lbl==770,:),data(lbl==771,:));            
    DPE{2} = reshape(FS,nf,16)';     
    
    % Automatic feature selection and classification
    tmpDPA = DPA{2}(:);
    [usl sortedDPA] = sort(tmpDPA,'descend');
    SelFeat = sortedDPA(1:5);
    cdata = data(:,SelFeat);
    [usl Err] = classify(cdata,cdata,lbl,'diaglinear');
    SAcc{2} = 100-100*Err;
else
    DPA{2} = [];
    DPS{2} = [];
    DPE{2} = [];
    SAcc{2} = [];
end

if(sum(ismember(ca,[771 769]))==2)
    % BFLH
    data = [D{2}; D{3}];
    data = reshape(data,size(data,1),nf*16);
    lbl = [L{2}; L{3}];
    %[COM,PWGR,V,vp,DISC]=cva_tun_opt(data,lbl);\
    FS = eegc3_fs2(data(lbl==771,:),data(lbl==769,:));            
    DPA{3} = reshape(FS,nf,16)';
    
    data = [D{2}(1:floor(size(D{2},1)/2),:,:); D{3}(1:floor(size(D{3},1)/2),:,:)];
    data = reshape(data,size(data,1),nf*16);
    lbl = [L{2}(1:floor(size(L{2},1)/2)); L{3}(1:floor(size(L{3},1)/2))];
    %[COM,PWGR,V,vp,DISC]=cva_tun_opt(data,lbl);
    FS = eegc3_fs2(data(lbl==771,:),data(lbl==769,:));                
    DPS{3} = reshape(FS,nf,16)';    
    
    data = [D{2}(floor(size(D{2},1)/2)+1:end,:,:); D{3}(floor(size(D{3},1)/2)+1:end,:,:)];
    data = reshape(data,size(data,1),nf*16);
    lbl = [L{2}(floor(size(L{2},1)/2)+1:end); L{3}(floor(size(L{3},1)/2)+1:end)];
    %[COM,PWGR,V,vp,DISC]=cva_tun_opt(data,lbl);
    FS = eegc3_fs2(data(lbl==771,:),data(lbl==769,:));                
    DPE{3} = reshape(FS,nf,16)';  
    
    % Automatic feature selection and classification
    tmpDPA = DPA{3}(:);
    [usl sortedDPA] = sort(tmpDPA,'descend');
    SelFeat = sortedDPA(1:5);
    cdata = data(:,SelFeat);
    [usl Err] = classify(cdata,cdata,lbl,'diaglinear');
    SAcc{3} = 100-100*Err;    
else
    DPA{3} = [];
    DPS{3} = [];
    DPE{3} = [];
    SAcc{3} = [];
end

% DP maps of classx vs baseline
if(sum(ismember(ca,770))==1)
    % RH
    data = [D{1}; BD];
    data = reshape(data,size(data,1),nf*16);
    lbl = [L{1}; BL];
    %[COM,PWGR,V,vp,DISC]=cva_tun_opt(data,lbl);
    FS = eegc3_fs2(data(lbl==770,:),data(lbl==-1,:));                
    MDPA{1} = reshape(FS,nf,16)';
    
    data = [D{1}(1:floor(size(D{1},1)/2),:,:); BD(1:floor(size(BD,1)/2),:,:)];
    data = reshape(data,size(data,1),nf*16);
    lbl = [L{1}(1:floor(size(L{1},1)/2)); BL(1:floor(size(BL,1)/2))];
    %[COM,PWGR,V,vp,DISC]=cva_tun_opt(data,lbl);
    FS = eegc3_fs2(data(lbl==770,:),data(lbl==-1,:));                    
    MDPS{1} = reshape(FS,nf,16)';    
    
    data = [D{1}(floor(size(D{1},1)/2)+1:end,:,:); BD(floor(size(BD,1)/2)+1:end,:,:)];
    data = reshape(data,size(data,1),nf*16);
    lbl = [L{1}(floor(size(L{1},1)/2)+1:end); BL(floor(size(BL,1)/2)+1:end)];
    %[COM,PWGR,V,vp,DISC]=cva_tun_opt(data,lbl);
    FS = eegc3_fs2(data(lbl==770,:),data(lbl==-1,:));                    
    MDPE{1} = reshape(FS,nf,16)';    
else
    MDPA{1} = [];
    MDPS{1} = [];
    MDPE{1} = [];
end

if(sum(ismember(ca,769))==1)
    % LH
    data = [D{2}; BD];
    data = reshape(data,size(data,1),nf*16);
    lbl = [L{2}; BL];
    %[COM,PWGR,V,vp,DISC]=cva_tun_opt(data,lbl);
    FS = eegc3_fs2(data(lbl==769,:),data(lbl==-1,:));                    
    MDPA{2} = reshape(FS,nf,16)';
    
    data = [D{2}(1:floor(size(D{2},1)/2),:,:); BD(1:floor(size(BD,1)/2),:,:)];
    data = reshape(data,size(data,1),nf*16);
    lbl = [L{2}(1:floor(size(L{2},1)/2)); BL(1:floor(size(BL,1)/2))];
    %[COM,PWGR,V,vp,DISC]=cva_tun_opt(data,lbl);
    FS = eegc3_fs2(data(lbl==769,:),data(lbl==-1,:));                        
    MDPS{2} = reshape(FS,nf,16)';    
    
    data = [D{2}(floor(size(D{2},1)/2)+1:end,:,:); BD(floor(size(BD,1)/2)+1:end,:,:)];
    data = reshape(data,size(data,1),nf*16);
    lbl = [L{2}(floor(size(L{2},1)/2)+1:end); BL(floor(size(BL,1)/2)+1:end)];
    %[COM,PWGR,V,vp,DISC]=cva_tun_opt(data,lbl);
    FS = eegc3_fs2(data(lbl==769,:),data(lbl==-1,:));                        
    MDPE{2} = reshape(FS,nf,16)';  
else
    MDPA{2} = [];
    MDPS{2} = [];
    MDPE{2} = [];
end

if(sum(ismember(ca,771))==1)
    % BF
    data = [D{3}; BD];
    data = reshape(data,size(data,1),nf*16);
    lbl = [L{3}; BL];
    %[COM,PWGR,V,vp,DISC]=cva_tun_opt(data,lbl);
    FS = eegc3_fs2(data(lbl==771,:),data(lbl==-1,:));                        
    MDPA{3} = reshape(FS,nf,16)';
    
    data = [D{3}(1:floor(size(D{3},1)/2),:,:); BD(1:floor(size(BD,1)/2),:,:)];
    data = reshape(data,size(data,1),nf*16);
    lbl = [L{3}(1:floor(size(L{3},1)/2)); BL(1:floor(size(BL,1)/2))];
    %[COM,PWGR,V,vp,DISC]=cva_tun_opt(data,lbl);
    FS = eegc3_fs2(data(lbl==771,:),data(lbl==-1,:));    
    MDPS{3} = reshape(FS,nf,16)';    
    
    data = [D{3}(floor(size(D{3},1)/2)+1:end,:,:); BD(floor(size(BD,1)/2)+1:end,:,:)];
    data = reshape(data,size(data,1),nf*16);
    lbl = [L{3}(floor(size(L{3},1)/2)+1:end); BL(floor(size(BL,1)/2)+1:end)];
    %[COM,PWGR,V,vp,DISC]=cva_tun_opt(data,lbl);
    FS = eegc3_fs2(data(lbl==771,:),data(lbl==-1,:));    
    MDPE{3} = reshape(FS,nf,16)';      
else
    MDPA{3} = [];
    MDPS{3} = [];
    MDPE{3} = [];
end