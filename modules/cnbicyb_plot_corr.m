SubList = {'AN14VE', 'MA25VE'};
NumSubjects = length(SubList);

datapath  = [pwd '/analysis/'];
figuredir = './figures/';

PatternLocationsId = {[4 9 14], [2 7 12 6 11 16]};
PatternLocationsLb = {'FCz-Cz-CPz', 'FC3-C3-CP3-CP4-C4-FC4'}; 

for sId = 1:NumSubjects
    
    csubject = SubList{sId};

    cracetimefilepath = [datapath '/' csubject '.time.race.mat'];
    cpadaccfilepath = [datapath '/' csubject '.padaccuracy.race.mat'];
    cdiscrfilepath = [datapath '/' csubject '.discriminancy.maps.mat'];
    ctopfilepath = [datapath '/' csubject '.timeonpad.race.mat'];
    cmetafilepath = [datapath '/' csubject '.metadata.mat'];
    
    cdatatime = load(cracetimefilepath);
    cdatapadacc = load(cpadaccfilepath);
    cdatadiscr = load(cdiscrfilepath);
    cdatatop = load(ctopfilepath);
    cdatameta = load(cmetafilepath);
    
    % Race times
    crt = cdatatime.time.run.values;
    crtDk = cdatatime.time.run.label.Dk;
    crtDl = mat2cell(cdatatime.time.run.label.Dl,ones(1,size(cdatatime.time.run.label.Dl,1)),8);
    
    % Command accuracies (Pad accuracy actually)
    cpa = cdatapadacc.pad.accuracy.run.values;
    cpaDk = cdatapadacc.pad.label.run.Dk;
    cpaDl = mat2cell(cdatapadacc.pad.label.session.Dl,ones(1,size(cdatapadacc.pad.label.session.Dl,1)),8);
    
    % Fisher scores
    cfs = cdatadiscr.discriminancy.run.fisherscore;
    cfsDk = cdatadiscr.discriminancy.run.label.Dk;
    cfsDl = mat2cell(cdatadiscr.discriminancy.run.label.Dl,ones(1,size(cdatadiscr.discriminancy.run.label.Dl,1)),8);
    cfsMk = cdatadiscr.discriminancy.run.label.Mk;
    
    % Time on pad
    ctopa = cdatatop.timeonpad.topall.values;
    ctoplbl = cdatatop.timeonpad.topall.labels;
    for rId=1:(length(ctopa)/18)
        for ctype=1:4
            tmplbl = ctoplbl((rId-1)*18+1:rId*18);
            tmptop = ctopa((rId-1)*18+1:rId*18);
            ctop(rId,ctype) = nanmean(tmptop(tmplbl==ctype));
        end
    end
        
    cfsDk = cdatadiscr.discriminancy.run.label.Dk;
    cfsDl = mat2cell(cdatadiscr.discriminancy.run.label.Dl,ones(1,size(cdatadiscr.discriminancy.run.label.Dl,1)),8);
    cfsMk = cdatadiscr.discriminancy.run.label.Mk;
    
    FreqGrid = cdatadiscr.discriminancy.freqs;
    BetaFreqs = 22:2:32;
    [~, SelBetaFreqIds] = intersect(FreqGrid, BetaFreqs);
    
    % Create the common attributes
    AllSessions = union(union(crtDl,cpaDl),cfsDl);
    
    ss = 0;
    crtF = {};
    cpaAllF = {};
    cpaSpinF = {};
    cpaJumpF = {};
    cfsMedialF = {};
    cfsLateralF = {};
    crtFA = [];
    cpaAllFA = [];
    cpaSpinFA = [];
    cpaJumpFA = [];
    cfsMedialFA = [];
    cfsLateralFA = [];    
    Num = [];
    for s=1:length(AllSessions)
        % Check if this session exists forall metrics, this is basically
        % for the fisher scores only
        ExistsIn = sum(strcmp(AllSessions{s},crtDl))+sum(strcmp(AllSessions{s},cpaDl))...
            +sum(strcmp(AllSessions{s},cfsDl));
        if(ExistsIn == 3)
            ss=ss+1;
            
            % Race Time
            crtF{ss} = crt(find(crtDk==find(strcmp(AllSessions{s},crtDl))));
           
            % Do not count slide
            cpaAllF{ss} = nanmean(cpa(1:2,find(cpaDk==find(strcmp(AllSessions{s},cpaDl)))));
            cpaSpinF{ss} = squeeze(cpa(1,find(cpaDk==find(strcmp(AllSessions{s},cpaDl)))));
            cpaJumpF{ss} = squeeze(cpa(2,find(cpaDk==find(strcmp(AllSessions{s},cpaDl)))));
            
            % Exclude also  runs that are not races (Mk<2)
            tmpfs = cfs(:,intersect(find(cfsMk>=2),find(cfsDk==find(strcmp(AllSessions{s},cfsDl)))),4);
            tmpfs = reshape(tmpfs,[23 16 size(tmpfs,2)]);
            % Alternative taskset [770 771]
            tmpfs2 = cfs(:,intersect(find(cfsMk>=2),find(cfsDk==find(strcmp(AllSessions{s},cfsDl)))),1);
            tmpfs2 = reshape(tmpfs2,[23 16 size(tmpfs2,2)]);
            
            cfsMedialF{ss} = squeeze(nanmean((squeeze(nanmean(tmpfs(SelBetaFreqIds,PatternLocationsId{1},:),1))),1));
            cfsLateralF{ss} = squeeze(nanmean((squeeze(nanmean(tmpfs(SelBetaFreqIds,PatternLocationsId{2},:),1))),1));
            
            cfsMedialF2{ss} = squeeze(nanmean((squeeze(nanmean(tmpfs2(SelBetaFreqIds,PatternLocationsId{1},:),1))),1));
            cfsLateralF2{ss} = squeeze(nanmean((squeeze(nanmean(tmpfs2(SelBetaFreqIds,PatternLocationsId{2},:),1))),1));
            
            Num(ss,:) = [size(crtF{ss},1) size(cpaAllF{ss},2) size(cfsMedialF{ss},2)];
            
            if(Num(ss,1) > Num(ss,2))
                % There is one qeird session for Eric which has onemore
                % race time more, without padacc and fs..., the second is manually deleted
                % so I should delete also here
                crtF{ss}(2)=[];
            end
            Num(ss,:) = [size(crtF{ss},1) size(cpaAllF{ss},2) size(cfsMedialF{ss},2)];
            
%             crtFA(ss) = nanmean(crtF{ss});
%             cpaAllFA(ss) = nanmean(cpaAllF{ss});
%             cpaSpinFA(ss) = nanmean(cpaSpinF{ss});
%             cpaJumpFA(ss) = nanmean(cpaJumpF{ss});
%             cfsMedialFA(ss) = nanmean(cfsMedialF{ss});
%             cfsLateralFA(ss) = nanmean(cfsLateralF{ss});

        end
    end

    crtC = cell2mat(crtF');
    
    ctopAllC = nanmean(ctop(:,1:3),2);
    ctopSpinC = ctop(:,1);
    ctopJumpC = ctop(:,2);
    
    cpaAllC = cell2mat(cpaAllF)';
    cpaSpinC = cell2mat(cpaSpinF)';
    cpaJumpC = cell2mat(cpaJumpF)';
    cfsMedialC = cell2mat(cfsMedialF)';
    cfsLateralC = cell2mat(cfsLateralF)';
    cfsMedialC2 = cell2mat(cfsMedialF2)';
    cfsLateralC2 = cell2mat(cfsLateralF2)';
    if(size(cfsLateralC,1)==size(cfsLateralC2,1))
        cfsMedialC = nansum([cfsMedialC cfsMedialC2],2);
        cfsLateralC = nansum([cfsLateralC cfsLateralC2],2);
    end
    
    cfsAllC = cfsMedialC+cfsLateralC;
    

    % Race Time vs Accuracy
    disp(SubList(sId));
    [r p] = corr(cpaAllC,crtC);
    disp(['RT vs Ovderall Accuracy: r=' num2str(r) ', p=' num2str(p)]);
    [r p] = corr(crtC(~isnan(cfsAllC)),cfsAllC(~isnan(cfsAllC)));
    
    
    disp(['Spin Accuracy vs Lateral Fisher Score: r=' num2str(r) ', p=' num2str(p)]);
    [r p] = corr(cpaJumpC(~isnan(cfsMedialC)),cfsMedialC(~isnan(cfsMedialC)));
    disp(['Jump Accuracy vs Medial Fisher Score: r=' num2str(r) ', p=' num2str(p)]);
    [r p] = corr(cpaAllC(~isnan(cfsAllC)),cfsAllC(~isnan(cfsAllC)));
    disp(['Overall Accuracy vs Overall Fisher Score: r=' num2str(r) ', p=' num2str(p)]);
    
    
    [r p] = corr(ctopSpinC(~isnan(cfsLateralC)),cfsLateralC(~isnan(cfsLateralC)));
    disp(['Spin TOP vs Lateral Fisher Score: r=' num2str(r) ', p=' num2str(p)]);
    [r p] = corr(ctopJumpC(~isnan(cfsMedialC)),cfsMedialC(~isnan(cfsMedialC)));
    disp(['Jump TOP vs Medial Fisher Score: r=' num2str(r) ', p=' num2str(p)]);
    [r p] = corr(ctopAllC(~isnan(cfsAllC)),cfsAllC(~isnan(cfsAllC)));
    disp(['Overall TOP vs Overall Fisher Score: r=' num2str(r) ', p=' num2str(p)]);
end
