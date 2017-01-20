function [AvgAcc AvgCM AvgCMReal] = cvk(afeats, alabels, atrials, K, NFeat)

% Accuracy with cross validation (balanced groups, random shuffling)
%cvind = crossvalind('Kfold',alabels,K);
cvind = crossvalindTrials(alabels,atrials,K);

for i=K:-1:1
    tstfoldind = find(cvind==i);
    trfoldind = setdiff([1:size(afeats,1)],find(cvind==i))';
   
    
    % Do feature selection on training set only
    for fr=1:size(afeats,2)
        for ch=1:size(afeats,3)
            cfeat = squeeze(afeats(:,fr,ch));
            clabels = alabels;
            FSmat(fr,ch) = fisherScore(cfeat,clabels); 
        end
    end

    [SortedFS,SortInd] = sort(FSmat(:),'descend');
    [i1,i2] = ind2sub(size(FSmat),SortInd);
    i1 = i1(1:NFeat);
    i2 = i2(1:NFeat);

    trdata = afeats;
    trdata(tstfoldind,:,:) = [];
    strdata = [];
    strlabels = alabels;
    strlabels(tstfoldind) = [];
    for s=1:size(trdata,1)
        for j=1:NFeat
            strdata(s,j) = trdata(s,i1(j),i2(j));
        end
    end
    
    
    ststdata = [];
    tstdata = afeats;
    tstdata(trfoldind,:,:)=[];
    ststlabels = alabels;
    ststlabels(trfoldind)=[];
    for s=1:size(tstdata,1)
        for j=1:NFeat
            ststdata(s,j) = tstdata(s,i1(j),i2(j));
        end
    end
    
    [~, err_tr(i)] = classify(strdata,strdata,strlabels);
    %[~, err_tr(i)] = classify(strdata,strdata,strlabels,'diaglinear');
    
    AccCVtr(i) = 100*(1-err_tr(i));
    fclass = classify(ststdata,strdata,strlabels);
    %fclass = classify(ststdata,strdata,strlabels,'diaglinear');
    [CMCVtst(i,:,:),~,AccCVtst(i), ~, CMRealCVtst(i,:,:) ] = confusion_matrix(fclass, ststlabels, 2);
end
AvgAcc = mean(AccCVtst);
AvgCM = squeeze(mean(CMCVtst,1));
AvgCMReal = squeeze(mean(CMRealCVtst,1));
