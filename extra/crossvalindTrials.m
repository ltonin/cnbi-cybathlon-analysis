function Indices = crossvalindTrials(alabels,atrials,K)

TrLbl=[];
Utrials = unique(atrials)'; 
for tr=Utrials
    tmptr = alabels(find(atrials==tr));
    TrLbl(tr) = tmptr(1);
end
cvind = crossvalind('Kfold',TrLbl,K);

Indices = alabels;
for tr=Utrials
    Indices(find(atrials==tr)) = cvind(tr);
end