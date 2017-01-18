function FSmat = fisherScore(featmat, labels)

ulbl = unique(labels);
M1 = squeeze(mean(featmat(find(labels==ulbl(1)),:,:),1));
S1 = squeeze(std(featmat(find(labels==ulbl(1)),:,:),1));

M2 = squeeze(mean(featmat(find(labels==ulbl(2)),:,:),1));
S2 = squeeze(std(featmat(find(labels==ulbl(2)),:,:),1));

FSmat = abs(M1-M2)./sqrt(S1.^2 + S2.^2);