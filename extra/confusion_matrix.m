function [ConfMatClass, ConfMatAll, Accuracy, Error, ConfMatReal] = confusion_matrix(TrueLabels, FoundLabels, NClass)

% function [ConfMatClass ConfMatAll Accuracy Error] =
% eegc3_confusion_matrix(TrueLabels, FoundLabels)
%
% Function to recover a confusion matric of classification problem given
% the groiund truth classes and the classes found by a classifier
%
% Inputs:
%
% TrueLabels: Ground truth data labels
%
% FoundLabels: Labels predicted by a classifier
%
% Outputs: 
%
% ConfMatClass: Confusion Matrx ClassNum x ClassNum in percentages per
% class
%
% ConfMatAll: Confusion Matrx ClassNum x ClassNum in overall percentages
% 
% Accuracy: Total accuracy % (diagonal of Confusion Matrix)
%
% Error: Total Error % (100 - Accuracy)

if(~isvector(TrueLabels) || ~isvector(FoundLabels))
    disp('[eegc3_confusion_matrix] Labels must be vectors');
    ConfMat = [];
    Accuracy = [];
    Error = [];
    return;
end

N1 = length(TrueLabels);
N2 = length(FoundLabels);

if(N1 ~= N2)
    disp('[eegc3_confusion_matrix] True and predicted labels must be the same size');
    ConfMat = [];
    Accuracy = [];
    Error = [];
    return;
else
    NSample = N1;
end

% Compute Confusion Matrix

ConfMat = zeros(NClass);
ConfMatAll = zeros(NClass);
for i=1:NSample
    ConfMat(TrueLabels(i),FoundLabels(i)) = ...
        ConfMat(TrueLabels(i),FoundLabels(i)) + 1;
end

ConfMatReal = ConfMat;
ConfMatAll = 100*ConfMat/sum(ConfMat(:));

ConfMatClass = nan(NClass);
for i=1:NClass
    if(sum(ConfMat(i,:))==0)
        ConfMatClass(i,:) = NaN;
    else
        ConfMatClass(i,:) = 100*ConfMat(i,:)/sum(ConfMat(i,:));
    end
end

Accuracy = nansum(diag((ConfMatAll)));

Error = 100 - Accuracy;