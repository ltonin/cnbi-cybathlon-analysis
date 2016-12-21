function [Dr, Drk] = cnbiproc_reshape_ts_bc(D, Dk)
% [Dr, Drk] = cnbiproc_reshape_ts_bc(D, Dk)
%
% The function reshape the given dataset D and the trial labels Dk:
%
% [Samples x Frequencies x Channels x Trials] -> [(Samples x Trials) x (Frequencies x Channels)]
%
% Input:
%   D   [sample x freqs x chans x trials]
%   Dk  [trials x 1]
%
% Last dimension (Trials) is not required.

    NumDims = ndims(D);  
    rshlabel = true;
    
    if nargin < 2
        rshlabel = false;
        Dk = [];
    end
    
    if ((NumDims < 3) || (NumDims > 4))
        error('chk:dim', 'Input matrix must have 3 or 4 dimensions to perform the reshaping');
    end
    
    NumSamples = size(D, 1);
    NumFreqs   = size(D, 2);
    NumChans   = size(D, 3);
    NumTrials  = size(D, 4);
    
    NumFeatures     = NumFreqs*NumChans;
    NumObservations = NumSamples*NumTrials;
    
    if NumDims == 4
        D = permute(D, [1 4 2 3]);
    end

    Dr = reshape(D, [NumObservations NumFeatures]);
   
    
    if rshlabel == true
        if NumTrials ~= length(Dk)
            error('chk:arg', 'size(D, 4) must be equal to length(Dk)');
        end
    
        Drk = reshape(repmat(Dk, [1 NumSamples])', [NumObservations 1]);
    end
end