function [INDEX, P] = cnbiproc_get_outlayers(DATA, NSTD)
% [INDEX, P] = cnbiproc_get_outlayers(DATA, NSTD)
% 
% The function identifies outlayers in the given DATA based on values
% greater than NINCLUSION*(standard deviation). If DATA is a matrix (2 dimensions), 
% the function computes outlayers across the second dimension.
%
% Input:
%   - DATA          Data matrix [samples x variables]
%   - NSTD          Number of standard deviations to be considered outlayer
%
%   Output:
%   - INDEX         Logical vector/matrix of the same size of DATA : true
%                   means outlayer 
%   - P             Percentage of outlayers with respect to original data
%                   size

    if ismatrix(DATA) == false
        error('chk:dim', 'DATA must have maximum 2 dimensions');
    end
    
    INDEX  = false(size(DATA));
    P = zeros(size(DATA, 2), 1);
    
    for vId = 1:size(DATA, 2)
        cdata = DATA(:, vId);
       
        csigma  = std(cdata);
        INDEX(:, vId) = abs(cdata) > NSTD*csigma;
        P(vId) = 100*sum(INDEX(:, vId))./length(cdata);
    end

    
end