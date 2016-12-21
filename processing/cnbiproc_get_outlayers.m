function INDEX = cnbiproc_get_outlayers(DATA, N)
% indexes = cnbiproc_get_outlayers(DATA, N)
% 
% The function identifies outlayers in the given DATA based on values
% greater than N*(standard deviation). If DATA is a matrix (2 dimensions), 
% the function computes outlayers across the second dimension.
%
% Input:
%   - DATA          Data matrix [samples x variables]
%   - N             Number of standard deviations
%
%   Output:
%   - INDEX         Logical vector/matrix of the same size of DATA with
%                   with outlayers

    if ismatrix(DATA) == false
        error('chk:dim', 'DATA must have maximum 2 dimensions');
    end
    
    INDEX  = false(size(DATA));
    
    for vId = 1:size(DATA, 2)
        cdata = DATA(:, vId);
       
        csigma  = std(cdata);
        INDEX(:, vId) = abs(cdata) >= N*csigma;
    end

end