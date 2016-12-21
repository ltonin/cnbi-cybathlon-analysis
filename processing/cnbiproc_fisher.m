function F = cnbiproc_fisher(P, Pk, rmsize)
% F = proc_fisher(P, Pk)
% 
% The function computes the fisher score for each pair channel-frequency.
%
% Input:
%   - P 	data matrix  [samples x frequencies x channels]
%   - Pk    label vector (Only two classes are allowed) [samples x 1]
%
% It returns a vector F of fisher score. The vector is in the format
% [(channelsxfrequencies) x 1]
    
    % Default value for rmsize (outlayer removal)
    if nargin == 2
        rmsize = [];
    end

    % Check number of classes
    Classes = unique(Pk);
    NumClasses = length(Classes);
    
    if NumClasses ~= 2
        error('chk:classes', 'Number of classes must be 2');
    end
    
    % Check dimensionality
    if isequal(size(P, 1), length(Pk)) == false
        error('chk:dim', 'First dimension of P and length of Pk must be the same');
    end
    
    
    % Reshaping data matrix [samples x (channels x frequencies)]
    rP = cnbiproc_reshape_ts_bc(P);
    
    % Identify outlayers
    outlayers = false(size(rP));
    if isempty(rmsize) == false
        outlayers = cnbiproc_get_outlayers(rP, rmsize);
    end
    
    
    F = zeros(size(rP, 2), 1);
    for fId = 1:size(rP, 2)
        coutlayers = outlayers(:, fId);
        cindex1 = Pk == Classes(1) & coutlayers == false;
        cindex2 = Pk == Classes(2) & coutlayers == false;
        
        m1 = mean(rP(cindex1, fId));
        s1 = std(rP(cindex1, fId));
        
        m2 = mean(rP(cindex2, fId));
        s2 = std(rP(cindex2, fId));
        F(fId) = abs(m2 - m1) ./ sqrt(s1.^2 + s2.^2);
    end
    

    

end