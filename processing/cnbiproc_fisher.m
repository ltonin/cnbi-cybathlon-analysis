function F = cnbiproc_fisher(P, Pk, NSTD)
% F = cnbiproc_fisher(P, Pk [, NSTD])
% 
% The function computes the fisher score for each pair channel-frequency.
% Nan values are not considered. Optional argument NSTD to remove outlayers 
% (number of standard deviation to consider a sample as outlayer). 
% If required, the outlayer removal is performed for each class, separately. 
%
% Input:
%   - P         data matrix  [samples x frequencies x channels]
%   - Pk        label vector (Only two classes are allowed) [samples x 1]
%   - NSTD      Outlayers removal. Number of standard deviation to 
%               consider a sample as outlayer 
%
% It returns a vector F of fisher score. The vector is in the format
% [(channelsxfrequencies) x 1]
%
% SEE ALSO: cnbiproc_get_outlayers
    
    % Default value for rmsize (outlayer removal)
    if nargin == 2
        NSTD = [];
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
    
    % Class indexes
    Index1 = Pk == Classes(1);
    Index2 = Pk == Classes(2);
    
    F = zeros(size(rP, 2), 1);
    for fId = 1:size(rP, 2)
        
        % Getting current data for each class and the given feature
        cdata1 = rP(Index1, fId);
        cdata2 = rP(Index2, fId);
        
        % If rmsize is provided, remove outlayers per class
        if isempty(NSTD) == false
            out1 = cnbiproc_get_outlayers(cdata1, NSTD);
            out2 = cnbiproc_get_outlayers(cdata2, NSTD);
            cdata1 = cdata1(out1 == false);
            cdata2 = cdata2(out2 == false);
        end
        
        % Computing mean and standard deviation for each class
        m1 = nanmean(cdata1);
        s1 = nanstd(cdata1);
        
        m2 = nanmean(cdata2);
        s2 = nanstd(cdata2);
        
        % Computing feature score for the given feature
        F(fId) = abs(m2 - m1) ./ sqrt(s1.^2 + s2.^2);
    end
    

    

end