function [ERSP, RSTR, ASTR] = cnbiproc_ersp(reference, activity, rmsize)
% [ERSP, RSTR, ASTR] = cnbiproc_ersp(reference, activity, rmsize)
%
% The function computes the ERSP (Event Related Spectral Perturbation) for
% each element k in the activity vector X_a with respect to the reference 
% distribution. 
% In detail, it computes MU_r and SIGMA_r for the reference distribution 
% and, thus, the z-score for each element X_a(k) in the activity vector, 
% according to the forumula:
%
%   ERSP(k) = (X_a(k) - MU_r)./SIGMA_r 
%
% By default, the function remove nan values and values greater than 3 
% times the standard deviation of two vectors.
%
% Input:
%   - reference         Vector referring to the reference period [samples x 1]
%   - activity          Vector referring to the activity period [samples x 1]
%   - rmsize            Optional argument to control outlayer removal. It
%                       represents the number of standard deviations to be
%                       set for outlayer removal. [Default: 3]
%
% Output:
%   - ERSP              Event Related Spectral Perturbation [samples x 1]
%   - RSTR              Structure with MU, SIGMA and percentage of removed 
%                       outlayers in the reference distribution
%   - ASTR              Structure with MU, SIGMA and percentage of removed 
%                       outlayers in the activity distribution

    if nargin == 2
        rmsize = 3;
    end

    nsamplesR = size(reference, 1);
    nsamplesA = size(activity, 1);
   
    % Check the nan values and the values greater than RMSIZE*SIGMA for
    % reference and activity distributions
    indexR  = isnan(reference) == false | abs(reference) < rmsize*nanstd(reference);
    removR = 1 - sum(indexR)./nsamplesR;
    
    indexA  = isnan(activity) == false  | abs(activity) < rmsize*nanstd(activity);
    removA = 1 - sum(indexA)./nsamplesA;
    
    % Raise a warning if the percentage of removed samples is greater than 30%
    if (removR > 0.3)
        warning('chk:rm', [num2str(100*removR) '% have been removed']);
    end
    
    if (removA > 0.3)
        warning('chk:rm', [num2str(100*removA) '% have been removed']);
    end
    
    % Compute MU and SIGMA for reference and activity distributions
    muR     = mean(reference(indexR));
    sigmaR  = std(reference(indexR));

    muA     = mean(activity(indexA));
    sigmaA  = std(activity(indexA));
    
    % Compute ERSP as z-score between each activity sample and the
    % reference distribution
    ERSP = (activity(indexA) - muR)./sigmaR;

    % Fill output structures
    RSTR.MU     = muR;
    RSTR.SIGMA  = sigmaR;
    RSTR.REJ    = 100*removR;
    
    ASTR.MU     = muA;
    ASTR.SIGMA  = sigmaA;
    ASTR.REJ    = 100*removA;


end