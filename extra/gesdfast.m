function [Rmax lambda r_out outlier] = gesdfast(x, r, sided, alpha, dist)
%
% The Generalized ESD procedure ('gesd'), tests for up to a prespecified
% number 'r' of outliers and it is specially recommended when testing for
% outliers among data coming from a normal distribution.
%
% Data in 'x' are organized so that columns are the time series and rows
% are the time intervals. All series contain the same number of
% observations.
%
% 'r' is the predefined number of outliers to check. It can be chosen
% somewhat higher than anticipated.
%
% 'sided' indicates if the critical values are from a two-sided
% distribution (default) or one-sided ('sided' = 1).
%
% 'alpha' specifies the significant level. By default 'alpha' = 0.05.
%
% The variable 'dist' indicates if the series are assumed to be normal
% (default) or log-normally distributed. It can take the values 0 (normal)
% or 1 (log-normal).
%
% [Rmax lambda r_out outlier outlier_num] = gesd(...) returns the following:
% Rmax          - indicates the maximum absolute z-score for iteration r_{i}.
%                 Iterations {1, ..., r} correspond to the rows.
% lambda        - critical values for each iteration {1, ..., r} (rows).
% r_out         - highest number of the iteration {1, ..., r} where Rmax >
%                 lambda, indicating the values to be considered outliers
% outlier_num   - matrix providing row and column numbers of the values in
%                 'x' considered as potential outliers.
%
% Created by Francisco Augusto Alcaraz Garcia
%            alcaraz_garcia@yahoo.com
%
% References:
%
% 1) B. Iglewicz; D.C. Hoaglin (1993). How to Detect and Handle Outliers.
% ASQC Basic References in Quality Control, vol. 16, Wisconsin, US.
%
% 2) B. Rosner (1983). Percentage Points for a Generalized ESD Many-Outlier
% Procedure. Technometrics 25(2), pp. 165-172.


% Check number of input arguements
if (nargin < 2) || (nargin > 5)
    error('Requires two to five input arguments.')
end

% Define default values
if nargin == 2,
    sided = 1;
    alpha = 0.05;
    dist = 0;
elseif nargin == 3,
    alpha = 0.05;
    dist = 0;
elseif nargin == 4,
    dist = 0;
end

% Normal transformation
if dist == 1,
    x = log(x);
end

if sided == 1,
   alpha = alpha/2;
end

% Check for validity of inputs
if ~isnumeric(x) || ~isreal(x),
    error('Input x must be a numeric array and x_date must be a string table.')
elseif alpha <= 0 || alpha >= 1,
    error('The confidence level must be between zero and one.')
end

[n, c] = size(x);
xr = x;
Rmax = zeros(r,c);
lambda = zeros(r,c);

row = cell(r,1);
col = cell(r,1);
for i = 1:r,
    nr = n - sum(isnan(xr));
    R = abs((xr - repmat(nanmean(xr),n,1))./repmat(nanstd(xr),n,1));
    Rmax(i,:) = nanmax(R);
    [i1, j1] = find(R == repmat(Rmax(i,:),n,1));
    row{i} = i1; % data structure
    col{i} = j1; % data structure
    xr(i1,j1) = NaN;
    
    p = 1 - alpha./nr;
    lambda(i,:) = tinv(p,nr-2).*(nr-1)./sqrt((nr-2+tinv(p,nr-2).^2).*nr);
end

pos = Rmax > lambda;

r_out = zeros(1,c);
outlier = [];
for i = 1:c,
    [i2, j2] = find(pos(:,i)==1);
    if ~isempty(max(i2)),
        r_out(1,i) = max(i2);
        for j=1:r_out(1,i),
            outlier = [outlier; row{j}(col{j}==i) ...
                repmat(i,sum(col{j}==i),1)];
        end
    end
end


