function [points] = gaussianPointPattern(nPoints, clusterRadius, clusterCenter, randSeed)
% Create a clustered point pattern composed of exactly one Gaussian
% cluster, where 'nPoints' are assigned locations about the parent point at
% 'clusterCenter' according to a Gaussian distribution with standard
% deviation 'clusterRadius' and mean 'clusterCenter'. The dimensionality of
% the resulting point pattern will match the dimensionality of
% 'clusterCenter'.
%
% INPUTS:
%
% nPoints - total number of points (sample size) composing the point
%           pattern, formatted as a positive integer.
%
% clusterRadius - cluster size, actualized as the standard deviation of the
%                 Guassian distribution used to place the points within the
%                 study area. Formatted as a numeric scalar.
%
% clusterCenter - Euclidean coordinates of the cluster's center, actualized
%                 as the mean of the Guassian distribution used to place
%                 the points within the study area. Formatted as a 1-by-d
%                 numeric matrix, where 'd' is the dimensionality of the
%                 point pattern.
%
% randSeed - optional input which sets the seed value for random number
%            generation, formatted as a positive integer. Given constant
%            values for all other input parameters, using the same seed
%            value will produce identical point patterns. 
% 
% OUTPUTS: 
%
% points - Euclidean coordinates of the resulting point pattern, formatted
%          as an n-by-d numeric matrix, where 'n' is the samples size
%          (equal to 'nPoints_total') and 'd' is the dimensionality of the
%          point pattern.
%
%
% AUTHORSHIP: 
%
% Author: Andrew M. Solitsz
% Contact: andysoltisz@gmail.com
%

    narginchk(3,4); % validation number of input arguments
    if nargin == 3
        % generate a random number seed if none provided
        randSeed = 'shuffle';
    end

    nDims_size = numel(clusterRadius);
    nDims_data = numel(clusterCenter); % determine dimensionality of point pattern (ie 2D vs 3D)
    if nDims_size ~= 1 && nDims_size ~= nDims_data
        error("Cluster size parameter must be a scalar or array with the same shape as cluster center.");
    end

    % place points in each cluster so they are normally distributed about
    % the cluster's center
    rng(randSeed, "twister"); s = rng; rng(s); % set random seed
    points = randn(nPoints, nDims_data) .* clusterRadius + clusterCenter; % sample point coordinates

end