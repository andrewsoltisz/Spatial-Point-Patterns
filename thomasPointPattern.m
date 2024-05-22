function [points, clusterIndex] = thomasPointPattern(studyArea, nPoints, nClusters, cluster_lambdaDaughter, cluster_meanRadius, cluster_stdevRadius, randSeed)
% Create a clustered point pattern following a Thomas point process, a kind
% of Neyman-Scott point process where cluster seeds ("parent" points) are
% distributed via a Poisson process, the number of ("daughter") points
% composing each cluster has a Poisson distribution, and the size of each
% cluster (radius) is normally distributed. The dimensionality of the
% resulting point pattern will match the dimensionality of 'studyArea'.
%
%
% INPUTS: 
%
% studyArea - bounds of the region that points will be placed in, formatted
%             as a d-by-2 matrix, where 'd' is the dimensionality of the
%             point pattern. This matrix defines a d-dimensional bounding
%             box, where row 'i' specifies the minimum (column 1) and
%             maximum (colum 2) bounds of the bounding box along dimension
%             'i'. For example, a cubic study area with side length 1 and
%             bottom left corner at the origin would be specified as 
%             [0,1; 0,1; 0,1]. 
%
% nPoints_total - total number of daughter points (sample size) to
%                 generate, formatted as a positive integer.
%
% clusterFraction - number of clusters ("parent" points), calculated as
%                   fraction of total points and formatted as a numeric
%                   scalar in the range (0,1].
%
% cluster_lambdaDaughter - lambda parameter of the Poisson distribution
%                          that will be sample to determine the number of
%                          "daughter" points assigned to each cluster,
%                          formatted as a positive numeric scalar.
%
% cluster_meanRadius - mean cluster size (radius), formatted as a numeric 
%                      scalar.
%
% cluster_stdevRadius - standard deviation of cluster size (radius),
%                       formatted as a numeric scalar.
%
% randSeed - optional input which sets the seed value for random number
%            generation, formatted as a positive integer. Given constant
%            values for all other input parameters, using the same seed
%            number will produce identical point patterns.
%
%
% OUTPUTS:
%
% points - Euclidean coordinates of the resulting point pattern formatted
%          as an n-by-d numeric matrix, where 'n' is the samples size
%          (equal to 'nPoints_total') and 'd' is the dimensionality of the
%          point pattern.
%
% clusterIndex - index of the cluster each point belongs to, formatted as
%                an n-by-1 matrix, where 'n' is the sample size.
%
%
% AUTHORSHIP: 
%
% Author: Andrew M. Solitsz
% Contact: andysoltisz@gmail.com
%

    % input validation
    narginchk(6,7); % validate number of input arguments
    if nargin == 6
        % generate a random number seed if none provided
        randSeed = 'shuffle';
    end

    % generate cluster centers
    cluster_centers = poissonPointPattern(studyArea, nClusters, randSeed); % "parent" point process

    % Determine number of daughter points for each cluster. Assign each
    % cluster a number of points 'ni', such that 'n' is Poisson distributed
    % and the sum of all 'ni' adds to the input 'nPoints'
    rng(randSeed, "twister"); s = rng; rng(s); % set random seed
    nPoints_iCluster = poissrnd(cluster_lambdaDaughter, nClusters, 1);
    if any(nPoints_iCluster == 0) 
        nPoints_iCluster = nPoints_iCluster + 1; % make sure all 'ni' > 0
    end
    nPoints_iCluster = nPoints_iCluster ./ sum(nPoints_iCluster); % normalize so sum==1
    nPoints_iCluster = round(nPoints_iCluster .* nPoints); % scale so that sum==total points 
    
    % fix error in total point count from rounding step above
    nPoints_error = sum(nPoints_iCluster) - nPoints;
    if nPoints_error ~= 0
        nPoints_absError = abs(nPoints_error);
        largestGlobalFix = floor(nPoints_absError / nClusters); % what is the largest nPoints we can take/add to all cluster?
        nPoints_left = nPoints_absError - (largestGlobalFix * nClusters); % how many nPoints are left to take/add?
        nPoints_fix = largestGlobalFix + [ones(nPoints_left,1); zeros(nClusters-nPoints_left,1)]; % put it all together
        nPoints_iCluster = nPoints_iCluster + (nPoints_fix .* -sign(nPoints_error)); % take/add points from clusters
    end

    % cluster sizes are normally distributed
    rng(randSeed, "twister"); s = rng; rng(s); % set random seed
    iCluster_radius = randn(nClusters,1) .* cluster_stdevRadius + cluster_meanRadius; % sample cluster sizes

    % place daughter points in each cluster so they are normally
    % distributed about the cluster's center
    nDims = size(studyArea, 1); % determine dimensionality of point pattern
    points = zeros(nPoints, nDims); % preallocate "daughter" point process
    clusterIndex = zeros(nPoints, 1); % preallocate
    iPoint = 1; % intial point index
    for iCluster = 1:nClusters
        nPoints_daughter = nPoints_iCluster(iCluster); % number of "daughter" points in this cluster
        fPoint = nPoints_daughter + iPoint - 1; % final point index
        pos = cluster_centers(iCluster, :); % grab position of this cluster's center 
        sz = iCluster_radius(iCluster); % grab size (radius) of this cluster
        points(iPoint:fPoint, :) = gaussianPointPattern(nPoints_daughter, sz, pos, randSeed); % generate normally distributed point pattern
        clusterIndex(iPoint:fPoint) = iCluster;
        iPoint = fPoint + 1; % update intial point index
    end    

end