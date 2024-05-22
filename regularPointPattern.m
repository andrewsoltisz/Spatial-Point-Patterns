function [points, spacing] = regularPointPattern(studyArea, nPoints, jitter_stdev)
% Create a regular point pattern with sample size of 'nPoints'. The
% dimensionality of the resulting point pattern will match the
% dimensionality of 'studyArea'. The total number of resulting points may
% differ slightly from 'nPoints' because the regular point pattern is
% realized as an N-dimensional grid of evenly spaced points, which means
% the number of points along each dimension will be an integer and the
% total number of resulting points will be the product of that set of
% integers. So the final sample size will be the closest number to
% 'nPoints' that can be achieved by a set of 'd' integers. Optionally, the
% resulting points can be jittered with Guassian positional error with mean
% of 0 and standard deviation equal to 'jitter_stdev'.
%
% INPUTS:
%
% studyArea - bounds of the region that points will be placed in, formatted
%             as a d-by-2 matrix, where 'd' is the dimensionality of the
%             point pattern. This matrix defines a d-dimensional bounding
%             box, where row 'i' specifies the minimum (column 1) and
%             maximum (colum 2) bounds of the bounding box along dimension
%             'i'. For example, a cubic study area with side length 1 would
%             be specified as [0,1; 0,1; 0,1].
%
% nPoints - total number of points (sample size) composing the point
%           pattern, formatted as a positive integer.
%
% jitter_stdev - standard deviation of the Guassian positional jitter
%                applied to the resulting regular point pattern, formatted
%                as a positive numeric scalar.
%
% 
% OUTPUTS: 
%
% points - Euclidean coordinates of the resulting point pattern formatted
%          as an n-by-d numeric matrix, where 'n' is the samples size
%          (equal to 'nPoints_total') and 'd' is the dimensionality of the
%          point pattern.
%
% spacing - the distance between nearboring points in the grid, before
%           jittering. Formatted as a positive scalar number.
%
%
% AUTHORSHIP: 
%
% Author: Andrew M. Solitsz
% Contact: andysoltisz@gmail.com
%

    % input validation
    if ~ismatrix(studyArea) || ~isnumeric(studyArea)
        error("The study area must be formatted as an n-by-2 numeric matrix.");
    end
    [nDims, nCols] = size(studyArea);
    if nCols ~= 2
        error("The study area must be formatted as an n-by-2 numeric matrix.");
    end
    if ~isscalar(nPoints) || ~isnumeric(nPoints)
        error("Point pattern sample size must be a positive integer.");
    end
    isint = @(y) (y - round(y)) == 0;
    if ~isint(nPoints)
        error("Point pattern sample size must be a positive integer.");
    end
    
    % parse study area
    studyArea_min = studyArea(:, 1);
    studyArea_max = studyArea(:, 2);
    studyArea_size = studyArea_max - studyArea_min;
    studyArea_volume = prod(studyArea_size);

    % calculate point pattern parameters
    spacing = (studyArea_volume / nPoints) .^ (1/nDims); % spacing between points in every dimension
    nPoints_dim = round(studyArea_size ./ spacing); % number of points along each dimension

    % calculate point coordinates along each dimension
    gridCoordinates = cell(1, nDims); % preallocate
    for iDim = 1:nDims
        gridCoordinates{iDim} = linspace(studyArea_min(iDim), studyArea_max(iDim), nPoints_dim(iDim));
    end

    % generate euclidean coordinates for every point in the pattern
    points = combvec(gridCoordinates{:})';

    if nargin == 3
        % input validation
        if jitter_stdev < 0  || ~isscalar(jitter_stdev)
            error("Jitter standard deviation must be a positive number.");
        end

        % apply Guassian positional jitter with mean=0 and user-defined
        % standard deviation
        nPoints_actual = size(points, 1);
        jitter = jitter_stdev .* randn(nPoints_actual, nDims);
        points = points + jitter;
    end

end