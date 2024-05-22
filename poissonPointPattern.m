function [points] = poissonPointPattern(studyArea, nPoints, randSeed)
% Create a Poisson point pattern with sample size of 'nPoints'. The
% dimensionality of the resulting point pattern will match the
% dimensionality of 'studyArea'.
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

    if nargin == 3
        % set random seed if one is provided. Use to create same random
        % data each function call, useful for replication and validation.
        rng(randSeed, "twister"); s = rng; rng(s); 
    end

    nDims = size(studyArea, 1);
    studyArea_min = studyArea(:, 1)';
    studyArea_max = studyArea(:, 2)';
    points = rand(nPoints, nDims) .* (studyArea_max - studyArea_min) + studyArea_min;

end
