function  [seeds] = im2seeds_V2(I,filtersize,viewResults,dilsize)
%IM2SEEDS produces point-like watershed seeds for a 2D image based on
%regional minima of the filtered image.
%
% Input:    I:              2D image.
%           filtersize:     Gaussian filter size. Must have the function
%                           "filterImage3DpaddedEdges.m"
%           viewResults:    If set to 1 an overlay of the seeds will be
%                           shown
%
% Output:   seeds = mask of points of minima in the image for seeded
%                   watershed.
%
% 12/18/12 Timothy Vanderleest
% 02/05/2020 Roopa Madhu (modified)



% filter image
[Igf] = filterImage3DpaddedEdges(I, 'Gauss', filtersize);

% find the regional minimum in the filtered image (the output are points).
IRM = imregionalmin(Igf);


% slightly dilate these points
SE2 = strel('disk',dilsize);
seeds = imdilate(IRM,SE2);


% if viewResults = 1 show overlay of the seeds
if nargin >2
    if viewResults
        imoverlay(I,seeds,[1 0 0]);
    end
end


