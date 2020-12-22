function [imageSegmented,imageBWlabel,seeds_out,mask_out] = wsSegmentSingleImageFunction_seasquirt2d(image,seeds_in,mask_in)
% perform segmentation on single image given the seeds and mask logical
% arrays, and the filtering parameters. The output is the segmented image,
% the label image, the output seeds and mask.

se = 2;
sm=2;
% impose minima for both the watershed seeds and the mask
imageMasked = imimposemin(image, seeds_in | mask_in);

% perform watershed, set all labeled regions to 1, convert to logical (BW)
imageSegmented = watershed(imageMasked,8);
imageSegmented(imageSegmented>0) = 1;
imageSegmented = logical(imageSegmented);

% convert BW to label image
imageBWlabel = bwlabel(imageSegmented,8);

% separate mask labels from seed regions
masklabels = unique(imageBWlabel(mask_in));
mask_out = false(size(seeds_in));
seeds_out = imageSegmented;
for i=1:length(masklabels)
    idx = imageBWlabel == masklabels(i);
    mask_out(idx) = true;
    seeds_out(idx) = false;
    % set mask regions of label matrix to -1
    imageBWlabel(idx) = -1;
end

% erode mask regions so that they don't overlay with interfaces in next frame
SE = strel('disk',se);
mask_out = imerode(mask_out,SE);

% erode seed regions so that they don't overlap with interfaces in next frame
SE = strel('disk',sm);
seeds_out = imerode(seeds_out,SE);

% Note that any areas that touch the edge of the image are set to value
% -1 to mark them as background - this will ensure in the subsequent
% analysis that these cells are not used for direct characterization of
% nodes, since they do not yield full neighborhood information
edgepixels1 = imageBWlabel(:,1:3);
edgepixels2 = imageBWlabel(1:3,:);
edgepixels3 = imageBWlabel(:,(end-2):end);
edgepixels4 = imageBWlabel((end-2):end,:);
%contiguous areas that touch the edge    
edgeAreas = unique( [edgepixels1(:);edgepixels2(:);edgepixels3(:);edgepixels4(:)] );
% set these areas to value -1
for ie=1:length(edgeAreas)
    cvalue = edgeAreas(ie);
    if cvalue>0
        imageBWlabel(imageBWlabel==cvalue) = -1;
    end
end
%In ImageBWlabel, -1--> mask; 0-->seg line
tmp_len = -1:1:length(unique(imageBWlabel))-2;
cellids = unique(imageBWlabel);
if any(cellids ~= tmp_len')
    imageBWlabel=imageBWlabel-1;
    imageBWlabel(imageBWlabel==-1)=0;
    imageBWlabel(imageBWlabel==-2)=-1;

end

end % of subfunction