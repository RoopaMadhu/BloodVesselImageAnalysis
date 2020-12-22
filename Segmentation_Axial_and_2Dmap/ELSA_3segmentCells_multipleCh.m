function [] = ELSA_3segmentCells_multipleCh(data,MovieNum,checkInterval,tvec,seedesize,maskesize)
% this function segments all images specified by the image file list and
% the time vector "tvec" using a seeded watershed technique. Before running
% this function the seeds must be initialized for the first frame. Because
% this function uses the seeds from the previous frame to segment each
% frame you must start tvec at 2 (so it will use the seeds from frame 1).
%
% Input:    data: structure array that contains lists of image files (such as created by the
%                   function ELSA_1loadImageList) and the source to the
%                   images for each movie.
%           MovieNum: Data structure array index that corresponds to the movie to be
%                     processed.
%           tvec: timevec specifying the frames to be analyzed
%                   (e.g. [1:20])
%
% Output:   no function output - results are written into the specified
%           directory
%

% 01/02/2019 Roopa Madhu (R.M)

% record original directory (and return to it at the end)
od = cd;

% filter size. If only performing spatial filtering filtersize should be a
% single value which is the Gaussian sigma, if spatial and temporal
% filtering include 3 values for the sigma in each dimension, e.g., 
%filtersize = [sigmaX sigmaY sigmaT].
filtersize =2;
%  filtersize=[6,6,2];


% extract list for the current movie from the data structure
list = data(MovieNum).ImageFileList;

movdir = data(MovieNum).Source;
cd(movdir);
cd ../Channel2;
data2 = load('dataimglist');
list2 = data2.data.ImageFileList;
movdir2 = cd;
% cd ../Channel1;
% data1 = load('dataimglist');
% list1 = data1.data.ImageFileList;
% movdir1 = cd;

cd(od);

if nargin>3
    if ~isempty(tvec)
        tlvec = tvec;
    end
    if tlvec(2)>tlvec(1)
        tstep = 1;  % time step is 1
    else
        tstep =-1;  % if iterating downward in time time step is -1
    end
end


% loop over all timepoints
for ti = 1:length(tlvec)
    
    t = tlvec(ti);
    
    % display current progress of processing
    fprintf('segmentation @ timepoint %04d',t);
    
    % load image from list and convert to scaled double and perform spatial
    % filtering
    if length(filtersize)==3
        [Itf] = temporalfilterV2(list,t,filtersize);
    else
        image = mat2gray(imread(list{t}) + imread(list2{t}));
        image = mat2gray(imread(list{t}));
        image = medfilt2(image);
        Itf = filterImage3DpaddedEdges(image,'Gauss',filtersize);
    end

    
    % load seeds and mask from last time point
    cd(data(MovieNum).Source);
    cd('SegmentationData');
    lframefoldername = sprintf('frame%04d',t-tstep);
    cd(lframefoldername);
    loadseeds = load('seeds.mat');
    lseeds = loadseeds.seeds;
    loadmask = load('mask.mat');
    lmask = loadmask.mask;
    [xm,ym] = find(lmask==0);
    idx_mask = convhull(xm,ym);
    lmask = poly2mask(ym(idx_mask),xm(idx_mask),size(lmask,1),size(lmask,2));
    lmask = ~imdilate(lmask,strel('disk',maskesize));
    [xs,ys]=find(lseeds);
    idx_convhull = convhull(xs,ys);
    lseeds = poly2mask(ys(idx_convhull),xs(idx_convhull),size(lseeds,1),size(lseeds,2));
    lseeds = imerode(lseeds,strel('disk',seedesize));
%     lmask = imerode(~lseeds,strel('disk',25));
    
    clear xs ys idx_convhull xm ym idx_mask;
    
    % at frames defined by the checkInterval make any necessary changes to
    % the seeds and mask
    if ~logical(mod(t,checkInterval))
        titlestr = ['Frame ',num2str(t)];
        [lseeds,lmask,eflag] = modifySeedsAndMaskV4(Itf,lseeds,lmask,titlestr);
%         lseeds = imerode(lseeds,strel('disk',2));
%         lmask = imerode(lmask,strel('disk',2));
%         
         
        % now save modifications to last time point
        seeds = lseeds;
        mask = lmask;
        save('seeds','seeds')
        save('mask','mask')
        
        % if exit flag is switched from within modifySeedsAndMaskV5, then
        % end function
        if eflag
            break
        end
  
    end 
    

    % move to current timepoint subfolder...
    cd ..
    cframefoldername = sprintf('frame%04d',t);   
    cd(cframefoldername);


    % perform watershed segmentation, which also automatically yields a
    % label matrix and updated seeds and mask for current time point
    [ImageSegment,ImageBWlabel,seeds,mask] = wsSegmentSingleImageV5(Itf,lseeds,lmask);
    
    % if there is an overlap between seeds and mask, then rectify it:
    
    seedsidx = find(seeds);
    maskidx = find(mask);
    [intelements,intidxseeds,intidxmask] = intersect(seedsidx,maskidx);
    if ~isempty(intelements)
        
        oldSarea = bwarea(lseeds);
        newSarea = bwarea(seeds);
        diffSarea = newSarea - oldSarea;
       
        if diffSarea <= -200
            seeds = convhull(seeds);
            mask(intidxmask) = 0;
       
        elseif diffSarea >= 200
            mask = convhull(mask);
            seeds(intidxseeds) = 0;
            
        end %of diffSarea if loop

    end %of intelements if loop
    
    Nzeros = numel(find(ImageSegment==0));
    Nones  = numel(find(ImageSegment));
    erosion=1;
    lmasknew=lmask;
    lseedsnew=lseeds;
    while erosion <=2
    if ( (Nzeros==numel(ImageSegment)) || (Nones==numel(ImageSegment)) )
        
        lmasknew = imerode(lmasknew,strel('disk',6));        
        [ImageSegment,ImageBWlabel,seeds,mask] = wsSegmentSingleImageV5(Itf,lseeds,lmasknew);
        erosion=erosion+1;
    else
        break
    end %of if loop
    end %of while loop
    
    if ( (Nzeros==numel(ImageSegment)) || (Nones==numel(ImageSegment)) )
        
        [ImageSegment,ImageBWlabel,seeds,mask] = wsSegmentSingleImageV5(Itf,loadseeds.seeds,loadmask.mask);
    end
    while (isempty(find(lseeds)) || isempty(find(lmask)))
             lseeds = imerode(lseeds,strel('disk',6));
             lmask = imerode(lmask,strel('disk',8));
             [ImageSegment,ImageBWlabel,seeds,mask] = wsSegmentSingleImageV5(Itf,lseeds,lmask);
     end
    % save results into the current results folder
    save('ImageSegment', 'ImageSegment');
    save('ImageBWlabel', 'ImageBWlabel');
    save('seeds','seeds');
    save('mask','mask');
  
       
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
    
end % of for t

fprintf('\n');
cd(od);
    
end % of function


%% =======================================================================
%
%                               subfunctions
%
%  =======================================================================

function [imageSegmented,imageBWlabel,seeds_out,mask_out] = wsSegmentSingleImageV5(image,seeds_in,mask_in)
% perform segmentation on single image given the seeds and mask logical
% arrays, and the filtering parameters. The output is the segmented image,
% the label image, the output seeds and mask.


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
SE = strel('disk',5);
mask_out = imerode(mask_out,SE);


% erode seed regions so that they don't overlap with interfaces in next frame
SE = strel('disk',25);
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


end % of subfunction