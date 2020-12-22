
function [seeds,mask,data] = segmeting2dimages(seeds,mask,I,centroids)


% 01/14/2019 Roopa Madhu

cd Channel2/;

% centSeeds = zeros(size(seeds));
% centroids = round(centroids);
% 
% for i = 1:length(centroids)    
% centSeeds(centroids(i,2),centroids(i,1)) = 1;
% end
% 
% seeds = imdilate(centSeeds,strel('disk',50));

I = NonUniformIllumCorrection(I,2);

[seeds,mask] = modifySeedsAndMaskV5_seasquirt(I,seeds,mask);

[ImageSegment,ImageBWlabel,seeds,mask] = wsSegmentSingleImageFunction_seasquirt2d(I,seeds,mask);

[seeds1,mask1] = modifySeedsAndMask_DAPIseeds(I,seeds,mask,centroids,1);

[ImageSegment,ImageBWlabel,seeds,mask] = wsSegmentSingleImageFunction_seasquirt2d(I,seeds1,mask1);

[seeds1,mask1] = modifySeedsAndMaskV5_seasquirt(I,seeds,mask);

[ImageSegment,ImageBWlabel,~,~] = wsSegmentSingleImageFunction_seasquirt2d(I,seeds1,mask1);

seeds = seeds1;
mask = mask1;

%As we are dealing with a single frame
%In ImageBWlabel, -1--> mask; 0-->seg line
tmp_len = -1:1:length(unique(ImageBWlabel))-2;
cellids = unique(ImageBWlabel);
if any(cellids ~= tmp_len')
    ImageBWlabel=ImageBWlabel-1;
    ImageBWlabel(ImageBWlabel==-1)=0;
    ImageBWlabel(ImageBWlabel==-2)=-1;

end

mkdir VesselSegmentationData;
cd VesselSegmentationData/;

mkdir SegmentationData/frame0001;

save SegmentationData/frame0001/seeds seeds;
save SegmentationData/frame0001/mask mask;
save SegmentationData/frame0001/ImageSegment ImageSegment;
save SegmentationData/frame0001/ImageBWlabel ImageBWlabel;

imwrite(mat2gray(I),'im_001_000.tif');

data(1).Source = cd;
a = dir('*.tif');
data(1).ImageFileList = string(strcat(a.folder,'/',a.name));
save dataimglist data;

ELSA_52extractNodes_V2_2(data(1).ImageFileList,1);
[matRes,trackingMatrixZT,trackingVectorLI]=gridAnalysis_interfaceOrientationNewZT6_SeaSquirt(data(1).Source,[],data(1).ImageFileList,'');
save trackingMats trackingMatrixZT trackingVectorLI matRes;


end %of the main function