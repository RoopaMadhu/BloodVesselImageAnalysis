


function [seeds,mask,centroids,PhGaussNorm,DAPInorm,imPh2] = segmentingBloodVessel2Dmap_V2(source,erosionsize)

%smooth the images
cd(source);
load Channel1/Data2Analyse.mat IntAlongVesselCapped;
imDAPI = IntAlongVesselCapped;
DAPInorm = NonUniformIllumCorrection(imDAPI,2);
DAPInorm = medfilt2(DAPInorm);

load Channel2/Data2Analyse.mat IntAlongVesselCapped;
imPh = IntAlongVesselCapped;
Ph1Gauss = filterImage3DpaddedEdges(imPh,'Gauss',10);
% imPh(imPh>100) = 100;
% imPh = medfilt2(imPh);
load Channel3/Data2Analyse.mat IntAlongVesselCapped;
imPh2 =  IntAlongVesselCapped;
% imPh2(imPh2<50)=0;
% imPh2 = medfilt2(imPh2);
Ph2Gauss = filterImage3DpaddedEdges(imPh2,'Gauss',10);
imPh2 = (Ph1Gauss + Ph2Gauss)/2;
% PhGauss = filterImage3DpaddedEdges(imPh,'Gauss',10);
PhGaussNorm = NonUniformIllumCorrection(Ph1Gauss,2);
PhGaussNorm = NonUniformIllumCorrection(PhGaussNorm,1);


DAPI_binary = imbinarize(DAPInorm,'adaptive','Sensitivity',0.1);
seeds=imopen(DAPI_binary,strel('disk',erosionsize));

CCInfo = bwconncomp(seeds);
seedprop = regionprops(CCInfo);

% propfilt = arrayfun( @(x) x.Area > 3200, seedprop);
% seeds_idx = CCInfo.PixelIdxList(propfilt);

% tmpimg = zeros(size(rimg.imDAPI));

centroids = cat(1,seedprop.Centroid);

% for i = 1:length(seeds_idx)
% tmpimg(seeds_idx{i}) = 1;
% end

mask = zeros(size(seeds));



end %of main function







