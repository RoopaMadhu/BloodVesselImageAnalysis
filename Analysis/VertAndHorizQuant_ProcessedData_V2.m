function [IntAve] = VertAndHorizQuant_ProcessedData_V2(MoviePathList,imgtype)

wd = cd;

for m = 1:length(MoviePathList)
    
    cd(MoviePathList{m});
    md = cd;
    %     load ExploratoryAnalysis/rimg_snorm.mat;
    load MovieInfo.mat;
    pxsz = LifImgInfo.Pxsizes.x;
    fiberlen = round(5/pxsz);
    if strcmp(imgtype,'ph')==1
        dilsizes(m) = round(1.5/pxsz);
    elseif strcmp(imgtype,'ecad')==1
        dilsizes(m) = round(5/pxsz);
    end
    
        npx_airydiscdia(m) = 0.21/pxsz; %0.21um is the airy disc dia
    filtersize(m) = npx_airydiscdia/6; %(1/3 of airy disc radius)

    %load the imagesegment matrix to delete the interfaces
    load ExploratoryAnalysis/VesselSegmentationData/SegmentationData/frame0001/ImageSegment.mat;
    load ExploratoryAnalysis/VesselSegmentationData/SegmentationData/frame0001/mask.mat;
    
    cd(md)
    
    
    
    [IntAve(m).MeanN] = VertAndHorizQuant_2_v3(ImageSegment,mask,dilsizes(m),fiberlen,imgtype,filtersize(m));
    
    cd(wd);
    clear rimg*;
end %of mth loop


end %of main function
