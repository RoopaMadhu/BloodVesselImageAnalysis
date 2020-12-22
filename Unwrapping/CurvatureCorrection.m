function [MetaSegData,MetaCorrectionData] = CurvatureCorrection(movielist)


for i = 1:length(movielist)
    
    cd(movielist{i});dd
    
    
    md =cd;
    
    
    %load ImageBWlabel
    
    if (exist('Channel2/VesselSegmentationData') == 7)
        cd Channel2/;        
    else
        if (exist('Channel3/VesselSegmentationData') == 7)
cd Channel3/;          
        else
            
cd Channel4/;          
        end
    end
            load VesselSegmentationData/SegmentationData/frame0001/ImageBWlabel.mat;
            load VesselSegmentationData/SegmentationData/frame0001/ImageSegment.mat;

               cd(md);

            MetaSegData(i).ImageBWlabel = ImageBWlabel;
            MetaSegData(i).ImageSegment = ImageSegment;
            MetaSegData(i).movpath = cd;
            
save ExploratoryAnalysis/RadiusMap ImageBWlabel ImageSegment -append;

cellid = unique(ImageBWlabel(:));
cellid(cellid<=0)=[];
cellsnext2mask = CellsNextToMask(ImageBWlabel);
cellid(cellsnext2mask) = [];

%load the weight maps

load ExploratoryAnalysis/RadiusMap lenmap2D_CorrectedCapped lenmap2D_RawCapped  len_of_arc_corrected len_of_arc;
load MovieInfo;

weightmap = (lenmap2D_CorrectedCapped - lenmap2D_RawCapped);
weightmapUncapped = (len_of_arc_corrected - len_of_arc);

for c = 1:length(cellid)
    
    MetaCorrectionData(i).AreaCorrection(c) = sum((lenmap2D_CorrectedCapped(ImageBWlabel==cellid(c)))) * (LifImgInfo.Pxsizes.x ^2);
    
end

AreaCorrection = MetaCorrectionData(i).AreaCorrection;
save ExploratoryAnalysis/RadiusMap.mat AreaCorrection weightmap* -append;

end %of for loop over movies

end %of main function
% imshow(weightmapUncapped,[])
% colormap(gca,hot);
% colorbar
% axis on
% xlabel('Circumferential Direction')
% ylabel('Axial Direction')
% title('Weightmap - Uncapped');
% SetFigureDefaults(gca);
% print -dpsc -r300 -painters WeightMaps.ps -append;
% 
% openfig('ExploratoryAnalysis/AxSec



