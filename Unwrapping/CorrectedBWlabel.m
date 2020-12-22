function [MetaData_Corrected] = CorrectedBWlabel(movielist)


% Load ImageBWlabel and weightmap and do poly2mask on each cell and find out
% the height and widht of each cell and also do regionprops


for i = 1:length(movielist)
    
    cd(movielist{i});
    
    
    md =cd;
    
    %load ImageBWlabel and lenmap2D_Capped - this is the correction map
    %where value at each pixel represents how much the actual length of the
    %arc is. The approximate length of the arc is 1. So the values in
    %lenmap2D_Capped are centered around 1.
    
    load ExploratoryAnalysis/RadiusMap.mat lenmap2D_CorrectedCapped;
    load ExploratoryAnalysis/2Dmaps ImageBWlabel;
    load MovieInfo;
%     load ExploratoryAnalysis/CellIDInfo;
    ImageBWlabel = ImageBWlabel.ImageBWlabel;
    
    pxsize = LifImgInfo.Pxsizes.x;
%     
%         cellsnext2mask = CellsNextToMask(ImageBWlabel);
%         cellID = unique(ImageBWlabel(:));
%         cellID(cellID<=0) = [];
%         cellID(cellsnext2mask) = [];
%     regpropstmp = regionprops(ImageBWlabel,'Area');
%         save the cellID for future use
%         save ExploratoryAnalysis/CellIDInfo cellID cellsnext2mask;
    
    for c = 1:length(cellID)
        
        [x,y] = find(ImageBWlabel == cellID(c));
        
        idxold = sub2ind(size(ImageBWlabel),x,y);
        [xnew,ynew] = rescaledata(x,y);
        
        szx = max(xnew) + 100;
        szy = max(ynew) + 100;
        
        singlecellmap = ones(szx,szy) * (-1);
        
        idxnew = sub2ind([szx,szy],xnew,ynew);
        
        singlecellmap(idxnew) = lenmap2D_CorrectedCapped(idxold);
        
        %convert the single cell map to binary bw label:
        singlecellbw = singlecellmap~=-1;
        
        Area_lenmap(c) = sum(singlecellmap(idxnew)) *(pxsize^2);
        
        % O/P of bwboundaries when used to create mask using poly2mask
        % doesnt create the same mask as the original due to the way
        % bwboundaries and poly2mask work.
        %Refer https://blogs.mathworks.com/steve/2014/03/27/comparing-the-geometries-of-bwboundaries-and-poly2mask/
        
        % Upsizing the bwlabel before using bwboundaries helps deal with
        % this discrepancy - each pixel is now represented by 3 pixels
        
        %         singlecellbw_resized = imresize(singlecellbw,3,'nearest');
        %         boundarypoints_resized = bwboundaries(singlecellbw_resized);
        %         boundarypoints_resized = cell2mat(boundarypoints_resized);
        
        % now scale back down to get the boundary points of the original
        % (not resized) single cell map:
        %         boundarypoints = [boundarypoints_resized(:,1)+1 , boundarypoints_resized(:,2)+1]./3;
        
        
        boundarypoints = bwboundaries(singlecellbw);
        boundarypoints = cell2mat(boundarypoints);
        
        boundarypoints_idx = sub2ind(size(singlecellmap),boundarypoints(:,1),boundarypoints(:,2));
        
        %taking cumulative sum along rows:
        singlecell_cumsum = cumsum(abs(singlecellmap),2);
        %because cumsum gives cumsum from beginning to the end of each row,
        %setting the elements outside the cell to be zeros:
        singlecell_cumsum(singlecellbw==0)=0;
        
        boundarypoints_mapval = singlecell_cumsum(boundarypoints_idx);
        
%         AreaCorrected(c) = polyarea(boundarypoints_mapval,boundarypoints(:,1)) * (pxsize^2);
        
        %generating a corrected bwlabel for each cell:
        
        singlecellbwlabel = poly2mask(boundarypoints_mapval,boundarypoints(:,1),szx,szy);
        
        %         imshowpair(singlecellbwlabel,ImageBWlabel==cellID(c))
        %         pause;
        
        ImgBWlabels(c).bwlabel = singlecellbwlabel;
        
        clear boundary*
        
        
    end %of looping over cells
    
    save ExploratoryAnalysis/CorrectedCellInfo ImgBWlabels -append;
    
%     CorrectedCellInfo_2D = findregprops(ImgBWlabels,pxsize);
%     
%     save ExploratoryAnalysis/CorrectedCellInfo CorrectedCellInfo_2D AreaCorrected Area_lenmap;
%     
%     MetaData_Corrected(i).CellInfo = CorrectedCellInfo_2D;
%     
%     clear ImgBWlabels CorrectedCellInfo_2D AreaCorrected Area_lenmap;
    
    
end %of for loop over movies

end %of main function








function [xnew,ynew] = rescaledata(x,y)

%% This function rescales data from 10 to range of old data;

%finding the range of data
xmin_old = min(x); xmax_old = max(x);
ymin_old = min(y); ymax_old = max(y);

xrange_old = (xmax_old-xmin_old);
yrange_old = (ymax_old - ymin_old);

%setting a new range
xmin_new = 100; xmax_new = xmin_new+xrange_old;
ymin_new = 100; ymax_new = ymin_new+yrange_old;

xrange_new =  xrange_old; yrange_new = yrange_old;

%rescaling x and y coordinates
xnew = ( (x - xmin_old) * (xrange_old / xrange_new) ) + xmin_new;
ynew = ( (y - ymin_old) * (yrange_old / yrange_new) ) + ymin_new;
end

function [CellInfo_2D] = findregprops(bwlabelstruct,pxsz)

for l = 1:length(bwlabelstruct)
    
    
    CellInfo_2D = regionprops('struct',bwlabelstruct(l).bwlabel,'Area','Perimeter','MajorAxisLength','MinorAxisLength','Eccentricity','Orientation');
    
        if length(CellInfo_2D) > 1
    
            [~,idx] = max([CellInfo_2D.Area]);
    
            CellInfo_2D = CellInfo_2D(idx);
        end
    
    % change Area into um^2

CellInfo_2Dnew.Area(l) = CellInfo_2D.Area*(pxsz^2);
    
    % change Perimeter into um
    CellInfo_2Dnew.Perimeter(l) = CellInfo_2D.Perimeter*pxsz;
    
    % change MajorAxisLength into um
    CellInfo_2Dnew.MajorAxisLength(l) = CellInfo_2D.MajorAxisLength*pxsz;
    
    % change MinorAxisLength into um
    CellInfo_2Dnew.MinorAxisLength(l) = CellInfo_2D.MinorAxisLength*pxsz;
    
    % change Orientation range such that positive x-axis is 0 deg and positive
    % y-axis is 90 deg
    tmpOr = CellInfo_2D.Orientation;
    tmpOridx = find(tmpOr<0);
    % tmpOr(tmpOridx)=abs(tmpOr(tmpOridx));
    tmpOr(tmpOridx) = tmpOr(tmpOridx)+180;
    CellInfo_2Dnew.Orientation(l) = tmpOr;
    
    CellInfo_2Dnew.Eccentricity(l) = CellInfo_2D.Eccentricity;
    
    [CellInfo_2Dnew.ht(l),CellInfo_2Dnew.wd(l)] = cellhtwd(bwlabelstruct(l).bwlabel);
    
    
    
end %end of for loop over bwlabels/cells

[CellInfo_2Dnew.ht] = [CellInfo_2Dnew.ht].*pxsz;
[CellInfo_2Dnew.wd] = [CellInfo_2Dnew.wd].*pxsz;

%I am doing this just because I want to use the name "CellInfo_2D"
%as the saved variable - there are functions which use this for
%downstream processing
CellInfo_2D = CellInfo_2Dnew;

end %of subfunction findregprops





