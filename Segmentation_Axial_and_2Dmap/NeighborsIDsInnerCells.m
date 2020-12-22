function [ CellInfo ] = NeighborsIDsInnerCells( wsimage )
%given a watershed segmentation image, calculate the number of neighbors
%for all the inner cells that aren't in contact with the background
% INPUT: ImageBWlabel_trackT =  watershed segmentation image; interfaces are
%                               zero, background regions are -1
% OUPUT: numNeiHist =   normalized histogram with number of neihjbors for
%                       inner cells
%        cellNStat =    vector with status (touching background or no) and
%                       number of neighbor for each cell, dimension is
%                       [number of cells x 2]
%
% DL 04/01/2014

% edited by RM 10/25/2017


% minimum and maximum index number in the watershed segmentation image
rmin = min(wsimage(:));
rmax = max(wsimage(:));

% initialize matrix with cell neighbor values
CellInfo.cellNStat = nan*zeros(rmax,2);

% loop over all index values
for i=1:rmax
    % make copy of segmentation image where you retain only the current
    % cell of interest
    imCopy = (wsimage==i);
    % if there is a cell with this index...
    if sum(imCopy(:))>0
        % ... make Euclidican distance map
        imED = bwdist(imCopy);
        % identify the pixel positions in the image that have a distance of
        % between 1-2 pixels from the edge of the cell of interest
        ringPos = find((imED>0)&(imED<4));
        % record the cell indices at those positions (i.e. identify which
        % cells touch the cell of interest)
        ringVal = wsimage(ringPos);
        % extract the unique value of those index positions (i.e. throw out
        % multiples)
        ringValU = unique(ringVal);
        % status = whether or not the background (with index = -1) is among
        % the group of touching regions
        status = double((ringValU(1)==-1));
        % record number of touching cells (excluding background)
        nnn = length(find(ringValU>0));
        CellInfo.cellNStat(i,:) = [status nnn];
        CellInfo.neighID{i} = ringValU(ringValU>0);
    end
end

% % identify inner cells (i.e. cells that are not in direct contact with the
% % background)
% innerPos = find(cellNStat(:,1)==0);
% % extract the number of neighbors for these inner cells
% numNei = cellNStat(innerPos,2);
% % make histogram
% edges = 1:1:10;
% numNeiHist = histcounts(numNei,edges,'Normalization','probability');

% bar([1:1:10],numNeiHist);
% xlabel('number of neighbors');
% ylabel('normalized frequency');


end

