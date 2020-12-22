function [xplot,yplot,zplot,idxmat,faces] = TranformationTo3D_V2()

%% GAME PLAN:

% 1.
% 2.
% 
% load Channel3/Data2analyse.mat Data2Analyse cap angbinsize;
% load Channel3/dataimglist data;

load Channel3/Data2analyse2.mat Data2Analyse cap angbinsize;
load Channel3/dataimglist data;

tstart = data(1).SegTimeInterval(1);
tend = data(1).SegTimeInterval(2);
%finding total length of ang bins before being capped (this is the approx.
%circumference of axial section)
len_orig = length(-180:angbinsize:180);


% load node positions
load('Channel2/VesselSegmentationData/SegmentationData/frame0001/dstruct_nodes.mat');
nodePos_2D = round(dstruct_nodes.positions(:,1:2));

AxSecNumVec = nodePos_2D(:,1); CappedAngbinNumVec = nodePos_2D(:,2);


idxvec = [len_orig-cap:len_orig,1:len_orig,1:cap+1];
angvec = -180:angbinsize:180;

%% finding the actual x,y,z coordinates of the whole 2D map

for i = tstart:tend
    x = Data2Analyse(i).smoothedxy(:,1);
    y = Data2Analyse(i).smoothedxy(:,2);
    
    [angs,~,~,~] = AngleSorting(x,y);
    %interpolating x and y at required ang intervals
    xv = mean([x(1),x(end)]);
    yv = mean([y(1),y(end)]);
    xactual(i,:) = interp1(angs,x,angvec,'linear',xv);
    yactual(i,:) = interp1(angs,y,angvec,'linear',yv);
    zactual(i,:) = repmat(i,length(angvec),1);
    clear x y xv yv angs;
end

xactual = xactual(tstart:tend,:);
yactual = yactual(tstart:tend,:);
zactual = zactual(tstart:tend,:);

xactual = xactual(:,idxvec);
yactual = yactual(:,idxvec);
zactual = zactual(:,idxvec);


%% Find the actual 3D coordinates of Vertices

nv = length(nodePos_2D);
for i = 1:nv
    Vert3D_x(i) = xactual(AxSecNumVec(i),CappedAngbinNumVec(i));
    Vert3D_y(i) = yactual(AxSecNumVec(i),CappedAngbinNumVec(i));
    Vert3D_z(i) = zactual(AxSecNumVec(i),CappedAngbinNumVec(i));
end

%% Find the centroid coordinates

load Channel2/VesselSegmentationData/SegmentationData/frame0001/ImageBWlabel.mat

load Channel2/VesselSegmentationData/SegmentationData/frame0001/dstruct_cellCentroids.mat
Cen_Pos = dstruct_cellCentroids.positions;
clear AxSecNumVec CappedAngbinNumVec;

ncells = length(Cen_Pos);

CellsNextToMask = [];
for i = 1:ncells
    dilatedCell = imdilate(ImageBWlabel==i,strel('disk',2));
    CellNeighbors = unique(ImageBWlabel(dilatedCell));
    if ~isempty(find(CellNeighbors==-1)) %checking if any of the neighbors are mask(labelled zeros)
        CellsNextToMask = [CellsNextToMask i];
    end
    clear dilatedCell CellNeighbors;
end

Cen_Pos(CellsNextToMask,:)=[];

AxSecNumVec = Cen_Pos(:,2);
CappedAngbinNumVec = Cen_Pos(:,1);

ncells = length(Cen_Pos);

for i = 1:ncells
    Cen_x(i) = xactual(round(AxSecNumVec(i)),round(CappedAngbinNumVec(i)));
    Cen_y(i) = yactual(round(AxSecNumVec(i)),round(CappedAngbinNumVec(i)));
    Cen_z(i) = zactual(round(AxSecNumVec(i)),round(CappedAngbinNumVec(i)));
end

% DuplicateNodeIdx = FindingDuplicateNodes_V2(Cen_x,Cen_y,Cen_z,Cen_Pos);

load('Channel2/VesselSegmentationData/SegmentationData/frame0001/dstruct_nodeNodeMat.mat');

NodeNodeMat = full(dstruct_nodeNodeMat.matrix);
load Channel2/VesselSegmentationData/SegmentationData/frame0001/dstruct_nodeCellMat.mat
NodeCellMat = full(dstruct_nodeCellMat.matrix);
NodeCellMat = NodeCellMat(:,1:end-1);
NodeCellMat(:,CellsNextToMask) = [];

%% Identifying duplicate nodes


%% Generating an idxmat with the node connections in pairs

idxmat = [];
for i = 1:length(NodeNodeMat)
    nodeConnection = find(NodeNodeMat(i,:));
    if ~isempty(nodeConnection)
        k=1;
        while k <= length(nodeConnection)
            idx_tmp(k,:) = [ i nodeConnection(k)];
            k=k+1;
        end
        idxmat = [idxmat;idx_tmp];
        
    end
    clear nodeConnection k idx_tmp;
end

for i = 1:size(idxmat,1)
    xplot(i,:) = [Vert3D_x(idxmat(i,1)),Vert3D_x(idxmat(i,2))];
    yplot(i,:) = [Vert3D_y(idxmat(i,1)),Vert3D_y(idxmat(i,2))];
    zplot(i,:) = [Vert3D_z(idxmat(i,1)),Vert3D_z(idxmat(i,2))];
end



%% Delete the duplicate cells
% cells with centroids less than 5px away from each other are considered
% duplicate and one of them are deleted

%% Not working!!
% c2 = combnk(1:ncells,2);
% 
% centdiff = Centmat(c2(:,1),:) - Centmat(c2(:,2),:);
% 
% diffid = nanmean(abs(centdiff),2);
% 
% combid = find(diffid<=5);
% 
% cid = c2(combid,:);

%%

faces = [];
cells2del=[];
for i = 1:ncells
    
    vid = find(NodeCellMat(:,i));
    t = NodeNodeMat(vid,:);
    t = t(:,vid);
    t = tril(t);
    [tmp1,tmp2]=(find(t));
    if (length(tmp1) ~= length(tmp2))
        cells2del = [cells2del,i];
    else
        tmp_faces(:,1) = vid(tmp1);
        tmp_faces(:,2) = vid(tmp2);
        tmp_faces(:,3) = nv+i; %for centroid of the cell
        faces = [faces;tmp_faces];
        clear vid tmp1 tmp2 tmp_faces t;
    end
    
end %of for loop
%%


Centmat = [Cen_z',Cen_y',Cen_x'];
Vertmat = [Vert3D_z',Vert3D_y',Vert3D_x'];
newVertmat = [Vertmat;Centmat];

%% figures

cmapnew=jet(ncells);
cmapnew = cmapnew(randperm(length(cmapnew)),:);
for i = 1:size(faces,1)
    cmapfaces(i,:) = cmapnew(faces(i,3)-nv,:);
end

% cmapfaces2 = cmapfaces(randperm(length(cmapfaces)),:);
% for i = 1:size(faces,1)
%     cmapfaces2(i,:) = cmapnew(faces(i,3)-nv,:);
% end

%% each cell is a different color - all the triangles in a cell are of the same color

mkdir ExploratoryAnalysis;
figure
h=trisurf(faces,newVertmat(:,1),newVertmat(:,2),newVertmat(:,3));
set(h,'FaceVertexCData',cmapfaces,'FaceAlpha',0.8,'EdgeAlpha',1)
axis equal
hold on
h2=plot3(zplot',yplot',xplot','LineWidth',2,'Color','k');
pause(1)
saveas(gcf,'ExploratoryAnalysis/Trisurf2.fig')
close(gcf);
% 
save ExploratoryAnalysis/Trisurfdata2 -regexp plot$ Cen* Vert* CellsNextToMask newVertmat faces nv ncells;

 
clearvars -except MovieList m wd;
end %of main function







