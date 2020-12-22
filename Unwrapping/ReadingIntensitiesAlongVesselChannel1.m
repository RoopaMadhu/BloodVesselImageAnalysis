function [] = ReadingIntensitiesAlongVesselChannel1(MainMoviePath,ChannelName)

%% This function unwraps the DAPI and ecad channels of the blood vessels.

% This function reads intensities from every pixel along the surface of the 
% DAPI channel of the blood vessel (obtained by watershed segmentation of 
% phalloidin channel) and interpolates intensities along the surface (DAPI
% channel at constant angular bins.
% The interpolated intensities are mapped onto a 2D plane with axial length
% of blood vessel along the columns of the 2D matrix and circumferential
% angular bins along the rows. Angular bins determined while unwrapping
% phalloidin are used here.

%% INPUT: path of the movie directory

%% OUTPUT: 2D maps are stored in the movie directory as .fig, .tif and .mat files

% 05/06/2017 Roopa Madhu

%% Beginning of the main function

% filter size. If only performing spatial filtering filtersize should be a
% single value which is the Gaussian sigma, if spatial and temporal
% filtering include 3 values for the sigma in each dimension, e.g., 
%filtersize = [sigmaX sigmaY sigmaT].

filtersize = 2;

cd(MainMoviePath);

cd(sprintf('%s',ChannelName));

md = cd;

%load the dataimglist and Data2Analyse workspaces from Channel2
load([num2str(MainMoviePath),'/Channel3/dataimglist.mat']);
load([num2str(MainMoviePath),'/Channel3/Data2Analyse.mat'],'Data2Analyse','angbinsize','cap');

tstart = data.SegTimeInterval(1);
tstop  = data.SegTimeInterval(2);

clear data;

% load img list from Channel1
load([num2str(md),'/dataimglist.mat']);
imglist = data(1).ImageFileList;

bins = -180:angbinsize:180;
IntAlongVessel=zeros(tstop,length(-180:angbinsize:180));


for t = tstart:tstop
    sprintf('Frame%d',t)
    rimg = imread(imglist{t});
    medfiltimg = medfilt2(rimg);
    smoothedImage = filterImage3DpaddedEdges(medfiltimg,'Gauss',filtersize);
    [imgsz1,imgsz2] = size(smoothedImage);
    xsmooth = Data2Analyse(t).smoothedxy(:,1);
    ysmooth = Data2Analyse(t).smoothedxy(:,2);
    [xc,yc]=meshgrid(1:imgsz1,1:imgsz2);
    intensity_smoothedSegline = interp2(xc,yc,double(smoothedImage)',xsmooth,ysmooth);
    
    [angtmp,~,xsmooth,ysmooth] = AngleSorting(xsmooth,ysmooth);
    %         [~,~,binidxtmp] = histcounts(angtmp,-180:angbinsize:180);
    v = mean([intensity_smoothedSegline(1),intensity_smoothedSegline(end)]);
    IntAlongVessel(t,:)=interp1(angtmp,intensity_smoothedSegline,bins,'linear',v);
    %         IntAlongVessel(t,:)=accumarray(binidxtmp,intensity_smoothedSegline,[length(bins),1],@mean);
    Data2Analysenew(t).intSegline = IntAlongVessel(t,:);
    
    % finding the intensity along eroded segmentation line so as to use
    % it to reduce noise
    xsmooth_e = Data2Analyse(t).erodedxy(:,1);
    ysmooth_e = Data2Analyse(t).erodedxy(:,2);
    
    intensity_erdSegline = interp2(xc,yc,double(smoothedImage)',xsmooth_e,ysmooth_e);
    [angtmp_eroded,~,~,~] = AngleSorting(xsmooth_e,ysmooth_e);
    
    
    %     [~,~,binidxtmp_eroded] = histcounts(angtmp_eroded,-180:angbinsize:180);
    %     IntAlongErdSeg(t,:)=accumarray(binidxtmp_eroded,intensity_erodedSegline,[length(bins),1],@mean);
    
    ve = mean([intensity_erdSegline(1),intensity_erdSegline(end)]);
    IntAlongErdSeg(t,:)=interp1(angtmp_eroded,intensity_erdSegline,bins,'linear',ve);
    
    
    Data2Analysenew(t).intEroded = IntAlongErdSeg(t,:);
    
    clear xs ys intensity_smoothedSegline center_smooth angtmp binidxtmp rimg smoothedimg xc yc;
end
clear Data2Analyse;
Data2Analyse = Data2Analysenew;
IntAlongVesselCapped = [IntAlongVessel(tstart:tstop,end-cap:end),IntAlongVessel(tstart:tstop,1:end),IntAlongVessel(tstart:tstop,1:cap+1)];
IntAlongErdSegCapped = [IntAlongErdSeg(tstart:tstop,end-cap:end),IntAlongErdSeg(tstart:tstop,1:end),IntAlongErdSeg(tstart:tstop,1:cap+1)];

IntAlongVesselCapped_norm = NonUniformbgCorrection(IntAlongVesselCapped);
IntAlongErdSegCapped_norm = NonUniformbgCorrection(IntAlongErdSegCapped);

imwrite((IntAlongVesselCapped_norm),'IntAlongVesselCapped.tif');
figure
imshow(IntAlongVesselCapped_norm,[]);
saveas(gcf,'IntAlongVesselCapped.fig');

xlabel('Angle','FontSize',16,'FontWeight','bold');
ylabel('Axial Distance','FontSize',16,'FontWeight','bold');
h=gcf;
MovingAxisLabels(h);
saveas(gcf,'IntAlongVesselCapped_axislabels.tif');
saveas(gcf,'IntAlongVesselCapped_axislabels.fig');
close(gcf);

imwrite(uint16(IntAlongErdSegCapped_norm),'IntAlongErdSegCapped.tif');
figure
imshow(IntAlongErdSegCapped_norm,[]);
saveas(gcf,'IntAlongErdSegCapped');

xlabel('Angle','FontSize',16,'FontWeight','bold');
ylabel('Axial Distance','FontSize',16,'FontWeight','bold');
h=gcf;
MovingAxisLabels(h);
saveas(gcf,'IntAlongErdSegCapped_axislabels.tif');
saveas(gcf,'IntAlongErdSegCapped_axislabels.fig');
close(gcf);

%saving data for each movie in respective folders
save('Data2Analyse.mat','Data2Analyse','IntAlong*','angbinsize','cap','-v7.3');
clear Data2Analyse IntAlong*;

end %of the main function


