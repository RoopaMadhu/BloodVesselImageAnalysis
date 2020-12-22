function [] = ReadingIntensitiesAlongVessel_Ver2(Source)

%% This function unwraps the phalloidin channel of the blood vessels.


% This function reads intensities from every pixel along the surface of the 
% blood vessel (obtained by watershed segmentation) and interpolates 
% intensities along the surface at constant angular bins.
% The interpolated intensities are mapped onto a 2D plane with axial length
% of blood vessel along the columns of the 2D matrix and circumferential
% angular bins along the rows. Angular bins are chosen so that the length
% of the arc subtended by the angular bins is approximately equal to the
% pixel size along the axial direction of blood vessel.

%% INPUT: path of the movie directory

%% OUTPUT: 2D maps are stored in the movie directory as .fig, .tif and .mat files

% 12/28/2018 Roopa Madhu

%% Beginning of the main function

% filter size. If only performing spatial filtering filtersize should be a
% single value which is the Gaussian sigma, if spatial and temporal
% filtering include 3 values for the sigma in each dimension, e.g., 
%filtersize = [sigmaX sigmaY sigmaT].

filtersize = 2;


cd(Source);
md=cd;

%load the pixel sizes information
load ../MovieInfo.mat;
pxsz = LifImgInfo.Pxsizes.x;
sgolayfsize = round(23/pxsz); %Savtisky-golay filtersize is chosen to be 23um
while mod(sgolayfsize,2)==0 %if filter size is even, make it odd
    sgolayfsize = sgolayfsize + 1;
end

% load the image list 
load([num2str(Source),'/dataimglist.mat']);
imlist=data.ImageFileList;
tstart = data.SegTimeInterval(1);
tstop  = data.SegTimeInterval(2);

% load Data2Analyse.mat;

cd('SegmentationData');


Data2Analyse = struct('originalxy',[],'smoothedxy',[],'intensitysmoothed',[],'circum',[]);

%% Determining the appropriate angular bin size

for t = tstart : tstop
    
    load(sprintf('frame%04d/ImageSegment.mat',t));
    
    %     finding segementation line coordinates
    [x,y]=find(ImageSegment==0);
    
    % finding the center of the axial section and angle of each point on the surface from the center
    [~,~,xsorted,ysorted] = AngleSorting(x,y);
    
    % smoothing and getting rid of sharp edges
    xwrap = [xsorted;xsorted];
    ywrap = [ysorted;ysorted];
    xwrapsmooth=sgolayfilt(xwrap,3,sgolayfsize); ywrapsmooth=sgolayfilt(ywrap,3,sgolayfsize);
    veclen = size(xsorted);
    ai = round(veclen/2);
    xsmooth = xwrapsmooth(ai:ai+veclen-1);
    ysmooth = ywrapsmooth(ai:ai+veclen-1);
    
    [~,~,xsmooth,ysmooth] = AngleSorting(xsmooth,ysmooth);
    
    % tmp_angbwpoints=atan2d((ys-center(2)),(xs-center(1))) - atan2d(([ys(2:end);ys(1)]-center(2)),([xs(2:end);xs(1)]-center(1)));
    tmp_Distbwpoints = sqrt(diff([xsmooth(1);xsmooth]).^2+diff([ysmooth(1);ysmooth]).^2);
    tmp_DistFromOrigin = cumsum(tmp_Distbwpoints); %dist between a starting point on the cell and all the other points on the cell
    Data2Analyse(t).circum = round(tmp_DistFromOrigin(end));
end

clear t;

angbinsize = 360/mean([Data2Analyse.circum]);
bins = -180:angbinsize:180;

IntAlongVessel=zeros(tstop,length(-180:angbinsize:180));

for t=tstart:tstop
    sprintf('Frame%d',t)
    %load the raw image
    fullimgname = imlist{t};
    rawimg=imread(fullimgname);
    [imgsz1,imgsz2] = size(rawimg);
    
    load(sprintf('frame%04d/ImageSegment.mat',t));
    %finding segmentation line coordinates
    [x,y]=find(ImageSegment==0);
    
    %% finding the center of the axial section and angle of each point on the surface from the center
    [ang1,idx1,xsorted,ysorted] = AngleSorting(x,y);
    
    %% smoothing and getting rid of sharp edges
    
    xwrap = [xsorted;xsorted];
    ywrap = [ysorted;ysorted];
    xwrapsmooth=sgolayfilt(xwrap,3,sgolayfsize); ywrapsmooth=sgolayfilt(ywrap,3,sgolayfsize);
    
    veclen = size(xsorted);
    ai = round(veclen/2);
    xsmooth = xwrapsmooth(ai:ai+veclen-1);
    ysmooth = ywrapsmooth(ai:ai+veclen-1);
    
    [angtmp,~,xsmooth,ysmooth] = AngleSorting(xsmooth,ysmooth);
    
    % finding the intensity along eroded segmentation line so as to use
    % it to reduce noise
    tmp_erodedSegline = imerode(poly2mask(ysmooth,xsmooth,imgsz1,imgsz2),strel('disk',4));
    [xe,ye]=find(boundarymask(tmp_erodedSegline));
    [~,~,xe,ye] = AngleSorting(xe,ye);
    xsmooth_e = sgolayfilt(xe,3,sgolayfsize); ysmooth_e = sgolayfilt(ye,3,sgolayfsize);
    
    [angtmp_e,~,xsmooth_e,ysmooth_e] = AngleSorting(xsmooth_e,ysmooth_e);
    
    
    % finding the intensity along dilated segmentation line so as to use
    % it to increase signal
    tmp_DilatedSegline = imdilate(poly2mask(ysmooth,xsmooth,imgsz1,imgsz2),strel('disk',1));
    [xd,yd]=find(boundarymask(tmp_DilatedSegline));
    [~,~,xd,yd] = AngleSorting(xd,yd);
    xsmooth_d = sgolayfilt(xd,3,sgolayfsize); ysmooth_d = sgolayfilt(yd,3,sgolayfsize);
    
    
    
    
    [angtmp_d,~,xsmooth_d,ysmooth_d] = AngleSorting(xsmooth_d,ysmooth_d);
    
    
    %Saving sgolayfilted values, unsgolayfilted values and idx of sorted values
    Data2Analyse(t).originalxy=[xsorted,ysorted];
    Data2Analyse(t).smoothedxy=[xsmooth,ysmooth];
    Data2Analyse(t).erodedxy=[xsmooth_e,ysmooth_e];
    Data2Analyse(t).dilatedxy=[xsmooth_d,ysmooth_d];
    
    %         Data2Analyse(i).idxorder=idx;
    
    %% extracting intensities along the interpolated line from spatially sgolayfilted image - sgolayfilting is done based on pixel size
    medfilteredimg = medfilt2(rawimg,[3,3]);
    smoothedImage = filterImage3DpaddedEdges(medfilteredimg,'Gauss',filtersize);
    
    [xc,yc]=meshgrid(1:size(smoothedImage,1),1:size(smoothedImage,2));
    
    intensity_smoothedSegline = interp2(xc,yc,double(smoothedImage)',xsmooth,ysmooth);
    intensity_erodedSegline = interp2(xc,yc,double(smoothedImage)',xsmooth_e,ysmooth_e);
    intensity_dilatedSegline = interp2(xc,yc,double(smoothedImage)',xsmooth_d,ysmooth_d);
    
    v = mean([intensity_smoothedSegline(1),intensity_smoothedSegline(end)]);
    ve = mean([intensity_erodedSegline(1),intensity_erodedSegline(end)]);
    vd = mean([intensity_dilatedSegline(1),intensity_dilatedSegline(end)]);
    
    tmp_intSeg=interp1(angtmp,intensity_smoothedSegline,bins,'linear',v);
    tmp_intErd=interp1(angtmp_e,intensity_erodedSegline,bins,'linear',ve);
    tmp_intDil=interp1(angtmp_d,intensity_dilatedSegline,bins,'linear',vd);
    
    
    Data2Analyse(t).intensitysmoothed = tmp_intSeg;
    Data2Analyse(t).intensityEroded = tmp_intErd;
    Data2Analyse(t).intensityDilated = tmp_intDil;
    
    IntAlongVessel(t,:) = tmp_intSeg-tmp_intErd+tmp_intDil;
    
    
    %             atmp2=accumarray(binidxtmp,intensity_smoothedSegline,[length(bins),1],@median);
    %             btmp2=accumarray(binidxtmp_eroded,intensity_erodedSegline,[length(bins),1],@median);
    %             IntAlongVessel2(t,:) = atmp2-btmp2;
    
    
    %                         for k = 1:length(-180:angbinsize:180)
    %
    %                              IntAlongVessel(t,k) = mean(intensity_smoothedSegline(binidxtmp==k)) - mean(intensity_erodedSegline(binidxtmp_eroded==k));
    %                         end
    
    
    %         clear the variables created for each frame
    clear -regexp ^x ^y ^dx ^dy ^new *tmp fullimgname  IntEroded ^intensity_ ^v ^tmp
    
end %of t-th loop
cd(md);%move to the movie directory to store the output as a workspace

%         IntAlongVessel(IntAlongVessel<0)=0;
% Capping the intensity vector and saving as .tif image so it can be segmented
cap = round(size(IntAlongVessel,2)/3);
IntAlongVesselCapped = [IntAlongVessel(tstart:tstop,end-cap:end),IntAlongVessel(tstart:tstop,1:end),IntAlongVessel(tstart:tstop,1:cap+1)];

IntAlongVesselCapped_norm = NonUniformbgCorrection(IntAlongVesselCapped);

imwrite(uint8(IntAlongVesselCapped_norm),'IntAlongVesselCapped.tif');
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

%saving data for each movie in respective folders
save('Data2Analyse.mat','Data2Analyse','IntAlong*','angbinsize', 'cap','-v7.3');
clear Data2Analyse IntAlongVessel IntAlongVesselCapped IntAlongVesselCapped_norm;

end %of the main function


