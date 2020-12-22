function [] = ReadingIntensitiesAlongVessel_ErDil(Source,imlist)

filtersize = 2;

% sgolayfsize = 311;

cd(Source);
md=cd;

load ../MovieInfo.mat;
pxsz = LifImgInfo.Pxsizes.x;
dilsize = round(1/pxsz);

sgolayfsize = round(23/pxsz);
while mod(sgolayfsize,2)==0 %if filter size is even, make it odd
    sgolayfsize = sgolayfsize + 1;
end

load([num2str(Source),'/dataimglist.mat']);
tstart = data.SegTimeInterval(1);
tstop  = data.SegTimeInterval(2);

load Data2Analyse.mat Data2Analyse angbinsize;

cd('SegmentationData');


bins = -180:angbinsize:180;

IntAlongVessel=zeros(tstop,length(-180:angbinsize:180));

for t=tstart:tstop
 sprintf('Frame%d',t)
    %load the raw image
    fullimgname = imlist{t};
    rawimg=imread(fullimgname);
    [imgsz1,imgsz2] = size(rawimg);
    
   [xsmooth,ysmooth] =  deal(Data2Analyse(t).smoothedxy(:,1),Data2Analyse(t).smoothedxy(:,2));
    
       medfilteredimg = medfilt2(rawimg,[3,3]);
    smoothedImage = filterImage3DpaddedEdges(medfilteredimg,'Gauss',filtersize);
    
    [xc,yc]=meshgrid(1:size(smoothedImage,1),1:size(smoothedImage,2));

   for i = 1:5 %looping to read intensities 1 to 5 microns above and below ph signal
   
       dsz = dilsize*i;
       
        % finding the intensity along eroded segmentation line so as to use
    % it to reduce noise
    tmp_erodedSegline = imerode(poly2mask(ysmooth,xsmooth,imgsz1,imgsz2),strel('disk',dsz));
    [xe,ye]=find(boundarymask(tmp_erodedSegline));
    [~,~,xe,ye] = AngleSorting(xe,ye);
    xsmooth_e = sgolayfilt(xe,3,sgolayfsize); ysmooth_e = sgolayfilt(ye,3,sgolayfsize);
    
    [angtmp_e,~,xsmooth_e,ysmooth_e] = AngleSorting(xsmooth_e,ysmooth_e);
    
    
    % finding the intensity along dilated segmentation line so as to use
    % it to increase signal
    tmp_DilatedSegline = imdilate(poly2mask(ysmooth,xsmooth,imgsz1,imgsz2),strel('disk',dsz));
    [xd,yd]=find(boundarymask(tmp_DilatedSegline));
    [~,~,xd,yd] = AngleSorting(xd,yd);
    xsmooth_d = sgolayfilt(xd,3,sgolayfsize); ysmooth_d = sgolayfilt(yd,3,sgolayfsize);
    
    
    
    
    [angtmp_d,~,xsmooth_d,ysmooth_d] = AngleSorting(xsmooth_d,ysmooth_d);
    
    
    %Saving sgolayfilted values, unsgolayfilted values and idx of sorted values
   
    Data2Analysen(t,i).erodedxy=[xsmooth_e,ysmooth_e];
    Data2Analysen(t,i).dilatedxy=[xsmooth_d,ysmooth_d];
    
    
        %% extracting intensities along the interpolated line from spatially sgolayfilted image - sgolayfilting is done based on pixel size
    
    intensity_erodedSegline = interp2(xc,yc,double(smoothedImage)',xsmooth_e,ysmooth_e);
    intensity_dilatedSegline = interp2(xc,yc,double(smoothedImage)',xsmooth_d,ysmooth_d);
    
    ve = mean([intensity_erodedSegline(1),intensity_erodedSegline(end)]);
    vd = mean([intensity_dilatedSegline(1),intensity_dilatedSegline(end)]);
    
    tmp_intErd=interp1(angtmp_e,intensity_erodedSegline,bins,'linear',ve);
    tmp_intDil=interp1(angtmp_d,intensity_dilatedSegline,bins,'linear',vd);
    
    
    Data2Analysen(t,i).intensityEroded = tmp_intErd;
    Data2Analysen(t,i).intensityDilated = tmp_intDil;
 
   end %of ith for loop
   
   
   
end %of t-th for loop

save Data2AnalyseErDil Data2Analysen;
end %of main function
   


