# BloodVesselImageAnalysis

##### Softwares used: MATLAB_R2019a

##### Prerequisites: 3D blood vessel imaging data in the form of .lif files; If the data is available in the form of .tif files, it should be arranged such that each .tif file corresponds to a single axial section of the blood vessel. 

##### Examples of use : Codes should be run in the order given below as output from one function is needed for the next function to run.

##### Loading lif files

1. ReadinglifFiles('')

1. Lif2tif(‘’,’’) 

At this point a folder (parent folder) is created for each blood vessel. Each folder has 3 subfolders - Channel1, Channel2 and Channel3 corresponding to DAPI, ecad and phalloidin channels of the images respectively. Axial sections of the blood vessel are stored as a series of .tif files in each of these subfolders. 
Metadata about the movie is stored in MovieInfo.mat which is located in the parent folder

##### Axial segmentation

Axial segmentation is done on the phalloidin stained images. To do this, move to Channel 3 run the following code:

1. data.Source = cd;
1. data.ImageFileList = LoadImageListMultipleMovies_V3(‘Channel3’);  %this function requires user input. Select first and last images of the image list.
1. data.SegTimeInterval = [tstart,tstop]; % tstart and tstop are the first and last axial sections user wants to segment.
1. save dataimglist data;
1. ELSA_3segmentCellsV4_convhullseeds_V2(data,1,checkInterval,tvec,seederosionsize,maskerosionsize); % checkInterval is an integer - refers to the interval when user wants to check the segmentation. If checkInterval is 100 then segmentation is shown and can be modified once every 100th frame - frame 1 to frame 99 undergo automated segmentation with the parameters inputted and are not open for manual intervention.

##### Unwrapping phalloidin channel

By this step user should still be in “Channel3” directory with the variable “data” in the workspace.

1. ReadingIntensitiesAlongVessel_Ver2(data.Source); % This line unwraps the phalloidin channel of the blood vessel.

##### Unwrapping e-cadherin channel

1. cd ../Channel2/
1. data.Source = cd;
1. data.ImageFileList = LoadImageListMultipleMovies_V3(‘Channel2’);  %this function requires user input. Select first and last images of the image list.
1. save dataimglist data;
1. ReadingIntensitiesAlongVesselChannel1(data.Source,’Channel2’)

##### Unwrapping DAPI channel

1. cd ../Channel1/
1. data.Source = cd;
1. data.ImageFileList = LoadImageListMultipleMovies_V3(‘Channel1’);  %this function requires user input. Select first and last images of the image list.
1. save dataimglist data;
1. ReadingIntensitiesAlongVesselChannel1(data.Source,’Channel1’);

##### Segmenting 2D map

1. cd ../Channel3/
2. [seeds,mask,centroids,PhGaussNorm,DAPInorm,imPh2] = segmentingBloodVessel2Dmap_V2(pathtoparentfolder,seederosionsize); % seederosionsize is an integer.
3. [seeds,mask,data] = segmeting2dimages(seeds,mask,PhGaussNorm,centroids);

##### Curvature correction

1. lengthmap_2d_V2(pathtoparentfolder) % parent folder is the main directory - check “loading lif files” section
1. RadiusCorrectionMap(pathtoparentfolder);
1. [MetaSegData,MetaCorrectionData] = CurvatureCorrection(pathtoparentfolder); %this function also stores data in the parent folder. 

###### Quantification
Quantification studies on control (C), BAPN treated (B) and FAKI treated (F) animals: 

1. [MetaC] = MetaAnalysis_CorrectedData_movlist(movlist) %movlist can be “pathtoparentfolder” if only one movie is being analyzed.
1. [MetaR] = MetaAnalysis_CorrectedData_movlist(movlist) %movlist can be “pathtoparentfolder” if only one movie is being analyzed.
1. [MetaF] = MetaAnalysis_CorrectedData_movlist(movlist) %movlist can be “pathtoparentfolder” if only one movie is being analyzed.
1. [MetadataTable] = ExtractPropFromMetaData_wt_ttest(MetaC,MetaR,MetaF)
1. FinalFigures_CBF(MetaC,MetaR,MetaF)
1. OriAnglePolarPlots_weighted(Metadata,movtype,movcolor); 
1. MetaInt = ExtractPhInt(movlist,ChannelName); %for extracting fiber information from a single channel.
1. [FiberInfo] = ExtractPhInfo_All(ctrlmovlist,bapnmovlist,fakmovlist); %for extracting fiber information from all channels.












