function[MetaData]=LoadImageListMultipleMovies(movienumseq,CommonMainFolderName,MovieName)
%% This function is a workflow for loading imagefile list of multiple movies
% This function presumes each movie has multiple subfolders referring
% to different channels
%% INPUT:
%movienumseq - a vector of numbers
%CommonMovieName - Only one entry
%MovieName - name of the channel. Only one entry

%% MAIN FUNCTION:
od=cd;

%% specify the first and last images in a stack for all the movies
for i=1:length(movienumseq)
    for j = 1:length(MovieName)
        %print the movie number being analysed
        sprintf('Specifying Images : %s_%d_%s',CommonMainFolderName,movienumseq(i),MovieName{j})
        
        cd(sprintf('%s%03d',CommonMainFolderName,movienumseq(i)));
        cd(MovieName{j}); %changing to the specified channel you want to load images from
        
        %% user specifies the first image
        
        %         [Files(i,j).fileNameI1, Files(i,j).filePathI1] = uigetfile('.tif','choose first image of first sub-stack (e.g. t000_i000)');
        
        %%complete file name of the first image
        %        Files(i,j).oriImageName1 = strcat(Files(i,j).filePathI1, Files(i,j).fileNameI1);

        %% SELECTING FIRST IMAGE IN THE LIST AS FIRST IMAGE
        
        tmp = dir('*.tif');
        
%         Files(i,j).fileNameI1 = tmp(1).name;
%         
%         Files(i,j).filePathI1 = tmp(1).folder;
        
        % when the folder contains new images apart from raw images in .tif format:
        
        Files(i,j).fileNameI1 = 'im_00001_000.tif';
        
        Files(i,j).filePathI1 = tmp(1).folder;
        
        % complete file name of the first image
        Files(i,j).oriImageName1 = strcat(Files(i,j).filePathI1,'/', Files(i,j).fileNameI1);
        
        %% load second z-stack for subsequent time point
        
          %%  user specifies the second image
       %        [Files(i,j).fileNameI2, Files(i,j).filePathI2] = uigetfile('.tif','choose first image of second sub-stack (e.g. t001_i000)');

        %%complete file name of the first image
        %        Files(i,j).oriImageName2 = strcat(Files(i,j).filePathI2, Files(i,j).fileNameI2);

         %% SELECTING LAST IMG FROM THE LIST AS SECOND IMAGE
%         Files(i,j).fileNameI2 = tmp(end).name;
%         
%         Files(i,j).filePathI2 = tmp(end).folder;

         % when the folder contains new images apart from raw images in .tif format:
 
        Files(i,j).fileNameI2 = 'im_05184_000.tif';
        
        Files(i,j).filePathI2 = tmp(end).folder;
        
                % complete file name of the first image
        Files(i,j).oriImageName2 = strcat(Files(i,j).filePathI2,'/', Files(i,j).fileNameI2);
        
        clear tmp
        %%
     
        cd(od);
    end %of jth loop
    
end % of ith loop


for i = 1:length(movienumseq)
    for j = 1:length(MovieName)
        
        %print the movie number being analysed
        sprintf('Loading Image list : Movie_%d',movienumseq(i))
        
        cd(sprintf('%s%03d',CommonMainFolderName,movienumseq(i)));
        cd(MovieName{j}); %changing to the specified channel you want to load images from
        data(1).ImageFileList=ELSA_1loadImageList_modified('','',Files(i,j));
        data(1).Source=cd;
        MetaData(i,j).ImageFileList = data(1).ImageFileList;
        MetaData(i,j).Source=cd;
        
        save dataimglist data;
        clear data;
        cd(od); %move back to working  directory
        
    end %of jth loop
    
end %end of ith loop

save(sprintf('MetaDataImgList_%s_%s',CommonMainFolderName,MovieName{j}), 'MetaData');
end %of main funciton


