    function [] = Lif2tif(LifList,SeriesList)
    
    %This program reads the mentioned z-stacks from the mentioned .lif workspaces and writes the .tif files for each of the
    %.lif workspaces
    
    %% This code can be run from various locations with following rules:
    
    % 1. This code can be run from outside the directory containing LifImgName workspaces
    
    % 2. The list of lifworkspaces (LifList) can be given as an :
    % a.) INPUT - in this case, each workspace can be from a different directory (as path names are given as input,it doesnt effect the
    % code) (if full path names of workspaces are not given, they are considered to be in the current working direcotry)
    %  OR
    % b.) SELECTED BY THE USER THROUGH A DIALOG BOX - in this case, all the workspaces should be in the same directory even if they are
    % outside the current working directory
    
    
    
    %% INPUTS:
    
    % LifList      - (O/P of ReadingLifFiles.m - .mat files)a list of full path names(if full path names are not given, they are considered to be in the current direcotry)
    %                         of .lif workspaces(names of workspace which are outputs of ReadingLifFiles.m)
    % SeriesList   - a list of names of the Image of series(should be 3D images) that you want to write into .tif files
    
    %% OUTPUTS:
    
    % Folders for movies are created and .tif files are written into appropriate movie folders.
    
    
    %% Main function:
    
    wd = cd; %save the current directory as working directory
    mode=0;
    %If the LifImgNames are not given as an input, then ask the user to select the files
    
    if isempty(LifList)
        
        [LifFileNames,LifPaths] = uigetfile('*.mat','Select the .mat workspace','MultiSelect','on');
        LifFileNames = cellstr(LifFileNames);
        LifPaths = cellstr(LifPaths);
        %getting rid of .mat part of filenames
        for i = 1 : length(LifFileNames)
            [LifPaths{i},LifFileNames{i} ]= fileparts(LifFileNames{i});
        end %of for loop
        
        % tmp_a = dir('*.mat');
        % for i = 1:length(tmp_a)
        % [LifPaths{i},LifFileNames{i}] = fileparts(tmp_a(i).name);
        % end
        
    else
        
        
        for i = 1 : length(LifList)
            [LifPaths{i},LifFileNames{i}] = fileparts(LifList{i});
        end %of for loop
        
        
    end %of if loop
    
    % Creating an LifImgPath array based on the location of workspaces
    % some of the LifPaths can be empty which means they are in the working
    % directory. In such cases replace the empty cells with path to the
    % current directory
    for i = 1 : length(LifPaths)
        if strcmp(LifPaths{i},'')
            LifPaths{i} = cd;
        end
    end %of for loop
    
    
    %%
    for i = 1 : length(LifFileNames)
        
        Workspace = load([LifPaths{i} filesep LifFileNames{i}]); % load the lif workspace
        [LifImgInfo] = LifInfo(Workspace.lif2mat); %find the names of zstacks stored in each field of the loaded lif workspace
        % if SeriesList is empty then convert all the series into tif files
        if isempty(SeriesList)
            SeriesList = LifImgInfo.Name;
            mode=1;
        end % of if loop
        
        
        for j = 1 : length(SeriesList)
            movdirName = [LifFileNames{i},'_',SeriesList{j}];
            movdirfullName = ([LifPaths{i} filesep movdirName]);
            dircheck = isdir(movdirfullName);
            % make a movie directory if it doesnt exist already
            if dircheck == 0
                mkdir(movdirfullName);
            end % of if loop
            
            cd(movdirfullName); %change to the movie directory. This is needed as imwrite is not taking paths for storing images
            Zstack_idx = find( strcmp(LifImgInfo.Name,SeriesList{j}) ); %find the index number of the zstack (that is to be written into .tif) in lif workspace
            
            Pxratio = LifImgInfo.Pxsizes(Zstack_idx).z / LifImgInfo.Pxsizes(Zstack_idx).y;
            
            % if these are dual channel images, create two folders and
            % write each channel into a seperate folder
            nChannels = length(Workspace.lif2mat(Zstack_idx).imgStruct.Image);
            if nChannels>1
                for c = 1:nChannels
                    
                    ChannelNames{c} = ['Channel',num2str(c)];
                    mkdir(ChannelNames{c});
                    cd(ChannelNames{c});
                    
                    Zstack_img = Workspace.lif2mat(Zstack_idx).imgStruct.Image{c};%retrieve the 3D img from lif workspace
                    dim_rimg = size(Zstack_img);
                    
                    %Write the 3D image into .tif files
                    for z = 1:size(Zstack_img,1)
                        tmp_rimg = squeeze(Zsta2dck_img(z,:,:));
                        tmp_resizedimg = imresize(tmp_rimg,[dim_rimg(2),dim_rimg(3)*Pxratio]);
                        imwrite(tmp_resizedimg,sprintf('im_%05d_000.tif',z));
                        clear tmp*;
                    end % of z-for loop
                    clear Zstack_img ; %clearing the 3D img before retrieving a new one
                    cd(movdirfullName); %move back to the main movie directory
                end %of c loop - channels
            end %end of if loop - checking in nChannels i s>1
            
            save MovieInfo LifImgInfo;
        end % of j loop - movies
        
        clear Workspace ; %clearing the workspace before loading a new one
        
        % clear the SeriesList if it is not given as an input. Just to be make sure it is not being appended.
        if mode ==1
            SeriesList = '';
        end
        
    end % of i-for loop
    
    cd(wd) %go back to the working directory
    
end %of main function
