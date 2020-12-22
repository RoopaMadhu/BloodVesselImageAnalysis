function [modseeds,modmask,eflag] = modifySeedsAndMaskV2_RM(I,seeds,mask,titlestr)
%MODIFYSEEDSANDMASK allows a user to easily make changes to watershed seeds
%and background mask.
%
% INPUTS:       I:      The image for verifying seeds and masks.
%               seeds:  Watershed seeds (2D binary).
%               mask:   Mask of background region or regions (2D binary).
%               titlestr: title to be displayed on the figure (optional).
%
% OUTPUTS:      modseeds: Modified seeds.
%               modmask:  Modified mask.
%               eflag: true if user selects the 'Quit' option. eflag used to
%               break from the parent function.
%
% 12/18/12 Timothy Vanderleest

% 03/02/18 Modified by Roopa Madhu

eflag = false;

% rescale image (if not already the case)
I = mat2gray(I);

% initialize the modified seeds and mask to the input seeds and mask
modseeds = seeds;
modmask = mask;

% % for all seeds on the border, make point-like
% L = bwlabel(modseeds);
% STATS = regionprops(L,'Centroid');
% 
% % get edge labels
% edgepixels1 = L(:,1);
% edgepixels2 = L(1,:);
% edgepixels3 = L(:,end);
% edgepixels4 = L(end,:);
% %contiguous areas that touch the edge    
% edgeAreas = unique( [edgepixels1(:);edgepixels2(:);edgepixels3(:);edgepixels4(:)] );
% edgeAreas(edgeAreas == 0) = [];
% 
% BW = false(size(seeds));
% for i = 1:length(edgeAreas)
%     iL = edgeAreas(i);
%     
%     centroid = round(STATS(iL).Centroid);
%     BW(centroid(2),centroid(1)) = true;
%     modseeds(L == iL) = false;
% end
% % modseeds = modseeds | imdilate(BW,strel('disk',8));



% initialize RGB image where 2nd and 3rd layers will have seed and mask
% information
imgrgb = zeros([size(I),3]);
I = filterImage3DpaddedEdges(I,'Gauss',1);
imgrgb(:,:,1) = imadjust(I);


unsatisfactory = true;
while unsatisfactory
    % update translucent RGB image with modseeds and modmask and display
    imgrgb(:,:,2) = 0.5*im2double(modseeds)+0.55*imadjust(I);
    imgrgb(:,:,3) = 0.5*im2double(modmask)+0.55*imadjust(I);
    imshow(imgrgb,'InitialMagnification',240)
    if nargin > 3
        title(titlestr,'FontSize',16)
    end
    

    % Choose an option for modifying current layer
    choice = menu('Modification:','Continue','Add Point Seeds',...
        'Add Polygon Seed','Add Polygon Mask','Remove Regions',...
        'Remove Polygon','Remove Spot','Erode Regions','Dilate Regions','Remove All','MaskInverted','Quit');
    
    switch choice
        case 1 % Save and Quit
            unsatisfactory = false;
            
        case 2 % Add Point Seeds (can continue adding until user hits enter). 
            title(gca,'Add Point Seeds','FontSize',16);
            [y, x] = ginput;
            x = round(x); y = round(y);
            BW = false(size(I));
            for k=1:length(x)
                BW(x(k),y(k)) = true;
            end
            % slightly dilate so seed is not single pixel
            BW = imdilate(BW,strel('disk',40));
            modseeds = modseeds | BW;
            
        case 3 % Add Polygon Seed
            BW = roipoly();
            title(gca,'Add Polygon Seed','FontSize',16);                
            modseeds = modseeds | BW;
            
        case 4 % Add Polygon Mask    
            BW = roipoly();
            title(gca,'Add Polygon Mask','FontSize',16);                
            modmask = modmask | BW;
        case 5  % Remove Region 
            title(gca,'Remove Region','FontSize',16);
            [y, x] = ginput;
            x = round(x); y = round(y);
            
            for k=1:length(x)
                if modmask(x(k),y(k))
                    L = bwlabel(modmask);
                    label = L(x(k),y(k));
                    idx = L == label;
                    modmask(idx) = false;
                elseif modseeds(x(k),y(k))
                    L = bwlabel(modseeds);
                    label = L(x(k),y(k));
                    idx = L == label;
                    modseeds(idx) = false;
                end
            end
            
        case 6 % Remove Polygon (removes polygon from both seeds and mask)
            BW = roipoly();
            title(gca,'Remove Polygon','FontSize',16);                
            modseeds(BW) = false;
            modmask(BW) = false;
            

            
        case 7 % Remove Spot Region
            imshow(imgrgb,'InitialMagnification',240); 
            title(gca,'Select Spot to Remove','FontSize',16);
            [y, x] = ginput(1);
            x = round(x); y = round(y);
            BW = false(size(seeds));
            BW(x,y) = true;
            BW = imdilate(BW,strel('disk',8));
            modseeds(BW) = false;
            modmask(BW) = false;

        case 8 % Erode Region
            imshow(imgrgb,'InitialMagnification',240); 
            title(gca,'Select Regions to Erode','FontSize',16);
            [y, x] = ginput;
            x = round(x); y = round(y);
            for k=1:length(x)
                if modmask(x(k),y(k)) % if mask region is clicked erode it only
                    L = bwlabel(modmask);
                    label = L(x(k),y(k));
                    BW = L == label;
                    modmask(BW) = false;
                    BWerode = imerode(BW,strel('disk',1));
                    modmask(BWerode) = true;
                elseif modseeds(x(k),y(k)) % if seed is clicked erode it only
                    L = bwlabel(modseeds);
                    label = L(x(k),y(k));
                    BW = L == label;
                    modseeds(BW) = false;
                    BWerode = imerode(BW,strel('disk',1));
                    modseeds(BWerode) = true;
                end
            end

        case 9 % Dilate Regions
            imshow(imgrgb,'InitialMagnification',240); 
            title(gca,'Select Regions to Dilate','FontSize',16);
            [y, x] = ginput;
            x = round(x); y = round(y);
            for k=1:length(x)
                if modmask(x(k),y(k)) % if mask region is clicked dilate it only
                    L = bwlabel(modmask);
                    label = L(x(k),y(k));
                    BW = L == label;
                    BWdilate = imdilate(BW,strel('disk',1));
                    modmask(BWdilate) = true;
                elseif modseeds(x(k),y(k)) % if seed region is clicked dilate it only
                    L = bwlabel(modseeds);
                    label = L(x(k),y(k));
                    BW = L == label;
                    BWdilate = imdilate(BW,strel('disk',1));
                    modseeds(BWdilate) = true;
                end
            end
        case 10 % Clear all seeds and mask
            modseeds = false(size(I));
            modmask = false(size(I));
        case 11 % Mask is the inverse of selected region
            BW = roipoly();
            title(gca,'Add Polygon Mask-inverted','FontSize',16);                
            modmask = modmask | ~BW;
        case 12 % Exit Function
            eflag = true;
            break
                

    end
end



close all
end

