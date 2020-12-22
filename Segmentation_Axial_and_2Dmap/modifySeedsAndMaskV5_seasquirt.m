
function [modseeds,modmask,ImageSegment,ImageBWlabel,eflag] = modifySeedsAndMaskV5_seasquirt(I,seeds,mask,Imseg,imbwl,titlestr)
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
% 11/08/18 Roopa Madhu

eflag = false;

% rescale image (if not already the case)
I = mat2gray(I);

% initialize the modified seeds and mask to the input seeds and mask
modseeds = seeds;
modmask = mask;

% initialize RGB image where 2nd and 3rd layers will have seed and mask
% information
imgrgb = zeros([size(I),3]);
I = imadjust(I);
imgrgb(:,:,1) = I;

unsatisfactory = true;
while unsatisfactory
    % update translucent RGB image with modseeds and modmask and display
    imgrgb(:,:,2) = 0.5*im2double(modseeds)+0.55*imadjust(I);
    imgrgb(:,:,3) = 0.55*im2double(modmask)+0.5*imadjust(I);
    imshow(imgrgb,'InitialMagnification',240)
    if nargin > 3
        title(titlestr,'FontSize',16)
    end
    % Choose an option for modifying current layer
    choice = menu('Modification:','Continue','Add Point Seeds',...
        'Add Line Seeds','Remove Line Seeds','Add Polygon Seed','Add Polygon Mask','Remove Regions',...
        'Remove Polygon','Remove Spot','Erode Regions','Dilate Regions','Clear All','MaskInverted',...
        'SegmentAndShow','TranslateSeeds','Translate Mask','Quit');
    
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
            BW = imdilate(BW,strel('disk',22));
            modseeds = modseeds | BW;
            
        case 3 % Add line seed
            for l = 1:1
                h = imline;
                %             position = h.getPosition();
                %             x1(l) = position(1,2);
                %             y1(l) = position(1,1);
                %             x2(l) = position(2,2);
                %             y2(l) = position(2,1);
                
                position = wait(h);
                x1(l) = position(1,2);
                y1(l) = position(1,1);
                x2(l) = position(2,2);
                y2(l) = position(2,1);
            end
            for l = 1:1
                BW = lineMask(size(I), [x1(l) y1(l)], [x2(l) y2(l)],14);
                modseeds = modseeds | BW;
            end
            
        case 4 % Remove line seed
            for l = 1:1
                h = imline;
                position = wait(h);
                x1(l) = position(1,2);
                y1(l) = position(1,1);
                x2(l) = position(2,2);
                y2(l) = position(2,1);
            end
            for l = 1:1
                BW = lineMask(size(I), [x1(l) y1(l)], [x2(l) y2(l)],16);
                modseeds(BW) = false;
                modmask(BW) = false;
            end
            
        case 5 % Add Polygon Seed
            for ps = 1:2
                
                B(ps).BW = roipoly();
                title(gca,'Add Polygon Seed','FontSize',16);
            end
            for ps = 1:2
                modseeds = modseeds | B(ps).BW;
            end
            
        case 6 % Add Polygon Mask
            BW = roipoly();
            title(gca,'Add Polygon Mask','FontSize',16);
            modmask = modmask | BW;
        case 7  % Remove Region
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
            
        case 8 % Remove Polygon (removes polygon from both seeds and mask)
            BW = roipoly();
            title(gca,'Remove Polygon','FontSize',16);
            modseeds(BW) = false;
            modmask(BW) = false;
            
            
            
        case 9 % Remove Spot Region
            imshow(imgrgb,'InitialMagnification',240);
            title(gca,'Select Spot to Remove','FontSize',16);
            [y, x] = ginput(1);
            x = round(x); y = round(y);
            BW = false(size(seeds));
            BW(x,y) = true;
            BW = imdilate(BW,strel('disk',4));
            modseeds(BW) = false;
            modmask(BW) = false;
            
        case 10 % Erode Region
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
                    BWerode = imerode(BW,strel('disk',20));
                    modmask(BWerode) = true;
                elseif modseeds(x(k),y(k)) % if seed is clicked erode it only
                    L = bwlabel(modseeds);
                    label = L(x(k),y(k));
                    BW = L == label;
                    modseeds(BW) = false;
                    BWerode = imerode(BW,strel('disk',10));
                    modseeds(BWerode) = true;
                end
            end
            
        case 11 % Dilate Regions
            imshow(imgrgb,'InitialMagnification',240);
            title(gca,'Select Regions to Dilate','FontSize',16);
            [y, x] = ginput;
            x = round(x); y = round(y);
            for k=1:length(x)
                if modmask(x(k),y(k)) % if mask region is clicked dilate it only
                    L = bwlabel(modmask);
                    label = L(x(k),y(k));
                    BW = L == label;
                    BWdilate = imdilate(BW,strel('disk',5));
                    modmask(BWdilate) = true;
                elseif modseeds(x(k),y(k)) % if seed region is clicked dilate it only
                    L = bwlabel(modseeds);
                    label = L(x(k),y(k));
                    BW = L == label;
                    BWdilate = imdilate(BW,strel('disk',5));
                    modseeds(BWdilate) = true;
                end
            end
            
        case 12 % Clear all seeds and mask
            modseeds = false(size(I));
            modmask = false(size(I));
        case 13 % Mask is the inverse of selected region
            BW = roipoly();
            title(gca,'Add Polygon Mask-inverted','FontSize',16);
            modmask = modmask | ~BW;
            
        case 14 % Segment and show results
            [ImageSegment,ImageBWlabel,modseeds,modmask] = wsSegmentSingleImageFunction_seasquirt2d(I,modseeds,modmask);
            figure
            imoverlay(mat2gray(I),imdilate(~ImageSegment,strel('disk',2)));
            pause
            
        case 15 %translate seeds
            imshow(imgrgb,'InitialMagnification',240);
            prompt = 'Enter the x and y components of the translate vector and erosion size: ';
            x = input(prompt);
            modseeds = imtranslate(imerode(modseeds,strel('disk',x(3))),[x(1),x(2)]);
            
            
        case 16 %translate mask
            
            imshow(imgrgb,'InitialMagnification',240);
            prompt = 'Enter the x and y components of the translate vector and erosion size: ';
            x = input(prompt);
            modmask = imtranslate(imerode(modmask,strel('disk',x(3))),[x(1),x(2)]);
            
        
            
            
        case 17 % Exit Function
            eflag = true;
            break
            
            
    end
end

if exist('ImageSegment') == 1 
    
else
    ImageSegment = imgseg;
    ImageBWlabel = imbwl;
end

close all
end
