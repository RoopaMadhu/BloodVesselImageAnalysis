function [modseeds,modmask,eflag] = modifySeedsAndMask_DAPIseeds(I,seeds,mask,centroids,titlestr)
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
% 12/20/2019 Roopa Madhu (modified)


eflag = false;

% rescale image (if not already the case)
I = mat2gray(I);

% initialize the modified seeds and mask to the input seeds and mask
modseeds = seeds;
modmask = mask;

% initialize RGB image where 2nd and 3rd layers will have seed and mask
% information
imgrgb = zeros([size(I),3]);
imgrgb(:,:,1) = I;

unsatisfactory = true;
while unsatisfactory
    % update translucent RGB image with modseeds and modmask and display
    imgrgb(:,:,2) = 0.45*im2double(modseeds)+0.55*I;
    imgrgb(:,:,3) = 0.45*im2double(modmask)+0.55*I;
%     imshow(imgrgb,'InitialMagnification',240)
    imshow(imgrgb);
    set(gca, 'position',[0.01 0.01  0.98  0.98]);
    hold on
    plot(centroids(:,1),centroids(:,2),'ro')
    if nargin > 4
        title(num2str(titlestr),'FontSize',16)
    end
    

    % Choose an option for modifying current layer
    choice = menu('Modification:','Continue','Add Point Seeds',...
        'Add Line Seeds','Remove Line Seeds','Add Polygon Seed','Add Polygon Mask','Remove Regions',...
        'Remove Polygon','Remove Spot','Erode Regions','Dilate Regions','RemoveEverything','Mask Inverted','Quit');
    
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
            BW = imdilate(BW,strel('disk',23));
            modseeds = modseeds | BW;
            
        case 3 % Add line seed
           for l = 1:5
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
           for l = 1:5
            BW = lineMask(size(I), [x1(l) y1(l)], [x2(l) y2(l)],8);
            modseeds = modseeds | BW;
           end
            
        case 4 % Remove line seed
            for l = 1:5
            h = imline;
            position = wait(h);
            x1(l) = position(1,2);
            y1(l) = position(1,1);
            x2(l) = position(2,2);
            y2(l) = position(2,1);
            end
            for l = 1:5
            BW = lineMask(size(I), [x1(l) y1(l)], [x2(l) y2(l)],10);
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
            BW = imdilate(BW,strel('disk',20));
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
                    BWerode = imerode(BW,strel('disk',1));
                    modmask(BWerode) = true;
                elseif modseeds(x(k),y(k)) % if seed is clicked erode it only
                    L = bwlabel(modseeds);
                    label = L(x(k),y(k));
                    BW = L == label;
                    modseeds(BW) = false;
                    BWerode = imerode(BW,strel('disk',6));
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

        case 14 % Exit Function
            eflag = true;
            break
                

    end
end



close all
end

