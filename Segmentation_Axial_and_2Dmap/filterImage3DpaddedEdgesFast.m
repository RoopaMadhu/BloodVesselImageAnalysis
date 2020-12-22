function [ImageMatrixFiltered] = filterImage3DpaddedEdgesFast(image, filtername, filtersize)
% filter image in 3 dimension; pad edges to prevent boundary effects
% INPUT:    image       = 2D or 3D image
%           filtername  = filter method (see below) 
%           filtersize  = filter size (see below)
%          
%           OPTIONS:
%           filtername = 'Gauss' Gaussian filter
%           filtersize = sigma of Gauss, can be either 
%           [sigma] or [sigmax sigmay sigmaz];
%
%           filtername = 'ra' rolling average
%           filtersize = width of rolling average (+-ra) 
%           [width] or [widthx widthy widthz];
%
%           NOTE 1: If image is 2-dimensional, the filtering will be
%           performed two-dimensional, as well. If the image is
%           3-dimensional, then there are two options: If filtersize is of
%           length 1, then the filtering will be performed in a
%           two-dimensional way on successive layers of the image; if
%           filtersize has length 3 [sigx sigz sigz], the filtering will be
%           3-dimensional, but it will be performed sequentially, i.e. the
%           2D filtering will be performed first, and the result will be
%           filtered in the z-dimension
%       
%           NOTE 2: The image matrix is padded before filtering to avoid 
%           edge effects; the padded pixel layer corresponds to filtersize 
%           in that dimension


%% determine filter size
if length(filtersize)==1
    fsx = filtersize; fsy = filtersize; fsz = 0;
    ftype = 2;
else
    fsx = filtersize(1); fsy = filtersize(2); fsz = filtersize(3);
    ftype = 3;
end

%% define filter

switch filtername
    case 'ra'
        fdimx = 2*fsx+1; 
        fdimy = 2*fsy+1; 
        fdimz = 2*fsz+1;
        
        cfilter     = ones(fdimx,fdimy,fdimz);
        cfilterNorm = cfilter/sum(cfilter(:));
        
    case 'Gauss'
        sigmax = fsx; sigmay = fsy; sigmaz = fsz;
        if fsz==0, sigmaz = inf; end
        
        fsxfull = ceil(3*fsx); 
        fsyfull = ceil(3*fsy);  
        fszfull = ceil(3*fsz); 
        fdimx   = 2*fsxfull+1; 
        fdimy   = 2*fsyfull+1; 
        fdimz   = 2*fszfull+1;
        
        % xy-filter
        cfilterXY     = zeros(fdimx,fdimy);
        for i = -fsxfull:fsxfull
            ex = exp( -(i^2) / (2*sigmax^2) );
            for j = -fsyfull:fsyfull
                ey = exp( -(j^2) / (2*sigmay^2) );
                cfilterXY(i+fsxfull+1, j+fsyfull+1) = ex * ey;
            end
        end
        cfilterNormXY     = cfilterXY/sum(cfilterXY(:)); 
        
        %z-filter
        cfilterZ     = zeros(fdimz,1);   
        for k= -fszfull:fszfull
            ez = exp( -(k^2) / (2*sigmaz^2) );
            cfilterZ(k+fszfull+1) = ez;
        end
            
        cfilterNormZ     = cfilterZ/sum(cfilterZ(:)); 
        
        
        padx = (fdimx-1)/2; 
        pady = (fdimy-1)/2; 
        padz = (fdimz-1)/2; 
        
end


%% construct padded image

% size of image
[sx,sy,sz] = size(image);
% size of padded image
sxp = sx + 2*padx;
syp = sy + 2*pady;
szp = sz + 2*padz;

ImageMatrixPadded = zeros(sxp,syp,szp);

% fill center
ImageMatrixPadded(padx+1:padx+sx, pady+1:pady+sy,padz+1:padz+sz) = image;

% pad x (rows) - pad top rows, then pad bottom rows
ImageMatrixPadded(1:padx,:,:) = repmat(ImageMatrixPadded(padx+1,:,:),[padx 1 1]);
ImageMatrixPadded(padx+sx+1:sxp,:,:) = repmat(ImageMatrixPadded(padx+sx,:,:),[padx 1 1]);

% pad y (columns) - pad left columns, then pad right columns
ImageMatrixPadded(:,1:pady,:) = repmat(ImageMatrixPadded(:,pady+1,:),[1 pady 1]);
ImageMatrixPadded(:,pady+sy+1:syp,:,:) = repmat(ImageMatrixPadded(:,pady+sy,:),[1 pady 1]);

% pad z (layers) - pad top layers, the pad bottom layers
ImageMatrixPadded(:,:,1:padz) = repmat(ImageMatrixPadded(:,:,padz+1),[1 1 padz]);
ImageMatrixPadded(:,:,padz+sz+1:szp) = repmat(ImageMatrixPadded(:,:,padz+sz),[1 1 padz]);


%% filter image

% if image is 3D, filter all z-layers individually

fprintf(' filter XY...');
for iz = 1:szp
    % filter layers individually
    cimage = ImageMatrixPadded(:,:,iz);
    cimage_filter = convn(cimage,cfilterNormXY);
    ImageMatrixPaddedFiltered(:,:,iz) = cimage_filter;
end
fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b');
% if 3D-filtering is desired, the subdequently also filter the third
% dimension

if (sz>1) 
    if (ftype==3)
        
        fprintf(' filter Z...');
        ImageMatrixPaddedFiltered2 = ImageMatrixPaddedFiltered;
        % loop over all central layers
        for iz = fszfull+1:szp-fszfull
            % filtering layers
            cimagemat = ImageMatrixPaddedFiltered(:,:,iz-fszfull:iz+fszfull);
            for iz2 = 1:size(cimagemat,3)
                cimagematF(:,:,iz2) = cfilterNormZ(iz2)*cimagemat(:,:,iz2);
            end
            ImageMatrixPaddedFiltered2(:,:,iz) = sum(cimagematF,3);
                
        end
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b');
       
       ImageMatrixPaddedFiltered = ImageMatrixPaddedFiltered2;
    end
    % fprintf('\n');
end

% crop out center; note that x-y padding is increased by the filterinf by
% another padx, pady (hence 2*padx in the cropping), whereas tjhe padding
% in the z-direction stays padz
ImageMatrixFiltered = ImageMatrixPaddedFiltered(2*padx+1:2*padx+sx, 2*pady+1:2*pady+sy, padz+1:padz+sz);

% figure; imshow(ImageMatrix(:,:,1),[]);
% figure; imshow(ImageMatrixFiltered(:,:,1),[]);

end

