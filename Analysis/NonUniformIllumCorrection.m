function [rimg_snorm] = NonUniformIllumCorrection(rimg,direction)
%% Correct for non-uniform illumination

if direction == 2
    
    for i=1:size(rimg,2) %if the rimg has vertical bright and dark stripes
        pstart = max(i,i-20);
        pend = min(size(rimg,2),i+20);
        icrop = rimg(:,pstart:pend);
        imin = nanmin(icrop(:));
        imax = nanmean(icrop(:));
        rimg_snorm(:,i) = (rimg(:,i)-imin)./(imax-imin);
    end %of for loop
    
elseif direction ==1
    
    for i=1:size(rimg,1)
        pstart = max(i,i-20);
        pend = min(size(rimg,1),i+20);
        icrop = rimg(pstart:pend,:);
        imin = nanmin(icrop(:));
        imax = nanmean(icrop(:));
        rimg_snorm(i,:) = (rimg(i,:)-imin)./(imax-imin);
    end %of for loop
    
end %of if loop

end %of the main function