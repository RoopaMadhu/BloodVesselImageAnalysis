function [seeds,mask] = InitializeSeedsAndMaskbloodVessel(newimg)
Itf = filterImage3DpaddedEdges(newimg,'Gauss',6);

[seeds] = im2seeds(Itf,12,0);
mask = false(size(seeds));


[seeds,mask] = modifySeedsAndMaskV5(Itf,seeds,mask);

end