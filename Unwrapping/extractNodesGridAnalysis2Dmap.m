
function [] = extractNodesGridAnalysis2Dmap(source)

%source is the path to channel2 where IntAlongVesselCapped.tif is stored.

% record original directory (and return to it at the end)
od = cd;

cd(source);

imglist = repmat({fullfile(source,'IntAlongVesselCapped.tif')},3,1);
copyfile('VesselSegmentationData/frame0001','VesselSegmentationData/frame0002');
copyfile('VesselSegmentationData/frame0001','VesselSegmentationData/frame0003');
ELSA_52extractNodes_V2(imglist,1:3);
[matRes,trackingMatrixZT,trackingVectorLI]=gridAnalysis_interfaceOrientationNewZT6(source,[],imglist,[1 3]);
saveas(gcf,'GridAnalysisOP');
close(gcf);
IV1=[trackingMatrixZT(:,4:8:end),trackingMatrixZT(:,3:8:end)];
IV2=[trackingMatrixZT(:,4:8:end)+trackingMatrixZT(:,6:8:end),trackingMatrixZT(:,3:8:end)+trackingMatrixZT(:,5:8:end)];
save trackingMats trackingMatrixZT trackingVectorLI matRes IV1 IV2;

cd(od);

end %of the main function
