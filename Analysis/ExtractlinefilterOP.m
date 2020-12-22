function [MetaOP] = ExtractlinefilterOP(movlist,workspacename,movvec,movtype,wd)

% k=[10:16,20:30];

if isempty(movvec)
    movvec = 1:length(movlist);
end
angvec = 1:12;
for i = movvec
    cd(movlist{i});
    
    load(['ExploratoryAnalysis/',workspacename]);
    
    %     MetaOP(i).sumabs = [OPimg.sum];
    %     if isfield(OPimg,'mean')
    %         MetaOP(i).meanabs = [OPimg.mean];
    %     end
    load ExploratoryAnalysis/VesselSegmentationData/SegmentationData/frame0001/ImageBWlabel;
    
    %     if length(OpenImgn)==37
    %         angvec = 1:3:37;
    %     elseif length(OpenImgn)==12
    %         angvec = 1:13;
    %         OpenImgn(13).img = OpenImgn(1).img;
    %
    %     end
    cells = (unique(ImageBWlabel(:)));
    cells(cells<=0)=[];
    
    k=1;
    for ang = angvec
        %         k=1;
        MetaOP(i).mean(k) = nanmean(abs(OPimg(k).OP(:)));
        for n = 1:length(cells)
            cellmat = OPimg(k).OP(ImageBWlabel==cells(n));
            MetaOP(i).cellphsum(n,k) = (nansum(cellmat(:)));
            MetaOP(i).cellphmean(n,k) = (nanmean(cellmat(:)));
            MetaOP(i).cellphmeanabs(n,k) = nanmean(abs(cellmat(:)));
            MetaOP(i).cellphmeanpos(n,k) = nanmean(nanmean(cellmat(cellmat>0)));
        end
        k=k+1;
        
    end
    
    MetaOP(i).cellphabsratio = MetaOP(i).cellphmeanabs(:,1)./MetaOP(i).cellphmeanabs(:,7);
    MetaOP(i).cellphposratio = MetaOP(i).cellphmeanpos(:,1)./MetaOP(i).cellphmeanpos(:,7);
    MetaOP(i).cellphratio = MetaOP(i).cellphmean(:,1)./MetaOP(i).cellphmean(:,7);
    
    
    
    %     MetaOP(i).maxproj = nanmax(cat(3,OPimg.OP),3);
    
    %         OPimg = loadedVar.OPimg;
    %load movie info and calculate filtersizes
    %     load MovieInfo.mat;
    %     pxsz = LifImgInfo.Pxsizes.x;
    %     npx_airydiscdia = 0.21/pxsz;
    %     filtersize =(npx_airydiscdia/6);
    %     dilsize = round(2.5/pxsz);
    %
    %         load ExploratoryAnalysis/VesselSegmentationData/SegmentationData/frame0001/ImageSegment;
    %     dilatedSeg = imdilate(~ImageSegment,strel('square',dilsize));
    %     angbinsize = OPimg.angbinsize; cap = OPimg.cap;
    %     bins = -180:angbinsize:180;
    %     AngCapped = [bins(end-cap:end),bins(1:end),bins(1:cap+1)];
    %     MetaOP(i).angbinsize = angbinsize;
    %     MetaOP(i).cap = cap;
    %     stripylim1 = find(ismembertol(AngCapped,-45,0.2,'DataScale',1),1);
    %     stripylim2 = find(ismembertol(AngCapped,45,0.2,'DataScale',1),1);
    %
    %     if stripylim2<stripylim1
    %         stripylim2 = find(ismembertol(AngCapped,45,0.2,'DataScale',1),1,'last');
    %     end
    %     stripystrip = [stripylim1:stripylim2];
    %
    %     nonstripylim1 = find(ismembertol(AngCapped,-135,0.2,'DataScale',1),1);
    %     nonstripylim2 = find(ismembertol(AngCapped,-45,0.2,'DataScale',1),1);
    %
    %     if nonstripylim2<nonstripylim1
    %         nonstripylim2 = find(ismembertol(AngCapped,-45,0.2,'DataScale',1),1,'last');
    %     end
    %     nonstripystrip = [nonstripylim1:nonstripylim2];
    %
    %     for j = angvec
    %         sstr = OPimg(j).OP(:,stripystrip);
    %         nsstr = OPimg(j).OP(:,nonstripystrip);
    %         MetaOP(i).stripysumabs(j) = (nansum(abs(sstr(:))));
    %         MetaOP(i).nonstripysumabs(j) = (nansum(abs(nsstr(:))));
    %         MetaOP(i).stripymeanabs(j) = nanmean(abs(sstr(:)));
    %         MetaOP(i).nonstripymeanabs(j) = nanmean(abs(nsstr(:)));
    %
    %         MetaOP(i).stripysum(j) = (nansum((sstr(:))));
    %         MetaOP(i).nonstripysum(j) = (nansum((nsstr(:))));
    %         MetaOP(i).stripymean(j) = nanmean((sstr(:)));
    %         MetaOP(i).nonstripymean(j) = nanmean((nsstr(:)));
    %
    %         MetaOP(i).stripysumpos(j) = nansum(nansum(sstr(sstr>0)));
    %         MetaOP(i).nonstripysumpos(j) = nansum(nansum(nsstr(nsstr>0)));
    %         MetaOP(i).stripymeanpos(j) = nanmean(nanmean(sstr(sstr>0)));
    %         MetaOP(i).nonstripymeanpos(j) = nanmean(nanmean(sstr(sstr>0)));
    %
    %         MetaOP(i).sumpos(j) = nansum(nansum(OPimg(j).OP(OPimg(j).OP>0)));
    %         MetaOP(i).meanpos(j) = nanmean(nanmean(OPimg(j).OP(OPimg(j).OP>0)));
    %
    %         MetaOP(i).sum(j) = nansum(nansum(OPimg(j).OP(:)));
    %         MetaOP(i).mean(j) = nanmean(nanmean(OPimg(j).OP(:)));
    %     end
    %
    %     stripydseg = dilatedSeg(:,stripystrip); nonstripydseg = dilatedSeg(:,nonstripystrip);
    %
    %         stripyop3dmat = cat(3,OPimg.stripy.OP);nonstripyop3dmat = cat(3,OPimg.nonstripy.OP);
    %
    %         stripyop3dmat(repmat(stripydseg,1,1,12))=nan; nonstripyop3dmat(repmat(nonstripydseg,1,1,12))=nan;
    %     %
    %         MetaOP(i).stripysum = nansum(abs(reshape(stripyop3dmat,[],size(stripyop3dmat,3))));
    %         MetaOP(i).nonstripysum = nansum(abs(reshape(nonstripyop3dmat,[],size(nonstripyop3dmat,3))));
    %     MetaOP(i).stripymaxproj2 = nanmax(stripyop3dmat,3);
    %     MetaOP(i).nonstripymaxproj2 = nanmax(nonstripyop3dmat,3);
    %
    %     MetaOP(i).stripysum = [OPimg.stripy.sum];
    %     MetaOP(i).nonstripysum = [OPimg.nonstripy.sum];
    %
    %     MetaOP(i).stripyrimg = filterImage3DpaddedEdges(OPimg.stripy(1).rimg,'Gauss',3*filtersize);
    %     MetaOP(i).nonstripyrimg = filterImage3DpaddedEdges(OPimg.stripy(1).rimg,'Gauss',3*filtersize);
    %
    %     MetaOP(i).stripymaxproj = nanmax(cat(3,OPimg.stripy.OP),3);
    %     MetaOP(i).nonstripymaxproj = nanmax(cat(3,OPimg.nonstripy.OP),3);
    %
    %     for j = 1:12
    %         OPimg.stripy(j)
    %         MetaOP(i).stripysum2(j) = nansum(nansum(OPimg.stripy(j).OP));
    %         MetaOP(i).nonstripysum2(j) = nansum(nansum(OPimg.nonstripy(j).OP));
    %     end
    
    
    % MetaOP(i).cap = OPimg.cap;
    % MetaOP(i).angbinsize = OPimg.angbinsize;
    % MetaOP(k(i)).OP = OPimg;
    % for j = 1:length(OPimg)
    % % MetaOP(i).max(j) = nanmax(OPimg(j).OP(:));
    % % MetaOP(i).min(j) = nanmin(OPimg(j).OP(:));
    %  MetaOP(i).kernel{j} = (OPimg(j).kernel);
    % end
    %% figures
    
    %     cd(wd)
    %     figure(1)
    %     subplot(1,2,1)
    %     imshow(filterImage3DpaddedEdges(OPimg(1).rimg,'Gauss',filtersize),[])
    %     title(['Movie#',num2str(i)]);
    %     subplot(1,2,2)
    %     imshow(MetaOP(i).maxproj,[])
    %     title(['Max.Proj.of OP'])
    %     print(gcf,'-painters',['Doc',movtype,'275.ps'],'-dpsc','-r300','-append','-bestfit');
    %     subplot(1,2,1)
    %     plot(0:15:180,[MetaOP(i).sum,MetaOP(i).sum(1)],'b-','LineWidth',2)
    %     ylim([0,inf])
    %     subplot(1,2,2)
    %     plot(0:15:180,[MetaOP(i).sum,MetaOP(i).sum(1)]./mean(MetaOP(i).sum),'b-','LineWidth',2)
    %     ylim([0,inf])
    %     print(gcf,'-painters',['Doc',movtype,'275.ps'],'-dpsc','-r300','-append','-bestfit');
    
    
end %of mooping of movies - for loop

end %of main function