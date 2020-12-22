function [MetamapC,MetadenC] = Fig5e(Metadata,movtype)


PrefOr = PrefrentialOri_V2(Metadata);

threshval = 0.2;
proprange = -15:30:195;
% find axially oriented vs isotropic vessles

axorid = (PrefOr(:,2) > threshval);
circorid = (PrefOr(:,1) > threshval);
isotorid = (PrefOr(:,1) <= threshval & PrefOr(:,2) <= threshval);


if ~isempty(find(circorid))
    
    [MetamapC,MedCurvC] = areadiafun_propinterp(Metadata(circorid),5);
    MetadenC = areadiafun_propinterpgraphs(MetamapC,proprange);
    circinfomean = nanmean(cat(3,MetadenC.den),3);
    circinfostd = nanstd(cat(3,MetadenC.den),3);
    circ_mean = [circinfomean(1,1),circinfomean(4,1);circinfomean(1,5),circinfomean(4,5)];
    circ_std = [circinfostd(1,1),circinfostd(4,1);circinfostd(1,5),circinfostd(4,5)];
    makeheatmap(MetamapC,MetadenC,movtype,'Circumferentially')
end

if ~isempty(find(axorid))
    
    [MetamapC,MedCurvC] = areadiafun_propinterp(Metadata(axorid),5);
    MetadenC = areadiafun_propinterpgraphs(MetamapC,proprange);
    meandist = nanmean(cat(3,MetadenC.den),3);
    %     [~,p1] = kstest2(meandist(:,1),rand(1,1000))
    %     [~,p2] = kstest2(meandist(:,end),rand(1,1000))
    makeheatmap(MetamapC,MetadenC,movtype,'Axially')
    circinfomean = nanmean(cat(3,MetadenC.den),3);
    circinfostd = nanstd(cat(3,MetadenC.den),3);
    circ_mean = [circinfomean(1,1),circinfomean(4,1);circinfomean(1,5),circinfomean(4,5)];
    circ_std = [circinfostd(1,1),circinfostd(4,1);circinfostd(1,5),circinfostd(4,5)];
    
end

if ~isempty(find(isotorid))
    
    [MetamapC,MedCurvC] = areadiafun_propinterp(Metadata(isotorid),5);
    MetadenC = areadiafun_propinterpgraphs(MetamapC,proprange);
    
    makeheatmap(MetamapC,MetadenC,movtype,'Isotropically');
    circinfomean = nanmean(cat(3,MetadenC.den),3);
    circinfostd = nanstd(cat(3,MetadenC.den),3);
    circ_mean = [circinfomean(1,1),circinfomean(4,1);circinfomean(1,5),circinfomean(4,5)];
    circ_std = [circinfostd(1,1),circinfostd(4,1);circinfostd(1,5),circinfostd(4,5)];
    
    
end



end

function [] = makeheatmap(MetamapC,MetadenC,movtype,subtype)
imshow_denheatmaps(nanmean(cat(3,MetadenC.den),3),[0 0.72]);
axis on
xlabel('AspectRatio')
ylabel('Orientation Angle')
title([subtype,' oriented - ',movtype]);
yticks(1:7)
yticklabels(0:30:180)
xticks(1:size(MetadenC(1).den,2))
medprop = nanmean(cat(1,MetamapC.Medianprop));
medprop2 = num2cell(round(medprop,2));
xticklabels(medprop2)
SetFigureDefaults(gca);
set(gcf,'Position',[665,73,490,549]);
subtypenew = subtype(1:3);
pause(1);
% saveas(gcf,sprintf('AspHeatmap_%s_%s',movtype,subtypenew));
% saveas(gcf,sprintf('AspHeatmap_%s_%s.tif',movtype,subtypenew));
% print( gcf, '-painters', sprintf('AspHeatmap_%s_%s.ai',movtype,subtypenew), '-dpsc','-r300');


end