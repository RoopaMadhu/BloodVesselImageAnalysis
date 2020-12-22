function [] = Fig4(MetaDataC,MetaDataR,movtype1,movtype2)

figure
    

%% 

Maj_C = cat(1,MetaDataC.MajAxlen);
Maj_R = cat(1,MetaDataR.MajAxlen);

Min_R = cat(1,MetaDataR.MinAxlen);
Min_C = cat(1,MetaDataC.MinAxlen);

Wd_C =cat(1,MetaDataC.CellWd);
Wd_R =cat(1,MetaDataR.CellWd);

Ht_C =cat(1,MetaDataC.CellHt);
Ht_R =cat(1,MetaDataR.CellHt);

%%
bar(1:2,[mean(Maj_C);mean(Min_C)])
ylim([0,25])

SetFigureDefaults(gca);
ylabel('Length (um)')
title('Major and Minor Axis lengths - CTRL')
xticklabels({'Major Axis','Minor Axis'});
hold on
errorbar([mean(Maj_C);mean(Min_C)],[std(Maj_C),std(Min_C)],'k','LineStyle','none','LineWidth',3);
saveas(gcf,'MajMinC')
print( gcf, '-painters', 'MajMinC.ai', '-dpsc','-r300')

close all

%%
bar(1:2,[mean(Maj_R);mean(Min_R)])
ylim([0,25])

SetFigureDefaults(gca);
ylabel('Length (um)')
title('Major and Minor Axis lengths - FAK')
xticklabels({'Major Axis','Minor Axis'});
hold on
errorbar([mean(Maj_R);mean(Min_R)],[std(Maj_R),std(Min_R)],'k','LineStyle','none','LineWidth',3);
saveas(gcf,'MajMinR')
print( gcf, '-painters', 'MajMinR.ai', '-dpsc','-r300')

close all

%%
bar(1:2,[mean(Ht_R);mean(Wd_R)])
ylim([0,25])
SetFigureDefaults(gca);
ylabel('Length (um)')
title('Height and widht of cells - FAK')
xticklabels({'Height','Width'});
hold on
errorbar([mean(Ht_R);mean(Wd_R)],[std(Ht_R),std(Wd_R)],'k','LineStyle','none','LineWidth',3);
saveas(gcf,'HtWdR')
print( gcf, '-painters', 'HtWdR.ai', '-dpsc','-r300')
close all
%%
bar(1:2,[mean(Ht_C);mean(Wd_C)])
ylim([0,25])
SetFigureDefaults(gca);
ylabel('Length (um)')
title('Height and widht of cells - CTRL')
xticklabels({'Height','Width'});
hold on
errorbar([mean(Ht_C);mean(Wd_C)],[std(Ht_C),std(Wd_C)],'k','LineStyle','none','LineWidth',3);
saveas(gcf,'HtWdC')
print( gcf, '-painters', 'HtWdC.ai', '-dpsc','-r300');
close all
%%
OrC = cat(1,MetaDataC.OriAng)';
OrR = cat(1,MetaDataR.OriAng)';

OrC = [OrC-180, OrC, OrC+180];
OrR = [OrR-180, OrR, OrR+180];

polarhistogram(deg2rad(OrC),deg2rad(0:15:180),'Normalization','Probability');
h = gca;
h.ThetaLim = [0 180];
h.RLim = [0 0.18];
set(gca,'FontSize',16)
set(gca,'FontWeight','bold')
set(gca,'LineWidth',1)
title('Orientation Angle of Cells - CTRL')
saveas(gcf,'OrC')
print( gcf, '-painters', 'OrC.ai', '-dpsc','-r300');

polarhistogram(deg2rad(OrR),deg2rad(0:15:180),'Normalization','Probability');
h = gca;
h.ThetaLim = [0 180];
h.RLim = [0 0.18];
set(gca,'FontSize',16)
set(gca,'FontWeight','bold')
set(gca,'LineWidth',1)
title('Orientation Angle of Cells - FAK')
saveas(gcf,'OrR')
print( gcf, '-painters', 'OrR.ai', '-dpsc','-r300');


hold off

end