function [] = Fig5e_boxplots_weighted(MetaC,MetaR,MetaF)



color_c = 'b';
color_b = 'r';
color_f = 'k';

[aspinfoR] = Fig5c(MetaR);
[aspinfoC] = Fig5c(MetaC);
[aspinfoF] = Fig5c(MetaF);


figure

wtdboxplot(aspinfoC.circasp',aspinfoC.circwts',color_c,0.75);
hold on
wtdboxplot(aspinfoC.axasp',aspinfoC.axwts',color_c,1.75);
hold on
wtdboxplot(aspinfoR.circasp',aspinfoR.circwts',color_b,1);
hold on
wtdboxplot(aspinfoR.axasp',aspinfoR.axwts',color_b,2);
hold on
wtdboxplot(aspinfoF.circasp',aspinfoF.circwts',color_f,1.25);
hold on
wtdboxplot(aspinfoF.axasp',aspinfoF.axwts',color_f,2.25);

% legend({'CTRL','BAPN','FAK'})
SetFigureDefaults(gca)
title('Axial Fractions Vs Aspect Ratio')
xticklabels({'Circumferential','Axial'});
ylabel('Aspect Ratio')
grid on
xlim([0.5 2.5])
xticks(1:2);
% ylim([0 0.9]);

set(findobj(gca,'type','line'),'linew',2)

   saveas(gcf,'Fig5eBox_wt.fig');
   print( gcf, '-painters', ['Fig5eBox_wt.ai'], '-dpsc','-r300')

end

function [] = wtdboxplot(inputvec,wtvec,inputcolor,posnvec)

import iosr.*

a = @iosr.statistics.boxPlot;
a(posnvec,inputvec,'weights',wtvec,'boxColor',inputcolor,'boxWidth',0.1,...
    'ShowScatter',true,'notch',false,'ScatterMarker','o');


end



    function [meanpropmaj,meanpropmin,wts] =wtdprop(Metadata)
wts = arrayfun(@(x) length([x.OriAng]),Metadata);
wts = wts./sum(wts);

meanpropmaj = arrayfun(@(x) nanmean(x.MajAxlen),Metadata);
meanpropmin = arrayfun(@(x) nanmean(x.MinAxlen),Metadata);

    end





