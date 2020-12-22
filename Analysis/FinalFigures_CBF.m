

function [] = FinalFigures_CBF(MetaC,MetaR,MetaF)
%%
figure

subplot(3,1,1)
histogram(cat(1,MetaC.CellArea2D),0:20:500,'FaceColor','b')
ylim([0 400])
ylabel('# cells')
SetFigureDefaults(gca)
title('CTRL')

subplot(3,1,2)
histogram(cat(1,MetaR.CellArea2D),0:20:500,'FaceColor','r')
ylim([0 400])
ylabel('# cells')
SetFigureDefaults(gca)
title('BAPN')

subplot(3,1,3)
histogram(cat(1,MetaF.CellArea2D),0:20:500,'FaceColor','k')
ylim([0 400])
ylabel('# cells')
SetFigureDefaults(gca)
title('FAK')
xlabel('CellArea(um^2)')

% % saveas(gcf,'CellArea');
 print( gcf, '-painters', 'Area.ai', '-dpsc','-r300');
% print( gcf, '-painters', 'FinalFigzoo.ps', '-dpsc','-r300','-append');
% 
% 
% %% 
% 
% figure
% 
subplot(3,1,1)
histogram(cat(1,MetaC.Dia),0:5:70,'FaceColor','b')
xlim([15 65])
ylim([0 10])
ylabel('# vessels')
SetFigureDefaults(gca)
title('CTRL')

subplot(3,1,2)
histogram(cat(1,MetaR.Dia),0:5:70,'FaceColor','r')
xlim([15 65])
ylim([0 10])

ylabel('# vessels')
SetFigureDefaults(gca)
title('BAPN')

subplot(3,1,3)
histogram(cat(1,MetaF.Dia),0:5:70,'FaceColor','k')
xlim([15 65])
ylim([0 10])

ylabel('# vessels')
SetFigureDefaults(gca)
title('FAK')
xlabel('Vessel Diameter(um)')

% saveas(gcf,'VesselDiameter');
% print( gcf, '-painters', 'OrC.pdf', '-dpsc','-r300');
print( gcf, '-painters', 'FinalFigzoo.ps', '-dpsc','-r300','-append');

% 
% %% 
figure

plot(mean(cat(1,MetaC.NeighNum)),'ro-','LineWidth',2)
hold on
plot(mean(cat(1,MetaR.NeighNum)),'b*-','LineWidth',2)
plot(mean(cat(1,MetaF.NeighNum)),'k.-','LineWidth',2)
ylim([0 0.8])

SetFigureDefaults(gca)
grid on
legend({'CTRL','BAPN','FAK'})
title('Neighbor number distribution')
xlabel('# neighbors')

ylabel('Probability')

% saveas(gcf,'Neighnum');
print( gcf, '-painters', 'FinalFigzoo.ps', '-dpsc','-r300','-append');
%%

figure

htc = (cat(1,MetaC.CellHt));
wdc = (cat(1,MetaC.CellWd));

htr = (cat(1,MetaR.CellHt));
wdr = (cat(1,MetaR.CellWd));

htf = (cat(1,MetaF.CellHt));
wdf = (cat(1,MetaF.CellWd));

nc = length(htc); nr = length(htr); nf = length(htf);

htec = nanstd(htc);
wdec = nanstd(wdc);

hter = nanstd(htr);
wder = nanstd(wdr);

htef = nanstd(htf);
wdef = nanstd(wdf);

htc = nanmean(htc); wdc = nanmean(wdc); htr = nanmean(htr); wdr = nanmean(wdr);
htf = nanmean(htf); wdf = nanmean(wdf);

b = bar([htc,htr,htf; wdc, wdr, wdf]);
b(1).FaceColor = 'b';
b(2).FaceColor = 'r';
b(3).FaceColor = 'k';
hold on
errorbar([0.78,1,1.22;1.78,2,2.22],[htc,htr,htf; wdc, wdr, wdf],[htec,hter,htef; wdec, wder, wdef],'k-','LineStyle','none','LineWidth',3)
ylim([0 25])

xticklabels({'Height','Width'})
ylabel('Length(um)')
SetFigureDefaults(gca)
legend({'CTRL','BAPN','FAK'})
title('HtvsWd')
grid on

% saveas(gcf,'HtvsWd');
% print( gcf, '-painters', 'HtvsWd.pdf', '-dpsc','-r300');

print( gcf, '-painters', 'FinalFigzoo.ps', '-dpsc','-r300','-append');
%%

figure

majc = mean(cat(1,MetaC.MajAxlen));
minc = mean(cat(1,MetaC.MinAxlen));

majr = mean(cat(1,MetaR.MajAxlen));
minr = mean(cat(1,MetaR.MinAxlen));

majf = mean(cat(1,MetaF.MajAxlen));
minf = mean(cat(1,MetaF.MinAxlen));

majec = nanstd(cat(1,MetaC.MajAxlen));
minec = nanstd(cat(1,MetaC.MinAxlen));

majer = nanstd(cat(1,MetaR.MajAxlen));
miner = nanstd(cat(1,MetaR.MinAxlen));

majef = nanstd(cat(1,MetaF.MajAxlen));
minef = nanstd(cat(1,MetaF.MinAxlen));


b = bar([majc,majr,majf; minc, minr, minf]);
b(1).FaceColor = 'b';
b(2).FaceColor = 'r';
b(3).FaceColor = 'k';
hold on
errorbar([0.78,1,1.22;1.78,2,2.22],[majc,majr,majf; minc, minr, minf],[majec,majer,majef; minec, miner, minef],'k-','LineStyle','none','LineWidth',3)
grid on
ylim([0 25])


xticklabels({'MajAxLen','MinAxLen'})
ylabel('Length(um)')
SetFigureDefaults(gca)
legend({'CTRL','BAPN','FAK'})
title('MajvsMin')
% saveas(gcf,'MajvsMin');
% print( gcf, '-painters', 'MajvsMin.pdf', '-dpsc','-r300');
print( gcf, '-painters', 'FinalFigzoo.ps', '-dpsc','-r300','-append');

%%
figure

OrC = cat(1,MetaC.OriAng)';
OrR = cat(1,MetaR.OriAng)';
OrF = cat(1,MetaF.OriAng)';

% 
% OrC = [OrC-180, OrC, OrC+180];
% OrR = [OrR-180, OrR, OrR+180];
% OrF = [OrF-180, OrF, OrF+180];
% 
subplot(1,3,1)

polarhistogram(deg2rad(OrC),deg2rad(0:15:180),'Normalization','Probability');
h = gca;
h.ThetaLim = [0 180];
h.RLim = [0 0.25];
set(gca,'FontSize',16)
set(gca,'FontWeight','bold')
set(gca,'LineWidth',1)
title('CTRL')
set(gca,'Position',[0.05,0.05,0.25,0.8])
% saveas(gcf,'OrC')
% print( gcf, '-painters', 'OrC.pdf', '-dpsc','-r300');
subplot(1,3,2)

polarhistogram(deg2rad(OrR),deg2rad(0:15:180),'Normalization','Probability');
h = gca;
h.ThetaLim = [0 180];
h.RLim = [0 0.25];
set(gca,'FontSize',16)
set(gca,'FontWeight','bold')
set(gca,'LineWidth',1)
title('BAPN')
set(gca,'Position',[0.37,0.05,0.25,0.8])
% saveas(gcf,'OrR')
% print( gcf, '-painters', 'OrR.pdf', '-dpsc','-r300');

subplot(1,3,3)

polarhistogram(deg2rad(OrF),deg2rad(0:15:180),'Normalization','Probability');
h = gca;
h.ThetaLim = [0 180];
h.RLim = [0 0.25];
set(gca,'FontSize',16)
set(gca,'FontWeight','bold')
set(gca,'LineWidth',1)
title('FAK')
set(gca,'Position',[0.7,0.05,0.25,0.8])
% saveas(gcf,'Or')
% print( gcf, '-painters', 'Or.pdf', '-dpsc','-r300');
print( gcf, '-painters', 'FinalFigzoo.ps', '-dpsc','-r300','-append');

figure
PC = Fig5a_V3(MetaC,'CTRL');
print( gcf, '-painters', 'FinalFigzoo.ps', '-dpsc','-r300','-append');

figure
PR = Fig5a_V3(MetaR,'BAPN');
print( gcf, '-painters', 'FinalFigzoo.ps', '-dpsc','-r300','-append');

figure
PF = Fig5a_V3(MetaF,'FAK');
print( gcf, '-painters', 'FinalFigzoo.ps', '-dpsc','-r300','-append');

end








