function [] = Fig3(Metadata1,Metadata2,movtype1,movtype2)

figure
%%
histogram(cat(1,Metadata1.CellHt),0:2:40);
SetFigureDefaults(gca);
xlabel('Cell Height(um)');
ylabel('#Cells');
% ylim([0 600]);
% xticks(0:5:40);

saveas(gcf,['CellHt_',movtype1]);
print( gcf, '-painters', ['CellHt_',movtype1,'.ai'], '-dpsc','-r300')
clf(gcf);

histogram(cat(1,Metadata2.CellHt),0:2:40);

SetFigureDefaults(gca);
xlabel('Cell Height(um)');
ylabel('#Cells');
% ylim([0 600]);
xticks(0:5:40);

saveas(gcf,['CellHt_',movtype2]);
print( gcf, '-painters', ['CellHt_',movtype2,'.ai'], '-dpsc','-r300')
clf(gcf);
%%
histogram(cat(1,Metadata1.CellWd),0:2:40);
SetFigureDefaults(gca);
xlabel('Cell Width(um)');
ylabel('#Cells');
% ylim([0 600]);
% xticks(0:5:40);

saveas(gcf,['CellWd_',movtype1]);
print( gcf, '-painters', ['CellWd_',movtype1,'.ai'], '-dpsc','-r300')
clf(gcf);


histogram(cat(1,Metadata2.CellWd),0:2:40);

SetFigureDefaults(gca);
xlabel('Cell Width(um)');
ylabel('#Cells');
% ylim([0 600]);
% xticks(0:5:40);

saveas(gcf,['CellWd_',movtype2]);
print( gcf, '-painters', ['CellWd_',movtype2,'.ai'], '-dpsc','-r300')
clf(gcf);
%%
histogram(cat(1,Metadata1.CellArea2D),0:50:500);
SetFigureDefaults(gca);
xlabel('Cell Area(um^2)');
ylabel('#Cells');

saveas(gcf,['CellArea_',movtype1]);
print( gcf, '-painters', ['CellArea_',movtype1,'.ai'], '-dpsc','-r300')

clf(gcf);

histogram(cat(1,Metadata2.CellArea2D),0:50:500);

SetFigureDefaults(gca);
xlabel('Cell Area(um^2)');
ylabel('#Cells');

saveas(gcf,['CellArea_',movtype2]);
print( gcf, '-painters', ['CellArea_',movtype2,'.ai'], '-dpsc','-r300')

clf(gcf);
%%
end




