function [MetaPhCell,caratio] = ExtractPhIntQuantCell(movlist,Metadata,cinput)
% figure
for i = 1:length(movlist)
    
    MetaPhCell(i) = load([movlist{i},'/ExploratoryAnalysis/PhQuantcell_5.mat'],'PhQuantCell');
    circint = cat(1,MetaPhCell(i).PhQuantCell(:,1));
    axint   = cat(1,MetaPhCell(i).PhQuantCell(:,4));
    
    caratio(i) = nanmean(axint./circint);
    clear circint axint;
    
end



%% plot graphs
Dia = cat(1,Metadata.Dia);

plot(Dia,caratio,'o','Color',cinput,'MarkerSize',10,'LineWidth',2)
grid on
SetFigureDefaults(gca)
xlabel('Diameter')
ylabel('Axial to Circ Ph Ratio');
title('Mean axial to circ ph ratio');
hold on
xlim([19,80]);

end