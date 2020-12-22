function [aspinfo] = Fig5c(Metadata,color1,movtype)

PrefOr = PrefrentialOri_V2(Metadata);

prefax = PrefOr(:,2);
meanasp = arrayfun(@(x) (nanmean([x.MajAxlen]./[x.MinAxlen])),Metadata);
stdasp = arrayfun(@(x) (nanstd([x.MajAxlen]./[x.MinAxlen])),Metadata);

wts = arrayfun(@(x) length([x.OriAng]),Metadata);
wts_circ = wts(prefax<0.3);
ncirc = length(wts_circ);
wts_circ = wts_circ./sum(wts_circ);
wts_ax = wts(prefax>=0.3);
wts_ax = wts_ax./sum(wts_ax);
nax = length(wts_ax);

meancircasp = sum(meanasp(prefax<0.3).*wts_circ);
meanaxasp = sum(meanasp(prefax>=0.3).*wts_ax);
    
stdcircasp = std(meanasp(prefax<0.3),wts_circ);
stdaxasp = std(meanasp(prefax>=0.3),wts_ax);

aspinfo.circasp = meanasp(prefax<0.3);
aspinfo.circwts = wts_circ;

aspinfo.axasp = meanasp(prefax>=0.3);
aspinfo.axwts = wts_ax;


pval1=ttest2OnValues([meancircasp,stdcircasp,ncirc],[meanaxasp,stdaxasp,nax])
% plot(PrefOr(:,2),arrayfun(@(x) (nanmean([x.MajAxlen]./[x.MinAxlen])),Metadata),'o','LineWidth',2,'Color',color1);
% 
% SetFigureDefaults(gca);
% grid on
% 
% 
% xlabel('%Axially Oriented Cells');
% ylabel('Mean Aspect Ratio')
% 
% title(['Axial Orientation vs Aspect Ratio - ',movtype])

% print( gcf, '-painters', ['Fig5c.ai'], '-dpsc','-r300')
% saveas(gcf,'Fig5c')

end