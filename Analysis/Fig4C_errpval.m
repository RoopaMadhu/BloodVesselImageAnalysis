
%Fig 4C calculations of orientation angles for each movie

for i = 1:length(MetaDataC)
OrChist(i,:) = histcounts([MetaDataC(i).OriAng-180,MetaDataC(i).OriAng,MetaDataC(i).OriAng+180],'BinLimits',[0,179.9],'BinEdges',0:15:180,'Normalization','probability');
end


for i = 1:length(MetaDataR)
OrRhist(i,:) = histcounts([MetaDataR(i).OriAng-180,MetaDataR(i).OriAng,MetaDataR(i).OriAng+180],'BinLimits',[0,179.9],'BinEdges',0:15:180,'Normalization','probability');
end

R_axial = mean(OrRhist(:,6)+OrRhist(:,7))
std(OrRhist(:,6)+OrRhist(:,7))
R_circum = mean(OrRhist(:,1)+OrRhist(:,end))
std(OrRhist(:,1)+OrRhist(:,end))