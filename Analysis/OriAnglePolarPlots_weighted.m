function [] = OriAnglePolarPlots_weighted(Metadata,movtype,movcolor)

%% INPUTS: 
% Metadata: Output from MetaAnalysis_CorrectedData_movlist function
% movtype: string describing the movie type; e.g.: 'CTRL'
% movcolor: color of the graphs. e.g.: 'r' 

wts = arrayfun(@(x) length([x.OriAng]),Metadata);
wts = wts./sum(wts);

figure
for i = 1:length(wts)
    
h = histcounts([Metadata(i).OriAng],(-15:30:195),'Normalization','Probability');
h(1) = h(1) + h(end);
h(end) = [];
hmat(i,:) = h;
clear h;
end

meanvec = sum(hmat.*wts');
errvec = std(hmat,wts');
% errvec = errvec./sqrt(length(wts));
meanvec = [meanvec,meanvec,meanvec(1)];
errvec = [errvec,errvec,errvec(1)];
angvec = deg2rad(0:30:360);
p=polarplot(angvec,meanvec,'-','Color',movcolor,'lineWidth',2);
title(movtype);
set(gca,'RLim',[0 0.65]);
hold on

for i = 1:length(meanvec)
    
    polarplot([angvec(i)*ones(1,3)],[meanvec(i)-errvec(i),meanvec(i),meanvec(i)+errvec(i)],....
        'k-')
    
end

SetFigureDefaults_polar(gca);
print( gcf, '-painters', ['OriAng_std',movtype,'.ai'], '-dpsc','-r300');
end