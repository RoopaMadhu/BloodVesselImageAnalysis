
function [] = PhIntPlots_weighted(FiberInfo,movtype,Metadata,id1,id2)

if strcmp(movtype,'FAK') == 1
    int = FiberInfo.phfak;
    int2 = FiberInfo.efak;
    
elseif strcmp(movtype,'CTRL') == 1
    int = FiberInfo.phctrl;
    int2 = FiberInfo.ectrl;
    
elseif strcmp(movtype,'BAPN') == 1
    int = FiberInfo.phbapn;
    int2 = FiberInfo.ebapn;
    
    
end
%% calculate the wts and number of movies for each case

wt = arrayfun(@(x) length([x.CellArea2D]),Metadata);

% id1 = find(arrayfun(@(x) ~isempty(x.MeanInt),int));
% id2 = find(arrayfun(@(x) ~isempty(x.MeanInt),int2));
wt1 = wt(id1); wt2 = wt(id2);
wt1 = wt1./sum(wt1); %normalizing the wts.
wt2 = wt2./sum(wt2); %normalizing the wts.

n1 = length(id1); n2 = length(id2);


%% unweighted mean
% int1 = cat(1,int.MeanInt);
% meanint = nanmean(int1);
% meanvec1 = [meanint,meanint(1)];
% meanvec1 = [meanvec1,meanvec1(2:end)];
%
% stdint1 = nanstd(int1);
% stdint1 = [stdint1,stdint1(1)];
% stdint1 = [stdint1,stdint1(2:end)];
% errvec = stdint1./sqrt(n1);
% angvec = deg2rad(0:30:360);
%
% int2 = cat(1,int2.MeanInt);
%
% meanint2 = nanmean(int2);
% meanvec2 = [meanint2,meanint2(1)];
% meanvec2 = [meanvec2,meanvec2(2:end)];
%
% stdint2 = nanstd(int2);
% stdint2 = [stdint2,stdint2(1)];
% stdint2 = [stdint2,stdint2(2:end)];
% errvec2 = stdint2./sqrt(n2);

%% weighted mean and std calculation


int1 = cat(1,int.MeanInt);
int1 = int1./mean(int1,2);
int1_wt = int1.*wt1';

meanint = nansum(int1_wt);
meanvec1 = [meanint,meanint(1)];
meanvec1 = [meanvec1,meanvec1(2:end)];

stdint1 = sqrt(var(int1,wt1'));
stdint1 = [stdint1,stdint1(1)];
stdint1 = [stdint1,stdint1(2:end)];
% errvec = stdint1./sqrt(n1);
errvec = stdint1;
angvec = deg2rad(0:15:360);

int2 = cat(1,int2.MeanInt);
int2 = int2./mean(int2,2);%normalizing
int2_wt = int2.*wt2';
meanint2 = nansum(int2_wt);
meanvec2 = [meanint2,meanint2(1)];
meanvec2 = [meanvec2,meanvec2(2:end)];

stdint2 = sqrt(var(int2,wt2'));
stdint2 = [stdint2,stdint2(1)];
stdint2 = [stdint2,stdint2(2:end)];
% errvec2 = stdint2./sqrt(n2);
errvec2 = stdint2;

figure
h=polarplot(angvec,meanvec1);
h.LineWidth = 5;
% h.ThetaLim = [0,180];
hold on

for i = 1:length(meanvec1)
    
    polarplot([angvec(i)*ones(1,3)],[meanvec1(i)-errvec(i),meanvec1(i),meanvec1(i)+errvec(i)],....
        'k-','lineWidth',2)
    
       polarplot([angvec(i)*ones(1,3)],[meanvec2(i)-errvec2(i),meanvec2(i),meanvec2(i)+errvec2(i)],....
        'r-','linewidth',2)
    
end

y = gca;
y.RLim = [0,1.5];


SetFigureDefaults_polar(gca);
h=polarplot(angvec,meanvec2,'k');
h.LineWidth = 3;


% h1=polarplot(deg2rad(0:5:360),meanint-stdint);
% h1.LineStyle = '--';
% h1.LineWidth = 3;
%
% h2=polarplot(deg2rad(0:5:360),meanint+stdint);
% h2.LineStyle = '--';
% h2.LineWidth = 3;

% y = gca;
% y.RLim = [0,0.12];
% y.FontWeight = 'bold';

title(movtype);

% saveas(gcf,sprintf('AspHeatmap_%s_%s',movtype,subtypenew));
% saveas(gcf,sprintf('AspHeatmap_%s_%s.tif',movtype,subtypenew));
print( gcf, '-painters', sprintf('PhPolarplotsNew_%s.ai',movtype), '-dpsc','-r300');



end

