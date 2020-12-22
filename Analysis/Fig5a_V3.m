
function [PrefOr] = Fig5a_V3(Metadata,movtype)
figure
Dia = [Metadata.Dia];
interprange = [min(Dia)+1:5:max(Dia)-1];

%% Fig5a using interp1 instead of interpolatePoints2Vec
PrefOr = PrefrentialOri_V2(Metadata);

preid = all(~isnan(PrefOr),2);
PrefOr = PrefOr(preid,:);
% sig = ((interprange(2) - interprange(1)) /2) ;
plot(Dia,PrefOr(:,1),'go','LineWidth',2)
hold on
plot(Dia,PrefOr(:,2),'ro','LineWidth',2)
plot(Dia,PrefOr(:,3),'ko','LineWidth',2)

%moving average
% windowSize = 2;
% b = (1/windowSize)*ones(1,windowSize);
% a = 1;

% interp_ax = interp1(Dia,PrefOr(:,1),interprange,'linear');
% plot(interprange,sgolayfilt(interp_ax,3,5),'g.-','LineWidth',2);

% plot(interprange,filter(b,a,interp_ax),'g.-','LineWidth',2);
% plot(interprange,interp_ax,'g.-','LineWidth',2);


% interp_circum = interp1(Dia,PrefOr(:,2),interprange,'linear');
% plot(interprange,sgolayfilt(interp_circum,3,5),'r.-','LineWidth',2);
% plot(interprange,filter(b,a,interp_circum),'r.-','LineWidth',2);
% plot(interprange,interp_circum,'r.-','LineWidth',2);

% interp_other = interp1(Dia,PrefOr(:,3),interprange,'linear');
% plot(interprange,sgolayfilt(interp_other,3,5),'k.-','LineWidth',2);
% plot(interprange,filter(b,a,interp_other),'k.-','LineWidth',2);
% plot(interprange,interp_other,'k.-','LineWidth',2);

SetFigureDefaults(gca);
grid on
ylim([0 0.83]);
xlim([19,80]);

legend({'Circum','Axial','Other','Circum','Axial','Other'})

xlabel('Diameter');
ylabel('Fraction of cells')

title(['Fractions Of Cells vs Orientaions - ',movtype])

print( gcf, '-painters', ['Fig5a',movtype,'new.ai'], '-dpsc','-r300')
% saveas(gcf,'Fig5aR')
% 
end %of main function




