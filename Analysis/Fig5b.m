function [] = Fig5b(Dia,Metaph,color1,movtype)

plot(Dia,[Metaph.avgcaratio],'o','Color',color1,'LineWidth',2)


SetFigureDefaults(gca);
grid on
ylim([0.5 2]);

legend({'CONTROL'})

xlabel('Vessel Diameter');
ylabel('Circumferential to axial phalloidin ratio')

title(['Diameter vs phalloidin ratio (circum to axial) - ',movtype])

end %of main function
