
function [F] = PhInt_weightedNumbers(FiberInfo,MetadataC,MetadataR,MetadataF)

    intf = FiberInfo.phfak;
    intfcad = FiberInfo.efak;
    
    intc = FiberInfo.phctrl;
    intccad = FiberInfo.ectrl;
    
    intb = FiberInfo.phbapn;
    intbcad = FiberInfo.ebapn;
    
    

%% calculate the wts and number of movies for each case

[wtc,wt1c,wt2c] = findweightvecs(MetadataC,intc,intccad);
[wtb,wt1b,wt2b] = findweightvecs(MetadataR,intb,intbcad);
[wtf,wt1f,wt2f] = findweightvecs(MetadataF,intf,intfcad);


n1c = length(wt1c); n2c = length(wt2c);
n1b = length(wt1b); n2b = length(wt2b);
n1f = length(wt1f); n2f = length(wt2f);



%% weighted mean and std calculation

[F.meanintc,F.stdintc,F.stderrc] = findmeanAndstd(intc,wt1c,n1c);
[F.meanintb,F.stdintb,F.stderrb] = findmeanAndstd(intb,wt1b,n1b);
[F.meanintf,F.stdintf,F.stderrf] = findmeanAndstd(intf,wt1f,n1f);

[F.meanintccad,F.stdintccad,F.stderrccad] = findmeanAndstd(intccad,wt2c,n2c);
[F.meanintbcad,F.stdintbcad,F.stderrbcad] = findmeanAndstd(intbcad,wt2b,n2b);
[F.meanintfcad,F.stdintfcad,F.stderrfcad] = findmeanAndstd(intfcad,wt2f,n2f);



% F.pcc=ttest2OnValues([F.meanintc(1),F.stdintc(1),n1c],[F.meanintc(4),F.stdintc(4),n1c]);
% F.pbb=ttest2OnValues([F.meanintb(1),F.stdintb(1),n1b],[F.meanintb(4),F.stdintb(4),n1b]);
% F.pff=ttest2OnValues([F.meanintf(1),F.stdintf(1),n1f],[F.meanintf(4),F.stdintf(4),n1f]);

F.pcb=....
    ttest2OnValues([F.meanintc(4),F.stdintc(4),n1c],[F.meanintb(4),F.stdintb(4),n1b]);
F.pbf=ttest2OnValues([F.meanintb(4),F.stdintb(4),n1b],[F.meanintf(4),F.stdintf(4),n1f]);
F.pcf=ttest2OnValues([F.meanintc(4),F.stdintc(4),n1c],[F.meanintf(4),F.stdintf(4),n1f]);
% 
cint = cat(1,intc.MeanInt);
bint = cat(1,intb.MeanInt);
fint = cat(1,intf.MeanInt);

F.ratioc = cint(:,1)./cint(:,4);
F.ratiob = bint(:,1)./bint(:,4);
F.ratiof = fint(:,1)./fint(:,4);

F.meanprefC = nansum(F.ratioc.*wt1c');
F.meanprefR = nansum(F.ratiob.*wt1b');
F.meanprefF = nansum(F.ratiof.*wt1f');

F.meanprefstdC = std(F.ratioc,wt1c');
F.meanprefstdR = std(F.ratiob,wt1b');
F.meanprefstdF = std(F.ratiof,wt1f');


[~,F.pcc]=ttestnew(F.ratioc,1);
[~,F.pbb]=ttestnew(F.ratiob,1);
[~,F.pff]=ttestnew(F.ratiof,1);

% 
% [~,F.pcc] = ttest2OnValues([F.meanintccad,
% [~,F.pbb] =
% [~,F.pff] =





save FiberInfoNumbers;


end

