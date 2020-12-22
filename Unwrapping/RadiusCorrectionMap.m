function [] = RadiusCorrectionMap(movielist)

wd = cd;

for i = 1:length(movielist)
    
    cd(movielist{i});
    
    
    md =cd;
    
    
    % loading dataimglist by moving to the ph channel as it contains ax
    % seg. section information - tstart and tstop of seg.
    
    phdirpath = findthePhChannel(md);
    
    if ~isempty(phdirpath)
        
        cd(phdirpath);
        
        load dataimglist;
        
        load ('Data2Analyse.mat','cap','angbinsize','Data2Analyse');
        
        
        
        
        angbins = -180:angbinsize:180;
        
        tstart = data(1).SegTimeInterval(1);
        tstop = data(1).SegTimeInterval(2);
        
        for t = tstart:tstop
            
            x = Data2Analyse(t).smoothedxy(:,1);
            y = Data2Analyse(t).smoothedxy(:,2);
            
            center=[mean([min(x),max(x)]),mean([min(y),max(y)])];
            
            [angtmp,~,xsorted,ysorted] = AngleSorting(x,y);
            
            vx = mean([xsorted(1),xsorted(end)]);
            vy = mean([ysorted(1),ysorted(end)]);
            
            
            xinterp=interp1(angtmp,xsorted,angbins,'linear',vx);
            yinterp=interp1(angtmp,ysorted,angbins,'linear',vy);
            
            RadiusMap(t,:) = (sqrt(dist2(center,[xinterp',yinterp'])));
            
            
        end
        
        
        cd(md);
        
        save ExploratoryAnalysis/RadiusMap RadiusMap -append;
        
        clear RadiusMap;
        
    else
        
        sprintf('No phalloidin channel available');
        
    end %of if loop
    
    
end

end %of main function




