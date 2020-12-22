function [] = lengthmap_2D_V2(movielist)

wd = cd;

for i = 1:length(movielist)
    
    cd(movielist{i});
    
    
    md =cd;
    
  
    
    % loading dataimglist by moving to the ph channel as it contains ax
    % seg. section information - tstart and tstop of seg.
    
   phdir = findthePhChannel(md);
   
   cd(phdir);
   
   load dataimglist;
   
    tstart = data(1).SegTimeInterval(1);
    tend = data(1).SegTimeInterval(2);
    
    
    load ('Data2Analyse.mat','cap','angbinsize','Data2Analyse');
    
    angbins = -180:angbinsize:180;
    
    for t = tstart:tend
        
        x = Data2Analyse(t).smoothedxy(:,1);
        y = Data2Analyse(t).smoothedxy(:,2);
        [angtmp,~,xsort,ysort] = AngleSorting(x,y);
        
        xinterp=interp1(angtmp,xsort,angbins,'linear','extrap');
        yinterp=interp1(angtmp,ysort,angbins,'linear','extrap');
        
        tmp_Distbwpoints = sqrt(diff([xinterp,xinterp(1)]).^2+diff([yinterp,yinterp(1)]).^2);
        tmp_DistFromOrigin = cumsum(tmp_Distbwpoints);
        CircInfo(t).cirum_ori = tmp_DistFromOrigin(end);
        
        len_of_arc(t,:) = tmp_Distbwpoints;
        
        angbinsize_finer = 0.01;
        angbins_finer = -180:angbinsize_finer:180;
        
        xinterp_finer=interp1(angtmp,xsort,angbins_finer,'linear','extrap');
        yinterp_finer=interp1(angtmp,ysort,angbins_finer,'linear','extrap');
        
        tmp_Distbwpoints = sqrt(diff([xinterp_finer,xinterp_finer(1)]).^2+diff([yinterp_finer,yinterp_finer(1)]).^2);
        len_of_arc_finer = tmp_Distbwpoints;
        tmp_DistFromOrigin = cumsum(tmp_Distbwpoints);
        CircInfo(t).cirum_finer = tmp_DistFromOrigin(end);
        
%     angdiff_endsfiner = angbins_finer(end)-abs(angbins_finer(1));
%         if angdiff_endsfiner == 0
%                     angdiffvec_finer = cumsum([diff(angbins_finer)]);
%                            tmp_Distbwpoints = sqrt(diff([xinterp]).^2+diff([yinterp]).^2);
% 
%         else
% 
%             
%         angdiffvec_finer = cumsum([diff(angbins_finer),angbins_finer(end)-abs(angbins_finer(1))]);
%         end
%         angdiffvec = cumsum([diff(angbins),angbins(end)-abs(angbins(1))]);

        % corrected length map
        interpdist=interp1(cumsum(diff([angbins_finer,angbins_finer(1)])),tmp_DistFromOrigin,cumsum(diff([angbins,angbins(1)])),'linear','extrap');
%         interpdist=interp1(angdiffvec_finer,tmp_DistFromOrigin,angdiffvec,'linear','extrap');

        len_of_arc_corrected(t,:) =  [interpdist(1),diff(interpdist)];
        CircInfo(t).arc_ori_circ = [interpdist(end),tmp_DistFromOrigin(end)];
        
        
        
    end
    
    lenmap2D_RawCapped = [len_of_arc(tstart:tend,end-cap:end),len_of_arc(tstart:tend,1:end),len_of_arc(tstart:tend,1:cap+1)];
%   lenmap2D_tcropped = len_of_arc(tstart:tend,:);
    lenmap2D_CorrectedCapped =  [len_of_arc_corrected(tstart:tend,end-cap:end),len_of_arc_corrected(tstart:tend,1:end),len_of_arc_corrected(tstart:tend,1:cap+1)];
    lenmap2D_CorrectedtCrop =  len_of_arc_corrected(tstart:tend,:);

    
    cd(md);
    
    save ExploratoryAnalysis/RadiusMap lenmap* len_of_arc* CircInfo;
    
    clear x* y* angtmp angbinsize_finer angbins_finer tmp* CircInfo len*
end%of for loop over movies

end %of main function
