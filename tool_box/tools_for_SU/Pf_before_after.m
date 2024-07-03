function [peak_bef com_bef peak_after com_after]   = Pf_before_after(pos, spks,Shocks_filt, min_lap, in_lap, sigma,Xedges)
    
    in_lap = in_lapA;
    min_lap = 5; 
    pos = pos_ave;
    spks=spks_ave; 
    
    shock = Shocks_filt(1); 
    
    %Get first and last shock lap
    %Most of the times the shock is given in a non-valid lap (removed form
    %inLap), so look for the shock nearest lap 
    %Fist shock
    fshock_lap=[];
    for indx = 1: size(inLap,1) 
        lap_ini = inLap(indx,1);
        lap_end = inLap(indx,2);
        diff_ini= abs(shock- lap_ini);
        diff_end= abs(shock- lap_end);
        fshock_lap=[fshock_lap;diff_ini,diff_end,indx];
    end 
    [~, iindx]=min(fshock_lap(:,1)); 
    [~, eindx]=min(fshock_lap(:,2)); 
    if iindx>eindx
        fshock_lap=eindx; 
    elseif iindx<eindx
        fshock_lap=iindx; 
    elseif iindx==eindx
        fshock_lap=iindx;
    else
        disp('No shock inLap')
    end
    
    %Last shock 
    shock_last = Shocks_filt(end);
    slast_lap=[];
    for indx = 1: size(inLap,1) 
        lap_ini = inLap(indx,1);
        lap_end = inLap(indx,2);
        diff_ini= abs(shock_last- lap_ini);
        diff_end= abs(shock_last- lap_end);
        slast_lap=[slast_lap;diff_ini,diff_end,indx];
    end 
    [~, iindx]=min(slast_lap(:,1)); 
    [~, eindx]=min(slast_lap(:,2)); 
    if iindx>eindx
        slast_lap=eindx; 
    elseif iindx<eindx
        slast_lap=iindx; 
    elseif iindx==eindx
        slast_lap=iindx;
    end
    
    
    %Check n laps before the shock 
    if shock_lap > min_lap
        %Compute firing map before
        spks_b = Restrict(spks,in_lap(1:shock_lap-1,:)); % restrict spikes to before shock
        pos_b = Restrict(pos,in_lap(1:shock_lap-1,:));
        [map_bef , stats_bef] = FiringCurve(pos_b , spks_b , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 1, 'minPeak' , 0.1);
       
        peak_bef = find(map_bef==stats_bef.peak); 
        com_bef=nan; 
        if sum(sum(~isnan(map_bef.fieldX)))>0
               pf_lim = stats_bef.fieldX; 
               pc_frmap = map_bef.rate; 
               field = pc_frmap(pf_lim(1,1):pf_lim(1,2));
                
              %Center of mass of the field 
               c = 1:size(field,2);                             
               com = sum(c .* field) / sum(field);
               com_bef = com + pf_lim(1,1); % back in general scale bins 
               clear pf_lim pc_frmap field c com
        end 
        
        %Compute the firing map after 
        %Choose the same # of laps from the end 
        spks_a = Restrict(spks,in_lap(end-shock_lap:end,:)); % restrict spikes to after shock
        pos_a = Restrict(pos,in_lap(end-shock_lap:end,:));
        [map_after , stats_after] = FiringCurve(pos_a , spks_a , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 1, 'minPeak' , 0.1);

        peak_after = find(map_after==stats_after.peak); 
        
        com_after=nan; 
        if sum(sum(~isnan(map_after.fieldX)))>0
               pf_lim = stats_after.fieldX; 
               pc_frmap = map_after.rate; 
               field = pc_frmap(pf_lim(1,1):pf_lim(1,2));
                
              %Center of mass of the field 
               c = 1:size(field,2);                             
               com = sum(c .* field) / sum(field);
               com_after = com + pf_lim(1,1); % back in general scale bins 
               clear pf_lim pc_frmap field c com
        end 
        
    else 
        disp('Not enough laps before shock')
        peak_bef=nan; com_bef=nan; peak_after=nan; com_after=nan;
        
    end     
    



end