function [between] = Between_pc(pos_ave,spks_ave,pos_rew,spks_rew,bin_size,sigma,Xedges)
%
% Between_pc - Randomly splits spike and position data of the two conditions (ave and rew) in two and randomly selects one group 
% of each condition and calcualtes remapping parameters(spatial correlation, firing rate change and rate overlap)between them. 
% Repeat 500 times  
%
% [between] = Between_pc(pos_ave,spks_ave,pos_rew,spks_rew,bin_size,sigma,Xedges)
%
% ---> INPUTS
% pos: matrix, timestamps in sec and position, one variable for each condition (ave and
% rew)
% spks: vector, spike time in sec, one variable for each condition (ave and
% rew)
% bin_size: bin size in sec (1 recommended). Time window to binnarized pos
%           and spike data
% sigma: int,FiringCurve parameter for smoothing rate maps
% Xedges: int,FiringCurve parameter to define the nBins of the rate map
%
% ---> OUTPUTS
% between: vector with the mean value of the 500 iterations for the 3
%          remapping parameters(spatial correlation, firing rate change and
%          rate overlap).
% 
%
%other functions:FiringCurve (FMAtoolbox), Bins_half
%
%Azul Silva, 2023
between = nan(500,3);

for c=1:500

    [groups_ave] = Bins_half(pos_ave,spks_ave,bin_size); 
    [groups_rew] = Bins_half(pos_rew,spks_rew,bin_size);
 
    %Choose randomly one of the groups for the comparison
    id = round(rand)*2; 
    if id==0 % not the best way
        id =1;
    end 
    
    time_x_ave = groups_ave{id,1};
    tspk_ave = groups_ave{id,2};
    
    time_x_rew = groups_rew{id,1};
    tspk_rew = groups_rew{id,2};
    
    %Calculate remapping parameters 
    [curveA] = FiringCurve(time_x_ave, tspk_ave , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4 , 'minPeak' , 0.2);
    [curveR] = FiringCurve(time_x_rew, tspk_rew , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4 , 'minPeak' , 0.2);
                    
    fr_A= nanmean(curveA.rate);
    fr_R= nanmean(curveR.rate);
                    
    %Fr change
    fr_change = abs((fr_A - fr_R)/(fr_A + fr_R));

    %Rate overlap
    if fr_A<=fr_R 
        overlap = fr_A/fr_R;
    else 
        overlap = fr_R/fr_A;
    end
    
    %Spatial  corr
    s = corrcoef(curveA.rate, curveR.rate);
    spatial = s(1,2);
    
    %Save 
    between(c,1)=spatial;
    between(c,2)=fr_change;
    between(c,3)=overlap;
    
   
end 
between = mean(between);
end