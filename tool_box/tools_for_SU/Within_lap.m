function [within_mean] = Within_lap(pos_tmp,spks_tmp,in_lap,sigma,Xedges)
%
% Within_lap - Randomly splits laps in two and calcualtes remapping parameters
% (spatial correlation, firing rate change, rate overlap and peak shift) within them. Repeat 500 times  
%
% ---> INPUTS
% pos_tmp: matrix, position timestamps in sec
% spks_tmp: vector, spike time in sec
% in_lap: 2 column matrix with beginning and end of each lap 
% sigma: int,FiringCurve parameter for smoothing rate maps
% Xedges: int,FiringCurve parameter to define the nBins of the rate map
%
% ---> OUTPUTS
% within_lap: vector with spatial correlation, fr_change, overlap, shift 
%other functions:FiringCurve (FMAtoolbox)
%
%Azul Silva, 2023

within = nan(500,4);

for c=1:500
    %Randomly splits laps 
    temp = in_lap(randperm(size(in_lap,1)),:); 
    half = ceil(size(in_lap,1)/2); 
   
    g1= temp(1:half,:);
    g2=temp(half+1:end,:);
   
    %Split spk and pos from each lap 
    spks_1 = Restrict(spks_tmp,g1);
    pos_1 = Restrict(pos_tmp,g1); 
    
    spks_2 = Restrict(spks_tmp,g2);
    pos_2 = Restrict(pos_tmp,g2); 
    
    %Calculate remapping parameters 
    [curve1,stats1] = FiringCurve(pos_1, spks_1 , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4 , 'minPeak' , 0.2);
    [curve2,stats2] = FiringCurve(pos_2, spks_2 , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4 , 'minPeak' , 0.2);
                    
    fr_1= nanmean(curve1.rate);
    fr_2= nanmean(curve2.rate);
                    
    %Fr change
    fr_change = abs((fr_1 - fr_2)/(fr_1 + fr_2));

    %Rate overlap
    if fr_1<=fr_2 
        overlap = fr_1/fr_2;
    else 
        overlap = fr_2/fr_1;
    end
    
    %Spatial  corr
    s = corrcoef(curve1.rate, curve2.rate);
    spatial = s(1,2);
    
    %Peak shift 
    shift = abs(stats1.x(1) - stats2.x(1)); 
    
%     Control plot 
    figure;hold on; 
    subplot(2,1,1);imagesc(curve1.rate), colormap 'jet'
    subplot(2,1,2);imagesc(curve2.rate), colormap 'jet'
    sgtitle('Firing curvs laps');
    
    %Save 
    within(c,1)=spatial;
    within(c,2)=fr_change;
    within(c,3)=overlap;
    within(c,4)=shift;
end 
    within_mean = nanmean(within);
end