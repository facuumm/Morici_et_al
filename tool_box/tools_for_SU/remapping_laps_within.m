function [within_corr_matrix_1] = remapping_laps_within(in_lap,pos_1,spks_1, sigma, Xedges,min_size, min_peak, min_time)
%
% 
% remapping_laps
% Computes spatial correlation across laps. 
% ---> INPUTS
%in_lap: 2 column matrix with beginning and end of each lap. Spks and pos will be restricted to these laps 
% pos_: matrix, position timestamps in sec
% spks_: vector, spike time in sec
% sigma: int,FiringCurve parameter for smoothing rate maps
% Xedges: int,FiringCurve parameter to define the nBins of the rate map
%
% ---> OUTPUTS
% 

%other functions:FiringCurve (FMAtoolbox)
%Azul Silva, 2024
    laps_map = [];
    for l=1:size(in_lap,1)
        %Keep only spk and pos from one
         spks = Restrict(spks_1, in_lap(l,:));
         pos = Restrict(pos_1, in_lap(l,:)); 
         [curve1] = FiringCurve(pos, spks , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , min_size , 'minPeak' , min_peak,'minTime',min_time);
         laps_map = [laps_map; curve1.rate];% saving stack of lap rate maps 
     end 
%         figure;imagesc([0:20:140], [1:1:size(laps_map,1)],laps_map), colormap 'jet'; title('Odd or Even');
     within_corr_matrix_1 = corrcoef(laps_map'); 
end 

