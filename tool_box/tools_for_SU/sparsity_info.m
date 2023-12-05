


function sparsity= sparsity_info(RateMap, OccMap)
% Sparsity calculation based on Skaags et al 1993 and Jung et al 1994
%
% sparsity = sparsity(RateMap, OccMap)
%
% --- INPUTS ---
% RateMap: vector, rate map calculated using FiringMap_LinearTrack
% OccMap: vector, occupancy Map caqlculated using FiringMap_LinearTrack
%
% --- OUTPUTS ---
% info: float, Sparsity value
% Based on Skaags et al 1993: Theta Phase Precession in Hippocampal Neuronal ...
% Based on Jung et al 1994: Comparison of Spatial Firing Characteristics..
% Morici Juan Facundo

OccMap(isnan(RateMap))=nan;    
OccMap=OccMap/nansum(OccMap(:));
sel=(OccMap>0);
% m = mean(RateMap);

if sum(sel)==0
    sparsity=0;
else
%     RateMap = RateMap./m;
    sparsity = (nansum(RateMap(sel).*OccMap(sel)))^2;
    sparsity = sparsity/(nansum(RateMap(sel).^2.*OccMap(sel)));
    
end
end
