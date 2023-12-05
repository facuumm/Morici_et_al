function info= spatial_info_rate(RateMap, OccMap)
% Spatial Information calculation based on Skaags et al 1993
%
% info = spatial_info_rate(RateMap, pos)
%
% --- INPUTS ---
% RateMap: vector, rate map calculated using FiringMap_LinearTrack
% OccMap: vector, occupancy Map caqlculated using FiringMap_LinearTrack
%
% --- OUTPUTS ---
% info: float, Spatial Information value
% Based on Skaags et al: Theta Phase Precession in Hippocampal Neuronal ...
% Silva Azul 2023 - edited by Morici Juan Facundo

OccMap(isnan(RateMap))=nan;    
OccMap=OccMap/nansum(OccMap(:));
sel=(OccMap>0);

if sum(sel)==0
    info=0;
else
    mr = nansum(RateMap(sel).*OccMap(sel));
    info=nansum(OccMap(sel).*RateMap(sel)/mr.*log2(RateMap(sel)/mr));
end
end
