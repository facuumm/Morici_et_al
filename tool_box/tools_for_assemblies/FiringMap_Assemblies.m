function [map] = FiringMap_Assemblies(patterns , cond , SpikeTrain , Is , th , pos , events , Nbins)
% This function contruct a 2D map crossing position and Activity peak for
% each assemblie.
%
% [map] = FiringMap_Assemblies(patterns , cond , SpikeTrain , Is , th , pos)
%
% --- Inputs ---
% patterns: matrix with the weigths for each cell in each assembly.
%           Structure: Single-Units x Assemblies (rows x column)
%
% cond: logical, to select which pattern will be used (1, include / 0, not include)
%       Row vector with the same elements as number of assemblies.
%       Example,
%       patterns:     A1    A2    A3
%                 SU1 0.6   0.1   0.1       ---> cond: 0   1   0
%                 SU2 0.3   0.6   0.1       In this case, only the second
%                 SU3 0.1   0.3   0.3       assembly (A2) will be selected.
%                 SU4 0.1   0.1   0.6
%
% SpikeTrain: matrix, First column, time bins, Rest columns Spike Counts
%             Example:   t    SU1   SU2   SU3   SU4
%                        0     3     5     1     3    This Spike Train was
%                       0.1    1     2     1     0    constructed using a
%                       0.2    3     5     1     3    0.1 sec time bin.
%                       ...   ...   ...   ...   ...
%
% Is: structure, it contains the periods of time (in sec) of interest.
%     Struture:
%                 Is.baseline     Should be a vector with the same length
%                 Is.aversive     of the SpikeTrain with logical values to
%                 Is.reward       include or not the time bin.
%                 Is.runaversive
%                 Is.runreward
%
%     Example:
%             is.aversive:
%                         b1  b2  b3  b4  b5  b6  b7  b8  b9
%                         0   0   1   1   1   1   0   0   0
%
% th: float/int, threshold value for peaks detection. Note that this
%     function zscored the assemblies activity.
%
% pos: matrix, it contains the timestamps and the X and Y postions.
%
% Nbins: int, number of spatial bins
%
% --- OUTPUT --- CHECKEAR ESTOOOOO
% map: matrix, It contains all the maps for each assemblie.
%
%     Example:                  Spatial Bins
%                       Bin1    Bin2    Bin3    ...   BinN
%             map1 --> Rate1   Rate2   Rate3   ...   RateN
%             map2 --> Rate1   Rate2   Rate3   ...   RateN
%             mapM --> Rate1   Rate2   Rate3   ...   RateN
%
% requirments:
%       assembly_activity.m from Lopes-dos-Santos et al 2013 (*, see below)
%       FiringMap from FMA toolbox
%
% Morici Juan Facundo 12/2023

bins = SpikeTrain(:,1);
dt = bins(2)-bins(1); % delta time
spks = SpikeTrain(:,2:end);

%calculation of assemblies activity
a = assembly_activity(patterns(:,cond) , spks');

%zscored of assemblies activity
a = zscore(a,1,2);

map = [];
for i = 1:size(a,1)
    %detection of peaks
    [pks,loc] = findpeaks(a(i,:),bins,'MinPeakHeight',th);
    
    % restriction of peaks and positions to events
    loc = Restrict(loc,events);
    pos = Restrict(pos,events);
    
    % normalization of positions (0 to 1)
    pos(:,2) = pos(:,2)-min(pos(:,2));
    pos(:,2) = pos(:,2)./max(pos(:,2));
    
    %firing curve calculation
    [m,stats] = FiringCurve(pos,loc,'nBins' , Nbins , 'smooth',2);
    
    %storing the constructed map
    map = [map ; m.rate];
    
end
end