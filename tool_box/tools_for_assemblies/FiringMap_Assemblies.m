function [map pc between within] = FiringMap_Assemblies(patterns , cond , SpikeTrain , th , pos , events , Nbins , Between , Within)
% This function contruct a 1D map crossing position and Activity peak for
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
% th: float/int, threshold value for peaks detection. Note that this
%     function zscored the assemblies activity.
%
% pos: matrix, it contains the timestamps and the X and Y postions.
%
% events: cell, it should contains events separated in both contidions you
%         wanna compare.
%         Required structure: events{1} ---> movement periods during cond1
%                             events{2} ---> movement periods during cond2
%
% Nbins: int, number of spatial bins
%
% Within: logical, if is True it will calculate spatial correlation, firing
%         rate change and rate overlap within conditions.
%
% Between: logical, if is true will calculate same parameters described
%          above but between conditions.
%
% --- OUTPUT ---
% map: structure, It contains all the maps for each assemblie. It contains
%      two different tags: 'cond1' and 'cond2'. 
%
%     Example:                  Spatial Bins
%                       Bin1    Bin2    Bin3    ...   BinN
%             map1 --> Rate1   Rate2   Rate3   ...   RateN
%             map2 --> Rate1   Rate2   Rate3   ...   RateN
%             mapM --> Rate1   Rate2   Rate3   ...   RateN
%
% pc: logical vector, contains 1 if the assemblie showed place tunning.
%     We defined a field larger that 4 spatial bins as a criteria.
%
% requirments:
%       assembly_activity from Lopes-dos-Santos et al 2013 (*, see below)
%       FiringMap from FMA toolbox
%       Between_pc and Within_pc from this toolbox, made by Azul SIlva
%
% Morici Juan Facundo 12/2023
%
% * Lopes-dos-Santos V, Ribeiro S, Tort AB. Detecting cell assemblies in
% large neuronal populations. J Neurosci Methods. 2013.
% doi: 10.1016/j.jneumeth.2013.04.010. PMID: 23639919.

bins = SpikeTrain(:,1);
dt = bins(2)-bins(1); % delta time
spks = SpikeTrain(:,2:end);

%calculation of assemblies activity
a = assembly_activity(patterns(:,cond) , spks');

%zscored of assemblies activity
a = zscore(a,1,2);

% Generation of structures to store outputs
map.cond1 = []; map.cond2 = [];
pc.cond1 = []; pc.cond2 = [];
within.cond1 = []; within.cond2 = [];
between = [];

for i = 1:size(a,1)
    %detection of peaks
    [pks,loc] = findpeaks(a(i,:),bins,'MinPeakHeight',th);
    
    for ii = 1:2
        % restriction of peaks and positions to events
        l = Restrict(loc,events{ii});
        p = Restrict(pos,events{ii});
        
        % normalization of positions (0 to 1)
        p(:,2) = p(:,2)-min(p(:,2));
        p(:,2) = p(:,2)./max(p(:,2));
        
        %firing curve calculation
        [m,stats] = FiringCurve(p,l,'nBins' , Nbins , 'smooth' , 2 , 'minPeak' , 0.2 , 'minSize' , 4);
        
        if ii == 1 %storing the constructed map
            map.cond1 = [map.cond1 ; m.rate];
        else
            map.cond2 = [map.cond2 ; m.rate];
        end
        
        %check if field is larger than 4 bins
        %if so, it will be considered as a place assemblie
        if ii == 1
            if not(isempty(stats.field))
                if sum(stats.field(:,:,1))>= 4
                    pc.cond1 = [pc.cond1 ; true];
                else
                    pc.cond1 = [pc.cond1 ; false];
                end
            else
                pc.cond1 = [pc.cond1 ; false];
            end
        else
            if not(isempty(stats.field))
                if sum(stats.field(:,:,1))>= 4
                    pc.cond2 = [pc.cond2 ; true];
                else
                    pc.cond2 = [pc.cond2 ; false];
                end
            else
                pc.cond2 = [pc.cond2 ; false];
            end
        end
        
        if Within
            loc1 = Restrict(loc,events{1});
            pos1 = Restrict(pos,events{1});
            % normalization of positions (0 to 1)
            pos1(:,2) = pos1(:,2)-min(pos1(:,2));
            pos1(:,2) = pos1(:,2)./max(pos1(:,2));
            
            loc2 = Restrict(loc,events{2});
            pos2 = Restrict(pos,events{2});
            % normalization of positions (0 to 1)
            pos2(:,2) = pos2(:,2)-min(pos2(:,2));
            pos2(:,2) = pos2(:,2)./max(pos2(:,2));
            
            if ii == 1 %storing the constructed map
                [w] = Within_pc(pos1,loc1,1,2,Nbins);
                within.cond1 = [within.cond1 ; w];
            else
                [w] = Within_pc(pos2,loc2,1,2,Nbins);
                within.cond2 = [within.cond2 ; w];
            end
            clear loc1 loc2 pos1 pos2 w
        end
        
    end
    
    if Between
        loc1 = Restrict(loc,events{1});
        pos1 = Restrict(pos,events{1});
        % normalization of positions (0 to 1)
        pos1(:,2) = pos1(:,2)-min(pos1(:,2));
        pos1(:,2) = pos1(:,2)./max(pos1(:,2));
        
        loc2 = Restrict(loc,events{2});
        pos2 = Restrict(pos,events{2});
        % normalization of positions (0 to 1)
        pos2(:,2) = pos2(:,2)-min(pos2(:,2));
        pos2(:,2) = pos2(:,2)./max(pos2(:,2));
            
        [b] = Between_pc(pos1,loc1,pos2,loc2,1,2,Nbins);
        between = [between ; b];
        clear loc1 loc2 pos1 pos2 b
    end
    clear loc pks
end

end