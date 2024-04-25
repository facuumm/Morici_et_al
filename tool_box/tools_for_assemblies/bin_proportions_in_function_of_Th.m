function [P] = bin_proportions_in_function_of_Th(patterns , cond , SpikeTrain, thresholded, numberSU)
% This function will use different thresholds (th: 1-15) to deteact peaks and
% will calculates the percentage of bins that have cross-members activation
%
% [R] = reactivation_strength(patterns , cond , SpikeTrain , is , th , type)
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
% thresholded: matrix containing which SU is considered as a member.
%
% numberSU: row vector, first value is the number of dHPC SUs and the second
%           one is the number of vHPC SUs.
%
% --- OUTPUT ---
% P: Matrix, each row is an assemblie, and each column represent the th
%    used. It can go from 1 to 15 std.
%
% requirments:
%       assembly_activity.m from Lopes-dos-Santos et al 2013 (*, see below)
%
% Morici Juan Facundo 04/2024

bins = SpikeTrain(:,1);
dt = bins(2)-bins(1); % delta time
spks = SpikeTrain(:,2:end);
thresholded = thresholded(:,cond);


a = assembly_activity(patterns(:,cond) , spks');
a = zscore(a,1,2);

P = [];
template = spks >0;
for i = 1:16
    tmp = [];
    for ii = 1:size(a,1)
        
        template2 = template .* thresholded(:,ii)';
        template2 = and(sum(template2(:,1:numberSU(1)),2)>0 , sum(template2(:,numberSU(1)+1:end),2)>0);
        [pks.all,loc.all] = findpeaks(a(ii,:),'MinPeakHeight',i-1);
        
        template3 = template2(loc.all);
        result = (sum(template3) / size(template3,1))*100;
        tmp = [tmp ; result]; clear template2 template3 result
    end
    P = [P , tmp]; clear tmp
end

end