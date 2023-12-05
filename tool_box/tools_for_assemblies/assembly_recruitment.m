function [P] = assembly_recruitment(patterns , cond , SpikeTrain , limits , events, th , normalization)
% This function calculate Reactivation Strength (van de Ven et al (2016)).
% Plot the activity sourrounding the events introduced.
%
% [P] = assembly_recruitment(patterns , cond , SpikeTrain , limits , events, th , normalization)
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
% limits: float, limits in seconds to calculate the graph
%
% events: vector, it contains the time events (in sec) of interest.
%
% th: float, define the threshold for peak detection
%
% normalization: logical, if you want to zscored the assemblies activity.
%
% --- OUTPUT ---
% P: vector storing the percentage of events that recruit assemblies.
%
%   Example,
%            A1   A2   A3
%            90   50   90
%           ...  ...  ...
%
% requirments:
%       assembly_activity.m from Lopes-dos-Santos et al 2013 (*, see below)
%
% Morici Juan Facundo 09/2023

bins = SpikeTrain(:,1);
spks = SpikeTrain(:,2:end);
dt = bins(2)-bins(1);
win = round(2/dt);

a = assembly_activity(patterns(:,cond) , spks');

if normalization % normalization
    a = zscore(a,1,2);
end

P = [];
for i = 1:size(a,1)
    
    tmp = 0;
    for ii = 1:size(events,1)
        [~ , t] = min(abs(bins-events(ii)));
        if sum(a(i,t-win : t+win)>th)>0
            tmp = tmp+1;
        end
        clear t
    end
    
    P = [P ; (tmp/size(events,1))*100];
    clear tmp
end


end