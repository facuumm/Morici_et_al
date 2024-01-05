function [R] = cumulative_activation_strength(patterns , cond , SpikeTrain , th , events , duration,baseline)
% This function find the peaks of the zscored assemblies strength and
% returns the cummulative of peaks in function of the time.
% If you include events, it will contatenate the events you introduce.
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
% th: float/int, threshold value for peaks detection. Note that this
%     function zscored the assemblies activity.
%
% events: two-elements matrix. [Begining End] of the segment of session I
%         am interested in analyse
%
% duration: int, determine the time bins duration (in sec).
%
% baseline: two-elements matrix. [Begining End] of the baseline segment I
%         	will use to calculate the fold change of reactivation.
%           If is empty, output will not be normalized
%
% --- OUTPUT ---
% R: matrix, contains the activation peaks cumsum in bins which length was
%    determined by 'duration' input.
%
% requirments:
%       assembly_activity.m from Lopes-dos-Santos et al 2013 (*, see below)
%
% Morici Juan Facundo 09/2023
% see: DOI: 10.1038/s41467-023-43939-z
% also see DOI:10.1016/j.jneumeth.2013.04.010 (*)

bins = SpikeTrain(:,1);
dt = bins(2)-bins(1); % delta time
spks = SpikeTrain(:,2:end);
prc = 90; % percentile to check if the Strenght is higher 
iterations = 100; % iterations to create the surrogated distribution

a = assembly_activity(patterns(:,cond) , spks');
a = zscore(a,1,2);
R = [];

for i = 1:size(a,1)
    
    [p,l] = findpeaks(a(i,:),bins,'MinPeakHeight',th);
    
    if not(isempty(events))
        loc = [];
        tmp = 0;
        for ii = 1:size(events,1)
            L = Restrict(l,events(ii,:));
            P = Restrict([l,p'],events(ii,:));
            if ii>1
                tmp = tmp + (events(ii,1) - events(ii-1,2));
                L = L - tmp;
            end
            loc = [loc ; L P(:,2)]; clear L P
        end
    end
    
    % 
    edges = [events(1,1) : duration : events(end,2)-tmp]; clear tmp
    [bincounts,ind]= histc(loc(:,1) , edges);
    
    result = zeros(1,size(edges, 2));
    for i = 1 : size(edges, 2)
        if not(isempty(loc(ind==i,2)))
            result(i) = result(i) + nanmean(loc(ind==i,2));
        else
            result(i) = result(i);
        end
    end
    
    if not(isempty(baseline))
        m = Restrict([l,p'],baseline);
        m = nanmean(m(:,2));
        result =  cumsum(result./m);
%         result = result - m;
    else
        result =  cumsum(result);
%         result = result;
    end
%     % binning peaks according to duration input
%     [t,b]=binspikes(loc,1/duration,[events(1,1) events(end,2)-tmp]); 
%     t = cumsum(t); %cumsum computation
%     
%     %store R output
    R = [R ; result]; clear t
    
end

end