function [R] = triggered_activity(patterns , cond , SpikeTrain , limits , events, graph , normalization , selection , baseline , Sizes)
% DO DOCUMENTATIOOOOOOOOOOOn!!!!!!!
% This function calculate Reactivation Strength (van de Ven et al (2016)).
% Plot the activity sourrounding the events introduced.
%
% [R] = triggered_activity(patterns , cond , SpikeTrain , limits , events, graph , normalization , selection , baseline)
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
% graph: logical, if 1 will give you a mean graph.
%
% normalization: logical, to define if assemblie activity will be zscored.
%
% selection: logical vector, it contains 1 or 0 to define which assemblie
%            will be considered
%
% Size: vector, should contains the number of SU per structure. useful if
%       we want to filter ripples according to the spiking activity across
%       strcutures.
%
% --- OUTPUT ---
% R: column vector storing all the Reactivation Strength values for sleep
%    and awake periods.
%
%    S1 ---- *C1*  ---- S2 ----  ---- C2 ---- S3
%
%    * Assemblies detected during this session. It could do the same if C2
%    is the condition of interest. Just be sure you are choosing correctly
%    the type and config inputs.
%
%    Structure:  mean(S1)-mean(S2)   mean(C1)   mean(C2)
%                         R1           A11        A21
%                         R2           A12        A22
%                         R3           A13        A23
%                         ...          ...        ...
%
% requirments:
%       assembly_activity.m from Lopes-dos-Santos et al 2013 (*, see below)
%
% Morici Juan Facundo 09/2023
% also see DOI:10.1016/j.jneumeth.2013.04.010 (*)

bins = SpikeTrain(:,1); % time bins
spks = SpikeTrain(:,2:end); 
dt = bins(2)-bins(1);
win = round(limits/dt);
check = (limits*2)/dt; % construct the maximun of points for each curve

a = assembly_activity(patterns(:,cond) , spks');

% if not(isempty(Sizes))
% %     criteria = InIntervals(bins,[events(:,1) events(:,3)]);
% %     criteria = and(criteria , and(sum(spks(:,1:Sizes(1)),2) , sum(spks(:,Sizes(1)+1:end),2)>0));
%     criteria = ToIntervals(bins,and(sum(spks(:,1:Sizes(1)),2)>0 , sum(spks(:,Sizes(1)+1:end),2)>0));
%     criteria = InIntervals(events(:,2) , criteria);
%     events = events(criteria,:);
% end


if not(isempty(selection))
    a = a(selection,:);
end

if normalization % normalization
    if not(isempty(baseline))
        baseline = InIntervals(bins,baseline);
        gain = nanmean(a(:,baseline),2);
        S = nanstd(a(:,baseline),1,2);
    else
        a = zscore(a,1,2);
    end
end


is = [];
for i = 1:size(events,1)
    if size(events,2)>1
        [h p] = min(abs(events(i,2)-bins));
    else
        [h p] = min(abs(events(i)-bins));
    end
    is = [is , p]; clear h p
end



R = [];
for i = 1:size(a,1)
    tmp = [];
    for ii = 1:size(events)
        t = is(ii);
        if and(t-win>1 , t+win<size(bins,1))
%           [f ff] = max(a(i,t-win : t+win));
%           t = t-win + ff; clear f ff
          tmp = [tmp ; a(i,t-win : t+win)];
        end
        clear t
    end
    if not(isempty(baseline))
        R = [R ; (nanmean(tmp)-gain(i))./S(i)];
    else
        R = [R ; nanmean(tmp)];
    end
    clear tmp
    
end

if logical(graph)
    figure,plot(R')
end

end