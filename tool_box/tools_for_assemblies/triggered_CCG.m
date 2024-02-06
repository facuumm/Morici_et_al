function [R , t] = triggered_CCG(patterns , cond , SpikeTrain , duration, binSize , events, th , Baseline, mode)
% This function construct a peri-event histogram using ripple peaks and
% assemblies activity peaks. If you wanna change the output normalization,
% please check details in PHIST function and change it in this code.
%
% [R , t] = triggered_CCG(patterns , cond , SpikeTrain , limits , events, th , Baseline)
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
% duration: int, it define in seconds the total duration of the CCG
%
% binSize: float, duration of binSize
%
% events: matrix, it contains Ripples timestamps (in sec).
%         Example:   Begining1    Peak1    End1
%                    Begining2    Peak2    End2
%                       ...        ...     ...
%
% Baseline: matrix, it contains begining and end of reference periods to
%           substract the events and determine the baseline epochs for Gain
%           normalization.
%           Example:   Begining1    End1
%                      Begining2    End2
%                         ...       ...
%
% mode: str, 'Gain', 'NormGain', 'Probability', FiringRate if is empry
%
% --- OUTPUT ---
% R: matrix, contains the PHIST constructed using the inputs
% t: vector, time vector
%
% requirments:
%       assembly_activity.m from Lopes-dos-Santos et al 2013 (*, see below)
%       PHIST from this toolbox
%       FMA toolbox
%
% Morici Juan Facundo 12/2023

bins = SpikeTrain(:,1);
spks = SpikeTrain(:,2:end);
if not(isempty(Baseline))
    baseline = SubtractIntervals(Baseline,[events(:,1)-0.05 events(:,3)]+0.05);
else
    baseline = [];
end
a = assembly_activity(patterns(:,cond) , spks');
a = zscore(a,1,2);

R = [];
for i = 1:size(a,1)
    [pks,loc] = findpeaks(a(i,:),bins,'MinPeakHeight',th);
    if size(events,2)>1
        [p , t] = PHIST(events(:,2),loc,baseline,duration,binSize,1,mode);
    else
        [p , t] = PHIST(loc,events,baseline,duration,binSize,1,mode);
    end
    R = [R , p];
    
end

end