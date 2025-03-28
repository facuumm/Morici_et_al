function [P] = assembly_peaks_detection(patterns ,SpikeTrain ,th)
% This function detects the activation peaks of different assemblies.
%
% [P] = reactivation_strengthV2(patterns ,SpikeTrain ,th , type)
%
% --- Inputs ---
% patterns: matrix with the weigths for each cell in each assembly.
%           Structure: Single-Units x Assemblies (rows x column)
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
%
% --- OUTPUT ---
% P: cell containin N column vectors for N assemblies storing all the
%    Activation peaks in seconds
%
% requirments:
%       assembly_activity.m from Lopes-dos-Santos et al 2013 (*, see below)
%
% Morici Juan Facundo 02/2024

% Organizing the variables
bins = SpikeTrain(:,1);
dt = bins(2)-bins(1); % delta time
spks = SpikeTrain(:,2:end);

% Assemblies activity calculation
a = assembly_activity(patterns , spks');
A = a;

% Normalization of the assemblies activity
a = zscore(a,1,2);

P = cell(1,size(a,1));
for i = 1:size(a,1)
    th = prctile(a(i,:),95);
    % using the time vector for plotting
    [pks,loc] = findpeaks(a(i,:),bins,'MinPeakHeight',th);

    P{i} = loc;
    clear pks loc
end

end