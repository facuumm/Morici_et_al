function [P R] = reactivation_strength_members(patterns , cond , SpikeTrain, th ,thresholded, numberSU , timestamps)
% This function will calculate the Reactivation Rate normalized by the
% presleep considering only bins that showed cross-firing of members
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

% variables for cumulative construction
w = 60; % bin size in sec
tempo = 2700; %maximal value in sec


bins = SpikeTrain(:,1);
dt = bins(2)-bins(1); % delta time
spks = SpikeTrain(:,2:end);
thresholded = thresholded(:,cond);


a = assembly_activity(patterns(:,cond) , spks');
a = zscore(a,1,2);

R = [];
P = [];
template = spks >0;

for ii = 1:size(a,1)
    
    template2 = template .* thresholded(:,ii)';
    template2 = and(sum(template2(:,1:numberSU(1)),2)>0 , sum(template2(:,numberSU(1)+1:end),2)>0);
    [pks.all,loc.all] = findpeaks(a(ii,:),bins,'MinPeakHeight',th);
    
    template2 = bins(template2);
    
    template2 = ismember(loc.all,template2);
    template3 = loc.all(template2);
    
    timePost = sum(timestamps.PostSleep(:,2) - timestamps.PostSleep(:,1))/60;
    timePre = sum(timestamps.PreSleep(:,2) - timestamps.PreSleep(:,1))/60;
    strength = ((size(Restrict(template3,timestamps.PostSleep),1) / timePost) / (size(Restrict(template3,timestamps.PreSleep),1) / timePre));
    
    R = [R ; strength];
    
    % concatenation of peaks
    temporal = [];
    iterator = 0;
    for ii = 1 : size(timestamps.PostSleep,1)
        t = Restrict(template3,timestamps.PostSleep(ii,:));
        t = (t - timestamps.PostSleep(ii,1));
        temporal = [temporal ; t + iterator];
        iterator = (timestamps.PostSleep(ii,2)-timestamps.PostSleep(ii,1))+iterator; % concatenation of time segments
        clear t
    end
    h = histcounts(temporal,'BinEdges',[0 : w : tempo]);
    P = [P , h'./(size(Restrict(template3,timestamps.PreSleep),1) / timePre)]; clear h
    
    
    clear template2 template3 result temporal iterator
end

end