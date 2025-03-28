function [surrogate , percentile] = surrogate_assembly_activity(Times,Periods)
% This function returns the surrogate of rates and the 99th percentile of
% that surrogate distribution.
%
% syntax = [surrogate , percentile] = surrogate_assembly_activity(Times,Periods)
%
% INPUTS
% Times = timestamps of each activation
% Periods = Periods of the session to calculate Rate
%
% OUTPUT
% surrogate, vector containing the rate distribution after shuffling Times1
%
% percentile = 99th percentile of surrogate
%
% other functions: Restrict and SubsSubtractIntervals (FMAtoolbox)
%
% Morci Juan Facundo 03/2025


surrogate = [];
for i = 1:500
    s = ShuffleSpks(Times);
    s = Restrict(s,Periods);
    s = length(s)/sum(Periods(:,2)-Periods(:,1));
    surrogate = [surrogate ; s];
end

percentile = prctile(surrogate,99);

end