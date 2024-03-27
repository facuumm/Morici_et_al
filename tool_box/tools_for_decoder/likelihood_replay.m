function [R H P] = bayesian_replay(RateMap, nSpks, time , dt , d)
% Replay detection and significance testing.
% This function calculates likelihood R between the replayed trajectory and
% a constant-velocity trajectory inferred with the initial and final
% positions. To determine if R is higher than expected by chance, this
% function constructs two distributions of Rs values by:
%       1) Circularly shifting the estimate at each time bin by a random distance 1000 times.
%       2) Randomly permuting tuning curves 1000 times.
%
% Method developed by Davidson et al 2009, doi:10.1016/j.neuron.2009.07.027
%
% USAGE
%       R = bayesian_replay(RateMap, nSpks, start, stop)
%
% --- INPUTS ---
% RateMap: Matrix, rate map for all the cells across the space
%          (rows: cell ids / columns: spatial bins)
%
% nSpks: column vector, number of spks ocurring during the event for each
%        unit in the RateMap. ** length(nSpks) == length(RateMap) **
%
% time: row vector, contain the start and end of the putative-replay event
%
% dt: float, time interval for bining.
%
% d: float, Distance margin of acceptance
%
% --- OUTPUTS ---
% R: float, likelihood value.
%
% H: row vector, output of both shuffling methods 
%
% P: row vector, probability of each spatial position of being
%              replayed in the event.
%
% Morici Juan Facundo, 09/06/2023


% Variables definition
x = RateMap;
p = nSpks;
Pr = [];
time = lin
for i = 1:size(time,2)-1
    t = time(i+1) - time(i);
    f = x .^ p;
    y = prod(f , 1);
    y = y .* exp(-t*sum(x));
    c = sum(y);
    Pr = [Pr , y ./ c];
end

clear x p t f y c
end