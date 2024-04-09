function probability = bayesian_decoder(RateMap, nSpks, start, stop)
% Probability calculation of being at different positions using spiking
% activity during an event based on Bayesian probability.
% Wirtsgafter & Wilson 2021: Bayesian Algorithmic Decoding of Acceleration
% and Speed Software (BADASS). https://doi.org/10.1016/j.simpa.2021.100125
%
% probability = bayesian_replay(RateMap, nSpks, start, stop)
%
% --- INPUTS ---
% RateMap: Matrix, rate map for all the cells across the space
%          (rows: cell ids / columns: spatial bins)
%
% Spks: column vector, number of spks ocurring during the event for each
%        unit in the RateMap. ** length(nSpks) == length(RateMap) **
%
% start: float, starting time stamp of the event (same unit as spikes)
%
%
% stop: float, ending time stamp of the event (same unit as spikes)
%
% --- OUTPUTS ---
% probability: row vector, probability of each spatial position of being
%              replayed in the event.
%       

% Morici Juan Facundo, 09/06/2023


% Variables definition
x = RateMap;
p = nSpks;
t = stop - start;


f = x .^ p;
y = prod(f , 1);
y = y .* exp(-t*sum(x));
c = sum(y);
probability = y ./ c;
clear x p t f y c
end