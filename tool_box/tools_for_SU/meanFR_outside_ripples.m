function [m] = meanFR_outside_ripples(ripplesTS , periods , spks , binsize)
% This function returns the mean firing rate outside of ripples.
%
% syntax: [m] = meanFR_outside_ripples(ripplesTS,limits,spks)
%
% --- INPUTS ---
% ripplesTS: matrix, contains the time stamps of ripples [start peak stop].
%            It should follow [events , TS].
%
% periods: matrix, it should contains the timestamps to substract the
%         ripples.
%
% spks: column vector, it contains the timestaps of the SU to evaluate.
%
% --- OUPUTS ----
% m: float, it contain the mean firing rate outside ripples.
%
%
% other functions coming from FMA toolbox
% Facundo Morici 08/2024

bufferedripples = [ripplesTS(:,1)-0.1 ripplesTS(:,3)+0.1];
bufferedripples = Restrict(bufferedripples, periods);
[baseline,ind]=SubtractIntervals(periods,bufferedripples,'strict','on');

totalbaselinetime = sum(baseline(:,2)-baseline(:,1));

% Constructing an spiketrain
[spikes,bins]=binspikes(spks(:,1),1/binsize,[baseline(1,1) baseline(end,2)]);
m = InIntervals(bins,baseline);
m = nanmean(spikes(m)./binsize);

% Counting the spikes
% m = Restrict(spks(:,1),baseline);
% m = size(m,1)/totalbaselinetime;

end