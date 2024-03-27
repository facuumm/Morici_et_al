function [SpksTrains , Bins] = spike_train_construction_counts(Spks, clusters, binSize, limits)
% Spike Trains matrix construction
%
% [SpksTrains , Bins , C] = spike_train_construction(Spks, clusters, type, binSize, limits, events, normalization, smooth)
%
% --- INPUTS ---
% Spks: Column vector, spikes times
%         (1st column: cluster id / 2nd column: time stamps)
%
% clusters: column vector, contains the cluster ids
%
% binSize: float, time window (sec) for Spike Train construction
%
% limits: [start stop] of the recording segment
%
% smooth: if true, convolution with gaussian kernel will be applied.
%         SD will be defined as binSize/sqrt(12)
%         Refs. Kruskal et am 2007 and Gido M van de Ven et al 2016
%
% --- OUTPUTS ---
% SpikeTrains: Spike Trains restricted to events periods.
%              rows: Time bins / Columns: clusters
%
% Bins: Time bins restricted to events
%
%C: column vector containing the id cluster
%
% Morici Juan Facundo, 11/06/2023
% Other funtions: binspikes, gaussfilt 

freq = 1/binSize;
[SpksTrains , Bins] = binspikes(Spks(ismember(Spks(:,1),clusters),2),freq,limits);

end
