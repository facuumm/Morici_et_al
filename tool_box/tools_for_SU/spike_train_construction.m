function [SpksTrains , Bins , C] = spike_train_construction(Spks, clusters, type, binSize, limits, events, normalization, smooth)
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
% type: column vector, int value. Same length of clusters
% (1 or 0, is or is not of cellular type of interest)
%
% binSize: float, time window (sec) for Spike Train construction
%
% limits: [start stop] of the recording segment
%
% events: matrix, [Star Stop] of each event 
%
% normalization: logic value, (True for zscore)
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
SpksTrains = [];
C = [];
for ii=1:length(clusters)
    cluster = clusters(ii,1);
    celltype = logical(type(type(:,1) == cluster,2));
    if celltype
        spks = Spks(Spks(:,1)==cluster,2);
        [tmp,bins]=binspikes(spks,freq,limits);
        
        if smooth
            try
                tmp = gaussfilt(bins,tmp,(binSize)/sqrt(12));
            catch
                tmp = tmp;
            end
        end
        
        SpksTrains = [SpksTrains , tmp];
        C = [C ; cluster];
        clear spkscelltype cluster
    end
    clear spks
    
end

if normalization
%     times = SubtractIntervals(limits,events);
%     S = mean(SpksTrains(InIntervals(bins,times),:),1);
%     SS = std(SpksTrains(InIntervals(bins,times),:),1);
%     SpksTrains = (SpksTrains - S) ./ SS;
    SpksTrains = zscore(SpksTrains);
    clear SS S
else
    SpksTrains = SpksTrains./binSize;
end

if isempty(events)
    events = limits;
end

SpksTrains = SpksTrains(InIntervals(bins,events),:);
Bins = bins(InIntervals(bins,events));
end
