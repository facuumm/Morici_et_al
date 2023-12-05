% Firing Map calculated the Firing curve for a linear track maze
%
% [RateMap,RateMapRaw, RateMapRaw_fil, OccMap,OccMapRawR, Nspikes, NspikesRawR, Xs]=firingMap_JFM(positionAndTime,Tspk_mov,dt ,Xedges, minTime)
%
% --- INPUTS ---
% pos: matrix (column1: time / column2: position in x axis)
% spks : column vector containing the spike times
% dt: float value sampling periode of the position information
% Xedges: number of bins in space
% minTime: float value to define minimum occupancy to consider or not the bin
%          Lisa Roux used 0.15 sec
%
% --- OUTPUTS ---
% RateMap: matrix containing the Rate curve filtered by minTime
% RateMapRaw: matrix containing raw map without filtering
% RateMapRaw_fil: RateMapRaw filtering the position
% OccMap: matrix containing the time per spatial bin, filtered
% OccMapRawR: matrix containing the time per spatial bin, without filtering
% Nspikes: matrix containing counts of spikes per spatial bin, filtered
% NspikesRawR: matrix containing counts of spikes per spatial bin, without
% filtering
% Xs: interpolation for each spike times to postion times
% X: Normalized bin positions

function [RateMap,RateMapRaw, RateMapRaw_fil, OccMap,OccMapRawR, Nspikes, NspikesRawR, Xs , X]=firingMap_JFM(pos,spks,dt ,Xedges, minTime)

% Space normalization (0 to 1)
% pos(:,2) = pos(:,2) - min(pos(:,2));
% pos(:,2) = pos(:,2) ./ max(pos(:,2));

% Occupancy in spatial bins (counts)
[N_bin,X] = histcounts(pos(:,2),Xedges);
X = X(2:end);
% Transformation to time (sec)
OccMap = N_bin*dt;

% Storage occupancy without filtering
OccMapRaw = OccMap;

% Storage Occupancy after filtering
miSmooth = 9;
OccMap = imgaussfilt(OccMap,2,'FilterSize',miSmooth,'Padding','circular'); 

% Template pointing bins out of time criteria
indx_mascara = (OccMapRaw<minTime) & (N_bin < 4);

% Find interpolated position of each spike:
Xs = interp1(pos(:,1),pos(:,2),spks);

% Count of interporalted spike times across space
[Nspikes] = histcounts(Xs,Xedges);

% Store raw counts curve
NspikesRaw = Nspikes;

% filtering Nspikes
Nspikes = imgaussfilt(Nspikes,2,'FilterSize',miSmooth,'Padding','replicate');

% Restriction of bins out of time criteria
Nspikes(indx_mascara) = NaN;
OccMap(indx_mascara) = NaN;
OccMapRawR = OccMapRaw;
OccMapRawR(indx_mascara)= NaN;
NspikesRawR = NspikesRaw;
NspikesRawR(indx_mascara)= NaN;

%Rate Map
RateMap = Nspikes./OccMap;

%Rate Map RAW
RateMapRaw = NspikesRawR./OccMapRawR;
%Filtro Rate Map RAW
RateMapRaw_fil = imgaussfilt(RateMapRaw,2,'FilterSize',miSmooth,'Padding','replicate');

end