%Firing Map
%spikes en c/bin filtrado/ tiempo en cada bin filtrado  como Roux.
%INPUT:positionAndTime (col1 tiempo, col2posx), Tspk_mov (time spk en movimiento)
%Xedges (areana size)
%OUTPUT:

function [RateMap,RateMapRaw, RateMapRaw_fil, OccMap,OccMapRawR, Nspikes, NspikesRawR, Xs]=firingMap_JFM(positionAndTime,Tspk_mov,dt ,Xedges, minTime)

% Creo matriz de XedgesxYedges con cantidad de veces que estuvo el animal en cada bin:
[N_bin,~, ~] = histcounts(positionAndTime(:,2),Xedges);

%Calculo el tiempo(sec)que paso en cada bin de la matriz
OccMap = N_bin*dt;

% %Pongo eps en los 0:
% OccMap(OccMap==0) = eps;

%Guardo Ocupancia sin filtro y sin restricciones:
OccMapRaw = OccMap;

%Aplico filtro
miSmooth = 9;
OccMap = imgaussfilt(OccMap,2,'FilterSize',miSmooth,'Padding','circular'); 

%Creo mascara con restricciones place map sobre crudos:
indx_mascara = (OccMapRaw<minTime) & (N_bin < 4);

%Find interpolated position of each spike:
Xs = interp1(positionAndTime(:,1),positionAndTime(:,2),Tspk_mov);

% Creo una matriz con cantidad de spikes por bines:
[Nspikes] = histcounts(Xs,Xedges);

%Guardo Nspikes sin filtro:
NspikesRaw = Nspikes;

%Aplico filtro a Nspikes
Nspikes = imgaussfilt(Nspikes,2,'FilterSize',miSmooth,'Padding','replicate');

%Aplico restriccion Nspikes y OccuMap
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