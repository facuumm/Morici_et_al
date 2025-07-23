
%Firing Map
% #spikes en c/bin filtrado/ tiempo en cada bin filtrado  como Roux. 
%INPUT:positionAndTime (col1 tiempo, col2posx, col3 posy), Tspk_mov (time spk en movimiento)
%Xedges e Yedges (areana size)
%OUTPUT:

function [RateMap,RateMapRaw, RateMapRaw_fil, OccMap,OccMapRawR, Nspikes, NspikesRawR, Xs, Ys] = FiringMap2D(positionAndTime,Tspk_mov, Xedges, Yedges, minTime)

    % Creo matriz de XedgesxYedges con cantidad de veces que estuvo el animal en cada bin:  
     [N_bin,~, ~] = histcounts2(positionAndTime(:,2),positionAndTime(:,3),Xedges,Yedges);
     
    %Calculo la Frec de muestreo de la posicion
     dt = (mean(diff(positionAndTime(:,1))));

    %Calculo el tiempo(sec)que paso en cada bin de la matriz 
     OccMap = N_bin*dt; 

     %Pongo eps en los 0:
     OccMap(OccMap==0) = eps;
 
     %Guardo Ocupancia sin filtro y sin restricciones:
     OccMapRaw = OccMap;

     %Aplico filtro 
     miSmooth = 9;
     OccMap = imgaussfilt(OccMap,2,'FilterSize',miSmooth,'Padding','circular'); %??

     %Creo mascara con restricciones place map sobre crudos: 
     indx_mascara = (OccMapRaw<minTime) & (N_bin < 4);

       %Find interpolated position of each spike:
        Xs = interp1(positionAndTime(:,1),positionAndTime(:,2),Tspk_mov);
        Ys = interp1(positionAndTime(:,1),positionAndTime(:,3),Tspk_mov);

        % Creo una matriz con cantidad de spikes por bines: 
         [Nspikes] = histcounts2(Xs,Ys,Xedges,Yedges);

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