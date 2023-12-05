function q = SkaggsRandom(spks, pos, Nspikes, OccMap, Xedges)
% This function asign random positions to the spikes 100 times, each time 
% skaggs was calculated and then 90-quantile was defined.
%
% --- Inputs ---
% spks: vector, spikes times
% pos: matrix, positions in x and their respective time stamps (column1: time / column2: position in x axis)
% Nspikes: matrix, counts of spikes per spatial bin
% OccMap: matrix, counts of bin occupancy
%
% --- OUTPUT ---
% q: float, 50-quantile from random skaags values.
% Silva Azul 2023 - edited by Morici Juan Facundo

info = zeros(100,1);
for S = 1:100
    a = min(spks) + rand*10;% 1° tiempo
    isi = diff(spks);
    isi_per = isi(randperm(length(isi)));%shuffle isi
    
    %genero tiempos random para los spikes que respetan el isi original
    spks_R = [a; isi_per];
    spks_R =cumsum(spks_R);
    
    %         %Interpolo posiciones para los spike time random
    %         XsPROBEper = interp1(pos(:,1),pos(:,2), spks_R);
    %
    %         %Creo una matriz con cantidad de spikes por bines:
    %         [NspikesPROBE_per] = histcounts(XsPROBEper,Xedges);
    %         %Aplico filtro a Nspikes
    %         NspikesPROBE_per = imgaussfilt(NspikesPROBE_per,2,'FilterSize',9,'Padding','replicate');
    %         %Aplico restricciones
    %         NspikesPROBE_per(isnan(Nspikes))= nan;
    %
    %         %Calculo rate
    %         RateMap_per = NspikesPROBE_per./OccMap;
    %
    % %         subplot(2,1,1);
    % %         imagesc(OccMap)
    % %         subplot(2,1,2);
    % %         imAlpha=ones(size(RateMap_per));
    % %         imAlpha(isnan(RateMap_per))=0;
    % %         imagesc(RateMap_per,'AlphaData',imAlpha)
    %
    %         info(S,1) = spatial_info_rate(RateMap_per,OccMap);%OccMap_per
    [curve , stats] = FiringCurve(pos , spks_R , 'smooth' , 2 , 'nBins' , Xedges , 'minSize' , 4 , 'minPeak' , 0.2);
    info(S,1) = stats.specificity;

end
%   mean = mean(info)
%   deviation = std(info)

    q = quantile(info, 0.50);
    
%     figure(88)
%     histogram(info,20)
end 