function shuffle = ShuffleSpks(spks)
% This function shuffle spikes 100 keeping the inter-spike intervals.
%
% --- Inputs ---
% spks: vector, spikes times
%
% --- OUTPUT ---
% shuffle: vector, shuffle spike times.
% Silva Azul 2023 - edited by Morici Juan Facundo
    a = min(spks) + rand*10;% 1° tiempo
    isi = diff(spks);
    isi_per = isi(randperm(length(isi)));%shuffle isi
    
    %genero tiempos random para los spikes que respetan el isi original
    spks_R = [a; isi_per];
    shuffle =cumsum(spks_R);

end 