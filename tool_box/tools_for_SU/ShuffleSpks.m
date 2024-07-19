function shuffle = ShuffleSpks(spks)
% This function shuffle spikes keeping the inter-spike intervals.
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
    
    %Create random spks keeping the original isi
    spks_R = [a; isi_per];
    shuffle =cumsum(spks_R);
end 