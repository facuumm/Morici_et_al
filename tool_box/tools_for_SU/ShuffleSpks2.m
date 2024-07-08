function shuffle = ShuffleSpks2(spks)
% This function shuffle spikes without keeping the inter-spike intervals.
%
% --- Inputs ---
% spks: vector, spikes times
%
% --- OUTPUT ---
% shuffle: vector, shuffle spike times.
% Morici Juan Facundo 07/2024
    
    b = mod(1:size(spks,1),2)*2-1;
    b = b(randperm(size(b,2)));
    shuffle = spks + (rand(size(spks)).*b');% 1° time

end 