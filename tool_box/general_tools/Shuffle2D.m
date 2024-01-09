function shuffle = Shuffle2D(matrix)
% This function shuffle the data keeping the inter-event intervals.
%
% --- Inputs ---
% matrix: matrix, first column contains the TimeStamps and second column
%         contains the other value.
%
% --- OUTPUT ---
% shuffle: matrix, contains the shuffled data.
% Morici Juan Facundo 01/2024


    a = min(matrix(:,1)) + rand*10;% 1° time stamp
    isi = diff(matrix(:,1));
    isi_per = isi(randperm(length(isi)));%shuffle isi
    
    % generation of final timestamps by add the first random time
    shuffle1 = [a; isi_per];
    shuffle1 = cumsum(shuffle1);
    
    % shuffleling the data in 2nd column
    shuffle2 = matrix(randperm(size(matrix,1)),2);
    
    
    %Storing data
    shuffle = [shuffle1 , shuffle2];
    
    
    clear shuffle1 schuffle2
end 