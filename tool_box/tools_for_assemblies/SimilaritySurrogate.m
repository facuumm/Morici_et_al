function p = SimilaritySurrogate(patterns1,patterns2,permutations)
% This function determine the 99th percentile of the Simminlarity Index
% form patterns1 and pattenrs2 after shuffling the weights of each pattern
% across neurons 100 times (Almeida-Filho et al., 2014).
%
% INPUTS
% patterns1 / patterns2: matrix, rows: cells weigth, column: assemblie.
% permutations/ int, number of times to perform shuffleing of weigths
%
% OUTPUT
% m: matrix storing 99th percentile of a surrogated.
%
% Morici Juan Facundo 09/2023
m = [];
for iii = 1 : permutations
mm = [];
for i  = 1 : size(patterns1,2)
    x = patterns1(:,i);
    x = x(randperm(length(patterns1)));
    row = [];
    for ii = 1 : size(patterns2,2)
        y = patterns2(:,ii);
        y = y(randperm(length(patterns2)));
        row = [row ; abs(inner(x,y))];
        clear y
    end
    clear x
    mm = [mm , row];
end
m = cat(3,m,mm);
clear mm
end
clear iii i ii

p = [];
for i  = 1 : size(patterns1,2)
    row = [];
    for ii = 1 : size(patterns2,2)
        row = [row ; prctile(m(ii,i,:),99.9)];
    end
    p = [p , row];
end
