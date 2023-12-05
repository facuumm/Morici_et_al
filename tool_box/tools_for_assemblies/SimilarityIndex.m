function [m,p] = SimilarityIndex(patterns1,patterns2,permutation)
% This function calculates the simmilarity Index between patterns stored
% in patterns1 and pattenrs2.
% This index is equal to the absolute value of the inner-product of their
% weight vectors (Almeida-Filho et al., 2014).
%
% INPUTS
% patterns1 / patterns2: matrix, rows: cells weigth, column: assemblie.
%
% OUTPUT
% m: matrix storing the similarity index between patterns from both inputs.
% p: matrix, 1 if the value is higher comapred with the 99th percentile
%
%
% other functions: SimilaritySurrogate
% Morci Juan Facundo 09/2023


m = [];
for i  = 1 : size(patterns1,2)
    x = patterns1(:,i);
    row = [];
    for ii = 1 : size(patterns2,2)
        y = patterns2(:,ii);
        row = [row ; abs(inner(x,y))];
        clear y
    end
    clear x
    m = [m , row];   
end
    p = SimilaritySurrogate(patterns1,patterns2,1000);
    p = m>p;
end