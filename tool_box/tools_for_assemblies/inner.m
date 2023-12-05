
function y = inner(a,b);
% This function compute the inner product of two vectors a and b.  
% y = inner(a,b)
% Input: The two vectors a and b, they must be same length
% Output: The value of the inner product of a and b.
% Morici Juan Facundo 17/08/2023
c=0;      % intialize the variable c
for i = 1 : length(a)   % start the loop
        c = c + a(i) * b(i);  % update c by the k-th product in inner product
end
y = c ;     % print value of c = inner product