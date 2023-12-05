function [p , t] = PHIST(Times1,Times2,baseline,d,b,sm,mode)
% This function construct a CCG using Times1 to lock the occurrence 
% of Times2
%
% INPUTS
% Times1 / Times2: Column vector. Times1 = reference / Times2 = referenced.
% baseline = TimeStamps [Begining End ; ... ... ; Begining End] of baseline
%            periods for Gain calculation. if mode is not 'Gain', introduce
%            an empty vector.
% d = int, duration fo the time window sourrounding reference.
% b = float, binsize for CCG.
% sm = int, smooth factor.
% mode = Str, 'Gain', 'NormGain' or empty for Firing Rate
%
% OUTPUT
% Raster Plot, in y-axis Times1 events, y x-axis time.
%
% Morci Juan Facundo 12/2023


[s,ids,groups] = CCGParameters(Times1,ones(length(Times1),1),Times2,ones(length(Times2),1)*2);
[ccg,T] = CCG(s,ids,'binSize',b,'duration',d,'smooth',sm,'mode','ccg');
ccg = (ccg(:,1,2)./length(Times1))./b;

if strcmp(mode,'Gain')
    bb =  length(Restrict(Times2,baseline));
    bb = bb/sum(baseline(:,2)-baseline(:,1));
    ccg = ccg ./bb;
end

if strcmp(mode,'NormGain')
    bb =  length(Restrict(Times2,baseline));
    bb = bb/sum(baseline(:,2)-baseline(:,1));
    ccg = ccg ./bb;
    ccg = ccg - min(ccg);
    ccg = ccg./max(ccg);
end

p = ccg;
t = T;
clear s ids groups 

end