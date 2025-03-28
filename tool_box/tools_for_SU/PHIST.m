function [p , t] = PHIST(Times1,Times2,d,b,sm)
% This function construct a CCG using Ripples timeStamps to lock the
% occurrence of SU activity.
%
% syntax: [p , t] = PHIST_Ripple_SU(Times1,Times2,baseline,d,b,sm,mode)
%
% INPUTS
% Times1: reference, needs to be ripple matrix as FMA toolbox export it.
%
% Times2: referenced.
%
% d = int, duration fo the time window sourrounding reference.
%
% b = float, binsize for CCG.
%
% sm = int, smooth factor.
%
% OUTPUT
% Raster Plot, in y-axis Times1 events, y x-axis time.
%
% Morci Juan Facundo 12/2023


[s,ids,groups] = CCGParameters(Times1,ones(length(Times1),1),Times2,ones(length(Times2),1)*2);
[ccg,T] = CCG(s,ids,'binSize',b,'duration',d,'smooth',sm,'mode','ccg');
ccg = (ccg(:,1,2)./length(Times1))./b;

% Zscore of density
% ccg = zscore(ccg(:,1,2)./sum(Times2));

p = ccg;
t = T;
clear s ids groups 

end