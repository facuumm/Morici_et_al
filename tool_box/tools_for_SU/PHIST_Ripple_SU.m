function [p , t] = PHIST_Ripple_SU(Times1,Times2,baseline,d,b,sm,mode)
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


[s,ids,groups] = CCGParameters(Times1(:,2),ones(length(Times1),1),Times2,ones(length(Times2),1)*2);
[ccg,T] = CCG(s,ids,'binSize',b,'duration',d,'smooth',sm,'mode','ccg');
ccg = (ccg(:,1,2)./length(Times1))./b;

if strcmp(mode,'Gain')
    baseline = SubtractIntervals(baseline,[Times1(:,1)-0.05 Times1(:,3)+0.05]);
    bb =  length(Restrict(Times2,baseline));
    bb = bb/sum(baseline(:,2)-baseline(:,1));
    ccg = ccg ./bb;
end

if strcmp(mode,'NormGain')
%     baseline = SubtractIntervals(baseline,[Times1(:,1)-0.05 Times1(:,3)+0.05]);
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