function [pInc pDec surp] = RippleModulation_assemblies(RipplesTS,Spks,Periods)
% This function determine if a SU is modulated by the ripples.
%
% syntax = [pInc pDec surp] = ThetaModulation(RipplesTS,Spks,cluster,Periods)
%
% INPUTS
% RipplesTS = ripples timestamps, same output from RippleDetection (FMAtoolbox)
% Spks = timestamps of each spike (1st column: cluster , 2nd: timestamp)
% cluster = vector containing the id clusters of interest
% Periods = Periods of sleep that are of your interest[begining end]
%
% OUTPUT
% pInc = probability to be up-modulated
% pDec = probability to be down-modulated
% suprise = positive if pInc > pDec
%
% other functions: poissonTest (Eran Stark)
%                  Restrict and SubsSubtractIntervals (FMAtoolbox)
%
% Morci Juan Facundo 11/2023

pInc = [];
pDec = [];
surp = [];

spks=Spks;
iii = Restrict(RipplesTS,Periods);

bufferedripples=[iii(:,1)-0.1 iii(:,3)+0.1];
[baseline,ind]=SubtractIntervals(Periods,bufferedripples,'strict','on');

totalbaselinetime=sum(baseline(:,2)-baseline(:,1));
baselinespikes=Restrict(spks(:,1),baseline);

% Restrict further to the third of ripples with the highest amplitude
totalrippletime=sum(iii(:,3)-iii(:,1));
ripplespikes=Restrict(spks(:,1),[iii(:,1) iii(:,3)]);
ncellbaselinespikes=length(baselinespikes);
ncellripplespikes=length(ripplespikes);

if ncellbaselinespikes~=0 & ncellripplespikes~=0
    [pInc pDec surp] = poissonTest(ncellbaselinespikes/totalbaselinetime,ncellripplespikes,totalrippletime);
else
    pInc=NaN;
    pDec=NaN;
    surp=NaN;
end


end

% 
% 
% function [ pInc pDec surp ] = poissonTest( baseRate, inCount, inTime )
% 
% if nargin < 3 || ~isequal( size( baseRate ), size( inCount ), size( inTime ) )
%     return
% end
% siz         = size( baseRate );
% baseRate    = baseRate( : );
% inCount     = inCount( : );
% inTime      = inTime( : );
% 
% lambdas     = baseRate .* inTime;
% pInc        = 1 - poisscdf( inCount - 1, lambdas );
% pDec        = poisscdf( inCount, lambdas );
% surp        = log10( ( pDec + eps )./ ( pInc + eps ) );
% 
% pInc        = reshape( pInc, siz );
% pDec        = reshape( pDec, siz );
% surp        = reshape( surp, siz );
% 
% end