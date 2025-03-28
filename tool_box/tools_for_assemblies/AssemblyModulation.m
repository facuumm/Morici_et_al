function [pInc pDec surp] = AssemblyModulation(Event,Times,Periods)
% This function determine if a SU is modulated by the ripples.
%
% syntax = [pInc pDec surp] = ThetaModulation(RipplesTS,Spks,cluster,Periods)
%
% INPUTS
% Event, matrix containing the start and end of the event to see if
%        modulates the assembly
% Times = timestamps of each activation
% Periods = Periods of the session to calculate de Basleine[begining end]
%
% OUTPUT
% pInc = probability to be up-modulated
% pDec = probability to be down-modulated
% suprise = positive if pInc > pDec
%
% other functions: poissonTest (Eran Stark)
%                  Restrict and SubsSubtractIntervals (FMAtoolbox)
%
% Morci Juan Facundo 03/2025

pInc = [];
pDec = [];
surp = [];


[baseline,ind]=SubtractIntervals(Periods,Event,'strict','on');

totalbaselinetime=sum(baseline(:,2)-baseline(:,1));
baselineTimes=Restrict(Times(:,1),baseline);

% Restrict further to the third of ripples with the highest amplitude
totaleventtime=sum(Event(:,2) - Event(:,1));
eventTimes=Restrict(Times,Event);
baselineTimes=length(baselineTimes);
eventTimes=length(eventTimes);

if baselineTimes~=0 & eventTimes~=0
    [pInc pDec surp] = poissonTest(baselineTimes/totalbaselinetime,eventTimes,totaleventtime);
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