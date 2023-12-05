function [p,theta,rbar,delta] = ThetaModulation(Spks,cluster,l,Periods)
% This function determine if a SU is phase-locked to an specific phase of the
% filtered lfp. Only SU that fires at leats 25 times will be included.
%
% syntax = [p,theta,rbar,delta] = ThetaModulation(Spks,cluster,lfp,Periods)
%
% INPUTS
% Spks = timestamps of each spike (1st column: cluster , 2nd: timestamp) 
% cluster = vector containing the id clusters of interest
% l = matric containing the timestamps (1st column) and lfp (2nd column)
% Periods = Periods of sleep that are of your interest[begining end]
%
% OUTPUT
% p = p-value for the probability that the data are drawn from a uniform
%     distribution rather than a unimodal distribution of unknown mean
%     direction .
%
% theta = non-weighted circular mean
%
% rbar = mean resultant length
%
% delta = circular dispersion
%
% other functions: rayleigth (Daniel Rizzuto)
%                  Restrict and SubsSubtractIntervals (FMAtoolbox)
% Morci Juan Facundo 11/2023

p = [];
theta = [];
rbar = [];
delta = [];

l(:,2) = angle(hilbert(l(:,2)));

for i = 1: size(cluster,1)
    spks = Restrict(Spks(Spks(:,1)==cluster(i),2),Periods);
    iii = Restrict(l,Periods);
    
    if not(isempty(spks)) & not(isempty(iii))
        [interpolated,discarded] = Interpolate(iii,spks);
        if length(interpolated)>25
            [p(i),theta(i),rbar(i),delta(i)] = rayleigh(interpolated(:,2));
        else
            p(i)=NaN;
            theta(i)=NaN;
            rbar(i)=NaN;
            delta(i)=NaN;
        end
    else
        p(i)=NaN;
        theta(i)=NaN;
        rbar(i)=NaN;
        delta(i)=NaN;
    end
    clear interpolate discarded spks iii
end
    
end
