function [average , time , surrogate] = triggered_average_assemblies_response(events,patterns,cond,spks,limits)
% This function construct a triggered-average of assemblies activity to the
% shocks. Also, it create a surrogate distribution and define if the value
% of response within the shock is higher than the 90th percentile. If so,
% it will be considered as shock-responsive assembliy.
%
% INPUTS
% events: column vector, it contains the begining of the shock (in sec).
%
% patterns: matrix, it contians the assemblies patterns. (Neuron x Assembies)
%
% cond: row vector, it contains 1 or 0 if you wanna include or not the
%       pattern. Same shape in column dimention as patterns input.
%
% spks: matrix, First column, time bins. From the second column: SpikeTrains.
%
% limits: vector, it contains the limits to contrais the Triggered-average.
%         E.g: [-1 2]
%
% OUTPUT
% average: matrix, it contains the average response.
%
% time: row vector, it contains the time axis for the CCG plot.
%
% surrogate: row vector, it contians information regarding to the
%            shock-responsiveness of each assemblie.
%
% other functions: assembly_activity from Lopes-dos-Santos et al (2013)
%                  10.1016/j.jneumeth.2013.04.010 
%
% Morci Juan Facundo 04/2024

% Definition of variables
bins = spks(:,1);
dt = bins(2)-bins(1);
limits_pos = [round(limits(1)/dt) round(limits(2)/dt)];
average = [];
time = [limits(1) : dt : limits(2)+dt];

% Assemblies Activity calculation
a = assembly_activity(patterns(:,cond) , spks(:,2:end)');
a = zscore(a,1,2);

% itaration across Assemblies and Events
for i = 1 : size(a,1)
    tmp = [];
    for ii = 1 : size(events,1)
        [~ , index] = min(abs(events(ii)-bins));
        tmp = [tmp , a(i,index+limits_pos(1) : index+limits_pos(2))'];
    end
    average = [average , nanmean(tmp,2)];
end

% Surrogate construction
surrogate = [];
for i = 1 : 100
    s = ShuffleSpks(events);
    for ii = 1 : size(a,1)
        tmp = [];
        for iii = 1 : size(s,1)
            [~ , index] = min(abs(s(iii)-bins));
            tmp = [tmp , a(ii,index+limits_pos(1) : index+limits_pos(2))'];
        end
        surrogate(:,ii,i) = nanmean(tmp,2);
    end
end

% Definition of shock-responsiveness
[~ , i] = min(abs(time-0));
[~ , ii] = min(abs(time-1));
surrogate = mean(surrogate(i:ii,:,:));
surrogate = prctile(surrogate,90,3);
m = nanmean(average(i:ii,:));

surrogate = m > surrogate;

clear s i ii m
end