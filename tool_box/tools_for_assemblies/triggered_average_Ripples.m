function [average , time] = triggered_average_Ripples(events,patterns,cond,spks,limits,components,counts)
% This function construct a triggered-average of assemblies activity to the
% Ripples. 
%
% Syntax: [average , time] = triggered_average_Ripples(events,patterns,cond,spks,limits)
%
% INPUTS
% events: column vector, it contains the event to tune the activity (in sec).
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
% components: string, to select what component of the projector matrix to
%             use. 
%                   'all' --> dorsal+ventral activations
%                   'dHPC' --> dorsal activations
%                   'vHPC' --> ventral activations
%                   'both' --> just joint activations
%
% counts: vector, contains the number of dorsal and ventral SU (eg. [10 5]).
%
% OUTPUT
% average: matrix, it contains the average response.
%
% time: row vector, it contains the time axis for the CCG plot.
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

% templates for deleting components from the projector
template = ones(size(patterns,1),1);
template(counts(1)+1:end) = 0;
template1 = ones(size(patterns,1),1);
template1(1:counts(1)) = 0;

% Assemblies Activity calculation
if strcmp(components,'all')
    a = assembly_activity(patterns(:,cond) , spks(:,2:end)');
elseif strcmp(components,'dHPC')
    a = assembly_activity_only_joint(patterns(:,cond) , spks(:,2:end)',template, template1,'dHPC');
elseif strcmp(components,'vHPC')
    a = assembly_activity_only_joint(patterns(:,cond) , spks(:,2:end)',template, template1,'vHPC');
elseif strcmp(components,'both')
    a = assembly_activity_only_joint(patterns(:,cond) , spks(:,2:end)',template, template1,'both');
end


% a = zscore(a,1,2);

% itaration across Assemblies and Events
for i = 1 : size(a,1)
%     a(i,:) = smoothdata(a(i,:),'gaussian','SmoothingFactor',0.75,'SamplePoints',bins); %smooth the data a bit
    tmp = [];
    for ii = 1 : size(events,1)
        [~ , index] = min(abs(events(ii)-bins));
        tmp = [tmp , a(i,index+limits_pos(1) : index+limits_pos(2))'];
    end
    average = [average , nanmean(tmp,2)];
end

end