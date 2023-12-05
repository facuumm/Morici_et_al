% Count of spikes occurring within an event
%
% count_spks(spks, clusters, start, stop)
%
% --- INPUTS ---
% spks: Matrix, first column must contain cluster id, second column times
%
% clusters: column vector, cluster ids
%
% start: float, starting time stamp of the event (same unit as spikes)
%
%
% stop: float, ending time stamp of the event (same unit as spikes)
%
% --- OUTPUTS ---
% c: column vector, number of spikes of each cluster occurring within the
%    event.
%
% Morici Juan Facundo, 09/06/2023


function c = count_spks(spks, clusters, start, stop)

c = [];

for i = 1:length(clusters)
    c = [c ; length(Restrict(spks(spks(:,1)==clusters(i),2),[start stop]))];
end

end