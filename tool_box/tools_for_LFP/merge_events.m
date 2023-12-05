% Merge time events if they occurs within a certain time window.
%
% [Merged] = merge_events(x, value)
%
% --- INPUTS ---
% x: Matrix, Time stamps of events to be merged
%    (1st column: Start / 2nd column: End).
%
% value: float, inter-replay period to defined which events would be merged
%
% --- OUTPUTS ---
% Merged: Matrix, Time stamps of marged events 
%         (1st column: Start / 2nd column: End).
%
%
% Morici Juan Facundo, 14/06/2023

function [Merged] = merge_events(x, value)


        % Merge Replays if inter-replay period is too short
        firstPass = x;
        secondPass = firstPass;
        iri = x(2:end,1) - x(1:end-1,2);
        toMerge = iri<value;
        while any(toMerge),
            % Get indices of first ripples in pairs to be merged
            rippleStart = strfind([0 toMerge'],[0 1])';
            % Incorporate second ripple into first in all pairs
            rippleEnd = rippleStart+1;
            secondPass(rippleStart,2) = secondPass(rippleEnd,2);
            % Remove second ripples and loop
            secondPass(rippleEnd,:) = [];
            iri = secondPass(2:end,1) - secondPass(1:end-1,2);
            toMerge = iri<value;
        end
        
        Merged = secondPass;
        
end