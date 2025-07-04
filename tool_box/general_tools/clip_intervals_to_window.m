function clipped = clip_intervals_to_window(events, window)
% CLIP_INTERVALS_TO_WINDOW Clips events to a specified time window.
%
%   clipped = CLIP_INTERVALS_TO_WINDOW(events, window)
%
%   Inputs:
%       events - Nx2 matrix of event intervals [start, end] in seconds
%       window - 1x2 vector [window_start, window_end]
%
%   Output:
%       clipped - Mx2 matrix of intervals clipped to the window

    % Extract window bounds
    window_start = window(1);
    window_end = window(2);

    % Initialize result
    clipped = [];

    for i = 1:size(events,1)
        start_time = events(i,1);
        end_time = events(i,2);

        % Skip events completely outside the window
        if end_time <= window_start || start_time >= window_end
            continue;
        end

        % Clip start and end to the window
        new_start = max(start_time, window_start);
        new_end   = min(end_time, window_end);

        % Append to output
        clipped = [clipped; new_start, new_end];
    end
end