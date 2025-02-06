function shortenedIntervals = shortenIntervals(intervals, maxDuration, fromLast)
    % intervals: An N x 2 matrix where each row is [start_time, end_time] in seconds.
    % maxDuration: The total maximum allowable duration in seconds.
    % fromLast: A boolean parameter (true or false) that decides from which direction to start modifying intervals.
    %           If false (default), start modifying from the first interval. 
    %           If true, start modifying from the last interval.

    % Calculate the total duration of all intervals.
    totalDuration = sum(intervals(:, 2) - intervals(:, 1));

    % If the total duration is already less than or equal to maxDuration, return the intervals as is
    if totalDuration <= maxDuration
        shortenedIntervals = intervals;
        return;
    end

    % Initialize variables to track the cumulative time
    cumulativeDuration = 0;
    N = size(intervals, 1);  % Number of intervals
    shortenedIntervals = []; % Initialize the empty list to store intervals within the criteria
    
    if fromLast
        % If starting from the last interval, loop backward
        for i = N:-1:1
            intervalDuration = intervals(i, 2) - intervals(i, 1);
            cumulativeDuration = cumulativeDuration + intervalDuration;
            
            % If cumulative time exceeds maxDuration, modify the current interval
            if cumulativeDuration > maxDuration
                excessDuration = cumulativeDuration - maxDuration;
                modifiedInterval = intervals(i, :);
                modifiedInterval(1) = intervals(i, 1) + excessDuration;  % Adjust the start time to fit
                shortenedIntervals = [modifiedInterval; shortenedIntervals]; % Add the modified interval
                break;  % Exit the loop since we no longer need to process earlier intervals
            else
                % Add the interval to the list if it doesn't exceed the maxDuration
                shortenedIntervals = [intervals(i, :); shortenedIntervals];
            end
        end
    else
        % If starting from the first interval, loop forward
        for i = 1:N
            intervalDuration = intervals(i, 2) - intervals(i, 1);
            cumulativeDuration = cumulativeDuration + intervalDuration;
            
            % If cumulative time exceeds maxDuration, modify the current interval
            if cumulativeDuration > maxDuration
                excessDuration = cumulativeDuration - maxDuration;
                modifiedInterval = intervals(i, :);
                modifiedInterval(2) = intervals(i, 2) - excessDuration;  % Adjust the end time to fit
                shortenedIntervals = [shortenedIntervals; modifiedInterval]; % Add the modified interval
                break;  % Exit the loop since we no longer need to process further intervals
            else
                % Add the interval to the list if it doesn't exceed the maxDuration
                shortenedIntervals = [shortenedIntervals; intervals(i, :)];
            end
        end
    end
end