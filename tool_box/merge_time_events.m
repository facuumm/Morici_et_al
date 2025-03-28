function merged_events = merge_time_events(events, periods)
    % merge_events - This function merges overlapping or adjacent time periods
    % and returns the events that occur within the merged periods.
    %
    % Syntax:
    %   merged_events = merge_events(events, periods)
    %
    % Inputs:
    %   events  - A vector containing event times (1xN array)
    %   periods - A matrix with each row representing a time period as [start, end] (Mx2 matrix)
    %
    % Outputs:
    %   merged_events - A vector containing the events that occur during the merged periods
    %
    % Example:
    %   events = [1, 2, 3, 7, 8, 9, 12, 14, 15];
    %   periods = [1 5; 3 6; 7 10; 12 16];
    %   merged_events = merge_events(events, periods);
    %   The output will be merged events that fall within the continuous periods.

    % Sort the periods based on the start time (ascending order)
    [sorted_periods, idx] = sortrows(periods, 1);
    
    % Initialize an empty array for storing merged periods
    merged_periods = [];
    
    % Start by considering the first period as the current merged period
    current_period = sorted_periods(1, :);
    
    % Iterate through the sorted periods and merge overlapping or adjacent ones
    for i = 2:size(sorted_periods, 1)
        % Check if the current period overlaps with or is adjacent to the next one
        if sorted_periods(i, 1) <= current_period(2)
            % If they overlap or are adjacent, extend the end time of the current period
            current_period(2) = max(current_period(2), sorted_periods(i, 2));
        else
            % If they don't overlap, add the current period to the merged periods
            merged_periods = [merged_periods; current_period];
            % Start a new period to consider for merging
            current_period = sorted_periods(i, :);
        end
    end
    
    % Add the last merged period to the list
    merged_periods = [merged_periods; current_period];
    
    % Now, associate the events with the merged periods
    merged_events = [];
    event_idx = 1;
    
    % Iterate through each merged period
    for i = 1:size(merged_periods, 1)
        % Extract the start and end times of the current merged period
        start_time = merged_periods(i, 1);
        end_time = merged_periods(i, 2);
        
        % For each event, check if it falls within the current merged period
        while event_idx <= length(events) && events(event_idx) <= end_time
            if events(event_idx) >= start_time
                % If the event occurs within the period, add it to the merged events list
                merged_events = [merged_events, events(event_idx)];
            end
            % Move to the next event
            event_idx = event_idx + 1;
        end
    end
end