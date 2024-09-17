function downsampled_t = downsampleTimeVector(t, original_fps, desired_interval)
    % downsampleTimeVector: Downsample a time vector to match a desired time interval
    %
    % Inputs:
    % t - The original time vector
    % original_fps - The original frame rate (frames per second)
    % desired_interval - The desired time interval in seconds
    %
    % Output:
    % downsampled_t - The downsampled time vector
    
    % Check inputs
    if isempty(t)
        error('Input time vector cannot be empty.');
    end
    if original_fps <= 0
        error('Original FPS must be a positive number.');
    end
    if desired_interval <= 0
        error('Desired interval must be a positive number.');
    end
    
    % Calculate the number of points in the original and downsampled vectors
    original_time_span = t(end) - t(1);
    num_points_original = length(t);
    
    % Determine the number of points for the desired interval
    num_points_desired = round(original_time_span / desired_interval) + 1;
    
    % Create the new downsampled time vector using linspace
    downsampled_t = linspace(t(1), t(end), num_points_desired);
    
    % Interpolate the original time vector to match the downsampled time points
    downsampled_t = interp1(t, t, downsampled_t, 'linear');
    
    % Output the downsampled time vector
end