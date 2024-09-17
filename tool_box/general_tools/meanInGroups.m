function means = meanInGroups(v, n)
    % Calculate the mean of a vector in groups of n points,
    % padding with NaN if necessary.
    %
    % Inputs:
    % v - The column vector
    % n - Number of points in each group
    %
    % Outputs:
    % means - A row vector of the mean of each group

    % Calculate the length of the vector
    len = length(v);

    % Check if padding is needed
    if mod(len, n) ~= 0
        % Calculate how many NaN values are needed to make the length a multiple of n
        pad_len = ceil(len / n) * n - len;
        
        % Pad the vector with NaN values
        padded_v = [v; NaN(pad_len, 1)];
    else
        % No padding needed, use the original vector
        padded_v = v;
    end

    % Reshape the vector into groups of n points
    reshaped_v = reshape(padded_v, n, []);

    % Calculate the mean for each group, ignoring NaN values
    means = mean(reshaped_v, 1, 'omitnan');
end