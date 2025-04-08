function plot_mean_with_jitter(vector1, vector2, varargin)
% PLOT_MEAN_WITH_JITTER Plot means and individual points with jitter for two vectors
%
%   plot_mean_with_jitter(vector1, vector2) plots the means as distinct markers
%   and individual points with horizontal jitter
%
%   Optional parameters:
%       'GroupNames'    : Cell array of group names (default: {'Group1','Group2'})
%       'PointColor'    : Color for individual points (default: [0.4 0.6 0.8])
%       'MeanColor'     : Color for mean markers (default: [0.9 0.2 0.2])
%       'PointSize'     : Size of individual points (default: 40)
%       'MeanSize'      : Size of mean markers (default: 100)
%       'JitterWidth'   : Width of jitter (default: 0.4)
%       'ShowErrorBars' : Show error bars (default: true)

% Parse optional parameters
p = inputParser;
addParameter(p, 'GroupNames', {'Group1','Group2'}, @iscell);
addParameter(p, 'PointColor', [0.4 0.6 0.8], @(x) isnumeric(x) && numel(x)==3);
addParameter(p, 'MeanColor', [0.9 0.2 0.2], @(x) isnumeric(x) && numel(x)==3);
addParameter(p, 'PointSize', 40, @isnumeric);
addParameter(p, 'MeanSize', 100, @isnumeric);
addParameter(p, 'JitterWidth', 0.4, @isnumeric);
addParameter(p, 'ShowErrorBars', true, @islogical);
parse(p, varargin{:});

% Prepare data
data = {vector1(:), vector2(:)};
n_groups = length(data);
x_pos = 1:n_groups;
means = cellfun(@mean, data);
stds = cellfun(@std, data);

% Create figure
figure;
hold on;

% Plot individual points with jitter (no outline)
for i = 1:n_groups
    % Add jitter to x-position
    jitter = p.Results.JitterWidth * (rand(size(data{i})) - 0.5);
    x = i + jitter;
    
    % Plot points without outline
    scatter(x, data{i}, p.Results.PointSize, p.Results.PointColor, ...
           'filled', 'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'none');
end

% Plot mean markers (larger and different color)
scatter(x_pos, means, p.Results.MeanSize, p.Results.MeanColor, ...
       'filled', 'Marker', 'd', 'MarkerEdgeColor', 'k', 'LineWidth', 1);

% Plot error bars if requested
if p.Results.ShowErrorBars
    errorbar(x_pos, means, stds, 'k.', 'LineWidth', 1.5, 'CapSize', 15);
end

% Customize axes
set(gca, 'XTick', x_pos, 'XTickLabel', p.Results.GroupNames, 'FontSize', 12);
ylabel('Value', 'FontSize', 12);
xlim([0.5 n_groups+0.5]);
box off;

% Add light grid
grid on;
set(gca, 'GridAlpha', 0.2);

hold off;
end