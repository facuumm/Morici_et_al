function plotPrePost(data)
% plotPrePost(data)
% data: N x 2 matrix, col1 = pre-sleep, col2 = post-sleep
%
% Plots jittered scatter of paired values, connects them with lines,
% and overlays group means.

    if size(data,2) ~= 2
        error('Input must be an N x 2 matrix');
    end

    N = size(data,1);
    jitterAmount = 0.1;

    % X positions with jitter
    x_pre  = 1 + (rand(N,1)-0.5)*jitterAmount;
    x_post = 2 + (rand(N,1)-0.5)*jitterAmount;

    % Scatter plot
%     figure; 
    hold on;
    scatter(x_pre, data(:,1), 60, 'filled')
    scatter(x_post, data(:,2), 60, 'filled')

    % Connect paired points
    for i = 1:N
        plot([x_pre(i), x_post(i)], [data(i,1), data(i,2)], 'k-', 'LineWidth', 0.8)
    end

    % Plot means
    mean_pre  = mean(data(:,1));
    mean_post = mean(data(:,2));
    plot([1,2], [mean_pre, mean_post], 'ro-', 'MarkerSize',8, 'LineWidth',2)

    % Axis / labels
    xlim([0.5 2.5])
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Pre-sleep','Post-sleep'})
    ylabel('Value')
    box off
end