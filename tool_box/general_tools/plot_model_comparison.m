function plot_model_comparison(group1_data, group2_data)
    % PLOT_MODEL_COMPARISON Visualizes paired data from two groups with connecting lines
    %
    % Inputs:
    %   group1_data: Nx2 matrix [conditionA, conditionB] for first group
    %   group2_data: Mx2 matrix [conditionA, conditionB] for second group
    %
    % Features:
    % - Paired visualization within each group
    % - Distinct colors/markers for each group
    % - Connecting lines between paired measurements
    % - Boxplots showing distribution
    % - Statistical comparison

    % Validate inputs
    if size(group1_data,2) ~= 2 || size(group2_data,2) ~= 2
        error('Input matrices must have 2 columns each');
    end
    
    % Create figure
    figure('Position', [100 100 900 600], 'Color', 'w');
    hold on;
    
    % Set visual parameters
    group_colors = [0.2 0.6 0.8;   % Group 1 color (blue)
                    0.8 0.4 0.2];   % Group 2 color (orange)
    marker_types = {'o', 's'};      % Circle and square markers
    line_colors = [0.5 0.5 0.5 0.3; % Group 1 line color (gray)
                   0.7 0.5 0.3 0.3]; % Group 2 line color (brown)
    marker_size = 80;
    
    % Plot Group 1 data
    n1 = size(group1_data,1);
    x1 = [ones(n1,1), 2*ones(n1,1)]; % X positions for conditions
    x1_jitter = x1 + (rand(size(x1))*0.2 - 0.1; % Add jitter
    
    % Plot connecting lines for group 1
    for i = 1:n1
        plot(x1_jitter(i,:), group1_data(i,:), ...
            'Color', line_colors(1,:), 'LineWidth', 1.2);
    end
    
    % Plot points for group 1
    scatter(x1_jitter(:,1), group1_data(:,1), marker_size, ...
           'Marker', marker_types{1}, ...
           'MarkerFaceColor', group_colors(1,:), ...
           'MarkerEdgeColor', 'k', ...
           'MarkerFaceAlpha', 0.7);
       
    scatter(x1_jitter(:,2), group1_data(:,2), marker_size, ...
           'Marker', marker_types{1}, ...
           'MarkerFaceColor', group_colors(1,:), ...
           'MarkerEdgeColor', 'k', ...
           'MarkerFaceAlpha', 0.7);
    
    % Plot Group 2 data (offset on x-axis)
    n2 = size(group2_data,1);
    x2 = [3*ones(n2,1), 4*ones(n2,1)]; % X positions for conditions
    x2_jitter = x2 + (rand(size(x2))*0.2 - 0.1; % Add jitter
    
    % Plot connecting lines for group 2
    for i = 1:n2
        plot(x2_jitter(i,:), group2_data(i,:), ...
            'Color', line_colors(2,:), 'LineWidth', 1.2);
    end
    
    % Plot points for group 2
    scatter(x2_jitter(:,1), group2_data(:,1), marker_size, ...
           'Marker', marker_types{2}, ...
           'MarkerFaceColor', group_colors(2,:), ...
           'MarkerEdgeColor', 'k', ...
           'MarkerFaceAlpha', 0.7);
       
    scatter(x2_jitter(:,2), group2_data(:,2), marker_size, ...
           'Marker', marker_types{2}, ...
           'MarkerFaceColor', group_colors(2,:), ...
           'MarkerEdgeColor', 'k', ...
           'MarkerFaceAlpha', 0.7);
    
    % Add boxplots (transparent)
    boxplot([group1_data(:,1); group1_data(:,2); group2_data(:,1); group2_data(:,2)], ...
           [ones(n1,1); 2*ones(n1,1); 3*ones(n2,1); 4*ones(n2,1)], ...
           'Colors', repmat(group_colors,2,1), ...
           'Widths', 0.4, ...
           'Symbol', '', ...
           'BoxStyle', 'outline');
    
    % Make boxes transparent
    h = findobj(gca, 'Tag', 'Box');
    for j = 1:length(h)
        patch(get(h(j), 'XData'), get(h(j), 'YData'), group_colors(ceil(j/2),:), ...
             'FaceAlpha', 0.2, 'EdgeColor', group_colors(ceil(j/2),:), 'LineWidth', 2);
    end
    
    % Add reference line
    yl = ylim;
    plot(xlim, [0 0], '--k', 'LineWidth', 1);
    
    % Calculate and display statistical significance
    [pval1, ~] = signrank(group1_data(:,1), group1_data(:,2));
    [pval2, ~] = signrank(group2_data(:,1), group2_data(:,2));
    
    % Group 1 significance
    if pval1 < 0.05
        sig_y1 = max(group1_data(:)) + 0.05*diff(ylim);
        line([1 2], [sig_y1 sig_y1], 'Color', 'k', 'LineWidth', 1.5);
        text(1.5, sig_y1 + 0.02*diff(ylim), get_pstar(pval1), ...
             'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
    end
    
    % Group 2 significance
    if pval2 < 0.05
        sig_y2 = max(group2_data(:)) + 0.05*diff(ylim);
        line([3 4], [sig_y2 sig_y2], 'Color', 'k', 'LineWidth', 1.5);
        text(3.5, sig_y2 + 0.02*diff(ylim), get_pstar(pval2), ...
             'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
    end
    
    % Formatting
    set(gca, 'XTick', [1.5 3.5], 'XTickLabel', {'Group 1', 'Group 2'}, ...
             'FontSize', 12, 'LineWidth', 1.5, 'Box', 'off');
    ylabel('Spearman''s ρ', 'FontSize', 14);
    title('Model Performance Comparison', 'FontSize', 16, 'FontWeight', 'bold');
    
    % Add legend
    legend({'Group 1 pairs', 'Group 1 condition A', 'Group 1 condition B', ...
            'Group 2 pairs', 'Group 2 condition A', 'Group 2 condition B'}, ...
           'Location', 'bestoutside', 'FontSize', 10);
    
    % Adjust axes
    all_data = [group1_data(:); group2_data(:)];
    y_range = range(all_data);
    ylim([min(all_data)-0.1*y_range, max(all_data)+0.2*y_range]);
    xlim([0.5 4.5]);
end

function pstar = get_pstar(pval)
    % Helper function for significance stars
    if pval < 0.001
        pstar = '***';
    elseif pval < 0.01
        pstar = '**';
    elseif pval < 0.05
        pstar = '*';
    else
        pstar = 'n.s.';
    end
end