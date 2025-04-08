function plot_adaptive_comparison(group1_data, group2_data)
    % PLOT_ADAPTIVE_COMPARISON Handles both single-group and two-group visualization
    %
    % Inputs:
    %   group1_data: NxC matrix (required)
    %   group2_data: MxC matrix (optional - leave empty [] for single group)
    %
    % Features:
    % - Works with 1 or 2 groups
    % - Automatically adjusts visualization
    % - Performs appropriate statistical tests
    % - Clean, publication-ready output

    %% Input handling
    if nargin < 2 || isempty(group2_data)
        single_group_mode = true;
        n_conditions = size(group1_data,2);
        conditions = cellstr(char(64+(1:n_conditions))'); % A, B, C, ...
    else
        single_group_mode = false;
        if size(group1_data,2) ~= size(group2_data,2)
            error('Groups must have same number of conditions');
        end
        n_conditions = size(group1_data,2);
        conditions = cellstr(char(64+(1:n_conditions))');
    end

    %% Visualization setup
    fig = figure('Position', [100 100 400+150*n_conditions, 600], 'Color', 'w');
    ax = axes('Parent', fig);
    hold(ax, 'on');
    
    % Visual parameters
    group_colors = [0.2 0.6 0.8;  % Group 1 (blue)
                    0.8 0.4 0.2];  % Group 2 (orange)
    marker_size = 80;
    jitter_amount = 0.15;
    box_width = 0.4;
    
    if single_group_mode
        %% Single group visualization
        n1 = size(group1_data,1);
        x_positions = 1:n_conditions;
        x_jitter = repmat(x_positions,n1,1) + (rand(n1,n_conditions)*jitter_amount*2 - jitter_amount);
        
        % Plot individual points
        scatter(ax, x_jitter(:), group1_data(:), marker_size, ...
               'MarkerFaceColor', group_colors(1,:), ...
               'MarkerEdgeColor', 'k', ...
               'MarkerFaceAlpha', 0.7);
        
        % Add boxplots
        for c = 1:n_conditions
            q = quantile(group1_data(:,c), [0.25 0.5 0.75]);
            iqr = q(3) - q(1);
            whisker_low = max(min(group1_data(:,c)), q(1) - 1.5*iqr);
            whisker_high = min(max(group1_data(:,c)), q(3) + 1.5*iqr);
            
            % Box
            patch(ax, [c-box_width/2 c+box_width/2 c+box_width/2 c-box_width/2], ...
                  [q(1) q(1) q(3) q(3)], group_colors(1,:), ...
                  'FaceAlpha', 0.2, 'EdgeColor', group_colors(1,:), 'LineWidth', 1.5);
            
            % Median line
            plot(ax, [c-box_width/2 c+box_width/2], [q(2) q(2)], ...
                 'Color', group_colors(1,:), 'LineWidth', 2);
            
            % Whiskers
            plot(ax, [c c], [q(1) whisker_low], '--', 'Color', group_colors(1,:), 'LineWidth', 1);
            plot(ax, [c c], [q(3) whisker_high], '--', 'Color', group_colors(1,:), 'LineWidth', 1);
            plot(ax, [c-box_width/4 c+box_width/4], [whisker_low whisker_low], '-', 'Color', group_colors(1,:), 'LineWidth', 1);
            plot(ax, [c-box_width/4 c+box_width/4], [whisker_high whisker_high], '-', 'Color', group_colors(1,:), 'LineWidth', 1);
        end
        
        % Formatting for single group
        set(ax, 'XTick', x_positions, 'XTickLabel', conditions, ...
                'FontSize', 12, 'LineWidth', 1.5, 'Box', 'off');
        title(ax, sprintf('Single Group (%d Conditions)', n_conditions), 'FontSize', 14);
        
    else
        %% Two-group visualization
        n1 = size(group1_data,1);
        n2 = size(group2_data,1);
        group1_x = 1:n_conditions;
        group2_x = (1:n_conditions) + n_conditions + 0.5;
        
        % Plot Group 1
        x1_jitter = repmat(group1_x,n1,1) + (rand(n1,n_conditions)*jitter_amount*2 - jitter_amount);
        scatter(ax, x1_jitter(:), group1_data(:), marker_size, ...
               'MarkerFaceColor', group_colors(1,:), ...
               'MarkerEdgeColor', 'k', ...
               'MarkerFaceAlpha', 0.7);
        
        % Plot Group 2
        x2_jitter = repmat(group2_x,n2,1) + (rand(n2,n_conditions)*jitter_amount*2 - jitter_amount);
        scatter(ax, x2_jitter(:), group2_data(:), marker_size, ...
               'MarkerFaceColor', group_colors(2,:), ...
               'MarkerEdgeColor', 'k', ...
               'MarkerFaceAlpha', 0.7);
        
        % Add boxplots for both groups
        for g = 1:2
            curr_data = {group1_data, group2_data}{g};
            curr_x = {group1_x, group2_x}{g};
            
            for c = 1:n_conditions
                q = quantile(curr_data(:,c), [0.25 0.5 0.75]);
                iqr = q(3) - q(1);
                whisker_low = max(min(curr_data(:,c)), q(1) - 1.5*iqr);
                whisker_high = min(max(curr_data(:,c)), q(3) + 1.5*iqr);
                
                % Box
                patch(ax, [curr_x(c)-box_width/2 curr_x(c)+box_width/2 curr_x(c)+box_width/2 curr_x(c)-box_width/2], ...
                      [q(1) q(1) q(3) q(3)], group_colors(g,:), ...
                      'FaceAlpha', 0.2, 'EdgeColor', group_colors(g,:), 'LineWidth', 1.5);
                
                % Median line
                plot(ax, [curr_x(c)-box_width/2 curr_x(c)+box_width/2], [q(2) q(2)], ...
                     'Color', group_colors(g,:), 'LineWidth', 2);
                
                % Whiskers
                plot(ax, [curr_x(c) curr_x(c)], [q(1) whisker_low], '--', 'Color', group_colors(g,:), 'LineWidth', 1);
                plot(ax, [curr_x(c) curr_x(c)], [q(3) whisker_high], '--', 'Color', group_colors(g,:), 'LineWidth', 1);
                plot(ax, [curr_x(c)-box_width/4 curr_x(c)+box_width/4], [whisker_low whisker_low], '-', 'Color', group_colors(g,:), 'LineWidth', 1);
                plot(ax, [curr_x(c)-box_width/4 curr_x(c)+box_width/4], [whisker_high whisker_high], '-', 'Color', group_colors(g,:), 'LineWidth', 1);
            end
        end
        
        % Formatting for two groups
        set(ax, 'XTick', [mean(group1_x) mean(group2_x)], ...
                'XTickLabel', {'Group 1', 'Group 2'}, ...
                'FontSize', 12, 'LineWidth', 1.5, 'Box', 'off');
        title(ax, sprintf('Two-Group Comparison (%d Conditions)', n_conditions), 'FontSize', 14);
        
        % Add condition labels
        text(ax, mean(group1_x), min(ylim(ax))-0.1*diff(ylim(ax)), ...
             ['Conditions: ' strjoin(conditions)], ...
             'HorizontalAlignment', 'center', 'FontSize', 11);
    end
    
    %% Common formatting
    ylabel(ax, 'Performance Metric', 'FontSize', 13);
    grid(ax, 'on');
    
    % Adjust axes
    if single_group_mode
        all_values = group1_data(:);
        xlim(ax, [0.5 n_conditions+0.5]);
    else
        all_values = [group1_data(:); group2_data(:)];
        xlim(ax, [0.5 max(group2_x)+0.5]);
    end
    y_range = range(all_values);
    ylim(ax, [min(all_values)-0.1*y_range, max(all_values)+0.15*y_range]);
    
    %% Add reference line at zero
    plot(ax, xlim(ax), [0 0], '--k', 'LineWidth', 1);
    
    %% Perform and display appropriate statistics
    if single_group_mode
        % Friedman test for single group with multiple conditions
        if n_conditions > 1
            p = friedman(group1_data, 1, 'off');
            fprintf('\nFriedman test (non-parametric repeated measures ANOVA):\n');
            fprintf('p = %.4f for %d conditions\n', p, n_conditions);
            
            if p < 0.05 && n_conditions > 2
                fprintf('Significant difference found - consider post-hoc tests\n');
            end
        end
    else
        % Two-way mixed ANOVA for two groups
        if n_conditions > 1
            try
                y = [group1_data(:); group2_data(:)];
                subject = [repmat((1:n1)',n_conditions,1); repmat((1:n2)'+n1,n_conditions,1)];
                group = [repmat({'Group1'},n1*n_conditions,1); repmat({'Group2'},n2*n_conditions,1)];
                condition = repmat([repmat(1:n_conditions,n1,1)'; repmat(1:n_conditions,n2,1)'],1);
                
                [p,tbl,stats] = anovan(y, {group, condition, subject}, ...
                    'model', 'interaction', 'random', 3, ...
                    'varnames', {'Group','Condition','Subject'}, 'display', 'off');
                
                fprintf('\nMixed ANOVA Results:\n');
                fprintf('Group effect: p = %.4f\n', p(1));
                fprintf('Condition effect: p = %.4f\n', p(2));
                fprintf('Group×Condition interaction: p = %.4f\n', p(3));
            catch ME
                fprintf('\nANOVA failed: %s\n', ME.message);
            end
        end
    end
end