function perform2WayANOVA(data, condition, subpop)
    % Ensure data, condition, and subpop have the same length
    if length(data) ~= length(condition) || length(data) ~= length(subpop)
        error('Data, condition, and subpopulation vectors must have the same length');
    end
    
    % Create a table with all factors
    tbl = table(data(:), condition(:), subpop(:), 'VariableNames', {'Data', 'Condition', 'Subpopulation'});
    
    % Perform the two-way repeated measures ANOVA
    % `condition` is the repeated measure factor, and `subpop` is the between-subjects factor
    [p, tbl_anova, stats] = anovan(tbl.Data, {tbl.Condition, tbl.Subpopulation}, ...
        'model', 'interaction', 'random', 1, 'varnames', {'Condition', 'Subpopulation'});
    
    % Display the ANOVA results
    disp('2-Way Repeated Measures ANOVA Results:');
    disp(tbl_anova);
    
    % Check if p has more than one value (because it may contain p-values for both main effects and the interaction)
    if numel(p) == 3
        disp(['p-value for Condition main effect: ', num2str(p(1))]);
        disp(['p-value for Subpopulation main effect: ', num2str(p(2))]);
        disp(['p-value for Interaction effect: ', num2str(p(3))]);
    else
        disp(['p-value for interaction effect: ', num2str(p)]);
    end
    
    % Plotting
    figure;
    hold on;
    
    % Define colors for conditions and subpopulations
    colors = [0 0 1; 1 0 0]; % Blue for condition 1, Red for condition 2
    markers = {'o', 's'};  % Circle for subpop 1, Square for subpop 2

    % Loop through unique subpopulations and conditions for plotting
    for sp = unique(subpop)'
        for cond = unique(condition)'
            % Extract data for the current subpopulation and condition
            data_idx = subpop == sp & condition == cond;
            current_data = data(data_idx);
            
            % Add jitter to x-axis to prevent overlap
            jitter = randn(size(current_data)) * 0.1;
            
            % Plot the scatter for the current combination (subpopulation + condition)
            scatter(jitter + cond + (sp-1)*0.2, current_data, ...
                'Marker', markers{sp}, 'MarkerFaceColor', colors(cond, :), 'MarkerEdgeColor', 'k');
        end
    end
    
    % Plot the mean for each condition and subpopulation combination (as a line)
    means = grpstats(data, {condition, subpop});
    for sp = unique(subpop)'
        for cond = unique(condition)'
            % Find the mean for the current combination
            mean_val = means.mean_Data(means.Subpopulation == sp & means.Condition == cond);
            plot(cond + (sp-1)*0.2, mean_val, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8, 'LineWidth', 2);
        end
    end
    
    % Labeling the plot
    xlabel('Condition');
    ylabel('Data');
    title('Scatter Plot of Data with Jitter, Means for Each Subpopulation and Condition');
    xticks([1, 2]);
    xticklabels({'Condition 1', 'Condition 2'});
    legend({'Subpopulation 1', 'Subpopulation 2'}, 'Location', 'best');
    hold off;
end