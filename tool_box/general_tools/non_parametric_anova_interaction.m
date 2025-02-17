function non_parametric_anova_interaction(data, subpopulation, condition)
    % This function performs a non-parametric unpaired 2-way ANOVA
    % data: a column vector with observations (numeric values)
    % subpopulation: a numeric vector where 1 = Aversive, 2 = Reward
    % condition: a numeric vector where 1 = Coordinated, 2 = Uncoordinated
    
    % Convert subpopulation and condition to categorical variables
    subpop = categorical(subpopulation, [1, 2], {'Aversive', 'Reward'});
    cond = categorical(condition, [1, 2], {'Coordinated', 'Uncoordinated'});

    % Create a new grouping variable that combines both subpopulation and condition
    group = strcat(cellstr(subpop), '_', cellstr(cond));

    % Perform Kruskal-Wallis test
    p = kruskalwallis(data, group, 'off'); % 'off' suppresses the boxplot
    
    % Display the Kruskal-Wallis Test p-value
    fprintf('Kruskal-Wallis Test p-value (for main effects and interaction): %.4f\n', p);
    
    % If you want to perform pairwise comparisons (post-hoc test):
    [p_values, ~, stats] = kruskalwallis(data, group, 'off');
    
    % Pairwise comparisons
    multcompare(stats);
    
    % To assess interaction:
    % Kruskal-Wallis tests for differences between the four groups.
    fprintf('Kruskal-Wallis Test results:\n');
    disp(stats);
end