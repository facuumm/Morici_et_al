function perform2WayANOVANonPaired(data, condition, subpop)
    % Ensure data, condition, and subpop have the same length
    if length(data) ~= length(condition) || length(data) ~= length(subpop) || length(condition) ~= length(subpop)
        error('Data, condition, and subpopulation vectors must have the same length');
    end

    % Display lengths of inputs to help debug
    disp(['Length of data: ', num2str(length(data))]);
    disp(['Length of condition: ', num2str(length(condition))]);
    disp(['Length of subpop: ', num2str(length(subpop))]);

    % Create a table with all factors
    tbl = table(data(:), condition(:), subpop(:), 'VariableNames', {'Data', 'Condition', 'Subpopulation'});

    % Map numeric values to string labels for better interpretation
    tbl.Condition = categorical(tbl.Condition, [1 2], {'Uncoordinated', 'Coordinated'});
    tbl.Subpopulation = categorical(tbl.Subpopulation, [1 2], {'Reward', 'Aversive'});

    % Perform the two-way non-paired ANOVA
    [p, tbl_anova, stats] = anovan(tbl.Data, {tbl.Condition, tbl.Subpopulation}, ...
        'model', 'interaction', 'varnames', {'Condition', 'Subpopulation'}, 'display', 'off');

    % Display the ANOVA results
    disp('2-Way Non-Paired ANOVA Results:');
    disp(tbl_anova);
    if numel(p) == 3
        disp(['p-value for Condition main effect: ', num2str(p(1))]);
        disp(['p-value for Subpopulation main effect: ', num2str(p(2))]);
        disp(['p-value for Interaction effect: ', num2str(p(3))]);
    else
        disp(['p-value for interaction effect: ', num2str(p)]);
    end

    % Perform non-paired comparisons between Aversive vs Reward within each condition (Coordinated and Uncoordinated)
    disp('Performing non-paired comparisons between Aversive vs Reward within each condition...');
    
    % Define the two comparisons you want:
    comparisons = { ...
        {'Coordinated', 'Aversive', 'Reward'}, ...
        {'Uncoordinated', 'Aversive', 'Reward'} ...
    };
    
    % Set alpha level and adjust for Bonferroni (for 2 comparisons)
    alpha = 0.05;
    adjusted_alpha = alpha / 2;  % Bonferroni adjustment for 2 comparisons (0.05 / 2 = 0.025)
    
    disp(['Bonferroni adjusted alpha: ', num2str(adjusted_alpha)]);
    
    for i = 1:length(comparisons)
        condition_name = comparisons{i}{1};
        subpop1 = comparisons{i}{2};
        subpop2 = comparisons{i}{3};
        
        % Filter data for the current condition and subpopulations
        condition_data = tbl(tbl.Condition == condition_name, :);
        subpop1_data = condition_data.Data(condition_data.Subpopulation == subpop1);
        subpop2_data = condition_data.Data(condition_data.Subpopulation == subpop2);
        
        % If both subpopulations have data, perform non-paired comparison (e.g., t-test)
        if ~isempty(subpop1_data) && ~isempty(subpop2_data)
            % Perform two-sample t-test (non-paired) for current condition
            [~, p_value] = ttest2(subpop1_data, subpop2_data, 'Tail','right');
            
            % Display results
            disp(['Comparison: ', char(subpop1), ' vs ', char(subpop2), ' (' char(condition_name) ')']);
            disp(['Unadjusted p-value: ', num2str(p_value)]);
            
            % Check if the p-value is below the Bonferroni adjusted threshold
            if p_value < adjusted_alpha
                disp('Result is statistically significant after Bonferroni correction.');
            else
                disp('Result is NOT statistically significant after Bonferroni correction.');
            end
        else
            disp(['Missing data for ', char(condition_name), ' - ', char(subpop1), ' vs ', char(subpop2)]);
        end
    end
end