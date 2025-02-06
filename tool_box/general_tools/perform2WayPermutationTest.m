function perform2WayPermutationTest(data, condition, subpop)
    % Check that data, condition, and subpop are of the same length
    if length(data) ~= length(condition) || length(data) ~= length(subpop)
        error('Data, condition, and subpopulation vectors must have the same length');
    end
    
    % Get unique conditions and subpopulations
    conditions = unique(condition);
    subpopulations = unique(subpop);
    
    % Number of permutations
    num_permutations = 10000;
    
    % Calculate observed statistic (e.g., interaction effect)
    observed_stat = calcInteractionStatistic(data, condition, subpop, conditions, subpopulations);
    
    % Permutation loop
    permuted_stats = zeros(num_permutations, 1);
    for i = 1:num_permutations
        % Randomly permute the condition labels
        permuted_condition = condition(randperm(length(condition)));
        
        % Calculate the statistic for the permuted data
        permuted_stats(i) = calcInteractionStatistic(data, permuted_condition, subpop, conditions, subpopulations);
    end
    
    % Calculate p-value based on permutation results
    p_value = mean(abs(permuted_stats) >= abs(observed_stat));
    
    % Display results
    disp(['Observed interaction statistic: ', num2str(observed_stat)]);
    disp(['Permutation-based p-value: ', num2str(p_value)]);
end

function stat = calcInteractionStatistic(data, condition, subpop, conditions, subpopulations)
    % Initialize the statistic (e.g., mean difference or other interaction measure)
    
    % Initialize an array to store the means for each condition/subpopulation pair
    means = NaN(length(subpopulations), length(conditions));
    
    % Loop through each subpopulation and condition
    for sp_idx = 1:length(subpopulations)
        for cond_idx = 1:length(conditions)
            % Get the indices for the current subpopulation and condition
            idx = (subpop == subpopulations(sp_idx)) & (condition == conditions(cond_idx));
            
            % Calculate the mean for this subpopulation/condition pair
            means(sp_idx, cond_idx) = mean(data(idx), 'omitnan');
        end
    end
    
    % Calculate the interaction term (this could be, for example, a mean difference or a specific comparison)
    % Here, we'll calculate the difference between the means for different conditions within each subpopulation
    stat = abs(means(1,1) - means(2,2));  % Example: interaction between subpop 1, condition 1 and subpop 2, condition 2
end