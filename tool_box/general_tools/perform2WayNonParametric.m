function perform2WayNonParametric(data, condition, subpop)
% Ensure data, condition, and subpop have the same length
if length(data) ~= length(condition) || length(data) ~= length(subpop)
    error('Data, condition, and subpopulation vectors must have the same length');
end

% Create a table with all factors
tbl = table(data(:), condition(:), subpop(:), 'VariableNames', {'Data', 'Condition', 'Subpopulation'});

% Perform the Kruskal-Wallis test for main effects and interaction
disp('Performing Kruskal-Wallis Test...');
% Test for main effect of Condition
[p_condition, tbl_condition] = kruskalwallis(tbl.Data, tbl.Condition, 'off');
% Test for main effect of Subpopulation
[p_subpop, tbl_subpop] = kruskalwallis(tbl.Data, tbl.Subpopulation, 'off');

disp(['p-value for Condition main effect: ', num2str(p_condition)]);
disp(['p-value for Subpopulation main effect: ', num2str(p_subpop)]);

% Test for interaction between Condition and Subpopulation
disp('Testing for interaction between Condition and Subpopulation...');
% For non-parametric interaction, apply Kruskal-Wallis on data grouped by both Condition and Subpopulation
grouped_data = grp2cell(tbl.Data, tbl.Condition, tbl.Subpopulation);
% Perform the Kruskal-Wallis test on these groups
[p_interaction, tbl_interaction] = kruskalwallis(grouped_data, [], 'off');
disp(['p-value for Interaction effect: ', num2str(p_interaction)]);

% If the interaction term is significant, perform pairwise comparisons
if p_interaction < 0.05
    disp('Performing multiple pairwise comparisons for interaction...');
    % For pairwise comparisons between conditions within each subpopulation
    unique_subpopulations = unique(tbl.Subpopulation);
    
    for i = 1:length(unique_subpopulations)
        % Filter data for the current subpopulation
        subpop_condition_data = tbl(tbl.Subpopulation == unique_subpopulations(i), :);
        
        % Perform pairwise comparisons for this subpopulation using Dunn's test or Mann-Whitney U test
        disp(['Performing pairwise comparisons for Subpopulation ', num2str(unique_subpopulations(i)), '...']);
        % Group the data by condition
        group_data = grp2cell(subpop_condition_data.Data, subpop_condition_data.Condition);
        % Perform Dunn's test (or Mann-Whitney U test for pairwise comparisons)
        p_pairwise = dunnTest(group_data);  % Dunn's test for pairwise comparisons
        disp(['Pairwise Comparison Results for Subpopulation ', num2str(unique_subpopulations(i)), ':']);
        disp(p_pairwise);
    end
else
    disp('No significant interaction effect found.');
end

end