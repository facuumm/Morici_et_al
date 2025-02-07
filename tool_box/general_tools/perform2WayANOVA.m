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

% Perform multiple comparisons for each condition
disp('Performing multiple comparisons across conditions...');
% Specify the conditions for pairwise comparisons (i.e., each condition against each other)
c = multcompare(stats, 'Dimension', 1); % Compare across 'Condition' dimension

% Display multiple comparison results
disp('Multiple Comparison Results for Condition:');
disp(c);

% Perform separate multiple comparisons within each Subpopulation group
unique_subpopulations = unique(tbl.Subpopulation);

for i = 1:length(unique_subpopulations)
    % Filter data for the current subpopulation
    subpop_condition_data = tbl(tbl.Subpopulation == unique_subpopulations(i), :);
    
    % Perform the ANOVA for the current subpopulation
    [p_subpop, tbl_anova_subpop, stats_subpop] = anovan(subpop_condition_data.Data, {subpop_condition_data.Condition}, ...
        'model', 'full', 'random', 1, 'varnames', {'Condition'});
    
    % Check if ANOVA is significant for this subpopulation before running multcompare
    if p_subpop(1) < 0.05
        disp(['Performing multiple comparisons for Subpopulation ', num2str(unique_subpopulations(i)), '...']);
        % Perform multiple comparisons for each condition within the current subpopulation
        c_subpop = multcompare(stats_subpop, 'Dimension', 1); % Compare across 'Condition' dimension
        disp(['Multiple Comparison Results for Subpopulation ', num2str(unique_subpopulations(i)), ':']);
        disp(c_subpop);
    else
        disp(['No significant differences in Subpopulation ', num2str(unique_subpopulations(i)), '.']);
    end
end
end