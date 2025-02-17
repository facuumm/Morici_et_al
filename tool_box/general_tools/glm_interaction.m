function glm_interaction(data, subpopulation, condition)
    % This function fits a GLM for a 2-way design and includes interaction between factors
    % data: the dependent variable (numeric values)
    % subpopulation: a numeric vector where 1 = Aversive, 2 = Reward
    % condition: a numeric vector where 1 = Coordinated, 2 = Uncoordinated

    % Convert subpopulation and condition into categorical variables
    subpop = categorical(subpopulation, [1, 2], {'Aversive', 'Reward'});
    cond = categorical(condition, [1, 2], {'Coordinated', 'Uncoordinated'});

    % Create the interaction term by combining subpopulation and condition
    interaction = strcat(cellstr(subpop), '_', cellstr(cond));

    % Prepare a table for GLM fitting
    dataTable = table(data, subpop, cond, 'VariableNames', {'Data', 'Subpopulation', 'Condition'});

    % Fit a GLM with a log link function (assuming the data is count-based)
    glm_model = fitglm(dataTable, 'Data ~ Subpopulation * Condition', 'Distribution', 'poisson', 'Link', 'log');

    % Display the model results
    disp(glm_model);

    % Check if the interaction term is significant
    % The p-value for the interaction term will help you assess the interaction effect
    p_values = glm_model.Coefficients.pValue;
    fprintf('p-value for Interaction (Subpopulation*Condition): %.4f\n', p_values(4));  % Assuming the interaction is in the 4th row
end