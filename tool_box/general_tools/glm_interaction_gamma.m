function glm_interaction_gamma(data, subpopulation, condition)
    % This function fits a GLM for a 2-way design with Gamma distribution and includes interaction between factors
    % data: the dependent variable (continuous and can be any range of values, but should not contain non-positive values)
    % subpopulation: a numeric vector where 1 = Aversive, 2 = Reward
    % condition: a numeric vector where 1 = Coordinated, 2 = Uncoordinated

    % Find the minimum value of the data
    min_data = min(data);
    
    % Shift the data to make sure all values are positive
    shift_constant = abs(min_data) + 1e-6;  % Adding a small value to avoid zero
    data_transformed = data + shift_constant;  % Shift all values to be positive

    % Debugging: Display the transformed data
    disp('Transformed Data:');
    disp(data_transformed);

    % Convert subpopulation and condition into categorical variables
    subpop = categorical(subpopulation, [1, 2], {'Reward', 'Aversive'});
    cond = categorical(condition, [1, 2], {'Uncoordinated', 'Coordinated'});

    % Prepare a table for GLM fitting
    dataTable = table(data_transformed, subpop, cond, 'VariableNames', {'Data', 'Subpopulation', 'Condition'});

    % Fit a GLM with a Gamma distribution and log link function
    glm_model = fitglm(dataTable, 'Data ~ Subpopulation * Condition', 'Distribution', 'gamma', 'Link', 'log');

    % Display the model results
    disp(glm_model);

    % Check if the interaction term is significant
    p_values = glm_model.Coefficients.pValue;
    fprintf('p-value for Interaction (Subpopulation*Condition): %.4f\n', p_values(4));  % Assuming the interaction is in the 4th row
end