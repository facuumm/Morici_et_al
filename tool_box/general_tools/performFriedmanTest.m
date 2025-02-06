function performFriedmanTest(data, condition, subpop)
    % Ensure that data, condition, and subpop are the same length
    if length(data) ~= length(condition) || length(data) ~= length(subpop)
        error('Data, condition, and subpopulation vectors must have the same length');
    end
    
    % Get unique subpopulations and conditions
    subpopulations = unique(subpop);
    conditions = unique(condition);
    
    % Initialize a structure to hold testData for each subpopulation
    testDataBySubpop = struct();
    
    % Loop through each subpopulation
    for sp_idx = 1:length(subpopulations)
        sp = subpopulations(sp_idx);  % Current subpopulation
        
        % Get the indices for the current subpopulation
        subpop_idx = (subpop == sp);
        
        % Extract data and conditions for the current subpopulation
        subpop_data = data(subpop_idx);
        subpop_conditions = condition(subpop_idx);
        
        % Initialize the testData matrix for this subpopulation
        num_subjects = length(unique(subpop(subpop == sp)));  % Unique subjects in this subpopulation
        num_conditions = length(conditions);  % Number of conditions (coordinated/uncoordinated)
        testData = NaN(num_subjects, num_conditions);  % Initialize with NaN values
        
        % Loop through each subject
        subject_idx = 1;
        for subj = unique(subpop(subpop_idx))  % Loop over unique subjects within this subpopulation
            % Get the indices for the current subject in the current subpopulation
            subject_data_idx = (subpop == sp) & (subpop == subj);
            
            % Extract data for the current subject
            subject_data = subpop_data(subject_data_idx);
            subject_conditions = subpop_conditions(subject_data_idx);
            
            % Debugging print statements
            disp(['Processing Subject: ', num2str(subj)]);
            disp(['Indices for subject ', num2str(subj), ': ', num2str(find(subject_data_idx))]);
            
            % Loop through each condition and assign the data for the subject
            for cond_idx = 1:num_conditions
                % Get the condition indices for the current condition
                condition_idx = subject_conditions == conditions(cond_idx);
                
                % Debugging print statements
                disp(['Processing Condition: ', num2str(conditions(cond_idx))]);
                disp(['Condition indices: ', num2str(find(condition_idx))]);
                
                % Ensure there are valid condition indices
                if sum(condition_idx) > 0
                    testData(subject_idx, cond_idx) = mean(subject_data(condition_idx));
                else
                    disp(['Warning: No data for subject ', num2str(subj), ' under condition ', num2str(conditions(cond_idx))]);
                end
            end
            subject_idx = subject_idx + 1;  % Move to the next subject
        end
        
        % Store the testData for this subpopulation in the structure
        testDataBySubpop.(sprintf('Subpop_%d', sp)) = testData;
    end
    
    % Perform the Friedman test for each subpopulation
    for sp_idx = 1:length(subpopulations)
        sp = subpopulations(sp_idx);
        
        % Get the test data for this subpopulation
        testData = testDataBySubpop.(sprintf('Subpop_%d', sp));
        
        % Display the testData matrix for each subpopulation
        disp(['Test data for subpopulation ', num2str(sp)]);
        disp(testData);
        
        % Ensure the testData matrix has at least 2 rows and columns
        if size(testData, 1) < 2 || size(testData, 2) < 2
            warning('Skipping Friedman test for subpopulation %d due to insufficient data', sp);
            continue;
        end
        
        % Perform the Friedman test for this subpopulation
        try
            [p, tbl, stats] = friedman(testData, num_conditions, 0.05);
            disp('Friedman Test Results:');
            disp(tbl);  % Display the ANOVA-like table
            disp(['p-value: ', num2str(p)]);
            
            % Post-hoc analysis if the result is significant
            if p < 0.05
                disp('Post-hoc analysis:');
                multcompare(stats);
            end
        catch ME
            % Handle any error that occurs during the Friedman test
            disp('Error during Friedman test:');
            disp(ME.message);
        end
    end
end