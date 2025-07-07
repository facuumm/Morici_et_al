function perform2WayPairedANOVA_fixed(data, condition, subpop, subjectID)
% perform2WayPairedANOVA_fixed performs a 2-way repeated measures ANOVA
% using fitrm and ranova, assuming the subject ID is provided.

    % Convert inputs to categorical
    condition = categorical(condition);
    subpop = categorical(subpop);
    subjectID = categorical(subjectID);

    % Make long table
    tbl = table(subjectID, condition, subpop, data, ...
                'VariableNames', {'Subject', 'Condition', 'Subpopulation', 'Data'});

    % Get condition levels
    condLevels = categories(condition);
    nConds = numel(condLevels);
    subjects = unique(subjectID);
    nSubjects = numel(subjects);

    % Preallocate
    dataMatrix = NaN(nSubjects, nConds);
    subpopVec = categorical(strings(nSubjects,1));

    % Fill matrix
    for i = 1:nSubjects
        subjData = tbl(tbl.Subject == subjects(i), :);
        subpopVec(i) = subjData.Subpopulation(1);
        for j = 1:nConds
            c = condLevels{j};
            val = subjData.Data(subjData.Condition == c);
            if numel(val) > 1
                warning('Multiple values for subject %s in condition %s. Averaging.', string(subjects(i)), c);
                val = mean(val);
            end
            if ~isempty(val)
                dataMatrix(i, j) = val;
            end
        end
    end

    % Wide table for fitrm
    condNames = strcat("Cond", condLevels);
    dataTbl = array2table(dataMatrix, 'VariableNames', condNames);
    dataTbl.Subpopulation = subpopVec;

    % Define within-subject design
    WithinDesign = table(condLevels, 'VariableNames', {'Condition'});

    % Fit repeated measures model
    rm = fitrm(dataTbl, sprintf('%s-%s ~ Subpopulation', condNames{1}, condNames{end}), ...
               'WithinDesign', WithinDesign);

    % Run ANOVA
    ranovatbl = ranova(rm);
    disp('--- Repeated Measures ANOVA ---');
    disp(ranovatbl);

    % Post-hoc
    try
        disp('--- Post-hoc: Condition Comparisons ---');
        c = multcompare(rm, 'Condition');
        disp(c);
    catch
        warning('Post-hoc comparison failed.');
    end
end