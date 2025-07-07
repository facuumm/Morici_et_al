function [pVals, postHocResults] = perform2WayWithinANOVA(y, factor1, factor2, subjectID)
    % Convert to categorical
    factor1 = categorical(factor1);
    factor2 = categorical(factor2);
    subjectID = categorical(subjectID);
    
    % Create a table with all variables
    dataTable = table(y, factor1, factor2, subjectID, ...
        'VariableNames', {'Y','Factor1','Factor2','Subject'});
    
    % Fit repeated measures model
    rm = fitrm(dataTable, 'Y ~ Factor1*Factor2', 'WithinDesign', dataTable(:,2:3));
    
    % Run repeated measures ANOVA
    ranovatbl = ranova(rm);
    
    % Extract p-values
    pVals.Factor1 = ranovatbl.pValueGG(strcmp(ranovatbl.Properties.RowNames, 'Factor1'));
    pVals.Factor2 = ranovatbl.pValueGG(strcmp(ranovatbl.Properties.RowNames, 'Factor2'));
    pVals.Interaction = ranovatbl.pValueGG(strcmp(ranovatbl.Properties.RowNames, 'Factor1:Factor2'));
    
    % Post-hoc tests (for significant effects)
    postHocResults = struct();
    if pVals.Factor1 < 0.05
        postHocResults.Factor1 = multcompare(rm, 'Factor1');
    end
    if pVals.Factor2 < 0.05
        postHocResults.Factor2 = multcompare(rm, 'Factor2');
    end
    if pVals.Interaction < 0.05
        postHocResults.Interaction = multcompare(rm, 'Factor1', 'By', 'Factor2');
    end
end