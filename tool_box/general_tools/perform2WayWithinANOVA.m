function perform2WayWithinANOVA(y, subpop, condition, subjectID)
% perform2WayWithinANOVA
% Performs a 2-way repeated-measures ANOVA with paired Bonferroni-corrected post-hoc tests.
% Automatically handles:
%   - Averaging multiple entries per subject-condition
%   - Removing subjects with missing conditions
%
% Inputs:
%   y         - data vector
%   subpop    - between-subjects factor (e.g., region)
%   condition - within-subjects factor (e.g., stimulus)
%   subjectID - subject identifier (repeats across conditions)
%
% All vectors must be the same length.

% Validate inputs
assert(numel(y) == numel(subpop) && numel(y) == numel(condition) && numel(y) == numel(subjectID), ...
    'All input vectors must have the same length.');

% Convert to columns
y         = y(:);
subpop    = subpop(:);
condition = condition(:);
subjectID = subjectID(:);

% Convert subjectID to categorical (needed for anovan)
subjectID_cat = categorical(subjectID);

% --- Run 2-way repeated-measures ANOVA ---
[p, tbl_anova, stats] = anovan(y, ...
    {condition, subpop, subjectID_cat}, ...
    'model', 'interaction', ...
    'random', 3, ...
    'varnames', {'Condition', 'Subpopulation', 'Subject'}, ...
    'display', 'off');

fprintf('\n--- 2-Way Repeated Measures ANOVA ---\n');
disp(tbl_anova);
fprintf('p(Condition):     %.4f\n', p(1));
fprintf('p(Subpopulation): %.4f\n', p(2));
fprintf('p(Interaction):   %.4f\n', p(3));

% --- Paired Post-hoc tests (Bonferroni) within each subpopulation ---
fprintf('\n--- Paired Post-Hoc Tests (Bonferroni corrected) ---\n');

uniqueSubpops = unique(subpop);
uniqueConds   = unique(condition);
nConds        = numel(uniqueConds);

for s = 1:numel(uniqueSubpops)
    currPop = uniqueSubpops(s);
    fprintf('\nSubpopulation: %d\n', currPop);

    % Extract data for this subpopulation
    mask      = subpop == currPop;
    y_pop     = y(mask);
    cond_pop  = condition(mask);
    subj_pop  = subjectID(mask);

    % Get unique subject list
    subjects = unique(subj_pop);

    % Initialize subject x condition matrix
    dataMatrix = NaN(numel(subjects), nConds);

    for i = 1:numel(subjects)
        for j = 1:nConds
            sel = subj_pop == subjects(i) & cond_pop == uniqueConds(j);
            vals = y_pop(sel);
            if isempty(vals)
                % leave as NaN
            elseif numel(vals) > 1
                warning('Averaging multiple entries for subj %d, cond %d (n=%d)', ...
                        subjects(i), uniqueConds(j), numel(vals));
                dataMatrix(i,j) = mean(vals);
            else
                dataMatrix(i,j) = vals;
            end
        end
    end

    % Drop subjects with missing data
    validRows = all(~isnan(dataMatrix), 2);
    nValid    = sum(validRows);

    if nValid < 2
        fprintf('❌ Not enough complete data for post-hoc comparisons in Subpop %d\n', currPop);
        continue;
    end

    dataMatrix = dataMatrix(validRows, :);

    % Paired t-tests with Bonferroni correction
    comparisons = nchoosek(1:nConds, 2);
    nTests = size(comparisons, 1);

    for c = 1:nTests
        i = comparisons(c,1);
        j = comparisons(c,2);
        [~, pval, ~, stats] = ttest(dataMatrix(:,i), dataMatrix(:,j));
        pBonf = min(pval * nTests, 1);  % Bonferroni correction
        fprintf('Cond %d vs %d: t(%d) = %.2f, p = %.4f (Bonf = %.4f)\n', ...
            uniqueConds(i), uniqueConds(j), stats.df, stats.tstat, pval, pBonf);
    end
end
end