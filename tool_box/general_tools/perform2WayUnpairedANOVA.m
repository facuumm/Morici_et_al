function perform2WayUnpairedANOVA(y, structure, condition)
% perform2WayPairedANOVA performs a two-way ANOVA with paired design factors,
% followed by specific post-hoc unpaired t-tests with Bonferroni correction.
%
% Inputs:
%   y          - Numeric vector of observations (dependent variable)
%   structure  - Numeric vector for "structure" factor: 1=dHPC, 2=vHPC
%   condition  - Numeric vector for "condition" factor: 1=shock, 2=decorrelated
%
% Example usage:
%   perform2WayPairedANOVA(y, structure, condition)
%
% The function performs:
%   - Two-way ANOVA with interaction between structure and condition
%   - Post-hoc tests for:
%       1) dHPC shock vs dHPC decorr
%       2) vHPC shock vs vHPC decorr
%       3) dHPC decorr vs vHPC decorr
%   - Bonferroni correction for the 3 post-hoc comparisons


if length(y) ~= length(structure) || length(y) ~= length(condition)
    error('Inputs y, structure, and condition must be the same length');
end

% Run 2-way ANOVA with interaction
[p, tbl, stats] = anovan(y, {structure, condition}, ...
    'model', 'interaction', 'varnames', {'Structure', 'Condition'});

fprintf('\nTwo-way ANOVA results:\n');
disp(tbl);

fprintf('p-values:\n');
fprintf('Structure: %.4f\n', p(1));
fprintf('Condition: %.4f\n', p(2));
fprintf('Interaction: %.4f\n\n', p(3));

% Extract groups for post-hoc tests
dHPC_shock = y(structure == 1 & condition == 1);
dHPC_decorr = y(structure == 1 & condition == 2);
vHPC_shock = y(structure == 2 & condition == 1);
vHPC_decorr = y(structure == 2 & condition == 2);

% Preallocate results
group1 = ["dHPC_shock", "vHPC_shock", "dHPC_decorr"];
group2 = ["dHPC_decorr", "vHPC_decorr", "vHPC_decorr"];
pvals = zeros(3,1);

% Perform unpaired t-tests for post-hoc comparisons
[~, pvals(1)] = ttest(dHPC_shock, dHPC_decorr);
[~, pvals(2)] = ttest(vHPC_shock, vHPC_decorr);
[~, pvals(3)] = ttest2(dHPC_decorr, vHPC_decorr);

% Bonferroni correction for 3 tests
bonferroni_pvals = min(pvals * 3, 1);

% Create results table
results = table(group1', group2', pvals, bonferroni_pvals, ...
    'VariableNames', {'Group1', 'Group2', 'pValue', 'Bonferroni_pValue'});

fprintf('Post-hoc pairwise comparisons (unpaired t-tests) with Bonferroni correction:\n');
disp(results);
end