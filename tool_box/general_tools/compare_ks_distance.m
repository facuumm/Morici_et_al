function p_value = compare_ks_distance(A, B, C, D, num_permutations)
% compare_ks_distance compares two pairs of distributions (A vs. B and C vs. D)
% using Kolmogorov-Smirnov (KS) distance and a permutation test.
%
% Inputs:
%   A               - (vector) Data from the first distribution (group A).
%   B               - (vector) Data from the second distribution (group B).
%   C               - (vector) Data from the third distribution (group C).
%   D               - (vector) Data from the fourth distribution (group D).
%   num_permutations - (scalar) Number of permutations for the test.
%
% Outputs:
%   p_value         - (scalar) p-value from permutation test assessing the
%                     difference between the two KS statistics.
%
% Example:
%   p_value = compare_ks_distance(A, B, C, D, 1000);

% Compute the KS statistics for both pairs
[~,~, ks_stat_AB] = kstest2(A, B);
[~,~, ks_stat_CD] = kstest2(C, D);

% Compute observed difference in KS statistics
observed_diff = ks_stat_AB - ks_stat_CD;

% Combine all samples for permutation
all_data = [A; B; C; D];
labels = [ones(size(A)); 2*ones(size(B)); 3*ones(size(C)); 4*ones(size(D))];

% Permutation test
perm_diffs = zeros(num_permutations, 1);
for i = 1:num_permutations
    shuffled_labels = labels(randperm(length(labels)));
    A_perm = all_data(shuffled_labels == 1);
    B_perm = all_data(shuffled_labels == 2);
    C_perm = all_data(shuffled_labels == 3);
    D_perm = all_data(shuffled_labels == 4);
    
    % Compute new KS statistics
    [~,~, ks_perm_AB] = kstest2(A_perm, B_perm);
    [~,~, ks_perm_CD] = kstest2(C_perm, D_perm);
    
    % Store the permuted KS difference
    perm_diffs(i) = ks_perm_AB - ks_perm_CD;
end

% Compute p-value (two-tailed test)
p_value = mean(abs(perm_diffs) >= abs(observed_diff));

end