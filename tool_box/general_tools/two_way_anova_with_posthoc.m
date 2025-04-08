function two_way_anova_with_posthoc(matrix1, matrix2)
    % Input:
    % matrix1: [n1 x 2] for dHPC -> column 1 = group (1–3), column 2 = value
    % matrix2: [n2 x 2] for vHPC -> same format

    % --- Combine matrices for 2-way ANOVA ---
    all_values = [matrix1(:,2); matrix2(:,2)];
    all_groups = [matrix1(:,1); matrix2(:,1)];
    all_regions = [ones(size(matrix1,1),1); 2*ones(size(matrix2,1),1)]; % 1 = dHPC, 2 = vHPC

    % Convert to table with categorical variables
    T = table(all_values, all_groups, all_regions, ...
              'VariableNames', {'Values', 'Group', 'Region'});
    T.Group = categorical(T.Group);
    T.Region = categorical(T.Region);

    % Run 2-way ANOVA with interaction
    disp('=== TWO-WAY ANOVA ===');
    [p, tbl, stats] = anovan(T.Values, {T.Group, T.Region}, ...
                             'model', 'interaction', ...
                             'varnames', {'Group', 'Region'}, ...
                             'display', 'on');

    % Post-hoc: Tukey for each region separately
    disp('=== Tukey Post-hoc Test for dHPC ===');
    do_tukey_posthoc(matrix1, 'dHPC');

    disp('=== Tukey Post-hoc Test for vHPC ===');
    do_tukey_posthoc(matrix2, 'vHPC');
    
    % Normality check
    check_normality(matrix1(:,2), 'dHPC');
    check_normality(matrix2(:,2), 'vHPC');

    % --- Plotting ---
    figure;
    hold on;
    colors = {[0.2 0.4 1], [1 0.2 0.2]}; % blue and red
    jitter_strength = 0.15;

    % dHPC
    for g = 1:3
        vals = matrix1(matrix1(:,1)==g, 2);
        x = g + (rand(size(vals)) - 0.5) * jitter_strength;
        scatter(x, vals, 30, 'filled', 'MarkerFaceColor', colors{1}, 'MarkerEdgeColor', 'none');
    end
    means_d = accumarray(matrix1(:,1), matrix1(:,2), [3,1], @mean);
    for g = 1:3
        plot([g-0.15 g+0.15], [means_d(g) means_d(g)], 'Color', 'k', 'LineWidth', 2);
    end

    % vHPC
    offset = 4;
    for g = 1:3
        vals = matrix2(matrix2(:,1)==g, 2);
        x = offset + g - 1 + (rand(size(vals)) - 0.5) * jitter_strength;
        scatter(x, vals, 30, 'filled', 'MarkerFaceColor', colors{2}, 'MarkerEdgeColor', 'none');
    end
    means_v = accumarray(matrix2(:,1), matrix2(:,2), [3,1], @mean);
    for g = 1:3
        x_pos = offset + g - 1;
        plot([x_pos-0.15 x_pos+0.15], [means_v(g) means_v(g)], 'Color', 'k', 'LineWidth', 2);
    end

    % Final plot tweaks
    xlim([0.5 6.5]);
    xticks([1 2 3 4 5 6]);
    xticklabels({'G1-dHPC', 'G2-dHPC', 'G3-dHPC', 'G1-vHPC', 'G2-vHPC', 'G3-vHPC'});
    ylabel('Values');
    title('Group comparison by region');
    legend({'dHPC', 'vHPC'}, 'Location', 'northwest');
    box on;
    hold off;
end

function do_tukey_posthoc(matrix, region_name)
    group = matrix(:,1);
    values = matrix(:,2);

    [~, ~, stats] = anova1(values, group, 'off');
    [c, ~, ~, gnames] = multcompare(stats, 'ctype', 'tukey-kramer', 'display', 'off');

    num_comparisons = size(c,1);
    bonf_p = min(c(:,6) * num_comparisons, 1); % Bonferroni

    disp(['Tukey HSD Results with Bonferroni Correction for ', region_name, ':']);
    disp('Comparison:     Mean Diff      Raw p-value     Bonferroni p');
    for i = 1:num_comparisons
        g1 = gnames{c(i,1)};
        g2 = gnames{c(i,2)};
        mean_diff = c(i,4);
        raw_p = c(i,6);
        corrected_p = bonf_p(i);
        fprintf('%s vs %s:    %.3f        p = %.4f     p_bonf = %.4f\n', ...
                g1, g2, mean_diff, raw_p, corrected_p);
    end
end

function check_normality(data, region_name)
    [h, p] = lillietest(data);
    if h == 0
        disp([region_name, ' data is normally distributed (p = ', num2str(p), ').']);
    else
        disp([region_name, ' data is NOT normally distributed (p = ', num2str(p), ').']);
    end
end