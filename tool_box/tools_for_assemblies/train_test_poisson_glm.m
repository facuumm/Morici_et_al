function [metrics_with, metrics_without, comparison] = train_test_poisson_glm(spike_data, shock_responsive, assembly_activity, train_window, test_window, varargin)
% Función mejorada con métricas robustas para comparación de modelos
% Opciones adicionales:
%   'n_permutations': Número de permutaciones (default: 1000)
%   'compute_ci': Calcular intervalos de confianza (default: true)

p = inputParser;
addParameter(p, 'n_permutations', 1000);
addParameter(p, 'compute_ci', true);
parse(p, varargin{:});

% --- 1. División train/test (original) ---
In = InIntervals(spike_data(:,1), train_window);
train_spikes = spike_data(In, 2:end); 
train_assembly = assembly_activity(InIntervals(assembly_activity(:,1), train_window));

In = InIntervals(spike_data(:,1), test_window);
test_spikes = spike_data(In, 2:end);
test_assembly = assembly_activity(InIntervals(assembly_activity(:,1), test_window));

% --- 2. Cálculo de características mejorado ---
% Medias y varianzas
features_train = compute_features(train_spikes, shock_responsive);
features_test_with = compute_features(test_spikes, shock_responsive);
features_test_without = compute_features(test_spikes, shock_responsive, 'exclude_responsive');

% --- 3. Modelado GLM ---
if sum(shock_responsive(:,2)) > 0
    model = fitglm(features_train, train_assembly, ...
        'Distribution', 'poisson', 'Link', 'log', 'Intercept', false);
    
    % Predicciones
    test_pred_with = predict(model, features_test_with);
    test_pred_without = predict(model, features_test_without);
    
    % --- 4. Métricas completas ---
    [metrics_with, metrics_without] = compute_basic_metrics(test_pred_with, test_pred_without, test_assembly);
    comparison = compute_advanced_comparison(test_pred_with, test_pred_without, test_assembly, p.Results);
else
    [metrics_with, metrics_without] = create_null_metrics();
    comparison = create_null_comparison();
end

% --- Funciones auxiliares ---
    function features = compute_features(spikes, responsive, mode)
        if nargin < 3, mode = 'all'; end
        
        % Medias
        mean_all = mean(spikes, 2);
        mean_resp = mean(spikes(:, responsive(:,2) == 1), 2);
        mean_nonresp = mean(spikes(:, responsive(:,2) == 0), 2);
        
        % Varianzas
        var_all = var(spikes, 0, 2);
        var_resp = var(spikes(:, responsive(:,2) == 1), 0, 2);
        var_nonresp = var(spikes(:, responsive(:,2) == 0), 0, 2);
        
        % Exclusión de neuronas shock-responsivas si se solicita
        if strcmp(mode, 'exclude_responsive')
            mean_resp(:) = 0;
            var_resp(:) = 0;
        end
        
        features = [mean_all, mean_resp, mean_nonresp, var_all, var_resp, var_nonresp, ones(size(mean_all))];
    end

    function [metrics_w, metrics_wo] = compute_basic_metrics(pred_w, pred_wo, actual)
        % Pearson
        [r_pearson_w, p_pearson_w] = corr(pred_w, actual);
        [r_pearson_wo, p_pearson_wo] = corr(pred_wo, actual);
        
        % Spearman
        [r_spearman_w, p_spearman_w] = corr(pred_w, actual, 'Type', 'Spearman');
        [r_spearman_wo, p_spearman_wo] = corr(pred_wo, actual, 'Type', 'Spearman');
        
        % MSE
        mse_w = mean((pred_w - actual).^2);
        mse_wo = mean((pred_wo - actual).^2);
        
        metrics_w = struct('Pearson_r', r_pearson_w, 'Pearson_p', p_pearson_w, ...
                         'Spearman_r', r_spearman_w, 'Spearman_p', p_spearman_w, ...
                         'MSE', mse_w);
        
        metrics_wo = struct('Pearson_r', r_pearson_wo, 'Pearson_p', p_pearson_wo, ...
                          'Spearman_r', r_spearman_wo, 'Spearman_p', p_spearman_wo, ...
                          'MSE', mse_wo);
    end

    function comp = compute_advanced_comparison(pred_w, pred_wo, actual, params)
        % 1. Test de permutación para MSE
        original_diff = mean((pred_wo - actual).^2) - mean((pred_w - actual).^2);
        perm_diffs = zeros(params.n_permutations, 1);
        
        parfor p = 1:params.n_permutations
            combined = [pred_w; pred_wo];
            shuffled = combined(randperm(length(combined)));
            diff_p = mean((shuffled(1:end/2) - actual).^2) - mean((shuffled(end/2+1:end) - actual).^2);
            perm_diffs(p) = diff_p;
        end
        pval_perm = mean(perm_diffs >= original_diff);
        
        % 2. Likelihood Ratio
        ll_w = sum(actual.*log(pred_w) - pred_w);
        ll_wo = sum(actual.*log(pred_wo) - pred_wo);
        llr = -2*(ll_wo - ll_w);
        pval_ll = 1 - chi2cdf(llr, 1);
        
        % 3. Intervalos de confianza Spearman (bootstrap)
        if params.compute_ci
            boots_w = bootstrp(1000, @(x,y) corr(x,y,'Type','Spearman'), pred_w, actual);
            boots_wo = bootstrp(1000, @(x,y) corr(x,y,'Type','Spearman'), pred_wo, actual);
            ci_w = prctile(boots_w, [2.5, 97.5]);
            ci_wo = prctile(boots_wo, [2.5, 97.5]);
            ci_overlap = ~(ci_wo(1) > ci_w(2) || ci_w(1) > ci_wo(2));
        else
            ci_w = [nan nan]; ci_wo = [nan nan]; ci_overlap = nan;
        end
        
        comp = struct(...
            'MSE_diff', original_diff, ...
            'MSE_ratio', mean((pred_wo - actual).^2)/mean((pred_w - actual).^2), ...
            'Permutation_pval', pval_perm, ...
            'Likelihood_ratio', llr, ...
            'Likelihood_pval', pval_ll, ...
            'Spearman_ci_with', ci_w, ...
            'Spearman_ci_without', ci_wo, ...
            'CIs_overlap', ci_overlap);
    end

    function [metrics_w, metrics_wo] = create_null_metrics()
        null = struct('r', nan, 'p', nan);
        metrics_w = struct('Pearson_r', nan, 'Pearson_p', nan, ...
                         'Spearman_r', nan, 'Spearman_p', nan, ...
                         'MSE', nan);
        metrics_wo = metrics_w;
    end

    function comp = create_null_comparison()
        comp = struct('MSE_diff', nan, 'MSE_ratio', nan, ...
                     'Permutation_pval', nan, ...
                     'Likelihood_ratio', nan, ...
                     'Likelihood_pval', nan, ...
                     'Spearman_ci_with', [nan nan], ...
                     'Spearman_ci_without', [nan nan], ...
                     'CIs_overlap', nan);
    end
end