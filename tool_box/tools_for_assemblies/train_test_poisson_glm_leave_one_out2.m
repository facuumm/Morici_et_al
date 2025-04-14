function [output_metrics] = train_test_poisson_glm_leave_one_out2(spike_data, shock_responsive, assembly_activity, train_window, test_window, region_boundaries, varargin)
% TRAIN_TEST_POISSON_GLM_LEAVE_ONE_OUT2
% Estimates neuron contributions to assembly activity using GLM and leave-one-out strategy.

    % ---------- Input parsing ----------
    p = inputParser;
    addRequired(p, 'spike_data', @isnumeric);
    addRequired(p, 'shock_responsive', @islogical);
    addRequired(p, 'assembly_activity', @isnumeric);
    addRequired(p, 'train_window', @isnumeric);
    addRequired(p, 'test_window', @isnumeric);
    addRequired(p, 'region_boundaries', @isnumeric);
    addParameter(p, 'n_permutations', 1000, @isnumeric);
    addParameter(p, 'compute_ci', true, @islogical);
    parse(p, spike_data, shock_responsive, assembly_activity, train_window, test_window, region_boundaries, varargin{:});

    % ---------- Neuron classification ----------
    is_member = shock_responsive(:,1) == 1;
    is_shock_resp = shock_responsive(:,2) == 1;

    dHPC_neurons = 1:region_boundaries(1);
    vHPC_neurons = (region_boundaries(1)+1):region_boundaries(2);

    neuron_groups = struct();
    neuron_groups.dHPC.shock_member = intersect(dHPC_neurons, find(is_member & is_shock_resp));
    neuron_groups.dHPC.member_only = intersect(dHPC_neurons, find(is_member & ~is_shock_resp));
    neuron_groups.dHPC.non_member  = intersect(dHPC_neurons, find(~is_member));

    neuron_groups.vHPC.shock_member = intersect(vHPC_neurons, find(is_member & is_shock_resp));
    neuron_groups.vHPC.member_only  = intersect(vHPC_neurons, find(is_member & ~is_shock_resp));
    neuron_groups.vHPC.non_member   = intersect(vHPC_neurons, find(~is_member));

    % ---------- Prepare training and test data ----------
    train_idx = InIntervals(spike_data(:,1), train_window);
    train_spikes_raw = spike_data(train_idx, 2:end);
    train_assembly = normalize_assembly(assembly_activity(InIntervals(assembly_activity(:,1), train_window), 2));

    test_idx = InIntervals(spike_data(:,1), test_window);
    test_spikes_raw = spike_data(test_idx, 2:end);
    true_assembly = normalize_assembly(assembly_activity(InIntervals(assembly_activity(:,1), test_window), 2));

    % ---------- Remove zero-variance neurons ----------
    zero_var_neurons = std(train_spikes_raw) == 0;
    train_spikes = train_spikes_raw(:, ~zero_var_neurons);
    test_spikes = test_spikes_raw(:, ~zero_var_neurons);

    retained_neuron_indices = find(~zero_var_neurons);

    % ---------- Fit full GLM ----------
    model = fitglm(train_spikes, train_assembly, 'Distribution','poisson', 'Link','log');

    % ---------- Evaluate full model ----------
    [output_metrics.full_model, ~] = evaluate_model(model, test_spikes, true_assembly);

    % ---------- Initialize ----------
    output_metrics.zero_var_neurons = find(zero_var_neurons);  % Save for reference
    regions = {'dHPC', 'vHPC'};
    group_types = {'shock_member', 'member_only', 'non_member'};

    % ---------- Loop through each region and group ----------
    for r = 1:length(regions)
        region = regions{r};
        for g = 1:length(group_types)
            group = group_types{g};

            all_neurons = neuron_groups.(region).(group);
            neurons_in_model = intersect(all_neurons, retained_neuron_indices);

            if isempty(neurons_in_model)
                output_metrics.(region).(group) = [];
                continue;
            end

            % Remap to new indices post zero-var removal
            mapped_indices = arrayfun(@(nid) find(retained_neuron_indices == nid, 1), neurons_in_model);

            for i = 1:length(mapped_indices)
                idx = mapped_indices(i);

                % Leave-one-out: remove neuron idx
                train_tmp = train_spikes;
                test_tmp = test_spikes;
                train_tmp(:, idx) = [];
                test_tmp(:, idx) = [];

                % Fit GLM without this neuron
                try
                    model_tmp = fitglm(train_tmp, train_assembly, 'Distribution','poisson', 'Link','log');
                    [metrics, ~] = evaluate_model(model_tmp, test_tmp, true_assembly);
                catch
                    metrics = struct('pearson', nan, 'spearman', nan, 'mse', nan, 'log_likelihood', nan);
                end

                output_metrics.(region).(group)(i) = struct( ...
                    'neuron_id', neurons_in_model(i), ...
                    'pearson', metrics.pearson, ...
                    'pearson_delta', delta(metrics.pearson, output_metrics.full_model.pearson), ...
                    'spearman', metrics.spearman, ...
                    'spearman_delta', delta(metrics.spearman, output_metrics.full_model.spearman), ...
                    'mse', metrics.mse, ...
                    'mse_delta', delta(output_metrics.full_model.mse, metrics.mse), ...
                    'log_likelihood', metrics.log_likelihood, ...
                    'sample_number', i);
            end
        end
    end

    % ---------- Helper functions ----------
    function [metrics, predictions] = evaluate_model(model, X, y)
        predictions = predict(model, X);
        valid = ~isnan(predictions) & ~isnan(y);

        if ~any(valid)
            metrics = struct('pearson', nan, 'spearman', nan, 'mse', nan, 'log_likelihood', nan);
            return;
        end

        pearson_r = corr(predictions(valid), y(valid));
        spearman_r = corr(predictions(valid), y(valid), 'Type','Spearman');
        mse = mean((predictions(valid) - y(valid)).^2);
        log_likelihood = sum(y(valid) .* log(predictions(valid)) - predictions(valid));

        metrics = struct('pearson', pearson_r, 'spearman', spearman_r, 'mse', mse, 'log_likelihood', log_likelihood);
    end

    function normalized = normalize_assembly(assembly)
        min_val = min(assembly);
        max_val = max(assembly);
        if max_val == min_val
            normalized = zeros(size(assembly));
        else
            normalized = (assembly - min_val) / (max_val - min_val);
        end
        normalized = max(0, min(1, normalized));
    end

    function d = delta(a, b)
        d = (abs(a) - abs(b));% / (abs(a) + abs(b));
    end
end