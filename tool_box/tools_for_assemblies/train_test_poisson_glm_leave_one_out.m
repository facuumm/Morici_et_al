function [output_metrics] = train_test_poisson_glm_leave_one_out(spike_data, shock_responsive, assembly_activity, train_window, test_window, region_boundaries,  varargin)
% TRAIN_TEST_POISSON_GLM_BALANCED Balanced comparison of neuron contributions
%
%   Inputs:
%       spike_data          : [time x neurons+1] matrix
%       shock_responsive    : [neurons x 2] matrix [membership, shock_response]
%       assembly_activity   : [time x 2] matrix
%       train_window        : [start, end] times for training
%       test_window         : [start, end] times for testing
%       region_boundaries   : [dHPC_end, vHPC_end] neuron indices
%
%   Optional Parameters:
%       'n_permutations'    : Number of permutations (default: 1000)
%       'compute_ci'        : Compute confidence intervals (default: true)

    % ==================== INPUT PARSING ====================
    p = inputParser;
    
    % Add required parameters
    addRequired(p, 'spike_data', @isnumeric);
    addRequired(p, 'shock_responsive', @islogical);
    addRequired(p, 'assembly_activity', @isnumeric);
    addRequired(p, 'train_window', @isnumeric);
    addRequired(p, 'test_window', @isnumeric);
    addRequired(p, 'region_boundaries', @isnumeric);
    
    % Add optional parameters with defaults
    addParameter(p, 'n_permutations', 1000, @isnumeric);
    addParameter(p, 'compute_ci', true, @islogical);
    
    parse(p, spike_data, shock_responsive, assembly_activity, train_window, test_window, region_boundaries, varargin{:});

    % ==================== NEURON CLASSIFICATION ====================
    is_member = shock_responsive(:,1) == 1;
    is_shock_resp = shock_responsive(:,2) == 1;
    
    % Define regions
    dHPC_neurons = 1:region_boundaries(1);
    vHPC_neurons = (region_boundaries(1)+1):region_boundaries(2);
    
    % Categorize neurons
    neuron_groups = struct();
    neuron_groups.dHPC.shock_member = intersect(dHPC_neurons, find(is_member & is_shock_resp));
    neuron_groups.dHPC.member_only = intersect(dHPC_neurons, find(is_member & ~is_shock_resp));
    neuron_groups.dHPC.non_member = intersect(dHPC_neurons, find(~is_member));
    
    neuron_groups.vHPC.shock_member = intersect(vHPC_neurons, find(is_member & is_shock_resp));
    neuron_groups.vHPC.member_only = intersect(vHPC_neurons, find(is_member & ~is_shock_resp));
    neuron_groups.vHPC.non_member = intersect(vHPC_neurons, find(~is_member));

    % ==================== DATA PREPARATION ====================
    % Training data
    train_idx = InIntervals(spike_data(:,1), train_window);
    train_spikes = spike_data(train_idx, 2:end);
    train_assembly = normalize_assembly(assembly_activity(InIntervals(assembly_activity(:,1), train_window), 2));
    
    % Testing data
    test_idx = InIntervals(spike_data(:,1), test_window);
    test_spikes = spike_data(test_idx, 2:end);
    true_assembly = normalize_assembly(assembly_activity(InIntervals(assembly_activity(:,1), test_window), 2));

    % ==================== MODEL TRAINING ====================
    X_train = [mean(train_spikes,2), var(train_spikes,0,2), ones(size(train_spikes,1),1)];
    model = fitglm(X_train, train_assembly, 'Distribution','poisson', 'Link','log');

    % ==================== BALANCED TESTING ====================
    % Get reference count (max shock-responsive members across regions)
    n_ref = max([length(neuron_groups.dHPC.shock_member), length(neuron_groups.vHPC.shock_member)]);
    
    % Initialize output structure
    output_metrics = struct();
    regions = {'dHPC', 'vHPC'};
    group_types = {'shock_member', 'member_only', 'non_member'};
    
    % Full model performance
    X_test_full = create_features(test_spikes);
    [output_metrics.full_model, ~] = evaluate_model(model, X_test_full, true_assembly);
    
    % Test each region and group
    for r = 1:length(regions)
        region = regions{r};
        
        for g = 1:length(group_types)
            group = group_types{g};
            neuron_ids = neuron_groups.(region).(group);
            
            if isempty(neuron_ids)
                output_metrics.(region).(group) = [];
                continue;
            end
            
            % Create balanced samples
            n_samples = min(n_ref, length(neuron_ids));
            sample_idx = repmat(1:length(neuron_ids), 1, ceil(n_ref/length(neuron_ids)));
            sample_idx = sample_idx(1:n_ref);
            
            % Test each sample
            for i = 1:n_ref
                neuron_idx = neuron_ids(sample_idx(i));
                
                temp_spikes = test_spikes;
                temp_spikes(:,neuron_idx) = [];
                
                X_test = create_features(temp_spikes);
                [metrics, ~] = evaluate_model(model, X_test, true_assembly);
                
                output_metrics.(region).(group)(i) = struct(...
                    'neuron_id', neuron_idx, ...
                    'pearson', metrics.pearson, ...
                    'pearson_delta', metrics.pearson - output_metrics.full_model.pearson, ...
                    'spearman', metrics.spearman, ...
                    'spearman_delta', metrics.spearman - output_metrics.full_model.spearman, ...
                    'mse', metrics.mse, ...
                    'mse_delta', metrics.mse - output_metrics.full_model.mse, ...
                    'sample_number', i);
            end
        end
    end

    % ==================== HELPER FUNCTIONS ====================
    function features = create_features(spike_matrix)
        features = [mean(spike_matrix,2), var(spike_matrix,0,2), ones(size(spike_matrix,1),1)];
    end

    function [metrics, predictions] = evaluate_model(model, X_test, y_true)
        predictions = predict(model, X_test);
        valid = ~isnan(predictions) & ~isnan(y_true);
        
        if ~any(valid)
            metrics = struct('pearson',nan, 'spearman',nan, 'mse',nan);
            return;
        end
        
        pearson_r = corr(predictions(valid), y_true(valid));
        spearman_r = corr(predictions(valid), y_true(valid), 'Type','Spearman');
        mse = mean((predictions(valid) - y_true(valid)).^2);
        
        metrics = struct('pearson',pearson_r, 'spearman',spearman_r, 'mse',mse);
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
end