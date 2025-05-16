function [member, shock_resp] = train_test_poisson_glm_member_shuffle(spike_data, shock_responsive, assembly_activity, train_window, test_window, region_boundaries, varargin)
% TRAIN_TEST_POISSON_GLM_MEMBER_SHUFFLE
% --------------------------------------------------------------
% Trains a Poisson GLM using the spike train of each neuron that is
% either a member or a shock-responsive member. Compares GLM prediction
% against shuffled spike trains using ShuffleSpks and spike_train_construction.
% --------------------------------------------------------------
% INPUTS:
%   spike_data         - structure with fields:
%       spks           - spike timestamps
%       clusters       - neuron IDs
%       cellulartype   - neuron type information (optional)
%       limits         - [start_time, end_time] for valid spike range
%   shock_responsive   - [N x 2] logical matrix: col 1 = is_member, col 2 = is_shock_resp
%   assembly_activity  - [T x 2] matrix: timestamps and assembly signal
%   train_window       - [1 x 2] vector: start and end time for training
%   test_window        - [1 x 2] vector: start and end time for testing
%   region_boundaries  - [1 x 2] vector: indices separating dHPC and vHPC neurons
% OPTIONAL:
%   'n_permutations'   - number of shuffles (default = 100)
% OUTPUT:
%   member, shock_resp - structures with fields:
%       dHPC, vHPC: each a matrix [neuron_id, pearson_real, mean_pearson_shuffle, gain]

    % ---------- Input parsing ----------
    p = inputParser;
    addRequired(p, 'spike_data', @isstruct);
    addRequired(p, 'shock_responsive', @islogical);
    addRequired(p, 'assembly_activity', @isnumeric);
    addRequired(p, 'train_window', @isnumeric);
    addRequired(p, 'test_window', @isnumeric);
    addRequired(p, 'region_boundaries', @isnumeric);
    addParameter(p, 'n_permutations', 100, @isnumeric);
    parse(p, spike_data, shock_responsive, assembly_activity, train_window, test_window, region_boundaries, varargin{:});

    n_permutations = p.Results.n_permutations;

    % ---------- SpikeTrains construction ----------
    [Spikes , bins , ~] = spike_train_construction(spike_data.spks, spike_data.clusters, spike_data.cellulartype, 0.025, spike_data.limits, [], false, true);

    % ---------- Neuron classification ----------
    is_member = shock_responsive(:,1) == 1;
    is_shock_resp = shock_responsive(:,2) & shock_responsive(:,1);

    dHPC_neurons = 1:region_boundaries(1);
    vHPC_neurons = (region_boundaries(1)+1):region_boundaries(2);

    % ---------- Prepare training and test data ----------
    train_idx = InIntervals(bins, train_window);
    train_assembly = normalize_assembly(assembly_activity(InIntervals(assembly_activity(:,1), train_window), 2));

    test_idx = InIntervals(bins, test_window);
    test_assembly = normalize_assembly(assembly_activity(InIntervals(assembly_activity(:,1), test_window), 2));

    % ---------- Initialize outputs ----------
    member.dHPC = [];
    member.vHPC = [];
    shock_resp.dHPC = [];
    shock_resp.vHPC = [];

    % ---------- Analysis function ----------
    function results = analyze_group(neuron_ids)
        results = zeros(length(neuron_ids), 10); % Adjusted size for all metrics
        
        for i = 1:length(neuron_ids)
            nid = neuron_ids(i);
            real_spk_times = spike_data.spks(spike_data.spks(:,1) == nid,2);
            spike_data_single.spks = real_spk_times;
            
            [Spikes_single, bins] = binspikes(spike_data_single.spks, 1/0.025, spike_data.limits);
            Spikes_single = gaussfilt(bins, Spikes_single, (0.025)/sqrt(12));
            
            X_train = Spikes_single(train_idx);
            X_test = Spikes_single(test_idx);
            
            if std(X_train) == 0 || std(X_test) == 0 || all(isnan(X_test))
                results(i,:) = [nid, NaN(1,9)];
                continue;
            end
            
            try
                model = fitglm(X_train, train_assembly, 'Distribution','poisson', 'Link','log');
                y_pred = predict(model, X_test);
                
                pearson_real = corr(y_pred, test_assembly, 'Rows','complete');
                spearman_real = corr(y_pred, test_assembly, 'Type','Spearman', 'Rows','complete');
                ll_real = sum(test_assembly .* log(y_pred) - y_pred);
                
            catch
                pearson_real = NaN;
                spearman_real = NaN;
                ll_real = NaN;
            end
            
            % Initialize shuffled metrics
            pearson_shuffled = nan(n_permutations,1);
            spearman_shuffled = nan(n_permutations,1);
            ll_shuffled = nan(n_permutations,1);
            
            for p = 1:n_permutations
                shuffled_spk_times = ShuffleSpks(real_spk_times);
                valid_spk_times = shuffled_spk_times(shuffled_spk_times >= spike_data.limits(1) & shuffled_spk_times <= spike_data.limits(2));
                if isempty(valid_spk_times)
                    continue;
                end
                
                spike_data_shuffled.spks = valid_spk_times;
                [Shuffled, bins] = binspikes(spike_data_shuffled.spks, 1/0.025, spike_data.limits);
                Shuffled = gaussfilt(bins, Shuffled, (0.025)/sqrt(12));
                shuffled_X_test = Shuffled(test_idx);
                
                try
                    shuffled_model = fitglm(X_train, train_assembly, 'Distribution','poisson', 'Link','log');
                    shuffled_pred = predict(shuffled_model, shuffled_X_test);
                    
                    pearson_shuffled(p) = corr(shuffled_pred, test_assembly, 'Rows','complete');
                    spearman_shuffled(p) = corr(shuffled_pred, test_assembly, 'Type','Spearman', 'Rows','complete');
                    ll_shuffled(p) = sum(test_assembly .* log(shuffled_pred) - shuffled_pred);
                    
                catch
                    continue;
                end
            end
            
            % Mean shuffled values
            mean_pearson = mean(pearson_shuffled, 'omitnan');
            mean_spearman = mean(spearman_shuffled, 'omitnan');
            mean_ll = mean(ll_shuffled, 'omitnan');
            
            % Gain calculation (add small epsilon to avoid div-by-zero)
            eps = 1e-8;
            pearson_gain = abs(pearson_real) / (abs(mean_pearson) + eps);
            spearman_gain = abs(spearman_real) / (abs(mean_spearman) + eps);
            ll_gain = ll_real / (mean_ll + eps);
            
            results(i,:) = [nid, pearson_real, mean_pearson, pearson_gain, ...
                spearman_real, mean_spearman, spearman_gain, ...
                ll_real, mean_ll, ll_gain];
        end
    end

    % ---------- Run analysis ----------
    member.dHPC = analyze_group(intersect(dHPC_neurons, find(is_member)));
    member.vHPC = analyze_group(intersect(vHPC_neurons, find(is_member)));
    shock_resp.dHPC = analyze_group(intersect(dHPC_neurons, find(is_shock_resp)));
    shock_resp.vHPC = analyze_group(intersect(vHPC_neurons, find(is_shock_resp)));

    % ---------- Helper ----------
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
