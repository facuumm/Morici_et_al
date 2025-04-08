function [Pred_with , Pred_without] = GLM_ShockSU(path)
% This function calculates the Pre and Post sleep assemblies rate.
%
% INPUTS
% path: cell, inside each cell you should put the path of each session you
%       want to analyze
%
% OUTPUT
% Pre, Post: Structure, it store the SU response sourrounding the ripple.
%
% T: vector, it contains time vector for visualization
%
% I: Structure, it contains the tag to determine if the SU is or not
%    shock responssive
%
% curve: Structure, it contains the response of the neuron sourrounding the
%        shock. It is zscored.
%
% Morci Juan Facundo 08/2024

% variables to use in the script
criteria_fr = 0;
criteria_type = 0; % criteria for celltype (0:pyr, 1:int, 2:all)
sm = 0;
bin = 0.025;

Pred_with.dHPC = [];   Pred_without.dHPC = [];
Pred_with.vHPC = [];   Pred_without.vHPC = [];


%% Main loop, to iterate across sessions
for tt = 1:length(path)
    %List of folders from the path
    files = dir(path{tt});
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
    clear files dirFlags
    
    for t = 1 : length(subFolders)-2
        disp(['-- Initiating analysis of folder #' , num2str(t) , ' from rat #',num2str(tt) , ' --'])
        session = [subFolders(t+2).folder,'\',subFolders(t+2).name];
        cd(session)
        
        %Loading TS of the sessions
        disp('Uploading session time stamps')
        load('session_organization.mat')
        load('behavioral_data.mat')
        
        %% load sleep states
        disp('Uploading sleep scoring')
        x = dir([cd,'\*-states.mat']);    states = load([cd,'\',x.name]);    states = states.states;
        REM.all = ToIntervals(states==5);    NREM.all = ToIntervals(states==3);
        clear x states
        
        baselineTS = baselineTS./1000;                  aversiveTS = aversiveTS./1000;                    rewardTS = rewardTS./1000;
        aversiveTS_run = aversiveTS_run./1000;          rewardTS_run = rewardTS_run./1000;
        
        NREM.baseline = Restrict(NREM.all,baselineTS);   NREM.aversive = Restrict(NREM.all,aversiveTS);   NREM.reward = Restrict(NREM.all,rewardTS);
        REM.baseline = Restrict(REM.all,baselineTS);     REM.aversive = Restrict(REM.all,aversiveTS);     REM.reward = Restrict(REM.all,rewardTS);
        
        %% Load ripples
        if exist('ripplesD_customized2.csv')
            ripplesD = table2array(readtable('ripplesD_customized2.csv'));
            RD = true;
        else
            RD = false;
        end
        
        if exist('ripplesV_customized2.csv')
            ripplesV = table2array(readtable('ripplesV_customized2.csv'));
            RV = true;
        else
            RV = false;
        end
        
        if and(RD,RV)
            RB = true;
            % coordination
            coordinated = [];
            coordinatedV = [];
            cooridnated_event = [];
            for i = 1:length(ripplesD)
                r = ripplesD(i,:);
                tmp = sum(and(ripplesV(:,2)>= r(1,2)-0.1, ripplesV(:,2)<= r(1,2)+0.1));
                if tmp>0
                    z = ripplesV(and(ripplesV(:,2)>= r(1,2)-0.1, ripplesV(:,2)<= r(1,2)+0.1),:);
                    coordinatedV = [coordinatedV ; z];
                    [p,indice] = min(abs(r(2)-z(:,2)));
                    coordinated = [coordinated ; r];
                    
                    %                     if r(2)<z(2) % keep only when dorsal happen first
                    % %                         cooridnated_event = [cooridnated_event ; r];
                    %                         coordinated = [coordinated ; r];
                    %                         coordinatedV = [coordinatedV ; z];
                    %                     end
                    
                    peak = min(min(z(indice,2),r(2)))+abs(z(indice,2)-r(2))/2;
                    low = min([r(1) , z(indice,1)]);
                    up = max([r(3) , z(indice,3)]);
                    cooridnated_event = [cooridnated_event ; low , peak , up];
                    
                    clear tmp2 tmp1 p indice z peak low up
                end
                clear r
            end
            clear x tmp i
            
            [C,IA,IC] = unique(coordinatedV(:,1));
            coordinatedV  = coordinatedV(IA,:); clear C IA IC
            
            % Store events time stamps
            % dRipples
            ripples.dHPC.baseline = Restrict(ripplesD , NREM.baseline);
            ripples.dHPC.reward = Restrict(ripplesD , NREM.reward);
            ripples.dHPC.aversive = Restrict(ripplesD , NREM.aversive);
            % vRipples
            ripples.vHPC.baseline = Restrict(ripplesV , NREM.baseline);
            ripples.vHPC.reward = Restrict(ripplesV , NREM.reward);
            ripples.vHPC.aversive = Restrict(ripplesV , NREM.aversive);
            % coordinated dRipples
            ripples.dHPC.coordinated.all = coordinated;
            ripples.dHPC.uncoordinated.all = ripplesD(not(ismember(ripplesD(:,2) , coordinated(:,2))),:);
            ripples.dHPC.coordinated.baseline = Restrict(coordinated , NREM.baseline);
            ripples.dHPC.coordinated.reward = Restrict(coordinated , NREM.reward);
            ripples.dHPC.coordinated.aversive = Restrict(coordinated , NREM.aversive);
            % coordinated vRipples
            ripples.vHPC.coordinated.all = coordinatedV;
            ripples.vHPC.uncoordinated.all = ripplesV(not(ismember(ripplesV(:,2) , coordinatedV(:,2))),:);
            ripples.vHPC.coordinated.baseline = Restrict(coordinatedV , NREM.baseline);
            ripples.vHPC.coordinated.reward = Restrict(coordinatedV , NREM.reward);
            ripples.vHPC.coordinated.aversive = Restrict(coordinatedV , NREM.aversive);
            %coordinated event
            cooridnated_event((cooridnated_event(:,3)-cooridnated_event(:,1)<0.04),:) = [];
            ripple_event.baseline = Restrict(cooridnated_event,baselineTS);
            ripple_event.reward = Restrict(cooridnated_event,rewardTS);
            ripple_event.aversive = Restrict(cooridnated_event,aversiveTS);
            ripple_event.all = cooridnated_event;
        else
            RB = false;
        end
        
        %% Spikes
        % Load Units
        disp('Uploading Spiking activity')
        cd 'Spikesorting'
        [clusters , numberD , numberV , spks , spks_dHPC , spks_vHPC , cellulartype] = load_SU_FM(cd,criteria_type,criteria_fr,aversiveTS_run,rewardTS_run);
        
        if isfile('dorsalventral_assemblies_aversiveVF.mat')
            iterator = true;
            disp('Loading Aversive template')
            [patterns , cond , Thresholded] = load_assemblies(cd , 'dorsalventral_assemblies_aversiveVF.mat', clusters, numberD , 'aversive');
            th = sum(Thresholded.aversive(:,cond.both),2)>=1;
            patterns = patterns.aversive(:,cond.both);
            Thresholded = Thresholded.aversive(:,cond.both);
        else
            th = zeros(size(clusters.all,1),1);
            iterator = false;
        end
        
        if exist('dHPC_responsivness_all.mat')
            load('dHPC_responsivness_all.mat')
            Responssive.dHPC = [];
            for i = 1 : size(clusters.dHPC,1)
                Responssive.dHPC = [Responssive.dHPC ; dHPC_resp.resp_ave(dHPC_resp.id == clusters.dHPC(i))==1];
            end
        end
        
        if exist('vHPC_responsivness_all.mat')
            load('vHPC_responsivness_all.mat')
            Responssive.vHPC = [];
            for i = 1 : size(clusters.vHPC,1)
                Responssive.vHPC = [Responssive.vHPC ; vHPC_resp.resp_ave(vHPC_resp.id == clusters.vHPC(i))==1];
            end
        end
        
        if and(RB , iterator)
            % SpikeTrains construction
            limits = [0 segments.Var1(end)/1000];
            events = [];
            [Spikes , bins , Clusters] = spike_train_construction([spks_dHPC;spks_vHPC], clusters.all, cellulartype, 0.025, limits, events, false, true);
%             clear limits events
            
            a = assembly_activity(patterns , Spikes');
            
            for i = 1 : size(patterns,2)
                template = logical([Thresholded(:,i) [Responssive.dHPC ;  zeros(size(Responssive.vHPC,1),1)]]);
%                 [metrics_with , metrics_without , comparison] = train_test_poisson_glm([bins' , Spikes], template , [bins' a(i,:)'], movement.aversive, [ripple_event.aversive(:,1) ripple_event.aversive(:,3)], 'n_permutations', 1000, 'compute_ci', true);
                [metrics_with , metrics_without , comparison] = train_test_poisson_glm([bins' , Spikes], template , [bins' a(i,:)'], movement.aversive, [ripple_event.aversive(:,1) ripple_event.aversive(:,3)]);
                Pred_with.dHPC = [Pred_with.dHPC ; metrics_with.Pearson_r , metrics_with.Spearman_r , metrics_with.MSE];
                Pred_without.dHPC = [Pred_without.dHPC ; metrics_without.Pearson_r , metrics_without.Spearman_r , metrics_without.MSE];
                clear test_pred_with test_pred_without train_models metrics_with metrics_without
                
                template = logical([Thresholded(:,i) [zeros(size(Responssive.dHPC,1),1) ;  Responssive.vHPC]]);
%                 [metrics_with , metrics_without , comparison] = train_test_poisson_glm([bins' , Spikes], template , [bins' a(i,:)'], movement.aversive, [ripple_event.aversive(:,1) ripple_event.aversive(:,3)], 'n_permutations', 1000, 'compute_ci', true);
                [metrics_with , metrics_without , comparison] = train_test_poisson_glm([bins' , Spikes], template , [bins' a(i,:)'], movement.aversive, [ripple_event.aversive(:,1) ripple_event.aversive(:,3)]);
                Pred_with.vHPC = [Pred_with.vHPC ; metrics_with.Pearson_r , metrics_with.Spearman_r , metrics_with.MSE];
                Pred_without.vHPC = [Pred_without.vHPC ; metrics_without.Pearson_r , metrics_without.Spearman_r , metrics_without.MSE];
                clear test_pred_with test_pred_without train_models metrics_with metrics_without
            end
        end
        disp(' ')
        clear A aversiveTS aversiveTS_run baselineTS rewardTS rewardTS_run
        clear behavior bins Cell_type_classification cellulartype cond
        clear is K Kinfo group_dHPC group_vHPC matrixC matrixD matrixV
        clear NREM REM WAKE segmentation tmp cond config
        clear spiketrains_dHPC spiketrains_vHPC opts MUA
        clear patterns Thresholded i  ii numberD numberV movement cross crossN
        clear Spikes bins Clusters Shocks_filt Rewards_filt config n_SU_D n_SU_V
        clear clusters coordinated coordinated_ripple_bursts coordinatedV
        clear cooridnated_event coordinatedV_refined coordinatedV_refined
        clear ripple_bursts ripple_event ripplesD ripplesV
        clear spks spks_dHPC spks_vHPC ripples cooridnated_event
        clear cooridnated_eventDV cooridnated_eventVD segments movement RD RV RB
    end
%     
% 
% y = (abs([Pred_with.dHPC(:,2)]) - abs([Pred_without.dHPC(:,2)]))./(abs([Pred_with.dHPC(:,2)]) + abs([Pred_without.dHPC(:,2)]));
% y1 = (abs([Pred_with.vHPC(:,2)]) - abs([Pred_without.vHPC(:,2)]))./(abs([Pred_with.vHPC(:,2)]) + abs([Pred_without.vHPC(:,2)]));
% grps = [ones(size(y,1),1) ; ones(size(y1,1),1)*2];
% figure,boxplot([y ; y1] , grps)
% 
% [h p] = signtest(y1)
% [h p] = signtest(y)
% [] = signrank(y,y1)
% 
% 
% 
% plot_model_comparison([y , y1])
%     
end



function plot_random_peak_window(spike_train, assembly_activity, color_split, responsive_neurons, special_neurons)
    % Parameters
    nNeurons = size(spike_train, 2);
    bin_size = 0.025; % 25ms bins
    fs = 1/bin_size; % 40 Hz sampling
    window_radius = 5; % 5 seconds before/after peak
    min_peak_dist = round(2/bin_size); % Minimum 2s between peaks (80 bins)
    
    % Validate marker vectors
    if length(responsive_neurons) ~= nNeurons || length(special_neurons) ~= nNeurons
        error('Marker vectors must be logical vectors of length %d', nNeurons);
    end
    
    % Find all significant peaks
    [pks, locs] = findpeaks(assembly_activity, ...
        'MinPeakHeight', mean(assembly_activity) + std(assembly_activity), ...
        'MinPeakDistance', min_peak_dist);
    
    % Randomly select one peak
    if isempty(locs)
        error('No significant peaks found in assembly activity');
    end
    selected_idx = randi(length(locs));
    peak_loc = locs(selected_idx);
    
    % Get time window (5s before to 5s after peak)
    time_window = max(1, peak_loc-round(fs*window_radius)) : ...
                 min(size(spike_train,1), peak_loc+round(fs*window_radius)-1);
    time_sec = (-window_radius:bin_size:(length(time_window)-1)*bin_size-window_radius)';
    
    % --- Plot 1: Neurons with Dual Markers ---
    figure('Position', [100 100 1200 600]);
    hold on;
    
    % Create custom y-tick labels with * and #
    neuron_labels = cell(1,nNeurons);
    for i = 1:nNeurons
        markers = '';
        if responsive_neurons(i), markers = [markers '*']; end
        if special_neurons(i), markers = [markers '#']; end
        neuron_labels{i} = sprintf('%d%s', i, markers);
    end
    
    % Plot neurons with style variations
    offsets = linspace(0, 10, nNeurons);
    for i = 1:nNeurons
        neuron_data = spike_train(time_window, i);
        if max(neuron_data) > min(neuron_data)
            norm_data = (neuron_data - min(neuron_data)) / (max(neuron_data) - min(neuron_data));
        else
            norm_data = zeros(size(neuron_data));
        end
        
        % Base color
        if i <= color_split
            line_color = 'k'; % Black for first group
        else
            line_color = 'r'; % Red for second group
        end
        
        % Line style customization
        if responsive_neurons(i) && special_neurons(i)
            line_style = '-'; % Both markers: dotted line
            line_width = 1.5;
        elseif responsive_neurons(i)
            line_style = '-'; % Responsive only: solid line
            line_width = 1.2;
        elseif special_neurons(i)
            line_style = '-'; % Special only: dashed line
            line_width = 1.0;
        else
            line_style = '-';
            line_width = 0.8;
        end
        
        plot(time_sec, norm_data + offsets(i), 'Color', line_color, ...
            'LineStyle', line_style, 'LineWidth', line_width);
    end
    
    % Mark peak time
    plot([0 0], [0 11], 'g--', 'LineWidth', 1.5);
    
    % Formatting
    xlabel('Time relative to peak (s)');
    ylabel('Neuron ID (*=responsive, #=special)');
    yticks(offsets);
    yticklabels(neuron_labels);
    title(sprintf('Neural Activity Around Random Peak (#%d of %d)\nBlack: %d | Red: %d | *: Responsive | #: Special', ...
        selected_idx, length(locs), color_split, nNeurons-color_split));
    grid on;
    xlim([-window_radius window_radius]);
    
    % --- Plot 2: Assembly Activity ---
    figure('Position', [100 100 1200 300]);
    plot(time_sec, assembly_activity(time_window), 'b', 'LineWidth', 2);
    hold on;
    plot([0 0], ylim, 'g--', 'LineWidth', 1.5);
    xlabel('Time relative to peak (s)');
    ylabel('Assembly Activity');
    title(sprintf('Assembly Activity (Peak at t=%.1fs)', peak_loc*bin_size));
    grid on;
    xlim([-window_radius window_radius]);
end