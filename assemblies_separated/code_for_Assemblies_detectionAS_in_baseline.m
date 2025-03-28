clear
clc
% close all

%% Parameters
path = {'E:\Rat126\Ephys\in_Pyr';'E:\Rat103\usable';'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat132\recordings\in_pyr';'E:\Rat165\in_pyr\'};%List of folders from the path

% path = {'\\Maryjackson\e\Rat127\Ephys\pyr';'\\Maryjackson\e\Rat128\Ephys\in_pyr\ready';'\\Maryjackson\e\Rat103\usable';'\\Maryjackson\e\Rat132\recordings\in_pyr'; '\\Maryjackson\e\Rat165\in_pyr'; ...
%     '\\Maryjackson\e\Rat126\Ephys\in_Pyr'};%List of folders from the path

% for SU
criteria_fr = 0; %criteria to include or not a SU into the analysis
criteria_n = [3 3]; % minimal number of neurons from each structure [vHPC dHPC]
criteria_type = 0; %criteria for celltype (0:pyr, 1:int, 2:all)
binSize = [0.015]; %for assemblie detection and activity strength

% Behavior
minimal_speed = 7; % minimal speed to detect quite periods
minimal_speed_time = 2; % minimal time to detect quite periods

percentages = [];

Number_of_assemblies = [];

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
        num_assemblies = [];
        
        disp(['-- Initiating analysis of folder #' , num2str(t) , ' from rat #',num2str(tt) , ' --'])
        session = [subFolders(t+2).folder,'\',subFolders(t+2).name];
        cd(session)
        
        %Loading TS of the sessions
        disp('Uploading session time stamps')
        load('session_organization.mat')
        
        %% Awake
        disp('Uploading behavioral data')
        load('behavioral_data.mat')
        
        %% Homecage data
        disp('Uploading sleep scoring')
        x = dir([cd,'\*-states.mat']);    states = load([cd,'\',x.name]);    states = states.states;
        REM.all = ToIntervals(states==5);    NREM.all = ToIntervals(states==3);    WAKE.all = ToIntervals(states==1);
        clear x states
        WAKE.baseline = Restrict(WAKE.all,baselineTS./1000);
        NREM.baseline = Restrict(NREM.all,baselineTS./1000);
        
        %% Spikes
        % Load Units
        disp('Uploading Spiking activity')
        cd 'Spikesorting'
        [clusters , numberD , numberV , spks , spks_dHPC , spks_vHPC , cellulartype] = load_SU_FM(cd,criteria_type,criteria_fr,aversiveTS_run./1000,rewardTS_run./1000);
        
        
        %% Assemblies detection
        if or(numberD > 3 , numberV > 3)
            for i = 1 : size(binSize,2)
                % --- Aversive ---
                disp('Lets go for the assemblies')
                %             if binSize(i) == 0.025
                %                 disp('Loading Aversive template')
                %                 load('dorsalventral_assemblies_aversive.mat')
                %             else
                disp('Detection of assemblies')
                % --- Options for assemblies detection ---
                opts.Patterns.method = 'ICA';
                opts.threshold.method= 'MarcenkoPastur';
                opts.Patterns.number_of_iterations= 500;
                opts.threshold.permutations_percentile = 0.9;
                opts.threshold.number_of_permutations = 500;
                opts.Patterns.number_of_iterations = 500;
                opts.Members.method = 'Sqrt';
                
                limits = [0 segments.Var1(end)]./1000;
                events = [];
                events = [NREM.baseline];
%                 events = sort([aversiveTS_run ; rewardTS_run]./1000);
%                 events = sort([movement.aversive ; movement.reward]);

                if and(not(isempty(clusters.dHPC)) , numberD > 3)
                    [SpksTrains.dHPC , Bins , Cluster.aversive.dHPC] = spike_train_construction(spks_dHPC, clusters.dHPC, cellulartype, binSize(i), limits, events, false,false);
                    [Th , pat , eig] = assembly_patternsJFM([SpksTrains.dHPC'],opts);
                    % save
                    Thresholded.dHPC = Th;
                    patterns.dHPC = pat;
                    clear cond Th pat
                    ND = size(patterns.dHPC,2);
                else
                    ND = 0;
                end
                
                if and(not(isempty(clusters.vHPC)) , numberV > 3)
                    [SpksTrains.vHPC , Bins , Cluster.aversive.vHPC] = spike_train_construction(spks_vHPC, clusters.vHPC, cellulartype, binSize(i), limits, events, false,false);
                    [Th , pat , eig] = assembly_patternsJFM([SpksTrains.vHPC'],opts);
                    % save
                    Thresholded.vHPC = Th;
                    patterns.vHPC = pat;
                    clear cond Th pat
                    NV = size(patterns.vHPC,2);
                else
                    NV = 0;
                end
                %             end
                num_assemblies = [num_assemblies , ND NV ]; clear ND NV
                
                
                save([cd,'\separated_assemblies_NREMBaseline.mat'],'Thresholded' , 'patterns' , 'criteria_fr' , 'criteria_n'),clear cond Th pat
                %             end
                
            end
            Number_of_assemblies = [Number_of_assemblies ; num_assemblies];
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
        clear cooridnated_eventDV cooridnated_eventVD segments events
    end
    clear num_assembliesA num_assembliesR
end



figure
bar([1 2] , [sum(Number_of_assemblies(:,1)) sum(Number_of_assemblies(:,2))])
