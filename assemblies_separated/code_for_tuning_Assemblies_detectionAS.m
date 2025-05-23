clear
clc
% close all

%% Parameters
path = {'E:\Rat126\Ephys\in_Pyr';'E:\Rat103\usable';'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat132\recordings\in_pyr';'E:\Rat165\in_pyr\'};%List of folders from the path

% for SU
criteria_fr = 0; %criteria to include or not a SU into the analysis
criteria_n = [3 3]; % minimal number of neurons from each structure [vHPC dHPC]
criteria_type = 0; %criteria for celltype (0:pyr, 1:int, 2:all)
binSize = [0.01]; %for qssemblie detection qnd qxctivity strength

% Behavior
minimal_speed = 7; % minimal speed to detect quite periods
minimal_speed_time = 2; % minimal time to detect quite periods

percentages = [];

Number_of_assemblies.aversive = [];
Number_of_assemblies.reward = [];

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
        num_assembliesR = [];
        num_assembliesA = [];
        percentagesA = [];
        
        disp(['-- Initiating analysis of folder #' , num2str(t) , ' from rat #',num2str(tt) , ' --'])
        session = [subFolders(t+2).folder,'\',subFolders(t+2).name];
        cd(session)
        
        %Loading TS of the sessions
        disp('Uploading session time stamps')
        load('session_organization.mat')
        
        %% Awake
        disp('Uploading behavioral data')
        load('behavioral_data.mat')

        %% Spikes
        % Load Units
        disp('Uploading Spiking activity')
        cd 'Spikesorting'
        spks = double([readNPY('spike_clusters.npy') readNPY('spike_times.npy')]);
        K = tdfread('cluster_group.tsv'); % import clusters ID and groups (MUA,Noise,Good)
        Kinfo = tdfread('cluster_info.tsv'); % import information of clusters
        K = [K.cluster_id(K.group(:,1) == 'g') , Kinfo.ch(K.group(:,1) == 'g'),]; % defining good clusters
        % Load neuronal classification
        load('Cell_type_classification')
        K = [K , Cell_type_classification(:,6:8)];
        group_dHPC = K(K(:,2) > 63,:);
        group_vHPC = K(K(:,2) <= 63,:);
        
        %Loop to select dorsal or ventral LFP and SU
        % z=1 --> dorsal
        % z=2 --> ventral
        for z = 1:2
            if z == 1
                spks_dHPC = spks(ismember(spks(:,1),group_dHPC(:,1)),:); %keep spks from good clusters
            else
                spks_vHPC = spks(ismember(spks(:,1),group_vHPC(:,1)),:); %keep spks from good clusters
            end
        end
        clear z
        spks_vHPC(:,2) = double(spks_vHPC(:,2))./20000;
        spks_dHPC(:,2) = double(spks_dHPC(:,2))./20000;
        spks(:,2) = double(spks(:,2))./20000;
        
        % Selection of celltype to analyze
        if criteria_type == 0 %pyr
            cellulartype = [K(:,1) , K(:,4)];
        elseif criteria_type == 1 % int
            cellulartype = [K(:,1) , not(K(:,4))];
        elseif criteria_type == 2 % all
            cellulartype = [K(:,1) , ones(length(K),1)];
        end
        
        %% Counting the Number f SU
        numberD = 0;
        clusters.all = [];
        clusters.dHPC = [];
        for ii=1:size(group_dHPC,1)
            cluster = group_dHPC(ii,1);
            Cellulartype = logical(cellulartype(cellulartype(:,1) == cluster,2));
            if Cellulartype
                a = length(Restrict(spks(spks(:,1)==cluster,2),aversiveTS_run./1000)) / ((aversiveTS_run(2)-aversiveTS_run(1))/1000);
                r = length(Restrict(spks(spks(:,1)==cluster,2),rewardTS_run./1000)) / ((rewardTS_run(2)-rewardTS_run(1))/1000);
                if or(a > criteria_fr ,  r > criteria_fr)
                    numberD = numberD+1;
                    clusters.all = [clusters.all ; cluster];
                    clusters.dHPC = [clusters.dHPC ; cluster];
                end
            end
            clear tmp cluster Cellulartype fr1 fr2 fr3 fr4 fr5 r a
        end
        
        numberV = 0;
        clusters.vHPC = [];
        for ii=1:size(group_vHPC,1)
            cluster = group_vHPC(ii,1);
            Cellulartype = logical(cellulartype(cellulartype(:,1) == cluster,2));
            if Cellulartype
                a = length(Restrict(spks(spks(:,1)==cluster,2),aversiveTS_run./1000)) / ((aversiveTS_run(2)-aversiveTS_run(1))/1000);
                r = length(Restrict(spks(spks(:,1)==cluster,2),rewardTS_run./1000)) / ((rewardTS_run(2)-rewardTS_run(1))/1000);
                if or(a > criteria_fr ,  r > criteria_fr)
                    numberV = numberV+1;
                    clusters.all = [clusters.all ; cluster];
                    clusters.vHPC = [clusters.vHPC ; cluster];
                end
            end
            clear tmp cluster Cellulartype fr1 fr2 fr3 fr4 fr5 r a
        end
        clear freq limits
        clear camara shock rightvalve leftvalve
        clear ejeX ejeY dX dY dX_int dY_int
        
        %% Assemblies detection
        if or(numberD > 3 , numberV > 3)
            for i = 1 : size(binSize,2)
            % --- Aversive ---
            disp('Lets go for the assemblies')
%             if binSize(i) == 0.025
%                 disp('Loading Aversive template')
%                 load('dorsalventral_assemblies_aversive.mat')
%             else
                disp('Detection of assemblies using Aversive template')
                % --- Options for assemblies detection ---
                opts.Patterns.method = 'ICA';
                opts.threshold.method= 'MarcenkoPastur';
                opts.Patterns.number_of_iterations= 500;
                opts.threshold.permutations_percentile = 0.9;
                opts.threshold.number_of_permutations = 500;
                opts.Patterns.number_of_iterations = 500;
                opts.Members.method = 'Sqrt';
                
                limits = aversiveTS_run./1000;
                events = [];
                events = movement.aversive;
                
                if not(isempty(clusters.dHPC))
                    [SpksTrains.dHPC.aversive , Bins.aversive , Cluster.aversive.dHPC] = spike_train_construction(spks_dHPC, clusters.dHPC, cellulartype, binSize(i), limits, events, false,false);
                    [Th , pat , eig] = assembly_patternsJFM([SpksTrains.dHPC.aversive'],opts);
                    % save
                    Thresholded.aversive.dHPC = Th;
                    patterns.aversive.dHPC = pat;
                    clear cond Th pat
                    ND = size(patterns.aversive.dHPC,2);
                else
                    ND = 0;
                end
                
                if not(isempty(clusters.vHPC))
                    [SpksTrains.vHPC.aversive , Bins.aversive , Cluster.aversive.vHPC] = spike_train_construction(spks_vHPC, clusters.vHPC, cellulartype, binSize(i), limits, events, false,false);
                    [Th , pat , eig] = assembly_patternsJFM([SpksTrains.vHPC.aversive'],opts);
                    % save
                    Thresholded.aversive.vHPC = Th;
                    patterns.aversive.vHPC = pat;
                    clear cond Th pat
                    NV = size(patterns.aversive.vHPC,2);
                else
                    NV = 0;
                end
%                 save([cd,'\dorsalventral_assemblies_aversiveVF.mat'],'Th' , 'pat' , 'eig' , 'criteria_fr' , 'criteria_n')
%             end
            num_assembliesA = [num_assembliesA , ND NV ]; clear ND NV
            
            % --- Reward ---
            disp('Loading Reward template')
%             if binSize(i) == 0.025
%                 load('dorsalventral_assemblies_reward.mat')
%             else
                disp('Detection of assemblies using Rewarded template')
                % --- Options for assemblies detection ---
                % --- Options for assemblies detection ---
                opts.Patterns.method = 'ICA';
                opts.threshold.method= 'MarcenkoPastur';
                opts.Patterns.number_of_iterations= 500;
                opts.threshold.permutations_percentile = 0.9;
                opts.threshold.number_of_permutations = 500;
                opts.Patterns.number_of_iterations = 500;
                opts.Members.method = 'Sqrt';
                
                limits = rewardTS_run./1000;
                events = [];
                events = movement.reward;
                
                if not(isempty(clusters.dHPC))
                    [SpksTrains.dHPC.reward , Bins.reward , Cluster.reward.dHPC] = spike_train_construction(spks_dHPC, clusters.dHPC, cellulartype, binSize(i), limits, events, false,false);
                    [Th , pat , eig] = assembly_patternsJFM([SpksTrains.dHPC.aversive'],opts);
                    % save
                    Thresholded.reward.dHPC = Th;
                    patterns.reward.dHPC = pat;
                    clear cond Th pat
                    ND = size(patterns.reward.dHPC,2);
                else
                    ND = 0;
                end
                
                if not(isempty(clusters.vHPC))
                    [SpksTrains.vHPC.reward , Bins.reward , Cluster.reward.vHPC] = spike_train_construction(spks_vHPC, clusters.vHPC, cellulartype, binSize(i), limits, events, false,false);
                    [Th , pat , eig] = assembly_patternsJFM([SpksTrains.vHPC.reward'],opts);
                    % save
                    Thresholded.reward.vHPC = Th;
                    patterns.reward.vHPC = pat;
                    clear cond Th pat
                    NV = size(patterns.reward.vHPC,2);
                else
                    NV = 0;
                end
                save([cd,'\assemblies_separated.mat'],'Th' , 'pat' , 'eig' , 'criteria_fr' , 'criteria_n')
%             end
            
            num_assembliesR = [num_assembliesR , ND NV]; clear ND NV
            
            %% Similarity Index Calculation
            if and(isfield(patterns.aversive,'dHPC') , isfield(patterns.reward,'dHPC'))
                [r.AR , p.AR] = SimilarityIndex(patterns.aversive.dHPC , patterns.reward.dHPC);
                A = sum(p.AR,1)>=1;
                R = sum(p.AR,2)>=1; R = R';
                AR.dHPC = [sum(sum(p.AR)) size(patterns.aversive.dHPC,2) size(patterns.reward.dHPC,2)];
            else
                AR.dHPC = [0 0 0];
            end
            
            if and(isfield(patterns.aversive,'vHPC') , isfield(patterns.reward,'vHPC'))
                [r.AR , p.AR] = SimilarityIndex(patterns.aversive.vHPC , patterns.reward.vHPC);
                A = sum(p.AR,1)>=1;
                R = sum(p.AR,2)>=1; R = R';
                AR.vHPC = [sum(sum(p.AR)) size(patterns.aversive.vHPC,2) size(patterns.reward.vHPC,2)];
            else
                AR.vHPC = [0 0 0];
            end            
            
            percentagesA = [percentagesA  , AR.dHPC  AR.vHPC ];
            clear A R r p
            clear patterns Thresholded cond 
            
            end
            Number_of_assemblies.aversive = [Number_of_assemblies.aversive ; num_assembliesA];
            Number_of_assemblies.reward = [Number_of_assemblies.reward ; num_assembliesR];
            percentages = [percentages ; percentagesA]; clear percentagesA
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
        clear cooridnated_eventDV cooridnated_eventVD segments
    end
        clear num_assembliesA num_assembliesR
    
end

figure
subplot(221),bar([1 2 3 4] , [sum(Number_of_assemblies.aversive(:,1)) sum(Number_of_assemblies.aversive(:,3)) sum(Number_of_assemblies.aversive(:,5)) sum(Number_of_assemblies.aversive(:,7))]),ylim([0 100])
subplot(222),bar([1 2 3 4] , [sum(Number_of_assemblies.aversive(:,2)) sum(Number_of_assemblies.aversive(:,4)) sum(Number_of_assemblies.aversive(:,6)) sum(Number_of_assemblies.aversive(:,8))]),ylim([0 150])

subplot(223),bar([1 2 3 4] , [sum(Number_of_assemblies.reward(:,1)) sum(Number_of_assemblies.reward(:,3)) sum(Number_of_assemblies.reward(:,5)) sum(Number_of_assemblies.reward(:,7))]),ylim([0 100])
subplot(224),bar([1 2 3 4] , [sum(Number_of_assemblies.reward(:,2)) sum(Number_of_assemblies.reward(:,4)) sum(Number_of_assemblies.reward(:,6)) sum(Number_of_assemblies.reward(:,8))]),ylim([0 150])


figure
subplot(121),bar([1 2 3 4] , [sum(percentages(:,1)) sum(percentages(:,7)) sum(percentages(:,13)) sum(percentages(:,19))]),ylim([0 20])
subplot(122),bar([1 2 3 4] , [sum(percentages(:,4)) sum(percentages(:,10)) sum(percentages(:,16)) sum(percentages(:,22))]),ylim([0 60])


ID1 = sum(percentages(:,1) + percentages(:,2) + percentages(:,3))./(Number_of_assemblies.aversive(:,1) + Number_of_assemblies.aversive(:,2) + Number_of_assemblies.aversive(:,3) + Number_of_assemblies.reward(:,1) + Number_of_assemblies.reward(:,2) +  Number_of_assemblies.reward(:,3));
ID2 = sum(percentages(:,4) + percentages(:,5) + percentages(:,6))./(Number_of_assemblies.aversive(:,4) + Number_of_assemblies.aversive(:,5) + Number_of_assemblies.aversive(:,6) + Number_of_assemblies.reward(:,4) + Number_of_assemblies.reward(:,5) +  Number_of_assemblies.reward(:,6));
ID3 = sum(percentages(:,7) + percentages(:,8) + percentages(:,9))./(Number_of_assemblies.aversive(:,7) + Number_of_assemblies.aversive(:,8) + Number_of_assemblies.aversive(:,9) + Number_of_assemblies.reward(:,7) + Number_of_assemblies.reward(:,8) +  Number_of_assemblies.reward(:,9));
ID4 = sum(percentages(:,10) + percentages(:,11) + percentages(:,12))./(Number_of_assemblies.aversive(:,10) + Number_of_assemblies.aversive(:,11) + Number_of_assemblies.aversive(:,10) + Number_of_assemblies.reward(:,10) + Number_of_assemblies.reward(:,11) +  Number_of_assemblies.reward(:,12));


IDs = [ID1 ; ID2 ; ID3 ; ID4];
grps = [ones(length(ID1),1) ; ones(length(ID2),1)*2 ; ones(length(ID3),1)*3 ; ones(length(ID4),1)*4];

figure,
boxplot(IDs,grps),hold on
kruskalwallis(IDs,grps)

figure
subplot(4,3,1),histogram(Number_of_assemblies.aversive(:,1),'BinEdges',[0:1:8]),ylim([0 30])
subplot(4,3,2),histogram(Number_of_assemblies.aversive(:,2),'BinEdges',[0:1:8]),ylim([0 30])
subplot(4,3,3),histogram(Number_of_assemblies.aversive(:,3),'BinEdges',[0:1:8]),ylim([0 30])

subplot(4,3,4),histogram(Number_of_assemblies.aversive(:,4),'BinEdges',[0:1:8]),ylim([0 30])
subplot(4,3,5),histogram(Number_of_assemblies.aversive(:,5),'BinEdges',[0:1:8]),ylim([0 30])
subplot(4,3,6),histogram(Number_of_assemblies.aversive(:,6),'BinEdges',[0:1:8]),ylim([0 30])

subplot(4,3,7),histogram(Number_of_assemblies.aversive(:,7),'BinEdges',[0:1:8]),ylim([0 30])
subplot(4,3,8),histogram(Number_of_assemblies.aversive(:,8),'BinEdges',[0:1:8]),ylim([0 30])
subplot(4,3,9),histogram(Number_of_assemblies.aversive(:,9),'BinEdges',[0:1:8]),ylim([0 30])

subplot(4,3,10),histogram(Number_of_assemblies.aversive(:,10),'BinEdges',[0:1:8]),ylim([0 30])
subplot(4,3,11),histogram(Number_of_assemblies.aversive(:,11),'BinEdges',[0:1:8]),ylim([0 30])
subplot(4,3,12),histogram(Number_of_assemblies.aversive(:,12),'BinEdges',[0:1:8]),ylim([0 30])