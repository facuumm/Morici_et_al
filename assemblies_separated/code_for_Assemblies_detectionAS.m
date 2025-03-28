clear
clc
% close all

%% Parameters
path = {'E:\Rat126\Ephys\in_Pyr';'E:\Rat103\usable';'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat132\recordings\in_pyr';'E:\Rat165\in_pyr\'};%List of folders from the path

% for SU
criteria_fr = 0; %criteria to include or not a SU into the analysis
criteria_n = [3 3]; % minimal number of neurons from each structure [vHPC dHPC]
criteria_type = 0; %criteria for celltype (0:pyr, 1:int, 2:all)
binSize = [0.015]; %for assemblie detection and activity strength

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
                
                if and(not(isempty(clusters.dHPC)) , numberD > 3)
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
                
                if and(not(isempty(clusters.vHPC)) , numberV > 3)
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
                
                if and(not(isempty(clusters.dHPC)) , numberD > 3)
                    [SpksTrains.dHPC.reward , Bins.reward , Cluster.reward.dHPC] = spike_train_construction(spks_dHPC, clusters.dHPC, cellulartype, binSize(i), limits, events, false,false);
                    [Th , pat , eig] = assembly_patternsJFM([SpksTrains.dHPC.reward'],opts);
                    % save
                    Thresholded.reward.dHPC = Th;
                    patterns.reward.dHPC = pat;
                    clear cond Th pat
                    ND = size(patterns.reward.dHPC,2);
                else
                    ND = 0;
                end
                
                if and(not(isempty(clusters.vHPC)) , numberV > 3)
                    [SpksTrains.vHPC.reward , Bins.reward , Cluster.reward.vHPC] = spike_train_construction(spks_vHPC, clusters.vHPC, cellulartype, binSize(i), limits, events, false,false);
                    [Th , pat , eig] = assembly_patternsJFM([SpksTrains.vHPC.reward'],opts);
                    % save
                    Thresholded.reward.vHPC = Th;
                    patterns.reward.vHPC = pat;
                    
                    NV = size(patterns.reward.vHPC,2);
                else
                    NV = 0;
                end
                save([cd,'\separated_assemblies.mat'],'Thresholded' , 'patterns' , 'criteria_fr' , 'criteria_n'),clear cond Th pat
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


%% Ven Graphs for similarity Index
% dHPC Assemblies
p1 = (sum(percentages(:,2))./sum(sum(percentages(:,2:3))))*100;
p2 = (sum(percentages(:,3))./sum(sum(percentages(:,2:3))))*100;

A = [ p1 p2 ];
I = (sum(percentages(:,1))./sum(sum(percentages(:,2:3))))*100;

subplot(121),venn(A,I), xlim([-7 11]), ylim([-6 6])

% vHPC Assemblies
p1 = (sum(percentages(:,5))./sum(sum(percentages(:,5:6))))*100;
p2 = (sum(percentages(:,6))./sum(sum(percentages(:,5:6))))*100;

A = [ p1 p2 ];
I = (sum(percentages(:,4))./sum(sum(percentages(:,5:6))))*100;

subplot(122),venn(A,I), xlim([-7 11]), ylim([-6 6])

%%

figure
subplot(121),bar([1 2] , [sum(Number_of_assemblies.aversive(:,1)) sum(Number_of_assemblies.aversive(:,2))])

subplot(122),bar([1 2] , [sum(Number_of_assemblies.reward(:,1)) sum(Number_of_assemblies.reward(:,2))])


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