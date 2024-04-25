clear
clc
% close all

%% Parameters
path = {'E:\Rat126\Ephys\in_Pyr';'E:\Rat103\usable';'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat132\recordings\in_pyr';'E:\Rat165\in_pyr\'};%List of folders from the path

% for SU
criteria_fr = 0; %criteria to include or not a SU into the analysis
criteria_n = [3 3]; % minimal number of neurons from each structure [vHPC dHPC]
criteria_type = 0; %criteria for celltype (0:pyr, 1:int, 2:all)
binSize = [0.005 0.01 0.025 0.05 ]; %for qssemblie detection qnd qxctivity strength

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
        x = dir([cd,'\*.cat.evt']);
        segments = readtable([cd,'\',x.name],'FileType','text');
        clear x
        % TimeStamps of begening and end of the sleep and awake trials
        % Reward and Aversive trials
        aversiveTS = [];
        aversiveTS_run = [];
        rewardTS = [];
        rewardTS_run = [];
        for y = 1 : height(segments)
            % Baseline sleep session TS detection
            if y == 1
                baselineTS(1,1) = segments.Var1(y);
            elseif y ==2
                baselineTS(1,2) = segments.Var1(y);
            end
            % Aversive sleep session TS detection
            if strcmp(segments.Var2{y},'aversive')
                if strcmp(segments.Var3{y},'End')
                    aversiveTS(1,1) = segments.Var1(y+1);
                    aversiveTS(1,2) = segments.Var1(y+2);
                    aversiveTS_run(1,1) = segments.Var1(y-1);
                    aversiveTS_run(1,2) = segments.Var1(y);
                    A = y;
                end
                % Rewarded sleep session TS detection
            elseif strcmp(segments.Var2{y},'reward')
                if strcmp(segments.Var3{y},'End')
                    rewardTS(1,1) = segments.Var1(y+1);
                    rewardTS(1,2) = segments.Var1(y+2);
                    rewardTS_run(1,1) = segments.Var1(y-1);
                    rewardTS_run(1,2) = segments.Var1(y);
                    R = y;
                end
            end
        end
        clear y A R
        
        % Defining what condition was first
        if aversiveTS_run(1) < rewardTS_run(1)
            config = 1;
        else
            config = 2;
        end
        
        %% Awake
        disp('Uploading digital imputs')
        % Load digitalin.mat
        load('digitalin.mat')
        
        % Behavioral calculations
        disp('Uploading DLC outputs')
        camara = ((camara(:,2)-camara(:,1))/2)+camara(:,1);
        % periods of movment during eacj condition
        if rewardTS_run(1) < aversiveTS_run(1)
            load('laps1.mat','posx','posy');
            [camaraR,~] = find((camara(:,1)-rewardTS_run(1)/1000)>0,1,'first'); %TimeStamp of the begining of aversive
            pos = [camara(camaraR : camaraR+length(posx)-1),posx,posy];
            %interpolation of dropped frames
            ejeX = pos(~isnan(pos(:,2)),1); dX = pos(~isnan(pos(:,2)),2); dX_int = interp1(ejeX , dX , pos(:,1));
            ejeY = pos(~isnan(pos(:,3)),1); dY = pos(~isnan(pos(:,3)),3); dY_int = interp1(ejeY , dY , pos(:,1));
            pos(:,2) =dX_int; pos(:,3) =dY_int;%saving corrected pos
            behavior.pos.reward = [pos];
            behavior.speed.reward = LinearVelocity(pos,0);
            behavior.speed.reward(:,2) = smoothdata(behavior.speed.reward(:,2),'gaussian',1,'SamplePoints',behavior.speed.reward(:,1));
            behavior.quiet.reward = QuietPeriods( behavior.speed.reward , minimal_speed , minimal_speed_time);
            clear pos camaraR posx posy
            load('laps2.mat','posx','posy');
            [camaraA,~] = find((camara(:,1)-aversiveTS_run(1)/1000)>0,1,'first'); %TimeStamp of the begining of aversive
            pos = [camara(camaraA : camaraA+length(posx)-1),posx,posy];
            %interpolation of dropped frames
            ejeX = pos(~isnan(pos(:,2)),1); dX = pos(~isnan(pos(:,2)),2); dX_int = interp1(ejeX , dX , pos(:,1));
            ejeY = pos(~isnan(pos(:,3)),1); dY = pos(~isnan(pos(:,3)),3); dY_int = interp1(ejeY , dY , pos(:,1));
            pos(:,2) =dX_int; pos(:,3) =dY_int;%saving corrected pos
            behavior.pos.aversive = [pos];
            behavior.speed.aversive = LinearVelocity(pos,0);
            behavior.speed.aversive(:,2) = smoothdata(behavior.speed.aversive(:,2),'gaussian',1,'SamplePoints',behavior.speed.aversive(:,1));
            behavior.quiet.aversive = QuietPeriods(behavior.speed.aversive , minimal_speed , minimal_speed_time);
            clear pos camaraR
        else
            load('laps2.mat','posx','posy');
            [camaraR,~] = find((camara(:,1)-rewardTS_run(1)/1000)>0,1,'first'); %TimeStamp of the begining of aversive
            %         [camaraR2,~] = find((camara(:,1)-rewardTS_run(2)/1000)<0,1,'last'); %TimeStamp of the ending of aversive
            pos = [camara(camaraR : camaraR+length(posx)-1),posx,posy];
            %interpolation of dropped frames
            ejeX = pos(~isnan(pos(:,2)),1); dX = pos(~isnan(pos(:,2)),2); dX_int = interp1(ejeX , dX , pos(:,1));
            ejeY = pos(~isnan(pos(:,3)),1); dY = pos(~isnan(pos(:,3)),3); dY_int = interp1(ejeY , dY , pos(:,1));
            pos(:,2) =dX_int; pos(:,3) =dY_int;%saving corrected pos
            behavior.pos.reward = [pos];
            behavior.speed.reward = LinearVelocity(pos,0);
            behavior.speed.reward(:,2) = smoothdata(behavior.speed.reward(:,2),'gaussian',1,'SamplePoints',behavior.speed.reward(:,1));
            behavior.quiet.reward = QuietPeriods( behavior.speed.reward , minimal_speed , minimal_speed_time);
            clear pos camaraR posx posy
            load('laps1.mat','posx','posy');
            [camaraA,~] = find((camara(:,1)-aversiveTS_run(1)/1000)>0,1,'first'); %TimeStamp of the begining of aversive
            pos = [camara(camaraA : camaraA+length(posx)-1),posx,posy];
            %interpolation of dropped frames
            ejeX = pos(~isnan(pos(:,2)),1); dX = pos(~isnan(pos(:,2)),2); dX_int = interp1(ejeX , dX , pos(:,1));
            ejeY = pos(~isnan(pos(:,3)),1); dY = pos(~isnan(pos(:,3)),3); dY_int = interp1(ejeY , dY , pos(:,1));
            pos(:,2) =dX_int; pos(:,3) =dY_int; %saving corrected pos
            behavior.pos.aversive = [pos];
            behavior.speed.aversive = LinearVelocity(pos,0);
            behavior.speed.aversive(:,2) = smoothdata(behavior.speed.aversive(:,2),'gaussian',1,'SamplePoints',behavior.speed.aversive(:,1));
            behavior.quiet.aversive = QuietPeriods(behavior.speed.aversive , minimal_speed , minimal_speed_time);
            clear pos camaraA posx posy
        end
        
        % Generation of no-movements periods
        % Reward
        start = behavior.speed.reward(1,1);   stop = behavior.speed.reward(end,1);
        movement.reward = InvertIntervals(behavior.quiet.reward , start , stop); %keep only those higher than criteria
%         movement.reward(movement.reward(:,2) - movement.reward(:,1) <1,:)=[]; %eliminate 1sec segments
        clear tmp start stop
        % Aversive
        start = behavior.speed.aversive(1,1);   stop = behavior.speed.aversive(end,1);
        movement.aversive = InvertIntervals(behavior.quiet.aversive , start , stop);%keep only those higher than criteria
%         movement.aversive(movement.aversive(:,2) - movement.aversive(:,1) <1,:)=[];
        clear tmp start stop

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
                [SpksTrains.all.aversive , Bins.aversive , Cluster.all.aversive] = spike_train_construction([spks_dHPC;spks_vHPC], clusters.all, cellulartype, binSize(i), limits, events, false,false);
                [Th , pat , eig] = assembly_patternsJFM([SpksTrains.all.aversive'],opts);
%                 save([cd,'\dorsalventral_assemblies_aversiveVF.mat'],'Th' , 'pat' , 'eig' , 'criteria_fr' , 'criteria_n')
%             end
            
            Thresholded.aversive.all = Th;
            patterns.all.aversive = pat;
            clear cond Th pat
            
            % Detection of members
            if not(isempty(Thresholded.aversive.all))
                if numberD>0
                    cond1 =  sum(Thresholded.aversive.all(1:size(clusters.dHPC,1),:),1)>0; %checking of dHPC SU
                    cond2 =  sum(Thresholded.aversive.all(size(clusters.dHPC,1)+1:end,:),1)>0; %checking of vHPC SU
                    cond.dHPC.aversive = and(cond1 , not(cond2));
                    cond.vHPC.aversive = and(cond2 , not(cond1));
                    cond.both.aversive = and(cond1 , cond2); clear cond1 cond2
                else
                    cond1 =  logical(zeros(1,size(Thresholded.aversive.all,2))); %checking of dHPC SU
                    cond2 =  logical(ones(1,size(Thresholded.aversive.all,2))); %checking of vHPC SU
                    cond.dHPC.aversive = and(cond1 , not(cond2));
                    cond.vHPC.aversive = and(cond2 , not(cond1));
                    cond.both.aversive = and(cond1 , cond2); clear cond1 cond2
                end
            else
                cond1 =  false; %checking of dHPC SU
                cond2 =  logical(0); %checking of vHPC SU
                cond.dHPC.aversive = and(cond1 , not(cond2));
                cond.vHPC.aversive = and(cond2 , not(cond1));
                cond.both.aversive = and(cond1 , cond2); clear cond1 cond2
            end
            num_assembliesA = [num_assembliesA , sum(cond.both.aversive) sum(cond.dHPC.aversive) sum(cond.vHPC.aversive)];
            
            % --- Reward ---
            disp('Loading Reward template')
%             if binSize(i) == 0.025
%                 load('dorsalventral_assemblies_reward.mat')
%             else
                disp('Detection of assemblies using Rewarded template')
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
                [SpksTrains.all.reward , Bins.reward , Cluster.all.reward] = spike_train_construction([spks_dHPC;spks_vHPC], clusters.all, cellulartype, binSize(i), limits, events, false,false);
                [Th , pat , eig] = assembly_patternsJFM([SpksTrains.all.reward'],opts);
%                 save([cd,'\dorsalventral_assemblies_rewardVF.mat'],'Th' , 'pat' , 'eig' , 'criteria_fr' , 'criteria_n')
%             end
            
            Thresholded.reward.all = Th;
            patterns.all.reward = pat;
            clear Th pat
            
            % Detection of members using
            if not(isempty(Thresholded.reward.all))
                if numberD>0
                    cond1 =  sum(Thresholded.reward.all(1:size(clusters.dHPC,1),:),1)>0; %checking of dHPC SU
                    cond2 =  sum(Thresholded.reward.all(size(clusters.dHPC,1)+1:end,:),1)>0; %checking of vHPC SU
                    cond.dHPC.reward = and(cond1 , not(cond2));
                    cond.vHPC.reward = and(cond2 , not(cond1));
                    cond.both.reward = and(cond1 , cond2); clear cond1 cond2
                else
                    cond1 =  logical(zeros(1,size(Thresholded.reward.all,2))); %checking of dHPC SU
                    cond2 =  logical(ones(1,size(Thresholded.reward.all,2))); %checking of vHPC SU
                    cond.dHPC.reward = and(cond1 , not(cond2));
                    cond.vHPC.reward = and(cond2 , not(cond1));
                    cond.both.reward = and(cond1 , cond2); clear cond1 cond2
                end
            else
                cond1 =  logical(0); %checking of dHPC SU
                cond2 =  logical(0); %checking of vHPC SU
                cond.dHPC.reward = and(cond1 , not(cond2));
                cond.vHPC.reward = and(cond2 , not(cond1));
                cond.both.reward = and(cond1 , cond2); clear cond1 cond2
            end
            num_assembliesR = [num_assembliesR , sum(cond.both.reward) sum(cond.dHPC.reward) sum(cond.vHPC.reward)];
            
            
            %% Similarity Index Calculation
            
            
            if and(not(isempty(patterns.all.aversive)) , not(isempty(patterns.all.reward)))
                [r.AR , p.AR] = SimilarityIndex(patterns.all.aversive , patterns.all.reward);
                A = sum(p.AR,1)>=1;
                R = sum(p.AR,2)>=1; R = R';
                
                AR.dHPC = cond.dHPC.aversive .* cond.dHPC.reward';
                AR.dHPC = and(AR.dHPC , p.AR);
                
                AR.vHPC = cond.vHPC.aversive .* cond.vHPC.reward';
                AR.vHPC = and(AR.vHPC , p.AR);
                
                AR.both = cond.both.aversive .* cond.both.reward';
                AR.both = and(AR.both , p.AR);
            else
                AR.dHPC = 0;
                AR.vHPC = 0;
                AR.both = 0;
            end
            
            percentagesA = [percentagesA  , sum(sum(AR.both)) , sum(sum(AR.dHPC)) , sum(sum(AR.vHPC))];
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
subplot(231),bar([1 2 3 4] , [sum(Number_of_assemblies.aversive(:,1)) sum(Number_of_assemblies.aversive(:,4)) sum(Number_of_assemblies.aversive(:,7)) sum(Number_of_assemblies.aversive(:,10))]),ylim([0 100])
subplot(232),bar([1 2 3 4] , [sum(Number_of_assemblies.aversive(:,2)) sum(Number_of_assemblies.aversive(:,5)) sum(Number_of_assemblies.aversive(:,8)) sum(Number_of_assemblies.aversive(:,11))]),ylim([0 150])
subplot(233),bar([1 2 3 4] , [sum(Number_of_assemblies.aversive(:,3)) sum(Number_of_assemblies.aversive(:,6)) sum(Number_of_assemblies.aversive(:,9)) sum(Number_of_assemblies.aversive(:,12))]),ylim([0 45])

subplot(234),bar([1 2 3 4] , [sum(Number_of_assemblies.reward(:,1)) sum(Number_of_assemblies.reward(:,4)) sum(Number_of_assemblies.reward(:,7)) sum(Number_of_assemblies.reward(:,10))]),ylim([0 100])
subplot(235),bar([1 2 3 4] , [sum(Number_of_assemblies.reward(:,2)) sum(Number_of_assemblies.reward(:,5)) sum(Number_of_assemblies.reward(:,8)) sum(Number_of_assemblies.reward(:,11))]),ylim([0 150])
subplot(236),bar([1 2 3 4] , [sum(Number_of_assemblies.reward(:,3)) sum(Number_of_assemblies.reward(:,6)) sum(Number_of_assemblies.reward(:,9)) sum(Number_of_assemblies.reward(:,12))]),ylim([0 45])

figure
subplot(131),bar([1 2 3 4] , [sum(percentages(:,1)) sum(percentages(:,4)) sum(percentages(:,7)) sum(percentages(:,10))]),ylim([0 20])
subplot(132),bar([1 2 3 4] , [sum(percentages(:,2)) sum(percentages(:,5)) sum(percentages(:,8)) sum(percentages(:,11))]),ylim([0 60])
subplot(133),bar([1 2 3 4] , [sum(percentages(:,3)) sum(percentages(:,6)) sum(percentages(:,9)) sum(percentages(:,12))]),ylim([0 25])


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