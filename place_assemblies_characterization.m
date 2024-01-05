clear
clc
close all

%% Parameters
path = {'E:\Rat103\usable';'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat132\recordings\in_pyr';'E:\Rat165\in_pyr\'};%List of folders from the path

%Sleep
time_criteria = 1; % minimal time to include a NREM epoch (in min)

% for SU
criteria_fr = 0; %criteria to include or not a SU into the analysis
criteria_n = [3 3]; % minimal number of neurons from each structure [vHPC dHPC]
criteria_type = 0; %criteria for celltype (0:pyr, 1:int, 2:all)
binSize = 0.025; %for qssemblie detection qnd qxctivity strength
n_SU_V = 0;
n_SU_D = 0;

win = 300; % time window for bin construction

% Behavior
minimal_speed = 7; % minimal speed to detect quite periods
minimal_speed_time = 2; % minimal time to detect quite periods

th = 5; % threshold for peak detection
s = true; % true if I wanna save maps at each folder

Number_of_assemblies.aversive = [];
Number_of_assemblies.reward = [];

map.bothA.aversive = []; map.bothA.reward = [];
map.bothR.aversive = []; map.bothR.reward = [];
map.dHPCA.aversive = []; map.dHPCA.reward = [];
map.dHPCR.aversive = []; map.dHPCR.reward = [];
map.vHPCA.aversive = []; map.vHPCA.reward = [];
map.vHPCR.aversive = []; map.vHPCR.reward = [];


Within.bothA.aversive = []; Within.bothA.reward = [];
Within.bothR.aversive = []; Within.bothR.reward = [];
Within.dHPCA.aversive = []; Within.dHPCA.reward = [];
Within.dHPCR.aversive = []; Within.dHPCR.reward = [];
Within.vHPCA.aversive = []; Within.vHPCA.reward = [];
Within.vHPCR.aversive = []; Within.vHPCR.reward = [];

Between.bothA = []; Between.bothR = [];
Between.dHPCA = []; Between.dHPCR = [];
Between.vHPCA = []; Between.vHPCR = [];

% Sacar el filtro que puse del FR en el counts de neuronas
%% Main loop, to iterate across sessions
for tt = 1:length(path)
    %List of folders from the path
    files = dir(path{tt});
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
    clear files dirFlags
    num_assembliesR = [];
    num_assembliesA = [];
    for t = 1 : length(subFolders)-2
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
        movement.reward(movement.reward(:,2) - movement.reward(:,1) <1,:)=[]; %eliminate 1sec segments
        clear tmp start stop
        % Aversive
        start = behavior.speed.aversive(1,1);   stop = behavior.speed.aversive(end,1);
        movement.aversive = InvertIntervals(behavior.quiet.aversive , start , stop);%keep only those higher than criteria
        movement.aversive(movement.aversive(:,2) - movement.aversive(:,1) <1,:)=[];
        clear tmp start stop
        
        %% load sleep states
        disp('Uploading sleep scoring')
        x = dir([cd,'\*-states.mat']);    states = load([cd,'\',x.name]);    states = states.states;
        REM.all = ToIntervals(states==5);    NREM.all = ToIntervals(states==3);    WAKE.all = ToIntervals(states==1);
        %         REM.all(REM.all(:,2)-REM.all(:,1)<60,:) = [];
        clear x states
        
        NREM.all(NREM.all(:,2)-NREM.all(:,1)<60*time_criteria,:)=[];
        NREM.baseline = Restrict(NREM.all,baselineTS./1000);
        NREM.aversive = Restrict(NREM.all,aversiveTS./1000);
        NREM.reward = Restrict(NREM.all,rewardTS./1000);
        
        REM.baseline = Restrict(REM.all,baselineTS./1000);
        REM.aversive = Restrict(REM.all,aversiveTS./1000);
        REM.reward = Restrict(REM.all,rewardTS./1000);
        
        WAKE.baseline = Restrict(WAKE.all,baselineTS./1000);
        WAKE.aversive = Restrict(WAKE.all,aversiveTS./1000);
        WAKE.reward = Restrict(WAKE.all,rewardTS./1000);

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
        K = [K , Cell_type_classification(:,6:7)];
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
            cellulartype = [K(:,1) , K(:,3)];
        elseif criteria_type == 1 % int
            cellulartype = [K(:,1) , K(:,4)];
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
        
        criteria_n = [3 3];
        
        %% Assemblies detection
        if and(numberD >= criteria_n(1),numberV >= criteria_n(2))
            disp('Lets go for the assemblies')
            % --- Options for assemblies detection ---
            opts.Patterns.method = 'ICA';
            opts.threshold.method= 'MarcenkoPastur';
            opts.Patterns.number_of_iterations= 500;
            opts.threshold.permutations_percentile = 0.9;
            opts.threshold.number_of_permutations = 500;
            opts.Patterns.number_of_iterations = 500;
            opts.Members.method = 'Sqrt';
            
            % --- Aversive ---
            disp('Loading Aversive template')
            if isfile('dorsalventral_assemblies_aversive3.mat')
                load('dorsalventral_assemblies_aversive3.mat')
                
                %             if not(exist('Th','var'))
                %                 disp('Detection of assemblies using Aversive template')
                %                 limits = aversiveTS_run./1000;
                %                 events = [];
                %                 events = movement.aversive;
                %                 [SpksTrains.all.aversive , Bins.aversive , Cluster.all.aversive] = spike_train_construction([spks_dHPC;spks_vHPC], clusters.all, cellulartype, binSize, limits, events, false,true);
                %                 [Th , pat] = assembly_patternsJFM([SpksTrains.all.aversive'],opts);
                %                 save([cd,'\dorsalventral_assemblies_aversive3.mat'],'Th' , 'pat' , 'criteria_fr' , 'criteria_n')
                %             end
                
                Thresholded.aversive.all = Th;
                patterns.all.aversive = pat;
                clear cond Th pat
                
                % Detection of members
                cond1 =  sum(Thresholded.aversive.all(1:size(clusters.dHPC,1),:),1)>0; %checking of dHPC SU
                cond2 =  sum(Thresholded.aversive.all(size(clusters.dHPC,1)+1:end,:),1)>0; %checking of vHPC SU
                cond.dHPC.aversive = and(cond1 , not(cond2));
                cond.vHPC.aversive = and(cond2 , not(cond1));
                cond.both.aversive = and(cond1 , cond2); clear cond1 cond2
                
                num_assembliesA = [num_assembliesA ; sum(cond.both.aversive) sum(cond.dHPC.aversive) sum(cond.vHPC.aversive)];
            end
            
            % --- Reward ---
            disp('Loading Reward template')
            if isfile('dorsalventral_assemblies_reward3.mat')
                load('dorsalventral_assemblies_reward3.mat')
                
                %             if not(exist('Th','var'))
                %                 disp('Detection of assemblies using Rewarded template')
                %                 limits = rewardTS_run./1000;
                %                 events = [];
                %                 events = movement.reward;
                %                 [SpksTrains.all.reward , Bins.reward , Cluster.all.reward] = spike_train_construction([spks_dHPC;spks_vHPC], clusters.all, cellulartype, binSize, limits, events, false,true);
                %                 [Th , pat] = assembly_patternsJFM([SpksTrains.all.reward'],opts);
                %                 save([cd,'\dorsalventral_assemblies_reward3.mat'],'Th' , 'pat' , 'criteria_fr' , 'criteria_n')
                %             end
                
                Thresholded.reward.all = Th;
                patterns.all.reward = pat;
                clear Th pat
                
                % Detection of members using
                cond1 =  sum(Thresholded.reward.all(1:size(clusters.dHPC,1),:),1)>0; %checking of dHPC SU
                cond2 =  sum(Thresholded.reward.all(size(clusters.dHPC,1)+1:end,:),1)>0; %checking of vHPC SU
                cond.dHPC.reward = and(cond1 , not(cond2));
                cond.vHPC.reward = and(cond2 , not(cond1));
                cond.both.reward = and(cond1 , cond2); clear cond1 cond2
                
                num_assembliesR = [num_assembliesR ; sum(cond.both.reward) sum(cond.dHPC.reward) sum(cond.vHPC.reward)];
            end
            
            [r.AR , p.AR] = SimilarityIndex(patterns.all.aversive , patterns.all.reward);
            A = sum(p.AR,1)>=1;
            R = sum(p.AR,2)>=1; R = R';
            
            AR.dHPC = cond.dHPC.aversive .* cond.dHPC.reward';
            AR.dHPC = and(AR.dHPC , p.AR);
            
            AR.vHPC = cond.vHPC.aversive .* cond.vHPC.reward';
            AR.vHPC = and(AR.vHPC , p.AR);
            
            AR.both = cond.both.aversive .* cond.both.reward';
            AR.both = and(AR.both , p.AR);
            
            
            clear A R r p
            
            % to save the clusters I used for further analysis
            %             save([cd,'\clusters_included_in_assemblies.mat'],'clusters')
            
            
            %% SpikeTrains construction
            limits = [0 segments.Var1(end)/1000];
            events = [];
            [Spikes , bins , Clusters] = spike_train_construction([spks_dHPC;spks_vHPC], clusters.all, cellulartype, binSize, limits, events, true, true);
            clear limits events
            
            %% Assemblies activation in the entier recording
            % Aversive
            if sum(cond.both.aversive)>=1
                pos = [behavior.pos.aversive(:,1:2) ; behavior.pos.reward(:,1:2)];
                [x xx] = sort(pos(:,1));
                pos = pos(xx,:);
                events = cell(2,1);
                events{1} = movement.aversive;
                events{2}  = movement.reward;
                
                [Maps pc between within] = FiringMap_Assemblies(patterns.all.aversive , cond.both.aversive , [bins' , Spikes] , th , pos , events , 60 , true , true);
                clear pos x xx
                
                if s
                    save([cd,'\Joint_Aversive_Assemblies_Maps.mat'],'Maps' , 'pc' , 'between' , 'within')
                end
                
                pc = or(pc.cond1 , pc.cond2); %logical to select maps to save
                
                map.bothA.aversive = [map.bothA.aversive ; Maps.cond1(pc,:)];
                map.bothA.reward = [map.bothA.reward ; Maps.cond2(pc,:)];
                
                Within.bothA.aversive = [Within.bothA.aversive ; within.cond1(pc,:)];
                Within.bothA.reward = [Within.bothA.reward ; within.cond2(pc,:)];
                
                Between.bothA = [Between.bothA ; between(pc,:)];
                
                clear within between pc Maps pos x xx events
            end
            
            if sum(cond.dHPC.aversive)>=1
                pos = [behavior.pos.aversive(:,1:2) ; behavior.pos.reward(:,1:2)];
                [x xx] = sort(pos(:,1));
                pos = pos(xx,:);
                events = cell(2,1);
                events{1} = movement.aversive;
                events{2}  = movement.reward;
                
                [Maps pc between within] = FiringMap_Assemblies(patterns.all.aversive , cond.dHPC.aversive , [bins' , Spikes] , th , pos , events , 60 , true , true);
                clear pos x xx
                
                if s
                    save([cd,'\dHPC_Aversive_Assemblies_Maps.mat'],'Maps' , 'pc' , 'between' , 'within')
                end
                
                pc = or(pc.cond1 , pc.cond2); %logical to select maps to save
                
                map.dHPCA.aversive = [map.dHPCA.aversive ; Maps.cond1(pc,:)];
                map.dHPCA.reward = [map.dHPCA.reward ; Maps.cond2(pc,:)];
                
                Within.dHPCA.aversive = [Within.dHPCA.aversive ; within.cond1(pc,:)];
                Within.dHPCA.reward = [Within.dHPCA.reward ; within.cond2(pc,:)];
                
                Between.dHPCA = [Between.dHPCA ; between(pc,:)];
                
                clear within between pc Maps pos x xx events
            end
            
            if sum(cond.vHPC.aversive)>=1
                pos = [behavior.pos.aversive(:,1:2) ; behavior.pos.reward(:,1:2)];
                [x xx] = sort(pos(:,1));
                pos = pos(xx,:);
                events = cell(2,1);
                events{1} = movement.aversive;
                events{2}  = movement.reward;
                
                [Maps pc between within] = FiringMap_Assemblies(patterns.all.aversive , cond.vHPC.aversive , [bins' , Spikes] , th , pos , events , 60 , true , true);
                clear pos x xx
                
                if s
                    save([cd,'\vHPC_Aversive_Assemblies_Maps.mat'],'Maps' , 'pc' , 'between' , 'within')
                end
                
                pc = or(pc.cond1 , pc.cond2); %logical to select maps to save
                
                map.vHPCA.aversive = [map.vHPCA.aversive ; Maps.cond1(pc,:)];
                map.vHPCA.reward = [map.vHPCA.reward ; Maps.cond2(pc,:)];
                
                Within.vHPCA.aversive = [Within.vHPCA.aversive ; within.cond1(pc,:)];
                Within.vHPCA.reward = [Within.vHPCA.reward ; within.cond2(pc,:)];
                
                Between.vHPCA = [Between.vHPCA ; between(pc,:)];
                
                clear within between pc Maps pos x xx events
            end
            
            % Reward
            if sum(cond.both.reward)>=1
                pos = [behavior.pos.aversive(:,1:2) ; behavior.pos.reward(:,1:2)];
                [x xx] = sort(pos(:,1));
                pos = pos(xx,:);
                events = cell(2,1);
                events{1} = movement.reward;
                events{2}  = movement.aversive;
                
                [Maps pc between within] = FiringMap_Assemblies(patterns.all.reward , cond.both.reward , [bins' , Spikes] , th , pos , events , 60 , true , true);
                clear pos x xx
                
                if s
                    save([cd,'\Joint_Reward_Assemblies_Maps.mat'],'Maps' , 'pc' , 'between' , 'within')
                end
                
                pc = or(pc.cond1 , pc.cond2); %logical to select maps to save
                
                map.bothR.reward = [map.bothR.reward ; Maps.cond1(pc,:)];
                map.bothR.aversive = [map.bothR.aversive ; Maps.cond2(pc,:)];
                
                Within.bothR.reward = [Within.bothR.reward ; within.cond1(pc,:)];
                Within.bothR.aversive = [Within.bothR.aversive ; within.cond2(pc,:)];
                
                Between.bothR = [Between.bothR ; between(pc,:)];
                
                clear within between pc Maps pos x xx events
            end
            
            if sum(cond.dHPC.reward)>=1
                pos = [behavior.pos.aversive(:,1:2) ; behavior.pos.reward(:,1:2)];
                [x xx] = sort(pos(:,1));
                pos = pos(xx,:);
                events = cell(2,1);
                events{1} = movement.reward;
                events{2}  = movement.aversive;
                
                [Maps pc between within] = FiringMap_Assemblies(patterns.all.reward , cond.dHPC.reward , [bins' , Spikes] , th , pos , events , 60 , true , true);
                clear pos x xx
                
                if s
                    save([cd,'\dHPC_Reward_Assemblies_Maps.mat'],'Maps' , 'pc' , 'between' , 'within')
                end
                
                pc = or(pc.cond1 , pc.cond2); %logical to select maps to save
                
                map.dHPCR.reward = [map.dHPCR.reward ; Maps.cond1(pc,:)];
                map.dHPCR.aversive = [map.dHPCR.aversive ; Maps.cond2(pc,:)];
                
                Within.dHPCR.reward = [Within.dHPCR.reward ; within.cond1(pc,:)];
                Within.dHPCR.aversive = [Within.dHPCR.aversive ; within.cond2(pc,:)];
                
                Between.dHPCR = [Between.dHPCR ; between(pc,:)];
                
                clear within between pc Maps pos x xx events
            end
            
            if sum(cond.vHPC.reward)>=1
                pos = [behavior.pos.aversive(:,1:2) ; behavior.pos.reward(:,1:2)];
                [x xx] = sort(pos(:,1));
                pos = pos(xx,:);
                events = cell(2,1);
                events{1} = movement.reward;
                events{2}  = movement.aversive;
                
                [Maps pc between within] = FiringMap_Assemblies(patterns.all.reward , cond.vHPC.reward , [bins' , Spikes] , th , pos , events , 60 , true , true);
                clear pos x xx
                
                if s
                    save([cd,'\vHPC_Reward_Assemblies_Maps.mat'],'Maps' , 'pc' , 'between' , 'within')
                end
                
                pc = or(pc.cond1 , pc.cond2); %logical to select maps to save
                
                map.vHPCR.reward = [map.vHPCR.reward ; Maps.cond1(pc,:)];
                map.vHPCR.aversive = [map.vHPCR.aversive ; Maps.cond2(pc,:)];
                
                Within.vHPCR.reward = [Within.vHPCR.reward ; within.cond1(pc,:)];
                Within.vHPCR.aversive = [Within.vHPCR.aversive ; within.cond2(pc,:)];
                
                Between.vHPCR = [Between.vHPCR ; between(pc,:)];
                
                clear within between pc Maps pos x xx events
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
    end
    
    Number_of_assemblies.aversive = [Number_of_assemblies.aversive ; sum(num_assembliesA)];
    Number_of_assemblies.reward = [Number_of_assemblies.reward ; sum(num_assembliesR)];
    clear num_assembliesA num_assembliesR
    
end

%% Plotting section
% The following section is used to plot the outputs of the script
% Please, be sure that you undersand the structure of the matrix you are
% intending to plot.

%% Both
figure
x = map.bothA.aversive - min(map.bothA.aversive,[],2);
x = x./max(x,[],2);

y = map.bothA.reward - min(map.bothA.reward,[],2);
y = y./max(y,[],2);

[c i] = max(x,[],2);
[c i] = sort(i);

subplot(121),imagesc(x(i,:))
subplot(122),imagesc(y(i,:))


figure
x = map.bothR.aversive - min(map.bothR.aversive,[],2);
x = x./max(x,[],2);

y = map.bothR.reward - min(map.bothR.reward,[],2);
y = y./max(y,[],2);

[c i] = max(y,[],2);
[c i] = sort(i);

subplot(121),imagesc(x(i,:))
subplot(122),imagesc(y(i,:))

%% dHPC
figure
x = map.dHPCA.aversive - min(map.dHPCA.aversive,[],2);
x = x./max(x,[],2);

y = map.dHPCA.reward - min(map.dHPCA.reward,[],2);
y = y./max(y,[],2);

[c i] = max(x,[],2);
[c i] = sort(i);

subplot(121),imagesc(x(i,:))
subplot(122),imagesc(y(i,:))


figure
x = map.dHPCR.aversive - min(map.dHPCR.aversive,[],2);
x = x./max(x,[],2);

y = map.dHPCR.reward - min(map.dHPCR.reward,[],2);
y = y./max(y,[],2);

[c i] = max(y,[],2);
[c i] = sort(i);

subplot(121),imagesc(x(i,:))
subplot(122),imagesc(y(i,:))

%% vHPC
figure
x = map.vHPCA.aversive - min(map.vHPCA.aversive,[],2);
x = x./max(x,[],2);

y = map.vHPCA.reward - min(map.vHPCA.reward,[],2);
y = y./max(y,[],2);

[c i] = max(x,[],2);
[c i] = sort(i);

subplot(121),imagesc(x(i,:))
subplot(122),imagesc(y(i,:))


figure
x = map.vHPCR.aversive - min(map.vHPCR.aversive,[],2);
x = x./max(x,[],2);

y = map.vHPCR.reward - min(map.vHPCR.reward,[],2);
y = y./max(y,[],2);

[c i] = max(y,[],2);
[c i] = sort(i);

subplot(121),imagesc(x(i,:))
subplot(122),imagesc(y(i,:))


%% Plot parameters
% Spatial Correlation
figure
subplot(3,2,1)
y = [Between.bothA(:,1) ; Within.bothA.reward(:,1) ; Within.bothA.aversive(:,1)];
x = [ones(length(Between.bothA(:,1)),1) ; ones(length(Within.bothA.reward(:,1)),1)*2 ; ones(length(Within.bothA.aversive(:,1)),1)*3];
scatter(x,y,'filled'),hold on,xlim([0 4]),ylim([-0.55 1.05])
scatter([1 2 3],[nanmedian(Between.bothA(:,1)) , nanmedian( Within.bothA.reward(:,1)), nanmedian( Within.bothA.aversive(:,1))],'filled')
% boxplot(y,x),ylim([-0.55 1.05])

[p,tbl,stats] = kruskalwallis(y,x);
c = multcompare(stats)


subplot(3,2,2)
y = [Between.bothR(:,1) ; Within.bothR.reward(:,1) ; Within.bothR.aversive(:,1)];
x = [ones(length(Between.bothR(:,1)),1) ; ones(length(Within.bothR.reward(:,1)),1)*2 ; ones(length(Within.bothR.aversive(:,1)),1)*3];
scatter(x,y,'filled'),hold on,xlim([0 4]),ylim([-0.55 1.05])
scatter([1 2 3],[nanmedian(Between.bothR(:,1)) , nanmedian( Within.bothR.reward(:,1)), nanmedian( Within.bothR.aversive(:,1))],'filled')
% boxplot(y,x),ylim([-0.55 1.05])

[p,tbl,stats] = kruskalwallis(y,x);
c = multcompare(stats)


subplot(3,2,3)
y = [Between.dHPCA(:,1) ; Within.dHPCA.reward(:,1) ; Within.dHPCA.aversive(:,1)];
x = [ones(length(Between.dHPCA(:,1)),1) ; ones(length(Within.dHPCA.reward(:,1)),1)*2 ; ones(length(Within.dHPCA.aversive(:,1)),1)*3];
scatter(x,y,'filled'),hold on,xlim([0 4]),ylim([-0.55 1.05])
scatter([1 2 3],[nanmedian(Between.dHPCA(:,1)) , nanmedian( Within.dHPCA.reward(:,1)), nanmedian( Within.dHPCA.aversive(:,1))],'filled')
% boxplot(y,x),ylim([-0.55 1.05])

[p,tbl,stats] = kruskalwallis(y,x);
c = multcompare(stats)

subplot(3,2,4)
y = [Between.dHPCR(:,1) ; Within.dHPCR.reward(:,1) ; Within.dHPCR.aversive(:,1)];
x = [ones(length(Between.dHPCR(:,1)),1) ; ones(length(Within.dHPCR.reward(:,1)),1)*2 ; ones(length(Within.dHPCR.aversive(:,1)),1)*3];
scatter(x,y,'filled'),hold on,xlim([0 4]),ylim([-0.55 1.05])
scatter([1 2 3],[nanmedian(Between.dHPCR(:,1)) , nanmedian( Within.dHPCR.reward(:,1)), nanmedian( Within.dHPCR.aversive(:,1))],'filled')
% boxplot(y,x),ylim([-0.55 1.05])

[p,tbl,stats] = kruskalwallis(y,x);
c = multcompare(stats)

subplot(3,2,5)
y = [Between.vHPCA(:,1) ; Within.vHPCA.reward(:,1) ; Within.vHPCA.aversive(:,1)];
x = [ones(length(Between.vHPCA(:,1)),1) ; ones(length(Within.vHPCA.reward(:,1)),1)*2 ; ones(length(Within.vHPCA.aversive(:,1)),1)*3];
scatter(x,y,'filled'),hold on,xlim([0 4]),ylim([-0.55 1.05])
scatter([1 2 3],[nanmedian(Between.vHPCA(:,1)) , nanmedian( Within.vHPCA.reward(:,1)), nanmedian( Within.vHPCA.aversive(:,1))],'filled')
% boxplot(y,x),ylim([-0.55 1.05])

[p,tbl,stats] = kruskalwallis(y,x);
c = multcompare(stats)

subplot(3,2,6)
y = [Between.vHPCR(:,1) ; Within.vHPCR.reward(:,1) ; Within.vHPCR.aversive(:,1)];
x = [ones(length(Between.vHPCR(:,1)),1) ; ones(length(Within.vHPCR.reward(:,1)),1)*2 ; ones(length(Within.vHPCR.aversive(:,1)),1)*3];
scatter(x,y,'filled'),hold on,xlim([0 4]),ylim([-0.55 1.05])
scatter([1 2 3],[nanmedian(Between.vHPCR(:,1)) , nanmedian( Within.vHPCR.reward(:,1)), nanmedian( Within.vHPCR.aversive(:,1))],'filled')
% boxplot(y,x),ylim([-0.55 1.05])

[p,tbl,stats] = kruskalwallis(y,x);
c = multcompare(stats)


%% Comparing Between groups across type of assemblies
figure
subplot(121)
y = [Between.bothA(:,1) ; Between.dHPCA(:,1) ;  Between.vHPCA(:,1)];
x = [ones(length(Between.bothA(:,1)),1) ; ones(length(Between.dHPCA(:,1)),1)*2 ; ones(length(Between.vHPCA(:,1)),1)*3];
scatter(x,y,'filled'),hold on,xlim([0 4]),ylim([-0.55 1.05])
scatter([1 2 3],[nanmedian(Between.bothA(:,1)) , nanmedian(Between.dHPCA(:,1)), nanmedian(Between.vHPCA(:,1))],'filled')

[p,tbl,stats] = kruskalwallis(y,x);
c = multcompare(stats)


subplot(122)
y = [Between.bothR(:,1) ; Between.dHPCR(:,1) ;  Between.vHPCR(:,1)];
x = [ones(length(Between.bothR(:,1)),1) ; ones(length(Between.dHPCR(:,1)),1)*2 ; ones(length(Between.vHPCR(:,1)),1)*3];
scatter(x,y,'filled'),hold on,xlim([0 4]),ylim([-0.55 1.05])
scatter([1 2 3],[nanmedian(Between.bothR(:,1)) , nanmedian(Between.dHPCR(:,1)), nanmedian(Between.vHPCR(:,1))],'filled')

[p,tbl,stats] = kruskalwallis(y,x);
c = multcompare(stats)