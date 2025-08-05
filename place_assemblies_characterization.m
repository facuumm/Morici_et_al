clear
clc
close all

%% Parameters
path = {'E:\Rat103\usable';'E:\Rat126\Ephys\in_Pyr';'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat132\recordings\in_pyr';'E:\Rat165\in_pyr\'};%List of folders from the path

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

Specificity.bothA = []; Specificity.bothR = [];
Specificity.dHPCA = []; Specificity.dHPCR = [];
Specificity.vHPCA = []; Specificity.vHPCR = [];

Iterator.bothA = []; Iterator.bothR = [];
Iterator.dHPCA = []; Iterator.dHPCR = [];
Iterator.vHPCA = []; Iterator.vHPCR = [];

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
        
        %% Loading TS of the sessions
        disp('Uploading session time stamps')
        x = dir([cd,'\*.cat.evt']);
        segments = readtable([cd,'\',x.name],'FileType','text');
        clear x
        
        %% TimeStamps of begening and end of the sleep and awake trials
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
        
        %% Load sleep states
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
        
        criteria_n = [3 3];
        if or(numberD > 3 , numberV > 3)
            %% --- Aversive ---
            disp('Lets go for the assemblies')
            if isfile('dorsalventral_assemblies_aversiveVF.mat')
                disp('Loading Aversive template')
                load('dorsalventral_assemblies_aversiveVF.mat')
            end
            
            Thresholded.aversive.all = Th;
            patterns.all.aversive = pat.*Th;
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
                cond1 =  logical(0); %checking of dHPC SU
                cond2 =  logical(0); %checking of vHPC SU
                cond.dHPC.aversive = and(cond1 , not(cond2));
                cond.vHPC.aversive = and(cond2 , not(cond1));
                cond.both.aversive = and(cond1 , cond2); clear cond1 cond2
            end
            num_assembliesA = [num_assembliesA ; sum(cond.both.aversive) sum(cond.dHPC.aversive) sum(cond.vHPC.aversive)];
            
            %% --- Reward ---
            disp('Loading Reward template')
            if isfile('dorsalventral_assemblies_rewardVF.mat')
                load('dorsalventral_assemblies_rewardVF.mat')
            end
            
            Thresholded.reward.all = Th;
            patterns.all.reward = pat.*Th;
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
            num_assembliesR = [num_assembliesR ; sum(cond.both.reward) sum(cond.dHPC.reward) sum(cond.vHPC.reward)];
            
            %% SpikeTrains construction
            limits = [0 segments.Var1(end)/1000];
            events = [];
            [Spikes , bins , Clusters] = spike_train_construction([spks_dHPC;spks_vHPC], clusters.all, cellulartype, binSize, limits, events, false, true);
            clear limits events
            
            %% Assemblies activation in the entier recording
            % Aversive
            if and(numberD >= criteria_n(1),numberV >= criteria_n(2))
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
                    
                    
                    
                    Pc = logical(pc.cond1); %logical to select maps to save
                    
                    Iterator.bothA = [Iterator.bothA ; Pc];
                    
                    Specificity.bothA = [Specificity.bothA ; pc.specificity.cond1];
                    
                    map.bothA.aversive = [map.bothA.aversive ; Maps.cond1];
                    map.bothA.reward = [map.bothA.reward ; Maps.cond2];
                    
                    Within.bothA.aversive = [Within.bothA.aversive ; within.cond1];
                    Within.bothA.reward = [Within.bothA.reward ; within.cond2];
                    
                    Between.bothA = [Between.bothA ; between];
                    
                    clear within between pc Maps pos x xx events Pc
                end
            end
            
            if sum(cond.dHPC.aversive)>=1
                pos = [behavior.pos.aversive(:,1:2) ; behavior.pos.reward(:,1:2)];
                [~, xx] = sort(pos(:,1));
                pos = pos(xx,:);
                events = cell(2,1);
                events{1} = movement.aversive;
                events{2}  = movement.reward;
                
                [Maps pc between within] = FiringMap_Assemblies(patterns.all.aversive , cond.dHPC.aversive , [bins' , Spikes] , th , pos , events , 60 , true , true);
                clear pos x xx
                
                if s
                    save([cd,'\dHPC_Aversive_Assemblies_Maps.mat'],'Maps' , 'pc' , 'between' , 'within')
                end
                
                Pc = logical(pc.cond1); %logical to select maps to save
                
                Iterator.dHPCA = [Iterator.dHPCA ; Pc];
                
                Specificity.dHPCA = [Specificity.dHPCA ; pc.specificity.cond1];
                
                map.dHPCA.aversive = [map.dHPCA.aversive ; Maps.cond1];
                map.dHPCA.reward = [map.dHPCA.reward ; Maps.cond2];
                
                Within.dHPCA.aversive = [Within.dHPCA.aversive ; within.cond1];
                Within.dHPCA.reward = [Within.dHPCA.reward ; within.cond2];
                
                Between.dHPCA = [Between.dHPCA ; between];
                
                clear within between pc Maps pos x xx events Pc
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
                
                Pc = logical(pc.cond1); %logical to select maps to save
                
                Iterator.vHPCA = [Iterator.vHPCA ; Pc];
                
                Specificity.vHPCA = [Specificity.vHPCA ; pc.specificity.cond1];
                
                map.vHPCA.aversive = [map.vHPCA.aversive ; Maps.cond1];
                map.vHPCA.reward = [map.vHPCA.reward ; Maps.cond2];
                
                Within.vHPCA.aversive = [Within.vHPCA.aversive ; within.cond1];
                Within.vHPCA.reward = [Within.vHPCA.reward ; within.cond2];
                
                Between.vHPCA = [Between.vHPCA ; between];
                
                clear within between pc Maps pos x xx events Pc
            end
            
            % Reward
            if and(numberD >= criteria_n(1),numberV >= criteria_n(2))
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
                    
                    Pc = logical(pc.cond2); %logical to select maps to save
                    
                    Iterator.bothR = [Iterator.bothR ; Pc];
                
                    Specificity.bothR = [Specificity.bothR ; pc.specificity.cond2];
                    
                    map.bothR.reward = [map.bothR.reward ; Maps.cond1];
                    map.bothR.aversive = [map.bothR.aversive ; Maps.cond2];
                    
                    Within.bothR.reward = [Within.bothR.reward ; within.cond1];
                    Within.bothR.aversive = [Within.bothR.aversive ; within.cond2];
                    
                    Between.bothR = [Between.bothR ; between];
                    
                    clear within between pc Maps pos x xx events Pc
                end
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
                
                Pc = logical(pc.cond2); %logical to select maps to save
                
                Iterator.dHPCR = [Iterator.dHPCR ; Pc];
                
                Specificity.dHPCR = [Specificity.dHPCR ; pc.specificity.cond2];
                
                map.dHPCR.reward = [map.dHPCR.reward ; Maps.cond1];
                map.dHPCR.aversive = [map.dHPCR.aversive ; Maps.cond2];
                
                Within.dHPCR.reward = [Within.dHPCR.reward ; within.cond1];
                Within.dHPCR.aversive = [Within.dHPCR.aversive ; within.cond2];
                
                Between.dHPCR = [Between.dHPCR ; between];
                
                clear within between pc Maps pos x xx events Pc
            end
            
            if sum(cond.vHPC.reward)>=1
                pos = [behavior.pos.aversive(:,1:2) ; behavior.pos.reward(:,1:2)];
                [x, xx] = sort(pos(:,1));
                pos = pos(xx,:);
                events = cell(2,1);
                events{1} = movement.reward;
                events{2}  = movement.aversive;
                
                [Maps pc between within] = FiringMap_Assemblies(patterns.all.reward , cond.vHPC.reward , [bins' , Spikes] , th , pos , events , 60 , true , true);
                clear pos x xx
                
                if s
                    save([cd,'\vHPC_Reward_Assemblies_Maps.mat'],'Maps' , 'pc' , 'between' , 'within')
                end
                
                Pc = logical(pc.cond2); %logical to select maps to save
                
                Iterator.vHPCR = [Iterator.vHPCR ; Pc];
                
                Specificity.vHPCR = [Specificity.vHPCR ; pc.specificity.cond2];
                
                map.vHPCR.reward = [map.vHPCR.reward ; Maps.cond1];
                map.vHPCR.aversive = [map.vHPCR.aversive ; Maps.cond2];
                
                Within.vHPCR.reward = [Within.vHPCR.reward ; within.cond1];
                Within.vHPCR.aversive = [Within.vHPCR.aversive ; within.cond2];
                
                Between.vHPCR = [Between.vHPCR ; between];
                
                clear within between pc Maps pos x xx events Pc
            end
        end
        disp(' ')
        
        %% Clear Work space
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

% Plot example of assemblies locked and non-locked
figure
% Both
subplot(3,2,1),
x = (map.bothA.aversive);
plot(x(33,:),'k'),hold on
plot(x(5,:),'r')

subplot(3,2,2),
y = map.bothR.reward;
plot(y(23,:),'k'),hold on
plot(y(18,:),'b')

%dHPC
subplot(3,2,3),
x = (map.dHPCA.aversive);
plot(x(130,:),'k'),hold on
plot(x(129,:),'r')

subplot(3,2,4),
x = (map.dHPCR.reward);
plot(x(145,:),'k'),hold on
plot(x(146,:),'b')

%vHPC
subplot(3,2,5),
x = (map.vHPCA.aversive);
plot(x(2,:),'k'),hold on
plot(x(11,:),'r')

subplot(3,2,6),
x = (map.vHPCR.reward);
plot(x(2,:),'k'),hold on
plot(x(25,:),'b')

% for plot Spatial Information
figure
subplot(131)
x = [ones(size(Specificity.bothA)) ; ones(size(Specificity.bothR))*2];
y = [Specificity.bothA ; Specificity.bothR];
scatter(x,y,'filled','jitter',0.1), xlim([0 3]),ylim([0 3]),hold on
scatter([1 2] , [nanmean(Specificity.bothA) nanmean(Specificity.bothR)],'filled')

subplot(132)
x = [ones(size(Specificity.dHPCA)) ; ones(size(Specificity.dHPCR))*2];
y = [Specificity.dHPCA ; Specificity.dHPCR];
scatter(x,y,'filled','jitter',0.1), xlim([0 3]),ylim([0 3]),hold on
scatter([1 2] , [nanmean(Specificity.dHPCA) nanmean(Specificity.dHPCR)],'filled')

subplot(133)
x = [ones(size(Specificity.vHPCA)) ; ones(size(Specificity.vHPCR))*2];
y = [Specificity.vHPCA ; Specificity.vHPCR];
scatter(x,y,'filled','jitter',0.1), xlim([0 3]),ylim([0 3]),hold on
scatter([1 2] , [nanmean(Specificity.vHPCA) nanmean(Specificity.vHPCR)],'filled')


% === Prepare data for two-way ANOVA ===

% Dependent variable: Spatial Specificity
y = [ ...
    Specificity.bothA ; Specificity.bothR ; ...
    Specificity.dHPCA ; Specificity.dHPCR ; ...
    Specificity.vHPCA ; Specificity.vHPCR ];

% Factor 1: Region
region = [ ...
    repmat({'Joint'}, length(Specificity.bothA), 1); ...
    repmat({'Joint'}, length(Specificity.bothR), 1); ...
    repmat({'dHPC'}, length(Specificity.dHPCA), 1); ...
    repmat({'dHPC'}, length(Specificity.dHPCR), 1); ...
    repmat({'vHPC'}, length(Specificity.vHPCA), 1); ...
    repmat({'vHPC'}, length(Specificity.vHPCR), 1) ];

% Factor 2: Condition
condition = [ ...
    repmat({'Aversive'}, length(Specificity.bothA), 1); ...
    repmat({'Reward'},   length(Specificity.bothR), 1); ...
    repmat({'Aversive'}, length(Specificity.dHPCA), 1); ...
    repmat({'Reward'},   length(Specificity.dHPCR), 1); ...
    repmat({'Aversive'}, length(Specificity.vHPCA), 1); ...
    repmat({'Reward'},   length(Specificity.vHPCR), 1) ];

% === Run two-way ANOVA with interaction ===
[p, tbl, stats] = anovan(y, {region, condition}, ...
    'model', 'interaction', ...
    'varnames', {'Region', 'Condition'});

% === Post-hoc comparisons: Bonferroni correction ===
[c, m, h, gnames] = multcompare(stats, ...
    'dimension', [1 2], ...
    'ctype', 'bonferroni');  % correct for multiple comparisons

% === Print only Reward vs Aversive comparisons within each region ===
comparisons = {
    'Joint & Aversive', 'Joint & Reward';
    'dHPC & Aversive',  'dHPC & Reward';
    'vHPC & Aversive',  'vHPC & Reward';
};

fprintf('\nPost-hoc comparisons (Reward vs Aversive within each region):\n')
for i = 1:size(comparisons,1)
    g1 = find(strcmp(gnames, comparisons{i,1}));
    g2 = find(strcmp(gnames, comparisons{i,2}));
    row = find((c(:,1)==g1 & c(:,2)==g2) | (c(:,1)==g2 & c(:,2)==g1));
    if ~isempty(row)
        fprintf('%s vs %s:\tp = %.4f\n', ...
            comparisons{i,1}, comparisons{i,2}, c(row,6));
    end
end







% for plot the percentages
figure,

% Define las variables en una celda para automatizar
fields = {'bothA', 'bothR', 'dHPCA', 'dHPCR', 'vHPCA', 'vHPCR'};

for i = 1:length(fields)
    data = Iterator.(fields{i});
    nMod = sum(data == 1);
    total = length(data); % asume que 0 es no modulado
    nNonMod = total - nMod;

    pMod = (nMod / total) * 100;
    pNonMod = 100 - pMod;

    % Crear etiquetas personalizadas
    labels = {
        sprintf('Modulated - %.0f%% (%d)', pMod, nMod), ...
        sprintf('Non-Modulated - %.0f%% (%d)', pNonMod, nNonMod)
    };

    subplot(3,2,i)
    pie([pMod, pNonMod], labels)
    title(fields{i}, 'Interpreter', 'none')
end