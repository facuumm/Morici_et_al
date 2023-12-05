clear
clc
% close all

%% Parameters
path = {'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat103\usable';'E:\Rat132\recordings\in_pyr'};%List of folders from the path

%Sleep
time_criteria = 1; % minimal time to include a NREM epoch (in min)

% for SU
criteria_fr = 0.0; %criteria to include or not a SU into the analysis
criteria_n = 3; % minimal number of neurons from each structure [vHPC dHPC]
criteria_type = 0; %criteria for celltype (0:pyr, 1:int, 2:all)
binSize = 0.02;

% Behavior
minimal_speed = 7; % minimal speed to detect quite periods
minimal_speed_time = 2; % minimal time to detect quite periods

n_SU_V = 0;
n_SU_D = 0;


EV.aversive.dvHPC = [];   EV.reward.dvHPC = [];
EV.aversive.dvHPCi = [];   EV.reward.dvHPCi = [];
FiringRate.dHPC.baseline = [];
FiringRate.dHPC.reward = [];
FiringRate.dHPC.aversive = [];
FiringRate.vHPC.baseline = [];
FiringRate.vHPC.reward = [];
FiringRate.vHPC.aversive = [];

number.assemblies.aversive = [];
number.assemblies.reward = [];

activity.assemblies.aversive.selective = [];
activity.assemblies.reward.selective = [];
activity.assemblies.aversive.all = [];
activity.assemblies.reward.all = [];
%% Main lo op, to iterate across sessions
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
        clear y
        
        %% load sleep states
        disp('Uploading sleep scoring')
        x = dir([cd,'\*-states.mat']);    states = load([cd,'\',x.name]);    states = states.states;
        REM.all = ToIntervals(states==5);    NREM.all = ToIntervals(states==3);    WAKE.all = ToIntervals(states==1);
        clear x states
        
        NREM.all(NREM.all(:,2)-NREM.all(:,1)<60*time_criteria,:)=[];
        NREM.baseline = Restrict(NREM.all,baselineTS./1000);
        NREM.aversive = Restrict(NREM.all,aversiveTS./1000);
        NREM.reward = Restrict(NREM.all,rewardTS./1000);
        
        REM.baseline = Restrict(REM.all,baselineTS./1000);
        REM.aversive = Restrict(REM.all,aversiveTS./1000);
        REM.reward = Restrict(REM.all,rewardTS./1000);
        
        %% load coordinated ripple bursts
        load('coordinated_ripple_bursts.mat')
        coordinated_ripple_bursts = [coordinated_ripple_bursts(:,1)  coordinated_ripple_bursts(:,3)];
        ripple_bursts.baseline = Restrict(coordinated_ripple_bursts,baselineTS./1000);
        ripple_bursts.reward = Restrict(coordinated_ripple_bursts,rewardTS./1000);
        ripple_bursts.aversive = Restrict(coordinated_ripple_bursts,aversiveTS./1000);
        ripple_bursts.all = coordinated_ripple_bursts;
        clear coordinated_ripple_bursts
        
        
        % Reduce the number of events to the minimal value across conditions
        m = min([length(ripple_bursts.baseline) , length(ripple_bursts.reward) , length(ripple_bursts.aversive)]);
        
        ripple_bursts.baseline = ripple_bursts.baseline(randperm(size(ripple_bursts.baseline,1)),:);
        ripple_bursts.baseline = ripple_bursts.baseline(1:m,:);
        ripple_bursts.baseline = sort(ripple_bursts.baseline);
        
        ripple_bursts.reward = ripple_bursts.reward(randperm(size(ripple_bursts.reward,1)),:);
        ripple_bursts.reward = ripple_bursts.reward(1:m,:);
        ripple_bursts.reward = sort(ripple_bursts.reward);
        
        ripple_bursts.aversive = ripple_bursts.aversive(randperm(size(ripple_bursts.aversive,1)),:);
        ripple_bursts.aversive = ripple_bursts.aversive(1:m,:);
        ripple_bursts.aversive = sort(ripple_bursts.aversive);
        
        
        %% Load ripples
        ripplesD = table2array(readtable('ripplesD_customized2.csv'));
        ripplesV = table2array(readtable('ripplesV_customized2.csv'));
        % coordination
        coordinated = [];
        coordinatedV = [];
        coordinatedV_refined = [];
        cooridnated_event = [];
        for i = 1:length(ripplesD)
            r = ripplesD(i,:);
            tmp = sum(and(ripplesV(:,2)>= r(1,2)-0.1, ripplesV(:,2)<= r(1,2)+0.1));
            if tmp>0
                z = ripplesV(and(ripplesV(:,2)>= r(1,2)-0.1, ripplesV(:,2)<= r(1,2)+0.1),:);
                coordinatedV = [coordinatedV ; z];
                [p,indice] = min(abs(r(2)-z(:,2)));
                coordinatedV_refined = [coordinatedV_refined ; z(indice,:)];
                coordinated = [coordinated ; r];
                
                cooridnated_event = [cooridnated_event ; min([r(1) , z(indice,1)]) , max([r(3) , z(indice,3)])];
                
                clear tmp2 tmp1 p indice z
            end
            clear r
        end
        clear x tmp i
        
        ripples.dHPC.baseline = Restrict(ripplesD , NREM.baseline);
        ripples.dHPC.reward = Restrict(ripplesD , NREM.reward);
        ripples.dHPC.aversive = Restrict(ripplesD , NREM.aversive);
        
        ripples.vHPC.baseline = Restrict(ripplesV , NREM.baseline);
        ripples.vHPC.reward = Restrict(ripplesV , NREM.reward);
        ripples.vHPC.aversive = Restrict(ripplesV , NREM.aversive);
        
        ripples.dHPC.coordinated.baseline = Restrict(coordinated , NREM.baseline);
        ripples.dHPC.coordinated.reward = Restrict(coordinated , NREM.reward);
        ripples.dHPC.coordinated.aversive = Restrict(coordinated , NREM.aversive);
        
        ripples.vHPC.coordinated.baseline = Restrict(coordinatedV_refined , NREM.baseline);
        ripples.vHPC.coordinated.reward = Restrict(coordinatedV_refined , NREM.reward);
        ripples.vHPC.coordinated.aversive = Restrict(coordinatedV_refined , NREM.aversive);
        ripple_event.baseline = Restrict(cooridnated_event,baselineTS./1000);
        ripple_event.reward = Restrict(cooridnated_event,rewardTS./1000);
        ripple_event.aversive = Restrict(cooridnated_event,aversiveTS./1000);
        ripple_event.all = cooridnated_event;
        
        % Reduce the number of events to the minimal value across conditions
        m = min([length(ripple_event.baseline) , length(ripple_event.reward) , length(ripple_event.aversive)]);
        
        ripple_event.baseline = ripple_event.baseline(randperm(size(ripple_event.baseline,1)),:);
        ripple_event.baseline = ripple_event.baseline(1:m,:);
        ripple_event.baseline = sort(ripple_event.baseline);
        
        ripple_event.reward = ripple_event.reward(randperm(size(ripple_event.reward,1)),:);
        ripple_event.reward = ripple_event.reward(1:m,:);
        ripple_event.reward = sort(ripple_event.reward);
        
        ripple_event.aversive = ripple_event.aversive(randperm(size(ripple_event.aversive,1)),:);
        ripple_event.aversive = ripple_event.aversive(1:m,:);
        ripple_event.aversive = sort(ripple_event.aversive);
        
        %% Awake
        disp('Uploading digital imputs')
        % Load digitalin.mat
        load('digitalin.mat')
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
        
        %% Spikes
        %Load Units
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
                n_SU_D = [n_SU_D ; length(group_dHPC)];
                spks_dHPC = spks(ismember(spks(:,1),group_dHPC(:,1)),:); %keep spks from good clusters
            else
                n_SU_V = [n_SU_V ; length(group_vHPC)];
                spks_vHPC = spks(ismember(spks(:,1),group_vHPC(:,1)),:); %keep spks from good clusters
            end
        end
        clear z
        spks_vHPC(:,2) = double(spks_vHPC(:,2))./20000;
        spks_dHPC(:,2) = double(spks_dHPC(:,2))./20000;
        % Selection of celltype to analyze
        if criteria_type == 0 %pyr
            cellulartype = [K(:,1) , K(:,3)];
        elseif criteria_type == 1 % int
            cellulartype = [K(:,1) , K(:,4)];
        elseif criteria_type == 2 % all
            cellulartype = [K(:,1) , ones(length(K),1)];
        end
        % Definition of the emotional transition
        if aversiveTS_run(1)<rewardTS_run(1)
            cond = 1; % 1 if the order was Aversive -> Reward
        else
            cond = 2;% 1 if the order was Reward -> Aversive
        end
        
        if and(not(isempty(group_vHPC)) , not(isempty(group_dHPC)))
            
            %% Neuronal Ensembles detection
            % --- Options for assemblies detection ---
            opts.Patterns.method = 'ICA';
            opts.threshold.method= 'MarcenkoPastur';
            opts.Patterns.number_of_iterations= 500;
            opts.threshold.permutations_percentile = 0.9;
            opts.threshold.number_of_permutations = 500;
            opts.Patterns.number_of_iterations = 500;
            %         opts.Members.method = 'Otsu';
            
            % --- Aversive ---
            limits = aversiveTS_run./1000;
            %         events = movement.aversive;
            events = [];
            [SpksTrains.dHPC.aversive , Bins.aversive , Cluster.dHPC.aversive] = spike_train_construction(spks_dHPC, group_dHPC(:,1), cellulartype, binSize, limits, events, false);
            [SpksTrains.vHPC.aversive , Bins.aversive , Cluster.vHPC.aversive] = spike_train_construction(spks_vHPC, group_vHPC(:,1), cellulartype, binSize, limits, events, false);
            
            [patterns.dHPC.aversive] = assembly_patterns([SpksTrains.dHPC.aversive'],opts);
            [patterns.vHPC.aversive] = assembly_patterns([SpksTrains.vHPC.aversive'],opts);
            %         [members.dHPC.aversive , patterns.dHPC.aversive] = assembly_patternsJFM([SpksTrains.dHPC.aversive'],opts);
            %         [members.vHPC.aversive , patterns.vHPC.aversive] = assembly_patternsJFM([SpksTrains.vHPC.aversive'],opts);
            
            % --- Reward ---
            limits = rewardTS_run./1000;
            %         events = movement.reward;
            events = [];
            [SpksTrains.dHPC.reward , Bins.reward , Cluster.dHPC.reward] = spike_train_construction(spks_dHPC, group_dHPC(:,1), cellulartype, binSize, limits, events, false);
            [SpksTrains.vHPC.reward , Bins.reward , Cluster.vHPC.reward] = spike_train_construction(spks_vHPC, group_vHPC(:,1), cellulartype, binSize, limits, events, false);
            
            [patterns.dHPC.reward] = assembly_patterns([SpksTrains.dHPC.reward'],opts);
            [patterns.vHPC.reward] = assembly_patterns([SpksTrains.vHPC.reward'],opts);
            %         [members.dHPC.reward , patterns.dHPC.reward] = assembly_patternsJFM([SpksTrains.dHPC.reward'],opts);
            %         [members.vHPC.reward , patterns.vHPC.reward] = assembly_patternsJFM([SpksTrains.vHPC.reward'],opts);
            
            if and(size(SpksTrains.dHPC.reward,2) > criteria_n , size(SpksTrains.vHPC.reward,2) > criteria_n )  % Checking of number of neurons per structure
                cond1 = and(not(isempty(patterns.dHPC.reward)) , not(isempty(patterns.dHPC.aversive)));
                cond2 = and(not(isempty(patterns.vHPC.reward)) , not(isempty(patterns.vHPC.aversive)));
                if and(cond1 , cond2)
                    %% SpikeTrains construction
                    limits = [0 segments.Var1(end)/1000];
                    events = [];
                    [SpksTrains.dHPC.all , bins , Cluster.dHPC.alls] = spike_train_construction(spks_dHPC, group_dHPC(:,1), cellulartype, binSize, limits, events, false);
                    [SpksTrains.vHPC.all , bins , Cluster.vHPC.all] = spike_train_construction(spks_vHPC, group_vHPC(:,1), cellulartype, binSize, limits, events, false);
                    
                    %% Assemblies activation in the entier recording
                    ActivityAssembles.dHPC.aversive = assembly_activity(patterns.dHPC.aversive , SpksTrains.dHPC.all');
                    ActivityAssembles.dHPC.reward = assembly_activity(patterns.dHPC.reward , SpksTrains.dHPC.all');
                    
                    ActivityAssembles.vHPC.aversive = assembly_activity(patterns.vHPC.aversive , SpksTrains.vHPC.all');
                    ActivityAssembles.vHPC.reward = assembly_activity(patterns.vHPC.reward , SpksTrains.vHPC.all');
                    
                    %% CCG using peaks of each assembly
                    % --- Peaks detection using the mean ---
                    peaks.dHPC.aversive = cell(1,size( ActivityAssembles.dHPC.aversive,1));
                    for i = 1 : size( ActivityAssembles.dHPC.aversive,1)
                        [p,pp] = findpeaks((ActivityAssembles.dHPC.aversive(i,:)));%,'Threshold',mean(ActivityAssembles.dHPC.aversive(i,:))+1*std(ActivityAssembles.dHPC.aversive(i,:)));
                        peaks.dHPC.aversive{i} = [pp' p'];
                    end
                    
                    peaks.dHPC.reward = cell(1,size( ActivityAssembles.dHPC.reward,1));
                    for i = 1 : size( ActivityAssembles.dHPC.reward,1)
                        [p,pp] = findpeaks((ActivityAssembles.dHPC.reward(i,:)));%,'Threshold',mean(ActivityAssembles.dHPC.reward(i,:))+1*std(ActivityAssembles.dHPC.reward(i,:)));
                        peaks.dHPC.reward{i} = [pp' p'];
                    end
                    
                    peaks.vHPC.aversive = cell(1,size( ActivityAssembles.vHPC.aversive,1));
                    for i = 1 : size( ActivityAssembles.vHPC.aversive,1)
                        [p,pp] = findpeaks((ActivityAssembles.vHPC.aversive(i,:)));%,'Threshold',mean(ActivityAssembles.vHPC.aversive(i,:))+1*std(ActivityAssembles.vHPC.aversive(i,:)));
                        peaks.vHPC.aversive{i} = [pp' p'];
                    end
                    
                    peaks.vHPC.reward = cell(1,size( ActivityAssembles.vHPC.reward,1));
                    for i = 1 : size( ActivityAssembles.vHPC.reward,1)
                        [p,pp] = findpeaks((ActivityAssembles.vHPC.reward(i,:)));%,'Threshold',mean(ActivityAssembles.vHPC.reward(i,:))+1*std(ActivityAssembles.vHPC.reward(i,:)));
                        peaks.vHPC.reward{i} = [pp' p'];
                    end
                                        
                    % --- CCG ---
                    tmp = [];
                    tmpA = [];
                    for i = 1 : size(ActivityAssembles.dHPC.aversive,1)
                        EventsD = bins(peaks.dHPC.aversive{i}(:,1));
                        for ii = 1 : size( ActivityAssembles.vHPC.aversive,1)
                            EventsV = bins(peaks.vHPC.aversive{ii}(:,1));
                            
                            % Baseline
                            x = Restrict(EventsD',ripple_event.baseline);
                            y = Restrict(EventsV',ripple_event.baseline);
                            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
                            [ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',0,'mode','ccg');
                            tmp1 = max(ccg(:,1,2)./sum(ccg(:,1,2)));
                            tmp1S = Threshold_CCG_Random(x , y , binSize, 2, 0, 'ccg', 100);
%                             tmp1 = max(ccg(:,1,2));
                            clear x y s ids groups ccg
                            
                            % Reward
                            x = Restrict(EventsD',ripple_event.reward);
                            y = Restrict(EventsV',ripple_event.reward);
                            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
                            [ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',0,'mode','ccg');
                            tmp2 = max(ccg(:,1,2)./sum(ccg(:,1,2)));
                            tmp2S = Threshold_CCG_Random(x , y , binSize, 2, 0, 'ccg', 100);
%                             tmp2 = max(ccg(:,1,2));
                            clear x y s ids groups ccg
                            
                            % Aversive
                            x = Restrict(EventsD',ripple_event.aversive);
                            y = Restrict(EventsV',ripple_event.aversive);
                            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
                            [ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',0,'mode','ccg');
                            tmp3 = max(ccg(:,1,2)./sum(ccg(:,1,2)));
                            tmp3S = Threshold_CCG_Random(x , y , binSize, 2, 0, 'ccg', 100);
%                             tmp3 = max(ccg(:,1,2));
                            clear x y s ids groups ccg                            
                            
                            cb = ttest(tmp1S,tmp1,'Tail','left');
                            cr = ttest(tmp2S,tmp2,'Tail','left');
                            ca = ttest(tmp3S,tmp3,'Tail','left');
                            
                            
                            if and(or(ca , cr) , not(cb))
                                tmp = [tmp ; tmp1 tmp2 tmp3];
                            end
                            
                            if and(and(ca , cr) , (cb))
                                tmpA = [tmpA ; tmp1 tmp2 tmp3];
                            end                            
                            
                            clear EventsV tmp1 tmp2 tmp3 cb ca cr
                        end
                        clear EventD
                    end
                    activity.assemblies.aversive.selective = [activity.assemblies.aversive.selective ; tmp]; clear tmp
                    activity.assemblies.aversive.all = [activity.assemblies.aversive.all ; tmpA]; clear tmpA

                    
                    tmp = [];
                    tmpA = [];
                    for i = 1 : size(ActivityAssembles.dHPC.reward,1)
                        EventsD = bins(peaks.dHPC.reward{i}(:,1));
                        for ii = 1 : size( ActivityAssembles.vHPC.reward,1)
                            EventsV = bins(peaks.vHPC.reward{ii}(:,1));
                            
                             % Baseline
                            x = Restrict(EventsD',ripple_event.baseline);
                            y = Restrict(EventsV',ripple_event.baseline);
                            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
                            [ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',0,'mode','ccg');
                            tmp1 = max(ccg(:,1,2)./sum(ccg(:,1,2)));
                            tmp1S = Threshold_CCG_Random(x , y , binSize, 2, 0, 'ccg', 100);
%                             tmp1 = max(ccg(:,1,2));
                            clear x y s ids groups ccg
                            
                            % Reward
                            x = Restrict(EventsD',ripple_event.reward);
                            y = Restrict(EventsV',ripple_event.reward);
                            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
                            [ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',0,'mode','ccg');
                            tmp2 = max(ccg(:,1,2)./sum(ccg(:,1,2)));
                            tmp2S = Threshold_CCG_Random(x , y , binSize, 2, 0, 'ccg', 100);
%                             tmp2 = max(ccg(:,1,2));
                            clear x y s ids groups ccg
                            
                            % Aversive
                            x = Restrict(EventsD',ripple_event.aversive);
                            y = Restrict(EventsV',ripple_event.aversive);
                            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
                            [ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',0,'mode','ccg');
                            tmp3 = max(ccg(:,1,2)./sum(ccg(:,1,2)));
                            tmp3S = Threshold_CCG_Random(x , y , binSize, 2, 0, 'ccg', 100);
%                             tmp3 = max(ccg(:,1,2));
                            clear x y s ids groups ccg
                            
                            cb = ttest(tmp1S,tmp1,'Tail','left');
                            cr = ttest(tmp2S,tmp2,'Tail','left');
                            ca = ttest(tmp3S,tmp3,'Tail','left');
                            
                            if and(or(ca , cr) , not(cb))
                                tmp = [tmp ; tmp1 tmp2 tmp3];
                            end
                            
                            if and(and(ca , cr) , (cb))
                                tmpA = [tmpA ; tmp1 tmp2 tmp3];
                            end                            
                            
                            clear EventsV tmp1 tmp2 tmp3 cb ca cr
                            
                        end
                        clear EventD
                    end
                    activity.assemblies.reward.selective = [activity.assemblies.reward.selective ; tmp]; clear tmp
                    activity.assemblies.reward.all = [activity.assemblies.reward.all ; tmpA]; clear tmpA
                end
            end
        end
        disp(['-- Analysis of folder #' , num2str(t) , ' from rat #',num2str(tt) , ' finished --'])
        disp(' ')
        clear spiketrains_dHPC_int spiketrains_dHPC_pyr spiketrains_vHPC_int spiketrains_vHPC_pyr
        clear aversiveTS aversiveTS_run rewardTS rewardTS_run behavior bins Cell_type_classification cellulartype
        clear group_dHPC group_vHPC is K Kinfo NREM REM WAKE ripplesD ripplesV segments
        clear spiketrains_dHPC spiketrains_vHPC ripple_bursts ripple_event baselineTS cond
        clear coordinated coordinatedV coordinatedV_refined movement spks
        clear spks_dHPC spks_vHPC patterns ActivityAssembles members opts
        clear SpksTrains cond1 cond2
    end
end

figure
subplot(121),boxplot(activity.assemblies.aversive.all)%,ylim([0.08 0.2])
[p ,t , stats] = kruskalwallis(activity.assemblies.aversive.all)
[c , m , h , nms] = multcompare(stats,'ctype','bonferroni')

subplot(122),boxplot(activity.assemblies.reward.all)%,ylim([0.08 0.2])
[p ,t , stats] = kruskalwallis(activity.assemblies.reward.all)
[c , m , h , nms] = multcompare(stats,'ctype','bonferroni')

figure
subplot(121),boxplot(activity.assemblies.aversive.selective)%,ylim([0.08 0.2])
[p ,t , stats] = kruskalwallis(activity.assemblies.aversive.selective)
[c , m , h , nms] = multcompare(stats,'ctype','bonferroni')

subplot(122),boxplot(activity.assemblies.reward.selective)%,ylim([0.08 0.2])
[p ,t , stats] = kruskalwallis(activity.assemblies.reward.selective)
[c , m , h , nms] = multcompare(stats,'ctype','bonferroni')