clear
clc
close all

%% Parameters
path = {'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat103\usable';'E:\Rat132\recordings\in_pyr'};%List of folders from the path

% for SU
criteria_fr = 0.02; %criteria to include or not a SU into the analysis
criteria_n = 6; % minimal number of neurons from each structure
criteria_type = 0; %criteria for celltype (0:pyr, 1:int, 2:all)
n_SU_V = [];
n_SU_D = [];
Xedges = 60; %number of bins for RateMap construction

binSize = 0.02; % bin size for replay events detection

% Behavior
minimal_speed = 5; % minimal speed to detect active periods
minimal_speed_time = 2; % minimal time to detect active periods

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
        quiet.reward = Restrict(behavior.quiet.reward , [start , stop]);
        movement.reward(movement.reward(:,2) - movement.reward(:,1) <1,:)=[];
        clear tmp start stop
        % Aversive
        start = behavior.speed.aversive(1,1);   stop = behavior.speed.aversive(end,1);
        movement.aversive = InvertIntervals(behavior.quiet.aversive , start , stop);%keep only those higher than criteria
        quiet.aversive = Restrict(behavior.quiet.aversive , [start , stop]);
        movement.aversive(movement.aversive(:,2) - movement.aversive(:,1) <1,:)=[];
        clear tmp start stop
        
        %% Shocks selection
        Shocks = Restrict(shock,aversiveTS_run ./1000);
        % Keep only the first shock of each TTL (first from 20)
        count = 1;
        deff = [];
        for i = 1:length(Shocks)
            if count == 1
                deff = [deff ; Shocks(i,1)];
            end
            if count ==20
                count = 0;
            end
            count = count + 1;
        end
        Shocks = deff;
        clear count deff shock i
        
        %Rewards selection
        [i ii] = sort([leftvalve(:,1) ; rightvalve(:,1)]);
        Reward = [leftvalve ; rightvalve];
        Reward = Restrict(Reward(ii,:),rewardTS_run ./1000);
        clear i ii
        
        %% load sleep states
        disp('Uploading sleep scoring')
        x = dir([cd,'\*-states.mat']);    states = load([cd,'\',x.name]);    states = states.states;
        REM.all = ToIntervals(states==5);    NREM.all = ToIntervals(states==3);    WAKE.all = ToIntervals(states==1);
        clear x states
        NREM.baseline = Restrict(NREM.all,baselineTS./1000);
        NREM.aversive = Restrict(NREM.all,aversiveTS./1000);
        NREM.reward = Restrict(NREM.all,rewardTS./1000);
        REM.baseline = Restrict(REM.all,baselineTS./1000);
        REM.aversive = Restrict(REM.all,aversiveTS./1000);
        REM.reward = Restrict(REM.all,rewardTS./1000);
        
        %% load coordinated ripple bursts
        load('coordinated_ripple_bursts.mat')
%         coordinated_ripple_bursts(coordinated_ripple_bursts(:,3) - coordinated_ripple_bursts(:,1) < 0.1,:) = [];
        ripple_bursts.baseline = Restrict(coordinated_ripple_bursts,baselineTS./1000);
        ripple_bursts.reward = Restrict(coordinated_ripple_bursts,rewardTS./1000);
        ripple_bursts.aversive = Restrict(coordinated_ripple_bursts,aversiveTS./1000);
        ripple_bursts.all = coordinated_ripple_bursts;
        clear coordinated_ripple_bursts
        
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
        clear cooridnated_event
        
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
        
        % keep only clusters that are the interested cellular type
        C = [];
        for ii=1:size(group_dHPC,1)
            cluster = group_dHPC(ii,1);
            celltype = logical(cellulartype(cellulartype(:,1) == cluster,2));
            if celltype
                C = [C ; cluster];
            end
        end
        group_dHPC = C;  clear C
        % keep only clusters that are the interested cellular type
        
        C = [];
        for ii=1:size(group_vHPC,1)
            cluster = group_vHPC(ii,1);
            celltype = logical(cellulartype(cellulartype(:,1) == cluster,2));
            if celltype
                C = [C ; cluster];
            end
        end
        group_vHPC = C;  clear C
        
        %% Neuronal Ensembles detection
        % --- Aversive ---
        limits = aversiveTS_run./1000;
        events = [];
        [SpksTrains.dHPC.aversive , Bins.aversive , Cluster.dHPC.aversive] = spike_train_construction(spks_dHPC, group_dHPC(:,1), cellulartype, binSize, limits, events, false);
        [SpksTrains.vHPC.aversive , Bins.aversive , Cluster.vHPC.aversive] = spike_train_construction(spks_vHPC, group_vHPC(:,1), cellulartype, binSize, limits, events, false);

        [assemblies.dHPC.aversive , patterns.dHPC.aversive] = assembly_patterns(SpksTrains.dHPC.aversive');
        [assemblies.vHPC.aversive , patterns.vHPC.aversive] = assembly_patterns(SpksTrains.vHPC.aversive');
        
        % Otsu's thresholding of the assemblie's members
        patterns.dHPC.aversive(:,sum(assemblies.dHPC.aversive)<2) = []; % if no SU pass through the threshold, those assemblies are delated
        assemblies.dHPC.aversive(:,sum(assemblies.dHPC.aversive)<2) = [];
        
        patterns.vHPC.aversive(:,sum(assemblies.vHPC.aversive)<2) = []; % if no SU pass through the threshold, those assemblies are delated
        assemblies.vHPC.aversive(:,sum(assemblies.vHPC.aversive)<2) = [];        

        % --- Reward ---
        limits = rewardTS_run./1000;
        events = [];
        [SpksTrains.dHPC.reward , Bins.reward , Cluster.dHPC.reward] = spike_train_construction(spks_dHPC, group_dHPC(:,1), cellulartype, 0.02, limits, events, false);
        [SpksTrains.vHPC.reward , Bins.reward , Cluster.vHPC.reward] = spike_train_construction(spks_vHPC, group_vHPC(:,1), cellulartype, 0.02, limits, events, false);

        [assemblies.dHPC.reward , patterns.dHPC.reward] = assembly_patterns(SpksTrains.dHPC.reward');
        [assemblies.vHPC.reward , patterns.vHPC.reward] = assembly_patterns(SpksTrains.vHPC.reward');
        
        patterns.dHPC.reward(:,sum(assemblies.dHPC.reward)<2) = []; % if no SU pass through the threshold, those assemblies are delated
        assemblies.dHPC.reward(:,sum(assemblies.dHPC.reward)<2) = [];
        
        patterns.vHPC.reward(:,sum(assemblies.vHPC.aversive)<2) = []; % if no SU pass through the threshold, those assemblies are delated
        assemblies.vHPC.reward(:,sum(assemblies.vHPC.aversive)<2) = [];      

        %% Assemblies activation in the entier recording
        limits = [0 segments.Var1(end)/1000];
        events = [];
        [SpksTrains.dHPC.all , Bins.all , Cluster.dHPC.alls] = spike_train_construction(spks_dHPC, group_dHPC(:,1), cellulartype, binSize, limits, events, false);
        [SpksTrains.vHPC.all , Bins.all , Cluster.vHPC.all] = spike_train_construction(spks_vHPC, group_vHPC(:,1), cellulartype, binSize, limits, events, false);
        
        ActivityAssembles.dHPC.aversive = assembly_activity(patterns.dHPC.aversive , SpksTrains.dHPC.all');
        ActivityAssembles.vHPC.aversive = assembly_activity(patterns.vHPC.aversive , SpksTrains.vHPC.all');
        
        ActivityAssembles.dHPC.reward = assembly_activity(patterns.dHPC.reward , SpksTrains.dHPC.all');
        ActivityAssembles.vHPC.reward = assembly_activity(patterns.vHPC.reward , SpksTrains.vHPC.all');
        

        %% Calculation of firing rate of each assemblie by using the Spikes of each member
        Spikes.dHPC.aversive = cell(2,size(assemblies.dHPC.aversive,2));
        for i = 1 : size(assemblies.dHPC.aversive,2)
            clusters = group_dHPC(logical(assemblies.dHPC.aversive(:,i)));
            Spikes.dHPC.aversive{1,i} = spks_dHPC(ismember(spks_dHPC(:,1) , clusters),2);
            [Spikes.dHPC.aversive{2,i},bins]=binspikes(Spikes.dHPC.aversive{1,i},1/binSize,limits);
            clear clusters
        end
        clear i
        
        Spikes.dHPC.reward = cell(2,size(assemblies.dHPC.reward,2));
        for i = 1 : size(assemblies.dHPC.reward,2)
            clusters = group_dHPC(logical(assemblies.dHPC.reward(:,i)));
            Spikes.dHPC.reward{1,i} = spks_dHPC(ismember(spks_dHPC(:,1) , clusters),2);
            [Spikes.dHPC.reward{2,i},bins]=binspikes(Spikes.dHPC.reward{1,i},1/binSize,limits);
            clear clusters
        end
        clear i 
        
        Spikes.vHPC.aversive = cell(2,size(assemblies.vHPC.aversive,2));
        for i = 1 : size(assemblies.vHPC.aversive,2)
            clusters = group_vHPC(logical(assemblies.vHPC.aversive(:,i)));
            Spikes.vHPC.aversive{1,i} = spks_vHPC(ismember(spks_vHPC(:,1) , clusters),2);
            [Spikes.vHPC.aversive{2,i},bins]=binspikes(Spikes.vHPC.aversive{1,i},1/binSize,limits);
            clear clusters
        end
        clear i
        
        Spikes.vHPC.reward = cell(2,size(assemblies.vHPC.reward,2));
        for i = 1 : size(assemblies.vHPC.reward,2)
            clusters = group_vHPC(logical(assemblies.vHPC.reward(:,i)));
            Spikes.vHPC.reward{1,i} = spks_vHPC(ismember(spks_vHPC(:,1) , clusters),2);
            [Spikes.vHPC.reward{2,i},bins]=binspikes(Spikes.vHPC.reward{1,i},1/binSize,limits);
            clear clusters
        end
        clear i 
        
        %% PHIST during ripple bursts
            for i = 1 : size(assemblies.vHPC.aversive,2)
                x = Restrict(Spikes.vHPC.aversive{1,i},NREM.baseline) ;
                y = ripple_bursts.baseline(:,2);
                [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
                [ccg,tttt] = CCG(s,ids,'binSize',binSize,'duration',1,'smooth',1,'mode','ccg'); %ccg calculation
                figure,plot(tttt,ccg(:,1,2)./sum(ccg(:,1,2)),'k'),hold on
                clear x y ccg tttt s ids groups

                
                x = Restrict(Spikes.vHPC.aversive{1,i},NREM.reward) ;
                y = ripple_bursts.reward(:,2);
                [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
                [ccg,tttt] = CCG(s,ids,'binSize',binSize,'duration',1,'smooth',1,'mode','ccg'); %ccg calculation
                plot(tttt,ccg(:,1,2)./sum(ccg(:,1,2)),'b'),hold on
                clear x y ccg tttt s ids groups
                
                x = Restrict(Spikes.vHPC.aversive{1,i},NREM.aversive) ;
                y = ripple_bursts.aversive(:,2);
                [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
                [ccg,tttt] = CCG(s,ids,'binSize',binSize,'duration',1,'smooth',1,'mode','ccg'); %ccg calculation
                plot(tttt,ccg(:,1,2)./sum(ccg(:,1,2)),'r'),hold on
                clear x y ccg tttt s ids groups
        end
                
            for i = 1 : size(assemblies.dHPC.aversive,2)
                x = Restrict(Spikes.dHPC.aversive{1,i},NREM.baseline) ;
                y = ripple_bursts.baseline(:,2);;
                [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
                [ccg,tttt] = CCG(s,ids,'binSize',binSize,'duration',1,'smooth',1,'mode','ccg'); %ccg calculation
                figure,plot(tttt,ccg(:,1,2)./sum(ccg(:,1,2)),'k'),hold on
                clear x y ccg tttt s ids groups

                
                x = Restrict(Spikes.dHPC.aversive{1,i},NREM.reward) ;
                y = ripple_bursts.reward(:,2);
                [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
                [ccg,tttt] = CCG(s,ids,'binSize',binSize,'duration',1,'smooth',1,'mode','ccg'); %ccg calculation
                plot(tttt,ccg(:,1,2)./sum(ccg(:,1,2)),'b'),hold on
                clear x y ccg tttt s ids groups
                
                x = Restrict(Spikes.dHPC.aversive{1,i},NREM.aversive) ;
                y = ripple_bursts.aversive(:,2);
                [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
                [ccg,tttt] = CCG(s,ids,'binSize',binSize,'duration',1,'smooth',1,'mode','ccg'); %ccg calculation
                plot(tttt,ccg(:,1,2)./sum(ccg(:,1,2)),'r'),hold on
                clear x y ccg tttt s ids groups
            end        
            
            for i = 1 : size(assemblies.vHPC.reward,2)
                x = Restrict(Spikes.vHPC.reward{1,i},NREM.baseline) ;
                y = ripple_bursts.baseline(:,2);
                [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
                [ccg,tttt] = CCG(s,ids,'binSize',binSize,'duration',1,'smooth',1,'mode','ccg'); %ccg calculation
                figure,plot(tttt,ccg(:,1,2)./sum(ccg(:,1,2)),'k'),hold on
                clear x y ccg tttt s ids groups

                
                x = Restrict(Spikes.vHPC.reward{1,i},NREM.reward) ;
                y = ripple_bursts.reward(:,2);
                [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
                [ccg,tttt] = CCG(s,ids,'binSize',binSize,'duration',1,'smooth',1,'mode','ccg'); %ccg calculation
                plot(tttt,ccg(:,1,2)./sum(ccg(:,1,2)),'b'),hold on
                clear x y ccg tttt s ids groups
                
                x = Restrict(Spikes.vHPC.reward{1,i},NREM.aversive) ;
                y = ripple_bursts.aversive(:,2);
                [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
                [ccg,tttt] = CCG(s,ids,'binSize',binSize,'duration',1,'smooth',1,'mode','ccg'); %ccg calculation
                plot(tttt,ccg(:,1,2)./sum(ccg(:,1,2)),'r'),hold on
                clear x y ccg tttt s ids groups
        end
                
            for i = 1 : size(assemblies.dHPC.reward,2)
                x = Restrict(Spikes.dHPC.reward{1,i},NREM.baseline) ;
                y = ripple_bursts.baseline(:,2);
                [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
                [ccg,tttt] = CCG(s,ids,'binSize',binSize,'duration',1,'smooth',1,'mode','ccg'); %ccg calculation
                figure,plot(tttt,ccg(:,1,2)./sum(ccg(:,1,2)),'k'),hold on
                clear x y ccg tttt s ids groups

                
                x = Restrict(Spikes.dHPC.reward{1,i},NREM.reward) ;
                y = ripple_bursts.reward(:,2);
                [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
                [ccg,tttt] = CCG(s,ids,'binSize',binSize,'duration',1,'smooth',1,'mode','ccg'); %ccg calculation
                plot(tttt,ccg(:,1,2)./sum(ccg(:,1,2)),'b'),hold on
                clear x y ccg tttt s ids groups
                
                x = Restrict(Spikes.dHPC.reward{1,i},NREM.aversive) ;
                y = ripple_bursts.aversive(:,2);
                [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
                [ccg,tttt] = CCG(s,ids,'binSize',binSize,'duration',1,'smooth',1,'mode','ccg'); %ccg calculation
                plot(tttt,ccg(:,1,2)./sum(ccg(:,1,2)),'r'),hold on
                clear x y ccg tttt s ids groups
            end   
            
        %% Xcorr of Assemblies activity
        
        id.burst.baseline = InIntervals(Bins.all,[ripple_bursts.baseline(:,1) ripple_bursts.baseline(:,3)]);
        id.burst.reward = InIntervals(Bins.all,[ripple_bursts.reward(:,1) ripple_bursts.reward(:,3)]);
        id.burst.aversive = InIntervals(Bins.all,[ripple_bursts.aversive(:,1) ripple_bursts.aversive(:,3)]);
        
        id.event.baseline = InIntervals(Bins.all,ripple_event.baseline);
        id.event.reward = InIntervals(Bins.all,ripple_event.reward);
        id.event.aversive = InIntervals(Bins.all,ripple_event.aversive);
        
        is.NREM.baseline = InIntervals(Bins.all,NREM.baseline);
        is.NREM.reward = InIntervals(Bins.all,NREM.reward);
        is.NREM.aversive = InIntervals(Bins.all,NREM.aversive);

        is.REM.baseline = InIntervals(Bins.all,REM.baseline);
        is.REM.reward = InIntervals(Bins.all,REM.reward);
        is.REM.aversive = InIntervals(Bins.all,REM.aversive);
        
        %% 
        
        i = InIntervals(Bins.all,[ripplesD(:,2)-0.1 ripplesD(:,2)+0.1]);
        [c , l] = xcorr(ActivityAssembles.dHPC.aversive(1,i) , ActivityAssembles.vHPC.aversive(2,i),25,'normalized');
        
        %% 
        for i = 1 : size(assemblies.dHPC.aversive,2)
%             a1 = zscore(ActivityAssembles.dHPC.aversive(i,:));
            for ii = 1 : size(assemblies.vHPC.aversive,2)
%                 a2 = zscore(ActivityAssembles.vHPC.aversive(ii,:));
                x = Restrict(Spikes.dHPC.aversive{1,i},ripple_event.baseline) ;
                y = Restrict(Spikes.vHPC.aversive{1,ii},ripple_event.baseline);
                [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
                [ccg,tttt] = CCG(s,ids,'binSize',binSize,'duration',1,'smooth',1,'mode','ccv'); %ccg calculation
                figure,plot(tttt,ccg(:,1,2),'k'),hold on
                clear x y ccg tttt s ids groups

                
                x = Restrict(Spikes.dHPC.aversive{1,i},ripple_event.reward) ;
                y = Restrict(Spikes.vHPC.aversive{1,ii},ripple_event.reward);
                [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
                [ccg,tttt] = CCG(s,ids,'binSize',binSize,'duration',1,'smooth',1,'mode','ccv'); %ccg calculation
                plot(tttt,ccg(:,1,2),'b'),hold on
                clear x y ccg tttt s ids groups
                
                x = Restrict(Spikes.dHPC.aversive{1,i},ripple_event.aversive) ;
                y = Restrict(Spikes.vHPC.aversive{1,ii},ripple_event.aversive);
                [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
                [ccg,tttt] = CCG(s,ids,'binSize',binSize,'duration',1,'smooth',1,'mode','ccv'); %ccg calculation
                plot(tttt,ccg(:,1,2),'r'),hold on
                clear x y ccg tttt s ids groups
            end
            clear a1
        end
                
        
        for i = 1 : size(assemblies.dHPC.reward,2)
%             a1 = zscore(ActivityAssembles.dHPC.aversive(i,:));
            for ii = 1 : size(assemblies.vHPC.reward,2)
%                 a2 = zscore(ActivityAssembles.vHPC.aversive(ii,:));
                x = Restrict(Spikes.dHPC.reward{1,i},ripple_event.baseline) ;
                y = Restrict(Spikes.vHPC.reward{1,ii},ripple_event.baseline);
                [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
                [ccg,tttt] = CCG(s,ids,'binSize',binSize,'duration',1,'smooth',1,'mode','ccg'); %ccg calculation
                figure,plot(tttt,ccg(:,1,2)./sum(ccg(:,1,2)),'k'),hold on
                clear x y ccg tttt s ids groups

                
                x = Restrict(Spikes.dHPC.reward{1,i},ripple_event.reward) ;
                y = Restrict(Spikes.vHPC.reward{1,ii},ripple_event.reward);
                [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
                [ccg,tttt] = CCG(s,ids,'binSize',binSize,'duration',1,'smooth',1,'mode','ccg'); %ccg calculation
                plot(tttt,ccg(:,1,2)./sum(ccg(:,1,2)),'b'),hold on
                clear x y ccg tttt s ids groups
                
                x = Restrict(Spikes.dHPC.reward{1,i},ripple_event.aversive) ;
                y = Restrict(Spikes.vHPC.reward{1,ii},ripple_event.aversive);
                [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
                [ccg,tttt] = CCG(s,ids,'binSize',binSize,'duration',1,'smooth',1,'mode','ccg'); %ccg calculation
                plot(tttt,ccg(:,1,2)./sum(ccg(:,1,2)),'r'),hold on
                clear x y ccg tttt s ids groups
            end
            clear a1
        end
        
    
        %% Ripple perievents histogram 
        % All events
        perievent.dHPC.aversive.all.baseline = [];
        perievent.dHPC.aversive.all.reward = [];
        perievent.dHPC.aversive.all.aversive = [];
        for ii = 1 : size(assemblies.dHPC.aversive,2)
            %Baseline        
            x = ripples.dHPC.baseline(:,2);
            y = Spikes.dHPC.aversive{1,ii};
            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,tttt] = CCG(s,ids,'binSize',binSize,'duration',1,'smooth',1,'mode','ccg'); %ccg calculation
            
            Ind = SubtractIntervals(baselineTS./1000,[x-0.05 x+0.05]); % mean FR for Gain calculation
            FR = sum(InIntervals(y,Ind)) / sum(Ind(:,2)-Ind(:,1));
            
            perievent.dHPC.aversive.all.baseline = [perievent.dHPC.aversive.all.baseline , ccg(:,1,2)./binSize./length(x)./FR];
            figure, plot(tttt,ccg(:,1,2)./binSize./length(x)./FR,'k'),hold on
            clear s ids groups ccg tttt FR Ind x y
            
            %Reward
            x = ripples.dHPC.reward(:,2);
            y = Spikes.dHPC.aversive{1,ii};
            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,tttt] = CCG(s,ids,'binSize',binSize,'duration',1,'smooth',1,'mode','ccg'); %ccg calculation
            
            Ind = SubtractIntervals(rewardTS./1000,[x-0.05 x+0.05]);
            FR = sum(InIntervals(y,Ind)) / sum(Ind(:,2)-Ind(:,1));
            
            perievent.dHPC.all.cooridnated.reward = [perievent.dHPC.aversive.all.reward , ccg(:,1,2)./binSize./length(x)./FR];
            plot(tttt,ccg(:,1,2)./binSize./length(x)./FR,'b'),hold on
            clear s ids groups ccg tttt FR Ind x y
            
            %Aversive
            x = ripples.dHPC.aversive(:,2);
            y = Spikes.dHPC.aversive{1,ii};
            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,tttt] = CCG(s,ids,'binSize',binSize,'duration',1,'smooth',1,'mode','ccg'); %ccg calculation
            
            
            Ind = SubtractIntervals(aversiveTS./1000,[x-0.05 x+0.05]);
            FR = sum(InIntervals(y,Ind)) / sum(Ind(:,2)-Ind(:,1));
            
            perievent.dHPC.aversive.all.aversive = [perievent.dHPC.aversive.all.reward , ccg(:,1,2)./binSize./length(x)./FR];
            plot(tttt,ccg(:,1,2)./binSize./length(x)./FR,'r'),hold on
            clear s ids groups ccg tttt FR Ind x y
        end

        perievent.vHPC.aversive.all.baseline = [];
        perievent.vHPC.aversive.all.reward = [];
        perievent.vHPC.aversive.all.aversive = [];
        for ii = 1 : size(assemblies.vHPC.aversive,2)
            %Baseline
            x = ripples.vHPC.baseline(:,2);
            y = Spikes.vHPC.aversive{1,ii};
            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,tttt] = CCG(s,ids,'binSize',binSize,'duration',1,'smooth',1,'mode','ccg'); %ccg calculation
            
            Ind = SubtractIntervals(baselineTS./1000,[x-0.05 x+0.05]); % mean FR for Gain calculation
            FR = sum(InIntervals(y,Ind)) / sum(Ind(:,2)-Ind(:,1));
                        
            perievent.vHPC.aversive.all.baseline = [perievent.vHPC.aversive.all.baseline , ccg(:,1,2)./binSize./length(x)./FR];
            figure, plot(tttt,ccg(:,1,2)./binSize./length(x)./FR,'k'),hold on
            clear s ids groups ccg tttt x y FR Ind
            
            %Reward
            x = ripples.vHPC.reward(:,2);
            y = Spikes.vHPC.aversive{1,ii};
            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,tttt] = CCG(s,ids,'binSize',binSize,'duration',1,'smooth',1,'mode','ccg'); %ccg calculation
            
            Ind = SubtractIntervals(rewardTS./1000,[x-0.05 x+0.05]); % mean FR for Gain calculation
            FR = sum(InIntervals(y,Ind)) / sum(Ind(:,2)-Ind(:,1));
                        
            perievent.vHPC.aversive.all.reward = [perievent.vHPC.aversive.all.reward , ccg(:,1,2)./binSize./length(x)./FR];
            plot(tttt,ccg(:,1,2)./binSize./length(x)./FR,'b'),hold on
            clear s ids groups ccg tttt x y FR Ind
            
            %Aversive
            x = ripples.vHPC.aversive(:,2);
            y = Spikes.vHPC.aversive{1,ii};
            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,tttt] = CCG(s,ids,'binSize',binSize,'duration',1,'smooth',1,'mode','ccg'); %ccg calculation
            
            Ind = SubtractIntervals(aversiveTS./1000,[x-0.05 x+0.05]); % mean FR for Gain calculation
            FR = sum(InIntervals(y,Ind)) / sum(Ind(:,2)-Ind(:,1));
                        
            perievent.vHPC.aversive.all.aversive = [perievent.vHPC.aversive.all.reward , ccg(:,1,2)./binSize./length(x)./FR];
            plot(tttt,ccg(:,1,2)./binSize./length(x)./FR,'r'),hold on
            clear s ids groups ccg tttt FR Ind x y
        end
        
        perievent.dHPC.reward.all.baseline = [];
        perievent.dHPC.reward.all.reward = [];
        perievent.dHPC.reward.all.aversive = [];
        for ii = 1 : size(assemblies.dHPC.reward,2)
            %Baseline
            x = ripples.dHPC.baseline(:,2);
            y = Spikes.dHPC.reward{1,ii};
            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,tttt] = CCG(s,ids,'binSize',binSize,'duration',1,'smooth',1,'mode','ccg'); %ccg calculation
            
            Ind = SubtractIntervals(baselineTS./1000,[x-0.05 x+0.05]); % mean FR for Gain calculation
            FR = sum(InIntervals(y,Ind)) / sum(Ind(:,2)-Ind(:,1));
            
            perievent.dHPC.reward.all.baseline = [perievent.dHPC.reward.all.baseline , ccg(:,1,2)./binSize./length(x)./FR];
            figure, plot(tttt,ccg(:,1,2)./binSize./length(x)./FR,'k'),hold on
            clear s ids groups ccg tttt x y FR Ind
            
            %Reward
            x = ripples.dHPC.reward(:,2);
            y = Spikes.dHPC.reward{1,ii};
            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,tttt] = CCG(s,ids,'binSize',binSize,'duration',1,'smooth',1,'mode','ccg'); %ccg calculation
            
            Ind = SubtractIntervals(rewardTS./1000,[x-0.05 x+0.05]); % mean FR for Gain calculation
            FR = sum(InIntervals(y,Ind)) / sum(Ind(:,2)-Ind(:,1));
            
            perievent.dHPC.reward.all.reward = [perievent.dHPC.reward.all.reward , ccg(:,1,2)./binSize./length(x)./FR];
            plot(tttt,ccg(:,1,2)./binSize./length(x)./FR,'b'),hold on
            clear s ids groups ccg tttt x y FR Ind
            
            %Aversive
            x = ripples.dHPC.aversive(:,2);
            y = Spikes.dHPC.reward{1,ii};
            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,tttt] = CCG(s,ids,'binSize',binSize,'duration',1,'smooth',1,'mode','ccg'); %ccg calculation
            
            Ind = SubtractIntervals(aversiveTS./1000,[x-0.05 x+0.05]); % mean FR for Gain calculation
            FR = sum(InIntervals(y,Ind)) / sum(Ind(:,2)-Ind(:,1));
            
            perievent.dHPC.reward.all.aversive = [perievent.dHPC.reward.all.reward , ccg(:,1,2)./binSize./length(x)./FR];
            plot(tttt,ccg(:,1,2)./binSize./length(x)./FR,'r'),hold on
            clear s ids groups ccg tttt x y FR Ind
        end

        perievent.vHPC.reward.all.baseline = [];
        perievent.vHPC.reward.all.reward = [];
        perievent.vHPC.reward.all.aversive = [];
        for ii = 1 : size(assemblies.vHPC.reward,2)
            x = ripples.vHPC.baseline(:,2);
            y = Spikes.vHPC.reward{1,ii};
            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,tttt] = CCG(s,ids,'binSize',binSize,'duration',1,'smooth',1,'mode','ccg'); %ccg calculation
            
            Ind = SubtractIntervals(baselineTS./1000,[x-0.05 x+0.05]); % mean FR for Gain calculation
            FR = sum(InIntervals(y,Ind)) / sum(Ind(:,2)-Ind(:,1));      
            
            perievent.vHPC.reward.all.baseline = [perievent.vHPC.reward.all.baseline , ccg(:,1,2)./binSize./length(x)./FR];
            figure, plot(tttt,ccg(:,1,2)./binSize./length(x)./FR,'k'),hold on
            clear s ids groups ccg tttt x y FR Ind
            
            %Reward
            x = ripples.vHPC.reward(:,2);
            y = Spikes.vHPC.reward{1,ii};
            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,tttt] = CCG(s,ids,'binSize',binSize,'duration',1,'smooth',1,'mode','ccg'); %ccg calculation
            
            Ind = SubtractIntervals(rewardTS./1000,[x-0.05 x+0.05]); % mean FR for Gain calculation
            FR = sum(InIntervals(y,Ind)) / sum(Ind(:,2)-Ind(:,1));             
            
            perievent.vHPC.reward.all.reward = [perievent.vHPC.reward.all.reward , ccg(:,1,2)./binSize./length(x)./FR];
            plot(tttt,ccg(:,1,2)./binSize./length(x)./FR,'b'),hold on
            clear s ids groups ccg tttt x y FR Ind
            
            %Aversive
            x = ripples.vHPC.aversive(:,2);
            y = Spikes.vHPC.reward{1,ii};
            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,tttt] = CCG(s,ids,'binSize',binSize,'duration',1,'smooth',1,'mode','ccg'); %ccg calculation
            
            Ind = SubtractIntervals(aversiveTS./1000,[x-0.05 x+0.05]); % mean FR for Gain calculation
            FR = sum(InIntervals(y,Ind)) / sum(Ind(:,2)-Ind(:,1)); 
            
            perievent.vHPC.reward.all.aversive = [perievent.vHPC.reward.all.reward , ccg(:,1,2)./binSize./length(x)./FR];
            plot(tttt,ccg(:,1,2)./binSize./length(x)./FR,'r'),hold on
            clear s ids groups ccg tttt x y FR Ind
        end
        
        % Coordinated events
        perievent.dHPC.aversive.cooridnated.baseline = [];
        perievent.dHPC.aversive.cooridnated.reward = [];
        perievent.dHPC.aversive.cooridnated.aversive = [];
        for ii = 1 : size(assemblies.dHPC.aversive,2)
            x = ripples.dHPC.coordinated.baseline(:,2);
            y = Spikes.dHPC.aversive{1,ii};
            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,tttt] = CCG(s,ids,'binSize',binSize,'duration',1,'smooth',1,'mode','ccg'); %ccg calculation
            perievent.dHPC.aversive.cooridnated.baseline = [perievent.dHPC.aversive.cooridnated.baseline , ccg(:,1,2)./binSize./length(x)];
            figure, plot(tttt,ccg(:,1,2)./binSize./length(x),'k'),hold on
            clear s ids groups ccg tttt
            
            x = ripples.dHPC.coordinated.reward(:,2);
            y = Spikes.dHPC.aversive{1,ii};
            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,tttt] = CCG(s,ids,'binSize',binSize,'duration',1,'smooth',1,'mode','ccg'); %ccg calculation
            perievent.dHPC.aversive.cooridnated.reward = [perievent.dHPC.aversive.cooridnated.reward , ccg(:,1,2)./binSize./length(x)];
            plot(tttt,ccg(:,1,2)./binSize./length(x),'b'),hold on
            clear s ids groups ccg tttt
            
            x = ripples.dHPC.coordinated.aversive(:,2);
            y = Spikes.dHPC.aversive{1,ii};
            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,tttt] = CCG(s,ids,'binSize',binSize,'duration',1,'smooth',1,'mode','ccg'); %ccg calculation
            perievent.dHPC.aversive.cooridnated.aversive = [perievent.dHPC.aversive.cooridnated.reward , ccg(:,1,2)./binSize./length(x)];
            plot(tttt,ccg(:,1,2)./binSize./length(x),'r'),hold on
            clear s ids groups ccg tttt
        end

        perievent.vHPC.aversive.cooridnated.baseline = [];
        perievent.vHPC.aversive.cooridnated.reward = [];
        perievent.vHPC.aversive.cooridnated.aversive = [];
        for ii = 1 : size(assemblies.vHPC.aversive,2)
            x = ripples.vHPC.coordinated.baseline(:,2);
            y = Spikes.vHPC.aversive{1,ii};
            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,tttt] = CCG(s,ids,'binSize',binSize,'duration',1,'smooth',1,'mode','ccg'); %ccg calculation
            perievent.vHPC.aversive.cooridnated.baseline = [perievent.vHPC.aversive.cooridnated.baseline , ccg(:,1,2)./binSize./length(x)];
            figure, plot(tttt,ccg(:,1,2)./binSize./length(x),'k'),hold on
            clear s ids groups ccg tttt
            
            x = ripples.vHPC.coordinated.reward(:,2);
            y = Spikes.vHPC.aversive{1,ii};
            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,tttt] = CCG(s,ids,'binSize',binSize,'duration',1,'smooth',1,'mode','ccg'); %ccg calculation
            perievent.vHPC.aversive.cooridnated.reward = [perievent.vHPC.aversive.cooridnated.reward , ccg(:,1,2)./binSize./length(x)];
            plot(tttt,ccg(:,1,2)./binSize./length(x),'b'),hold on
            clear s ids groups ccg tttt
            
            x = ripples.vHPC.coordinated.aversive(:,2);
            y = Spikes.vHPC.aversive{1,ii};
            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,tttt] = CCG(s,ids,'binSize',binSize,'duration',1,'smooth',1,'mode','ccg'); %ccg calculation
            perievent.vHPC.aversive.cooridnated.aversive = [perievent.vHPC.aversive.cooridnated.reward , ccg(:,1,2)./binSize./length(x)];
            plot(tttt,ccg(:,1,2)./binSize./length(x),'r'),hold on
            clear s ids groups ccg tttt
        end
        
        perievent.dHPC.reward.cooridnated.baseline = [];
        perievent.dHPC.reward.cooridnated.reward = [];
        perievent.dHPC.reward.cooridnated.aversive = [];
        for ii = 1 : size(assemblies.dHPC.reward,2)
            x = ripples.dHPC.coordinated.baseline(:,2);
            y = Spikes.dHPC.reward{1,ii};
            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,tttt] = CCG(s,ids,'binSize',binSize,'duration',1,'smooth',1,'mode','ccg'); %ccg calculation
            perievent.dHPC.reward.cooridnated.baseline = [perievent.dHPC.reward.cooridnated.baseline , ccg(:,1,2)./binSize./length(x)];
            figure, plot(tttt,ccg(:,1,2)./binSize./length(x),'k'),hold on
            clear s ids groups ccg tttt
            
            x = ripples.dHPC.coordinated.reward(:,2);
            y = Spikes.dHPC.reward{1,ii};
            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,tttt] = CCG(s,ids,'binSize',binSize,'duration',1,'smooth',1,'mode','ccg'); %ccg calculation
            perievent.dHPC.reward.cooridnated.reward = [perievent.dHPC.reward.cooridnated.reward , ccg(:,1,2)./binSize./length(x)];
            plot(tttt,ccg(:,1,2)./binSize./length(x),'b'),hold on
            clear s ids groups ccg tttt
            
            x = ripples.dHPC.coordinated.aversive(:,2);
            y = Spikes.dHPC.reward{1,ii};
            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,tttt] = CCG(s,ids,'binSize',binSize,'duration',1,'smooth',1,'mode','ccg'); %ccg calculation
            perievent.dHPC.reward.cooridnated.aversive = [perievent.dHPC.reward.cooridnated.reward , ccg(:,1,2)./binSize./length(x)];
            plot(tttt,ccg(:,1,2)./binSize./length(x),'r'),hold on
            clear s ids groups ccg tttt
        end

        perievent.vHPC.reward.cooridnated.baseline = [];
        perievent.vHPC.reward.cooridnated.reward = [];
        perievent.vHPC.reward.cooridnated.aversive = [];
        for ii = 1 : size(assemblies.vHPC.reward,2)
            x = ripples.vHPC.coordinated.baseline(:,2);
            y = Spikes.vHPC.reward{1,ii};
            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,tttt] = CCG(s,ids,'binSize',binSize,'duration',1,'smooth',1,'mode','ccg'); %ccg calculation
            perievent.vHPC.reward.cooridnated.baseline = [perievent.vHPC.reward.cooridnated.baseline , ccg(:,1,2)./binSize./length(x)];
            figure, plot(tttt,ccg(:,1,2)./binSize./length(x),'k'),hold on
            clear s ids groups ccg tttt
            
            x = ripples.vHPC.coordinated.reward(:,2);
            y = Spikes.vHPC.reward{1,ii};
            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,tttt] = CCG(s,ids,'binSize',binSize,'duration',1,'smooth',1,'mode','ccg'); %ccg calculation
            perievent.vHPC.reward.cooridnated.reward = [perievent.vHPC.reward.cooridnated.reward , ccg(:,1,2)./binSize./length(x)];
            plot(tttt,ccg(:,1,2)./binSize./length(x),'b'),hold on
            clear s ids groups ccg tttt
            
            x = ripples.vHPC.coordinated.aversive(:,2);
            y = Spikes.vHPC.reward{1,ii};
            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,tttt] = CCG(s,ids,'binSize',binSize,'duration',1,'smooth',1,'mode','ccg'); %ccg calculation
            perievent.vHPC.reward.cooridnated.aversive = [perievent.vHPC.reward.cooridnated.reward , ccg(:,1,2)./binSize./length(x)];
            plot(tttt,ccg(:,1,2)./binSize./length(x),'r'),hold on
            clear s ids groups ccg tttt
        end
                
        t
    end
    tt
end