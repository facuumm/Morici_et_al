clear
clc
close all

%% Parameters
path = {'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat103\usable';'E:\Rat132\recordings\in_pyr'};%List of folders from the path

% for SU
criteria_fr = 0.01; %criteria to include or not a SU into the analysis
criteria_n = 6; % minimal number of neurons from each structure
criteria_type = 0; %criteria for celltype (0:pyr, 1:int, 2:all)
n_SU_V = [];
n_SU_D = [];
Xedges = 60; %number of bins for RateMap construction

binSize = 0.001; % bin size for replay events detection

% Behavior
minimal_speed = 5; % minimal speed to detect active periods
minimal_speed_time = 2; % minimal time to detect active periods

Replay.candidate.dHPC.baseline=[]; Replay.candidate.dHPC.reward=[]; Replay.candidate.dHPC.aversive=[];
Replay.candidate.vHPC.baseline=[]; Replay.candidate.vHPC.reward=[]; Replay.candidate.vHPC.aversive=[];

Replay.selected.dHPC.baseline=[]; Replay.selected.dHPC.reward=[]; Replay.selected.dHPC.aversive=[];
Replay.selected.vHPC.baseline=[]; Replay.selected.vHPC.reward=[]; Replay.selected.vHPC.aversive=[];

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
        session = [subFolders(t+2).folder,'\',subFolders(t+2).name];
        cd(session)

        %Loading TS of the sessions
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
        %Shocks selection
        Shocks_filt = Restrict(shock,aversiveTS_run ./1000);
        % Keep only the first shock of each TTL (first from 20)
        count = 1;
        deff = [];
        for i = 1:length(Shocks_filt)
            if count == 1
                deff = [deff ; Shocks_filt(i,1)];
            end
            if count ==20
                count = 0;
            end
            count = count + 1;
        end
        Shocks_filt = deff;
        clear count deff shock i
        %Rewards selection
        Rewards_filt = Restrict([leftvalve ; rightvalve],rewardTS_run ./1000);
        disp('Uploading DLC outputs')
        camara = ((camara(:,2)-camara(:,1))/2)+camara(:,1);
        % periods of movment during eacj condition
        if rewardTS_run(1) < aversiveTS_run(1)
            load('laps1.mat','posx','posy');
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
        movement.reward(movement.reward(:,2) - movement.reward(:,1) <1,:)=[];
        clear tmp start stop
        % Aversive
        start = behavior.speed.aversive(1,1);   stop = behavior.speed.aversive(end,1);
        movement.aversive = InvertIntervals(behavior.quiet.aversive , start , stop);%keep only those higher than criteria
        movement.aversive(movement.aversive(:,2) - movement.aversive(:,1) <1,:)=[];
        clear tmp start stop
        
        % Definition of the place where the shocks were admninistrated
        position_shocks = [];
        for i = 1 : length(Shocks_filt)
            [~ , ii] = min(abs(behavior.pos.aversive(:,1)-Shocks_filt(i)));
            position_shocks = [position_shocks ; behavior.pos.aversive(ii,2)]; 
            clear ii
        end
        clear i
        position_shocks = position_shocks-min(behavior.pos.aversive(:,2));
        position_shocks = position_shocks./max(behavior.pos.aversive(:,2));

        %% load sleep states
        x = dir([cd,'\*-states.mat']);    states = load([cd,'\',x.name]);    states = states.states;
        
        REM = ToIntervals(states==5);    NREM = ToIntervals(states==3);    WAKE = ToIntervals(states==1);
        clear x states
        %% load coordinated ripple bursts
        load('coordinated_ripple_bursts.mat')
        ripple_bursts.baseline = Restrict(coordinated_ripple_bursts,baselineTS./1000);
        ripple_bursts.reward = Restrict(coordinated_ripple_bursts,rewardTS./1000);
        ripple_bursts.aversive = Restrict(coordinated_ripple_bursts,aversiveTS./1000);
        clear coordinated_ripple_bursts
        
        %% Load ripples
        ripplesD = table2array(readtable('ripplesD_customized2.csv'));
        ripplesV = table2array(readtable('ripplesV_customized2.csv'));
        % coordination
        coordinated = [];
        coordinatedV = [];
        coordinatedV_refined = [];
        for i = 1:length(ripplesD)
            r = ripplesD(i,:);
            tmp = sum(and(ripplesV(:,2)>= r(1,2)-0.1, ripplesV(:,2)<= r(1,2)+0.1));
            if tmp>0
                z = ripplesV(and(ripplesV(:,2)>= r(1,2)-0.1, ripplesV(:,2)<= r(1,2)+0.1),:);
                coordinatedV = [coordinatedV ; z];
                [p,indice] = min(abs(r(2)-z(:,2)));
                coordinatedV_refined = [coordinatedV_refined ; z(indice,:)];
                coordinated = [coordinated ; r];
                clear tmp2 tmp1 p indice z
            end
            clear r
        end
        clear x tmp i
        
        %% Spikes
        %Load Units
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
        
        %% Firing maps calculation
        maps_dHPC_A = [];        maps_dHPC_R = [];
        peak_dHPC_A = [] ;       peak_dHPC_R = [];
        field_dHPC_A = [];       field_dHPC_R = [];
        specificity_dHPC_A = []; specificity_dHPC_R = [];
        PC.dHPC.aversive=[]; PC.dHPC.reward=[];
        cluster_dHPC = [];
        for ii=1:size(group_dHPC,1)
            cluster = group_dHPC(ii,1);
            celltype = logical(cellulartype(cellulartype(:,1) == cluster,2));
            if celltype
                spks = spks_dHPC(spks_dHPC(:,1)==cluster,2);
                b = [movement.aversive ; movement.reward];
                tmp = Restrict(spks,b);
                m = length(tmp) / sum(b(:,2) - b(:,1));
%                 if and(m > criteria_fr , ~isnan(m))
                if ~isnan(m)
                    % --- Aversive ---
                    spks_tmp = Restrict(spks , movement.aversive);
                    pos_tmp = Restrict(behavior.pos.aversive(:,1:2) , movement.aversive);
%                     pos_tmp(or(pos_tmp(:,2)<25 , pos_tmp(:,2)>175),2) = nan; %Eliminate platform periods
                    pos_tmp(:,2) = pos_tmp(:,2)-min(pos_tmp(:,2)); pos_tmp(:,2) = pos_tmp(:,2)/max(pos_tmp(:,2)); %normalization of position
                    
                    %Firing curve construction
                    [curve , stats] = FiringCurve(pos_tmp , spks_tmp , 'smooth' , 2 , 'nBins' , Xedges , 'minSize' , 2 , 'minPeak' , 0.2);
                    maps_dHPC_A = [maps_dHPC_A ; curve.rate];
                    peak_dHPC_A = [peak_dHPC_A ; stats.peak(1)] ;
                    field_dHPC_A = [field_dHPC_A ; stats.fieldX(1,:)];
                    
                    %Quantile for place-cell definition
                    [~,~,curve1,OccMap,~,Nspikes,~,~,x] = FiringMap_LinearTrack(pos_tmp, spks_tmp, 1/30, Xedges, 0);
%                     maps_dHPC_A = [maps_dHPC_A ; curve];
%                     s = spatial_info_rate(curve, OccMap);
                    q = SkaggsRandom(spks_tmp, pos_tmp, Nspikes, OccMap, Xedges);  
                    sp = sparsity_info(curve1,OccMap);
                    specificity_dHPC_A = [specificity_dHPC_A ; stats.specificity , q , sp , cond];
                    
                    % Store of place-field location and cluster of PCs
                    if stats.specificity > 0.25
                        [ff,f] = max(curve.rate);
                        fff = [1/Xedges:1/Xedges:1];
                        PC.dHPC.aversive=[PC.dHPC.aversive ; cluster , ff , fff(f)];
                        clear f ff fff
                    end
                    clear s q curve OccMap Nspikes spks_tmp pos_tmp m curve stats sp curve1

                    % --- Reward ---
                    spks_tmp = Restrict(spks , movement.reward);
                    pos_tmp = Restrict(behavior.pos.reward(:,1:2) , movement.reward);
%                     pos_tmp(or(pos_tmp(:,2)<30 , pos_tmp(:,2)>150),2) = nan; %Eliminate platform periods
                    pos_tmp(:,2) = pos_tmp(:,2)-min(pos_tmp(:,2)); pos_tmp(:,2) = pos_tmp(:,2)/max(pos_tmp(:,2)); %normalization of position
                    
                    [curve , stats] = FiringCurve(pos_tmp , spks_tmp , 'smooth' , 2 , 'nBins' , Xedges , 'minSize' , 2 , 'minPeak' , 0.2);
                    maps_dHPC_R = [maps_dHPC_R ; curve.rate];
                    peak_dHPC_R = [peak_dHPC_R ; stats.peak(1)] ;
                    field_dHPC_R = [field_dHPC_R ; stats.fieldX(1,:)];

                    [~,~,curve1,OccMap,~,Nspikes,~,~,x] = FiringMap_LinearTrack(pos_tmp, spks_tmp, 1/30, Xedges, 0);
%                     maps_dHPC_R = [maps_dHPC_R ; curve];
%                     s = spatial_info_rate(curve, OccMap);
                    q = SkaggsRandom(spks_tmp, pos_tmp, Nspikes, OccMap, Xedges);   
                    sp = sparsity_info(curve1,OccMap);
                    specificity_dHPC_R = [specificity_dHPC_R ; stats.specificity , q , sp , cond];
                    
                    
                    % Store of place-field location and cluster of PCs
                    if stats.specificity > 0.25
                        [ff,f] = max(curve.rate);
                        fff = [1/Xedges:1/Xedges:1];
                        PC.dHPC.reward=[PC.dHPC.reward ; cluster , ff , fff(f)];
                        clear f ff fff
                    end
                    clear s q curve OccMap Nspikes spks_tmp pos_tmp m curve stats sp curve1
                    
                    cluster_dHPC = [cluster_dHPC ; cluster];
                    
                    condition1 = or(isnan(specificity_dHPC_R(end,1)) , isnan(specificity_dHPC_A(end,1))); %store only those that contain spatial information
                    condition2 = or(mean(maps_dHPC_R(end,:))<criteria_fr , mean(maps_dHPC_A(end,:))<criteria_fr); %store only those that fire at least 10 times
                    if or(condition1 , condition2)
                            specificity_dHPC_R(end,:) = [];
                            peak_dHPC_R(end,:)        = [];
                            field_dHPC_R(end,:)       = [];
                            maps_dHPC_R(end,:)        = [];
                            
                            specificity_dHPC_A(end,:) = [];
                            peak_dHPC_A(end,:)        = [];
                            field_dHPC_A(end,:)       = [];
                            maps_dHPC_A(end,:)        = [];
                            
                            cluster_dHPC(end,:) = [];
                    end
                    clear condition1 condition2
                end
            end
            clear celltype tmp b cluster
        end
        
        
        maps_vHPC_A = [];        maps_vHPC_R = [];
        peak_vHPC_A = [] ;       peak_vHPC_R = [];
        field_vHPC_A = [];       field_vHPC_R = [];
        specificity_vHPC_A = []; specificity_vHPC_R = [];
        PC.vHPC.aversive=[];     PC.vHPC.reward=[];
        cluster_vHPC = [];
        for ii=1:size(group_vHPC,1)
            cluster = group_vHPC(ii,1);
            celltype = logical(cellulartype(cellulartype(:,1) == cluster,2));
            if celltype
                spks = spks_vHPC(spks_vHPC(:,1)==cluster,2);
                b = [movement.aversive ; movement.reward];
                tmp = Restrict(spks,b);
                m = length(tmp) / sum(b(:,2) - b(:,1));
                %                 if and(m > criteria_fr , ~isnan(m))
                if ~isnan(m)
                    % --- Aversive ---
                    spks_tmp = Restrict(spks , movement.aversive);
                    pos_tmp = Restrict(behavior.pos.aversive(:,1:2) , movement.aversive);
                    %                     pos_tmp(or(pos_tmp(:,2)<30 , pos_tmp(:,2)>150),2) = nan; %Eliminate platform periods
                    pos_tmp(:,2) = pos_tmp(:,2)-min(pos_tmp(:,2)); pos_tmp(:,2) = pos_tmp(:,2)/max(pos_tmp(:,2)); %normalization of position
                    
                    %Firing curve construction
                    [curve , stats] = FiringCurve(pos_tmp , spks_tmp , 'smooth' , 2 , 'nBins' , Xedges , 'minSize' , 4 , 'minPeak' , 0.2);
                    maps_vHPC_A = [maps_vHPC_A ; curve.rate];
                    peak_vHPC_A = [peak_vHPC_A ; stats.peak(1)] ;
                    field_vHPC_A = [field_vHPC_A ; stats.fieldX(1,:)];
                    
                    [~,~,curve1,OccMap,~,Nspikes,~,~,x] = FiringMap_LinearTrack(pos_tmp, spks_tmp, 1/30, Xedges, 0);
                    %                     maps_vHPC_A = [maps_vHPC_A ; curve];
                    %                     s = spatial_info_rate(curve, OccMap);
                    q = SkaggsRandom(spks_tmp, pos_tmp, Nspikes, OccMap, Xedges);
                    sp = sparsity_info(curve1,OccMap);
                    specificity_vHPC_A = [specificity_vHPC_A ; stats.specificity , q , sp , cond];
                    
                    % Store of place-field location and cluster of PCs
                    if stats.specificity > 0.25
                        [ff,f] = max(curve.rate);
                        fff = [1/Xedges:1/Xedges:1];
                        PC.vHPC.aversive=[PC.vHPC.aversive ; cluster , ff , fff(f)];
                        clear f ff fff
                    end
                    clear s q curve OccMap Nspikes spks_tmp pos_tmp m curve stats sp curve1
                    
                    % --- Reward ---
                    spks_tmp = Restrict(spks , movement.reward);
                    pos_tmp = Restrict(behavior.pos.reward(:,1:2) , movement.reward);
                    %                     pos_tmp(or(pos_tmp(:,2)<30 , pos_tmp(:,2)>150),2) = nan; %Eliminate platform periods
                    pos_tmp(:,2) = pos_tmp(:,2)-min(pos_tmp(:,2)); pos_tmp(:,2) = pos_tmp(:,2)/max(pos_tmp(:,2)); %normalization of position
                    
                    %Firing curve construction
                    [curve , stats] = FiringCurve(pos_tmp , spks_tmp , 'smooth' , 2 , 'nBins' , Xedges , 'minSize' , 4 , 'minPeak' , 0.2);
                    maps_vHPC_R = [maps_vHPC_R ; curve.rate];
                    peak_vHPC_R = [peak_vHPC_R ; stats.peak(1)] ;
                    field_vHPC_R = [field_vHPC_R ; stats.fieldX(1,:)];
                    
                    [~,~,curve1,OccMap,~,Nspikes,~,~,x] = FiringMap_LinearTrack(pos_tmp, spks_tmp, 1/30, Xedges, 0);
                    %                     maps_vHPC_R = [maps_vHPC_R ; curve];
                    %                     s = spatial_info_rate(curve, OccMap);
                    q = SkaggsRandom(spks_tmp, pos_tmp, Nspikes, OccMap, Xedges);
                    sp = sparsity_info(curve1,OccMap);
                    specificity_vHPC_R = [specificity_vHPC_R ; stats.specificity , q , sp , cond];
                    
                    % Store of place-field location and cluster of PCs
                    if stats.specificity > 0.25
                        [ff,f] = max(curve.rate);
                        fff = [1/Xedges:1/Xedges:1];
                        PC.vHPC.reward=[PC.vHPC.reward ; cluster , ff , fff(f)];
                        clear f ff fff
                    end
                    clear s q curve OccMap Nspikes spks_tmp pos_tmp m curve stats sp curve1
                    
                    cluster_vHPC = [cluster_vHPC ; cluster];
                    
                    condition1 = or(isnan(specificity_vHPC_R(end,1)) , isnan(specificity_vHPC_A(end,1)));%store only those that contain spatial information
                    condition2 = or(mean(maps_vHPC_R(end,:))<criteria_fr , mean(maps_vHPC_A(end,:))<criteria_fr); %store only those that fire at least 10 times
                    if or(condition1 , condition2)
                        specificity_vHPC_R(end,:) = [];
                        peak_vHPC_R(end,:)        = [];
                        field_vHPC_R(end,:)       = [];
                        maps_vHPC_R(end,:)        = [];
                        
                        specificity_vHPC_A(end,:) = [];
                        peak_vHPC_A(end,:)        = [];
                        field_vHPC_A(end,:)       = [];
                        maps_vHPC_A(end,:)        = [];
                        
                        cluster_vHPC(end,:) = [];
                        
                    end
                    clear condition1 condition2
                end
            end
            clear celltype tmp b cluster
        end
        
        %         subplot(1,2,1), imagesc([1 : size(maps_vHPC_A,2)],[1 : size(maps_vHPC_A,1)]',maps_vHPC_A),xlabel('Position'),title('Aversive'); axis xy;
        %         subplot(1,2,2), imagesc([1 : size(maps_vHPC_R,2)],[1 : size(maps_vHPC_R,1)]',maps_vHPC_R),xlabel('Position'),title('Reward'); axis xy;
        
        %% Firing maps calculation
        PHIST.dHPC.aversive = [];        PHIST.dHPC.reward = [];
        clusterdHPC.aversive = [];      clusterdHPC.reward = [];
        for ii=1:size(group_dHPC,1)
            cluster = group_dHPC(ii,1);
            celltype = logical(cellulartype(cellulartype(:,1) == cluster,2));
            if celltype
                spks = spks_dHPC(spks_dHPC(:,1)==cluster,2);
                b = [Shocks_filt(:,1)-2 Shocks_filt(:,1)+2];
                cond1 = Restrict(spks,b);
                %                 base.aversive = SubtractIntervals(aversiveTS_run./1000 , b);
                
                b = [Rewards_filt(:,1)-2 Rewards_filt(:,1)+2];
                cond2 = Restrict(spks,b);
                %                 base.reward = SubtractIntervals(rewardTS_run./1000 , b);
                
%                 %Poisson
%                 base = InvertIntervals([ripplesD(:,1)-0.1, ripplesD(:,3)+0.1],NREM(:,1) , NREM(:,2));
%                 Base = Restrict(spks,base);
% 
%                 totalrippletime = sum(ripplesD(:,3)-ripplesD(:,1));
%                 ripplespikes = Restrict(spks,[ripplesD(:,1) ripplesD(:,3)]);
%                 nripplespikes = size(ripplespikes,1);
%                 
%                 ncellbaselinespikes = length(Base);
%                 ncellripplespikes = length(ripplespikes);
%                 totalbaselinetime = sum(base(:,2)-base(:,1));
%                 if ncellbaselinespikes~=0 & ncellripplespikes~=0
%                     [pInc pDec surp] = poissonTest(ncellbaselinespikes/totalbaselinetime,ncellripplespikes,totalrippletime);
%                 else
%                     pInc = NaN;
%                     pDec = NaN;
%                     surp = NaN;
%                 end
%                 
                
                if length(cond1)>5  % --- Aversive ---
                    x = Restrict(spks , [Shocks_filt(:,1)-2 Shocks_filt(:,1)+2]);
                    y = [Shocks_filt(:,1)];
                    [times,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
                    [ccg,TimeVector1] = CCG(times,groups,'binsize',0.05,'duration',4,'smooth',2);
                    ccg = ccg(:,1,2)./length(y)./0.05;
                    PHIST.dHPC.aversive = [PHIST.dHPC.aversive ; ccg'];
                    clear ccg x y
                    clusterdHPC.aversive = [clusterdHPC.aversive ; cluster];
                end
                if length(cond2)>5 % --- Reward ---
                    x = Restrict(spks , [Rewards_filt(:,1)-2 Rewards_filt(:,1)+2]);
                    y = [Rewards_filt(:,1)];
                    [times,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
                    [ccg,TimeVector1] = CCG(times,groups,'binsize',0.05,'duration',4,'smooth',2);
                    ccg = ccg(:,1,2)./length(y)./0.05;
                    PHIST.dHPC.reward = [PHIST.dHPC.reward ; ccg'];
                    clear ccg x y
                    
                    
                    clusterdHPC.reward = [clusterdHPC.reward ; cluster];
                end
            end
            
            clear celltype tmp b cluster pInc pDec surp
            clear base Base totalrippletime ripplespikes ncellbaselinespikes ncellripplespikes totalbaselinetime

        end
        clear base
        
        PHIST.vHPC.aversive = [];        PHIST.vHPC.reward = [];
        clustervHPC.aversive = [];      clustervHPC.reward = [];
        for ii=1:size(group_vHPC,1)
            cluster = group_vHPC(ii,1);
            celltype = logical(cellulartype(cellulartype(:,1) == cluster,2));
            if celltype
                spks = spks_vHPC(spks_vHPC(:,1)==cluster,2);
                b = [Shocks_filt(:,1)-2 Shocks_filt(:,1)+2];
                cond1 = Restrict(spks,b);
                
                b = [Rewards_filt(:,1)-2 Rewards_filt(:,1)+2];
                cond2 = Restrict(spks,b);
                
%                 %Poisson
%                 base = InvertIntervals([ripplesV(:,1)-0.1, ripplesV(:,3)+0.1],NREM(:,1) , NREM(:,2));
%                 Base = Restrict(spks,base);
%                 
%                 totalrippletime = sum(ripplesD(:,3)-ripplesD(:,1));
%                 ripplespikes = Restrict(spks,[ripplesV(:,1) ripplesV(:,3)]);
%                 nripplespikes = size(ripplespikes,1);
%                 
%                 ncellbaselinespikes = length(Base);
%                 ncellripplespikes = length(ripplespikes);
%                 totalbaselinetime = sum(base(:,2)-base(:,1));
%                 if ncellbaselinespikes~=0 & ncellripplespikes~=0
%                     [pInc pDec surp] = poissonTest(ncellbaselinespikes/totalbaselinetime,ncellripplespikes,totalrippletime);
%                 else
%                     pInc = NaN;
%                     pDec = NaN;
%                     surp = NaN;
%                 end                
                 
                if length(cond1)>5 % --- Aversive ---
                    
                    x = Restrict(spks , [Shocks_filt(:,1)-2 Shocks_filt(:,1)+2]);
                    y = [Shocks_filt(:,1)];
                    [times,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
                    [ccg,TimeVector1] = CCG(times,groups,'binsize',0.05,'duration',4,'smooth',2);
                    ccg = ccg(:,1,2)./length(y)./0.05;
                    PHIST.vHPC.aversive = [PHIST.vHPC.aversive ; ccg'];
                    clear ccg x y
                    
                    clustervHPC.aversive = [clustervHPC.aversive ; cluster];
                    
                end
                
                if length(cond2)>5 % --- Reward ---
                    x = Restrict(spks , [Rewards_filt(:,1)-2 Rewards_filt(:,1)+2]);
                    y = [Rewards_filt(:,1)];
                    [times,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
                    [ccg,TimeVector1] = CCG(times,groups,'binsize',0.05,'duration',4,'smooth',2);
                    ccg = ccg(:,1,2)./length(y)./0.05;
                    PHIST.vHPC.reward = [PHIST.vHPC.reward ; ccg'];
                    clear ccg x y
                    
                    clustervHPC.reward = [clustervHPC.reward ; cluster];
                    
                end
            end
            clear celltype tmp b cluster base Base pInc pDec surp
            clear base Base totalrippletime ripplespikes ncellbaselinespikes ncellripplespikes totalbaselinetime
        end
        
    
%         %% Firing maps related to each event
%         % Detection of space close to the event
%         pos_normalized = [behavior.pos.aversive(:,1),behavior.pos.aversive(:,2)];
%         pos_normalized(:,2) =  pos_normalized(:,2) - min( pos_normalized(:,2));
%         pos_normalized(:,2) =  pos_normalized(:,2) ./ max( pos_normalized(:,2));
%         laps.aversive = [];
%         for i = 1:length(Shocks_filt)
%             origin = Shocks_filt(i);
%             [~,ind] = min(abs(pos_normalized(:,1) - origin));
%             lap =  Restrict(pos_normalized , [origin-20 origin + 20]);
%             lap(:,2) = abs(lap(:,2) - pos_normalized(ind,2));%.*-1;
%             speed = Restrict(behavior.speed.aversive , [origin - 20 origin + 20]);
%                         
%             [~,ind] = max(lap(:,2));
%             
%             lap = lap(1:ind,:); speed = speed(1:ind,:);
%             
%             x = QuietPeriods(speed,minimal_speed,minimal_speed_time);
%             x = SubtractIntervals([origin-20 origin + 20],x);
%             x = x(x(:,2)-x(:,1)>1,:);
%             x = Restrict(lap,x);
%             laps.aversive = [laps.aversive ; x];
%             
% %             plot(x(:,1) , x(:,2),'*')
%             clear lap ind x
%         end
%         clear pos_normalized i 
%         
%         pos_normalized = [behavior.pos.reward(:,1),behavior.pos.reward(:,2)];
%         pos_normalized(:,2) =  pos_normalized(:,2) - min( pos_normalized(:,2));
%         pos_normalized(:,2) =  pos_normalized(:,2) ./ max( pos_normalized(:,2));
%         laps.reward = [];
%         Rewards_filt = sort(Rewards_filt);
%         for i = 2:length(Rewards_filt)
%             start = Rewards_filt(i-1,2);
%             stop = Rewards_filt(i,1);
%             if stop-start>1
%                 if and(start>pos_normalized(1,1) , stop<pos_normalized(end,1))
%                     [~,ind] = min(abs(pos_normalized(:,1) - start));
%                     speed = Restrict(behavior.speed.reward,[start stop]);
%                     speed = QuietPeriods(speed,minimal_speed,minimal_speed_time);
%                     [x , ~] = SubtractIntervals([start,stop],speed);
%                     lap =  Restrict(pos_normalized , x);
%                     lap(:,2) = abs((lap(:,2) - pos_normalized(ind,2)));
%                     
%                     %                 figure,plot(lap(:,1),lap(:,2)),hold on
%                     
%                     
%                     laps.reward = [laps.reward ; lap];
%                     
%                     %                 plot(lap(ind:ind2,1) , lap(ind:ind2,2),'*')
%                     clear lap ind x
%                 end
%             end
%         end
%         clear pos_normalized i
        
%         % Calculation of maps normalized to the event
%         PHIST.dHPC.normalizado.aversive = [];        PHIST.dHPC.normalizado.reward = [];
%         clusterdHPC.normalizado.aversive = [];      clusterdHPC.normalizado.reward = [];
%         for ii=1:length(group_dHPC)
%             cluster = group_dHPC(ii,1);
%             celltype = logical(cellulartype(cellulartype(:,1) == cluster,2));
%             if celltype
%                 spks = spks_dHPC(spks_dHPC(:,1)==cluster,2);
%                 b = [Shocks_filt(:,1)-2 Shocks_filt(:,1)+2];
%                 cond1 = Restrict(spks,b);
%                 
%                 b = [Rewards_filt(:,1)-2 Rewards_filt(:,1)+2];
%                 cond2 = Restrict(spks,b);
%                 
%                 if length(cond1)>5
%                     % --- Aversive ---
%                     spks_tmp = Restrict(spks ,[min(laps.aversive(:,1)) max(laps.aversive(:,1))]);
%                     [curve , stats] = FiringCurve(laps.aversive , spks_tmp , 'smooth' , 2 , 'nBins' , Xedges , 'minSize' , 2 , 'minPeak' , 0.2);
%                     
%                     PHIST.dHPC.normalizado.aversive = [PHIST.dHPC.normalizado.aversive ; curve.rate];
%                     clear curve stats
%                     clusterdHPC.normalizado.aversive = [clusterdHPC.normalizado.aversive ; cluster];
%                 end
%                 if length(cond2)>5
%                     % --- Reward ---
%                     spks_tmp = Restrict(spks ,[min(laps.reward(:,1)) max(laps.reward(:,1))]);
%                     [curve , stats] = FiringCurve(laps.reward , spks_tmp , 'smooth' , 2 , 'nBins' , Xedges , 'minSize' , 2 , 'minPeak' , 0.2);
%                     
%                     PHIST.dHPC.normalizado.reward = [PHIST.dHPC.normalizado.reward ; curve.rate];
%                     clear curve stats
%                     clusterdHPC.normalizado.reward = [clusterdHPC.normalizado.reward ; cluster];
%                 end
%             end
%             
%             clear celltype tmp b cluster
%         end
%         
%         PHIST.vHPC.normalizado.aversive = [];        PHIST.vHPC.normalizado.reward = [];
%         clustervHPC.normalizado.aversive = [];      clustervHPC.normalizado.reward = [];
%         for ii=1:length(group_vHPC)
%             cluster = group_vHPC(ii,1);
%             celltype = logical(cellulartype(cellulartype(:,1) == cluster,2));
%             if celltype
%                 spks = spks_vHPC(spks_vHPC(:,1)==cluster,2);
%                 b = [Shocks_filt(:,1)-2 Shocks_filt(:,1)+2];
%                 cond1 = Restrict(spks,b);
%                 
%                 b = [Rewards_filt(:,1)-2 Rewards_filt(:,1)+2];
%                 cond2 = Restrict(spks,b);
%                 
%                 if length(cond1)>5
%                     % --- Aversive ---
%                     spks_tmp = Restrict(spks ,[min(laps.aversive(:,1)) max(laps.aversive(:,1))]);
%                     [curve , stats] = FiringCurve(laps.aversive , spks_tmp , 'smooth' , 2 , 'nBins' , Xedges , 'minSize' , 2 , 'minPeak' , 0.2);
%                     
%                     PHIST.vHPC.normalizado.aversive = [PHIST.vHPC.normalizado.aversive ; curve.rate];
%                     clear curve stats
%                     clustervHPC.normalizado.aversive = [clustervHPC.normalizado.aversive ; cluster];
%                 end
%                 if length(cond2)>5
%                     % --- Reward ---
%                     spks_tmp = Restrict(spks ,[min(laps.reward(:,1)) max(laps.reward(:,1))]);
%                     [curve , stats] = FiringCurve(laps.reward , spks_tmp , 'smooth' , 2 , 'nBins' , Xedges , 'minSize' , 2 , 'minPeak' , 0.2);
%                     
%                     PHIST.vHPC.normalizado.reward = [PHIST.vHPC.normalizado.reward ; curve.rate];
%                     clear curve stats
%                     clustervHPC.normalizado.reward = [clustervHPC.normalizado.reward ; cluster];
%                 end
%             end
%             
%             clear celltype tmp b cluster
%         end


        %% Replay detection
        freq = 1/binSize;
        limits = [0 segments.Var1(end)/1000];
        spks = [];
        %Replay events in dHPC
        for ii=1:size(group_dHPC,1)
            cluster = group_dHPC(ii,1);
            celltype = logical(cellulartype(cellulartype(:,1) == cluster,2));
            if celltype
                spks = [spks;spks_dHPC(spks_dHPC(:,1)==cluster,2)];
            end
            clear celltype
        end
        [MUA.dHPC,bins]=binspikes(sort(spks,'ascend'),freq,limits);clear spks
        w = gausswin(30, 1);
        MUA.dHPC = filter(w , 1 ,MUA.dHPC);
        replay.dHPC = MUA.dHPC>mean(MUA.dHPC)+std(MUA.dHPC)*3;
        replay.dHPC = ToIntervals(bins',replay.dHPC);
        replay.dHPC = replay.dHPC(replay.dHPC(:,2)-replay.dHPC(:,1)>0.03,:);
        
        [replay.dHPC] = merge_events(replay.dHPC, 0.04);
        
        %filtering replay events by amount of PCs
        count = [];
        for i = 1 :  size(cluster_dHPC,1)
            cluster = cluster_dHPC(i,1);
            spks = spks_dHPC(spks_dHPC(:,1)==cluster,2);
            [status,interval,index] = InIntervals(spks,replay.dHPC);
            interval = unique(interval);
            count = [count ; interval(interval~=0)];
            clear spks cluster status interval index
        end
        [gc,grps] = groupcounts(count);
        replay.dHPC = replay.dHPC(grps(gc>length(cluster_dHPC)*0.15),:);
        
        
        %Replay events in vHPC
        spks = [];
        for ii=1:size(group_vHPC,1)
            cluster = group_vHPC(ii,1);
            celltype = logical(cellulartype(cellulartype(:,1) == cluster,2));
            if celltype
                spks = [spks;spks_vHPC(spks_vHPC(:,1)==cluster,2)];
            end
            clear celltype
        end
        [MUA.vHPC,bins]=binspikes(sort(spks,'ascend'),freq,limits);clear spks
        w = gausswin(30, 1);
        MUA.vHPC = filter(w , 1 ,MUA.vHPC);
        replay.vHPC = MUA.vHPC>mean(MUA.vHPC)+std(MUA.vHPC)*2;
        replay.vHPC = ToIntervals(bins',replay.vHPC);
        replay.vHPC = replay.vHPC(replay.vHPC(:,2)-replay.vHPC(:,1)>0.03,:);
        
        [replay.vHPC] = merge_events(replay.vHPC, 0.04);
        
        %filtering replay events by amount of PCs
        count = [];
        for i = 1 :  size(cluster_vHPC,1)
            cluster = cluster_vHPC(i,1);
            spks = spks_vHPC(spks_vHPC(:,1)==cluster,2);
            [status,interval,index] = InIntervals(spks,replay.vHPC);
            interval = unique(interval);
            count = [count ; interval(interval~=0)];
            clear spks cluster status interval index
        end
        [gc,grps] = groupcounts(count);
        replay.vHPC = replay.vHPC(grps(gc>length(cluster_vHPC)*0.15),:);
        
        %% Save Replay percentage
        % --- Baseline ---
        x = length(Restrict(Restrict(replay.dHPC , [ripplesD(:,2)-0.1 ripplesD(:,2)+0.1]),baselineTS./1000));
        y = baselineTS(2)/1000 - baselineTS(1)/1000;
        Replay.selected.dHPC.baseline=[Replay.selected.dHPC.baseline ; x/y]; clear x y
        % --- Reward ---
        x = length(Restrict(Restrict(replay.dHPC , [ripplesD(:,2)-0.1 ripplesD(:,2)+0.1]),rewardTS./1000));
        y = rewardTS(2)/1000 - rewardTS(1)/1000;
        Replay.selected.dHPC.reward=[Replay.selected.dHPC.reward ; x/y]; clear x y
        % --- Aversive ---
        x = length(Restrict(Restrict(replay.dHPC , [ripplesD(:,2)-0.1 ripplesD(:,2)+0.1]),aversiveTS./1000));
        y = aversiveTS(2)/1000 - aversiveTS(1)/1000;
        Replay.selected.dHPC.aversive=[Replay.selected.dHPC.reward ; x/y]; clear x y
        
        % --- Baseline ---
        x = length(Restrict(Restrict(replay.vHPC , [ripplesV(:,2)-0.1 ripplesV(:,2)+0.1]),baselineTS./1000));
        y = baselineTS(2)/1000 - baselineTS(1)/1000;
        Replay.selected.vHPC.baseline=[Replay.selected.vHPC.baseline ; x/y]; clear x y
        % --- Reward ---
        x = length(Restrict(Restrict(replay.vHPC , [ripplesV(:,2)-0.1 ripplesV(:,2)+0.1]),rewardTS./1000));
        y = rewardTS(2)/1000 - rewardTS(1)/1000;
        Replay.selected.vHPC.reward=[Replay.selected.vHPC.reward ; x/y]; clear x y
        % --- Aversive ---
        x = length(Restrict(Restrict(replay.vHPC , [ripplesV(:,2)-0.1 ripplesV(:,2)+0.1]),aversiveTS./1000));
        y = aversiveTS(2)/1000 - aversiveTS(1)/1000;
        Replay.selected.vHPC.aversive=[Replay.selected.vHPC.reward ; x/y]; clear x y        
        
        %% Bayesian Decoding across conditions
        % Position
        decoded_positions.dHPC.aversive = [];
        decoded_positions.vHPC.aversive = [];
        decoded_positions.dHPC.reward = [];
        decoded_positions.vHPC.reward = [];
        
   
        % --- Aversive ---
        tmp.vHPC = Restrict(Restrict(replay.vHPC , [ripplesV(:,2)-0.1 ripplesV(:,2)+0.1]),aversiveTS./1000);
        tmp.dHPC = Restrict(Restrict(replay.dHPC , [ripplesV(:,2)-0.1 ripplesV(:,2)+0.1]),aversiveTS./1000);
        
        % for ventral hippocampus
        realReplay.vHPC.aversive.forward = [];
        realReplay.vHPC.aversive.backward = [];
        for i = 1 : length(tmp.vHPC)
            % Decoding using dHPC SUs
            start = tmp.vHPC(i,1);
            stop = tmp.vHPC(i,2);
            bin = ((stop-start)/2) + start;
            bin = [start : 0.02 : stop];
            tmp1 = [];
            for ii = 1 : length(bin)-1
            nSpks = count_spks(spks_vHPC, cluster_vHPC, bin(ii), bin(ii+1));
            probability = bayesian_replay(PHIST.vHPC.normalizado.aversive, nSpks, bin(ii), bin(ii+1));
            
            [o p] = max(probability);
            
%             tmp1 = [tmp1 ; o p ii];
            tmp1 = [tmp1 , o p ii];

            clear nSpks probability o p
            end
            c = corrcoef(tmp1(:,3) , tmp1(:,2));
%             plot(tmp1(:,3) , tmp1(:,2),'*')
            
            shuffle = [];
            for ii = 1:1000
               p = corrcoef(tmp1(randperm(length(tmp1)),3) , tmp1(:,2));
                shuffle = [shuffle ; p(1,2)];
            end
            
            if not(isnan(c(1,2))) %if correlation exist
                if sign(c(1,2))<0 %for reverse replay
                    p = sum(shuffle<c(1,2))/1000;
                    if p<0.05
                        realReplay.vHPC.aversive.backward = [realReplay.vHPC.aversive.backward ; start stop];
                        figure,plot(tmp1(:,3) , tmp1(:,2),'*')
                    end
                else %for replay
                    p = sum(shuffle>c(1,2))/1000;
                    
                    if p<0.05
                        realReplay.vHPC.aversive.forward = [realReplay.vHPC.aversive.forward ; start stop];
                        plot(tmp1(:,3) , tmp1(:,2),'*')
                    end
                    
                end
                
%                 if p<0.05 %if
%                     figure, histogram(shuffle,100)
%                     xline(c(1,2))
%                 end
            end
            clear probability nSpks start stop p c
        end

        % for dorsal hippocampus
        realReplay.dHPC.aversive.forward = [];
        realReplay.dHPC.aversive.backward = [];
        for i = 1 : length(tmp.dHPC)
            % Decoding using dHPC SUs
            start = tmp.dHPC(i,1);
            stop = tmp.dHPC(i,2);
            bin = ((stop-start)/2) + start;
            bin = [bin-0.1 : 0.02 : bin+0.1];
            tmp1 = [];
            for ii = 1 : length(bin)-1
            nSpks = count_spks(spks_dHPC, clusterdHPC.normalizado.aversive, bin(ii), bin(ii+1));
            probability = bayesian_replay(PHIST.dHPC.normalizado.aversive, nSpks, bin(ii), bin(ii+1));
            
            [o p] = max(probability);
            
            tmp1 = [tmp1 ; o p ii];
            clear nSpks probability o p
            end
            c = corrcoef(tmp1(:,3) , tmp1(:,2));
%             plot(tmp1(:,3) , tmp1(:,2),'*')
            
            shuffle = [];
            for ii = 1:1000
               p = corrcoef(tmp1(randperm(length(tmp1)),3) , tmp1(:,2));
                shuffle = [shuffle ; p(1,2)];
            end
            
            if not(isnan(c(1,2))) %if correlation exist
                if sign(c(1,2))<0 %for reverse replay
                    p = sum(shuffle<c(1,2))/1000;
                    if p<0.05
                        realReplay.dHPC.aversive.backward = [realReplay.dHPC.aversive.backward ; start stop];
                        figure,plot(tmp1(:,3) , tmp1(:,2),'*')
                    end
                else %for replay
                    p = sum(shuffle>c(1,2))/1000;
                    
                    if p<0.05
                        realReplay.dHPC.aversive.forward = [realReplay.dHPC.aversive.forward ; start stop];
                        plot(tmp1(:,3) , tmp1(:,2),'*')
                    end
                    
                end
                
%                 if p<0.05 %if
%                     figure, histogram(shuffle,100)
%                     xline(c(1,2))
%                 end
            end
            clear probability nSpks start stop p c
        end        
  
        
%         
%       % --- Reward ---
%         tmp.vHPC = Restrict(replay.vHPC , [ripple_bursts.reward(:,1) ripple_bursts.reward(:,3)]);
%         tmp.dHPC = Restrict(replay.dHPC , [ripple_bursts.reward(:,1) ripple_bursts.reward(:,3)]);
%         
%         % for ventral hippocampus
%         realReplay.vHPC.reward.forward = [];
%         realReplay.vHPC.reward.backward = [];
%         for i = 1 : length(tmp.vHPC)
%             % Decoding using dHPC SUs
%             start = tmp.vHPC(i,1);
%             stop = tmp.vHPC(i,2);
%             bin = ((stop-start)/2) + start;
%             bin = [bin-0.1 : 0.02 : bin+0.1];
%             tmp1 = [];
%             for ii = 1 : length(bin)-1
%             nSpks = count_spks(spks_vHPC, clustervHPC.normalizado.aversive, bin(ii), bin(ii+1));
%             probability = bayesian_replay(PHIST.vHPC.normalizado.aversive, nSpks, bin(ii), bin(ii+1));
%             
%             [o p] = max(probability);
%             
%             tmp1 = [tmp1 ; o p ii];
%             clear nSpks probability o p
%             end
%             c = corrcoef(tmp1(:,3) , tmp1(:,2));
% %             plot(tmp1(:,3) , tmp1(:,2),'*')
%             
%             shuffle = [];
%             for ii = 1:1000
%                p = corrcoef(tmp1(randperm(length(tmp1)),3) , tmp1(:,2));
%                 shuffle = [shuffle ; p(1,2)];
%             end
%             
%             if not(isnan(c(1,2))) %if correlation exist
%                 if sign(c(1,2))<0 %for reverse replay
%                     p = sum(shuffle<c(1,2))/1000;
%                     if p<0.05
%                         realReplay.vHPC.reward.backward = [realReplay.vHPC.reward.backward ; start stop];
%                         figure,plot(tmp1(:,3) , tmp1(:,2),'*')
%                     end
%                 else %for replay
%                     p = sum(shuffle>c(1,2))/1000;
%                     
%                     if p<0.05
%                         realReplay.vHPC.reward.forward = [realReplay.vHPC.reward.forward ; start stop];
%                         plot(tmp1(:,3) , tmp1(:,2),'*')
%                     end
%                     
%                 end
%                 
% %                 if p<0.05 %if
% %                     figure, histogram(shuffle,100)
% %                     xline(c(1,2))
% %                 end
%             end
%             clear probability nSpks start stop p c
%         end
% 
%         
%         % for dorsal hippocampus
%         realReplay.dHPC.reward.forward = [];
%         realReplay.dHPC.reward.backward = [];
%         for i = 1 : length(tmp.dHPC)
%             % Decoding using dHPC SUs
%             start = tmp.dHPC(i,1);
%             stop = tmp.dHPC(i,2);
%             bin = ((stop-start)/2) + start;
%             bin = [bin-0.1 : 0.02 : bin+0.1];
%             tmp1 = [];
%             for ii = 1 : length(bin)-1
%             nSpks = count_spks(spks_dHPC, clusterdHPC.normalizado.aversive, bin(ii), bin(ii+1));
%             probability = bayesian_replay(PHIST.dHPC.normalizado.aversive, nSpks, bin(ii), bin(ii+1));
%             
%             [o p] = max(probability);
%             
%             tmp1 = [tmp1 ; o p ii];
%             clear nSpks probability o p
%             end
%             c = corrcoef(tmp1(:,3) , tmp1(:,2));
% %             plot(tmp1(:,3) , tmp1(:,2),'*')
%             
%             shuffle = [];
%             for ii = 1:1000
%                p = corrcoef(tmp1(randperm(length(tmp1)),3) , tmp1(:,2));
%                 shuffle = [shuffle ; p(1,2)];
%             end
%             
%             if not(isnan(c(1,2))) %if correlation exist
%                 if sign(c(1,2))<0 %for reverse replay
%                     p = sum(shuffle<c(1,2))/1000;
%                     if p<0.05
%                         realReplay.dHPC.reward.backward = [realReplay.dHPC.reward.backward ; start stop];
%                         figure,plot(tmp1(:,3) , tmp1(:,2),'*')
%                     end
%                 else %for replay
%                     p = sum(shuffle>c(1,2))/1000;
%                     
%                     if p<0.05
%                         realReplay.dHPC.reward.forward = [realReplay.dHPC.reward.forward ; start stop];
%                         plot(tmp1(:,3) , tmp1(:,2),'*')
%                     end
%                     
%                 end
%                 
% %                 if p<0.05 %if
% %                     figure, histogram(shuffle,100)
% %                     xline(c(1,2))
% %                 end
%             end
%             clear probability nSpks start stop p c
%         end       
%         
%         
%         
%         
        %% Saving data
        if and(tt == 1 , t == 1)
            maps.dHPC.aversive = maps_dHPC_A;
            maps.vHPC.aversive = maps_vHPC_A;
            maps.dHPC.reward   = maps_dHPC_R;
            maps.vHPC.reward   = maps_vHPC_R;
            
            peak.dHPC.aversive = peak_dHPC_A;
            peak.vHPC.aversive = peak_vHPC_A;
            peak.dHPC.reward   = peak_dHPC_R;
            peak.vHPC.reward   = peak_vHPC_R;
            
            field.dHPC.aversive = field_dHPC_A;
            field.vHPC.aversive = field_vHPC_A;
            field.dHPC.reward   = field_dHPC_R;
            field.vHPC.reward   = field_vHPC_R;
            
            specificity.dHPC.aversive = specificity_dHPC_A;
            specificity.vHPC.aversive = specificity_vHPC_A;
            specificity.dHPC.reward   = specificity_dHPC_R;
            specificity.vHPC.reward   = specificity_vHPC_R;
            
        else
            maps.dHPC.aversive = [maps.dHPC.aversive ; maps_dHPC_A];
            maps.vHPC.aversive = [maps.vHPC.aversive ; maps_vHPC_A];
            maps.dHPC.reward   = [maps.dHPC.reward ; maps_dHPC_R];
            maps.vHPC.reward   = [maps.vHPC.reward ; maps_vHPC_R];
            
            
            peak.dHPC.aversive = [peak.dHPC.aversive ; peak_dHPC_A];
            peak.vHPC.aversive = [peak.vHPC.aversive ; peak_vHPC_A];
            peak.dHPC.reward   = [peak.dHPC.reward ; peak_dHPC_R];
            peak.vHPC.reward   = [peak.vHPC.reward ; peak_vHPC_R];
            
            field.dHPC.aversive = [field.dHPC.aversive ; field_dHPC_A];
            field.vHPC.aversive = [field.vHPC.aversive ; field_vHPC_A];
            field.dHPC.reward   = [field.dHPC.reward ; field_dHPC_R];
            field.vHPC.reward   = [field.vHPC.reward ; field_vHPC_R];
            
            specificity.dHPC.aversive = [specificity.dHPC.aversive ; specificity_dHPC_A];
            specificity.vHPC.aversive = [specificity.vHPC.aversive ; specificity_vHPC_A];
            specificity.dHPC.reward   = [specificity.dHPC.reward ; specificity_dHPC_R];
            specificity.vHPC.reward   = [specificity.vHPC.reward ; specificity_vHPC_R];            
        end
        clear specificity_dHPC_A specificity_dHPC_R specificity_vHPC_A specificity_vHPC_R
        clear field_dHPC_A field_dHPC_R field_vHPC_A field_vHPC_R
        clear peak_dHPC_A peak_dHPC_R peak_vHPC_A peak_vHPC_R
        clear maps_dHPC_A maps_dHPC_R maps_vHPC_A maps_vHPC_R
        clear aversiveTS aversiveTS_run rewardTS rewardTS_run baselineTS
        clear behavior camara Cell_type_classification cellulartype cluster
        clear group_dHPC group_vHPC K Kinfo leftvalve rightvalve movement R
        clear Rewards_filt Shocks_filt spks_dHPC spks_vHPC

        t
    end
    tt
end

%% Selection of cells according to their specificity

%Plot of specificity
data = [mean(specificity.dHPC.reward) , mean(specificity.dHPC.aversive)]
subplot(121),bar([1 2],data), hold on
semR = std(specificity.dHPC.reward)/sqrt(length(specificity.dHPC.reward));
semA = std(specificity.dHPC.aversive)/sqrt(length(specificity.dHPC.aversive));
errhigh = [data(1)+semR , data(2)+semA];
errlow = [data(1)-semR , data(2)-semA];

err = errorbar([1 2] ,data, errlow , errhigh)
err.Color = [0 0 0];                            


data = [mean(specificity.vHPC.reward) , mean(specificity.vHPC.aversive)]
subplot(122),bar([1 2],data), hold on
semR = std(specificity.vHPC.reward)/sqrt(length(specificity.vHPC.reward));
semA = std(specificity.vHPC.aversive)/sqrt(length(specificity.vHPC.aversive));
errhigh = [data(1)+semR , data(2)+semA];
errlow = [data(1)-semR , data(2)-semA];

err = errorbar([1 2] ,data, errlow , errhigh)
err.Color = [0 0 0];                            

%% Detection of place cells
% PC.dHPC = or(specificity.dHPC.aversive(:,1)>specificity.dHPC.aversive(:,2) , specificity.dHPC.reward(:,1)>specificity.dHPC.reward(:,2));
% PC.vHPC = or(specificity.vHPC.aversive(:,1)>specificity.vHPC.aversive(:,2) , specificity.vHPC.reward(:,1)>specificity.vHPC.reward(:,2));

PC.dHPC = or(specificity.dHPC.aversive(:,1)>0.25 , specificity.dHPC.reward(:,1)>0.25);
PC.vHPC = or(specificity.vHPC.aversive(:,1)>0.25, specificity.vHPC.reward(:,1)>0.25);

figure,
% detection of place field
data = maps.dHPC.aversive(PC.dHPC,:);
dataA = data;
data = data - min(data,[],2);
data = data ./max(data,[],2);
f = [];
for i = 1:size(data)
    [~,ind]=max(data(i,:));
    f = [f ; ind];
    clear ind
end
[~,i] = sort(f,'descend');
data = data(i,:);
subplot(1,2,1), imagesc(x,[1 : size(data,1)]',data),xlabel('Position'),title('Aversive'); axis xy;
data = maps.dHPC.reward(PC.dHPC,:);
dataR = data;
data = data - min(data,[],2);
data = data ./max(data,[],2);
data = data(i,:);
subplot(1,2,2), imagesc(x,[1 : size(data,1)]',data),xlabel('Position'),title('Reward'); axis xy;
 clear data

correlation_dHPC=[];
for i = 1: size(maps.dHPC.reward)
    c = corrcoef(maps.dHPC.reward(i,:),maps.dHPC.aversive(i,:));
    if isnan(c(1,2))
        break
    end
    correlation_dHPC=[correlation_dHPC ; c(1,2)];
end
clear dataR dataA




figure,
% detection of place field
data = maps.vHPC.aversive(PC.vHPC,:);
dataA = data;
data = data - min(data,[],2);
data = data ./max(data,[],2);
f = [];
for i = 1:size(data)
    [~,ind]=max(data(i,:));
    f = [f ; ind];
    clear ind
end


[~,i] = sort(f,'descend');
data = data(i,:);
subplot(1,2,1), imagesc(x,[1 : size(data,1)]',data),xlabel('Position'),title('Aversive'); axis xy;
data = maps.vHPC.reward(PC.vHPC,:);
dataR = data;
data = data - min(data,[],2);
data = data ./max(data,[],2);
data = data(i,:);
subplot(1,2,2), imagesc(x,[1 : size(data,1)]',data),xlabel('Position'),title('Reward'); axis xy;
 clear data

correlation_vHPC=[];
for i = 1: size(maps.vHPC.reward)
    c = corrcoef(maps.vHPC.reward(i,:),maps.vHPC.aversive(i,:));

    correlation_vHPC=[correlation_vHPC ; c(1,2)];
end
clear dataR dataA

figure
y = [correlation_dHPC(PC.dHPC) ; correlation_vHPC(PC.vHPC)];
grp = [ones(length(correlation_dHPC(PC.dHPC)),1) ; ones(length(correlation_vHPC(PC.vHPC)),1)*2];
boxplot(y,grp),ylim([-1 1])


figure,
subplot(121),boxplot(specificity.dHPC.reward(PC.dHPC,1) , specificity.dHPC.reward(PC.dHPC,3)),ylim([0 3])
subplot(122),boxplot(specificity.dHPC.aversive(PC.dHPC,1) , specificity.dHPC.aversive(PC.dHPC,3)),ylim([0 3])

figure,
subplot(121),boxplot(specificity.vHPC.reward(PC.vHPC,1) , specificity.vHPC.reward(PC.vHPC,3)),ylim([0 3])
subplot(122),boxplot(specificity.vHPC.aversive(PC.vHPC,1) , specificity.vHPC.aversive(PC.vHPC,3)),ylim([0 3])

figure
y = correlation_dHPC(and(PC.dHPC,specificity.dHPC.aversive(:,3)==1));
yy = correlation_dHPC(and(PC.dHPC,specificity.dHPC.aversive(:,3)==2));
grp = [ones(length(y),1) ; ones(length(yy),1)*2];

boxplot([y ; yy] , grp),ylim([-1 1])
ttest2(y , yy)


figure
y = correlation_vHPC(and(PC.vHPC,specificity.vHPC.aversive(:,3)==1))
yy = correlation_vHPC(and(PC.vHPC,specificity.vHPC.aversive(:,3)==2))
grp = [ones(length(y),1) ; ones(length(yy),1)*2];
ttest2(y , yy)


boxplot([y ; yy] , grp),ylim([-1 1])








%% Sparsity
figure,subplot(121)
y = specificity.dHPC.aversive(and(PC.dHPC,specificity.dHPC.aversive(:,4)==2),3);
yy = specificity.dHPC.reward(and(PC.dHPC,specificity.dHPC.reward(:,4)==2),3);
grp = [ones(length(y),1) ; ones(length(yy),1)*2];
ttest2(y , yy)

boxplot([y ; yy] , grp),ylim([0 1])


subplot(122)
y = specificity.vHPC.aversive(and(PC.vHPC,specificity.vHPC.aversive(:,4)==2),3);
yy = specificity.vHPC.reward(and(PC.vHPC,specificity.vHPC.reward(:,4)==2),3);
grp = [ones(length(y),1) ; ones(length(yy),1)*2];
ttest2(y , yy)

boxplot([y ; yy] , grp),ylim([0 1])




%% Sparsity
figure,subplot(121)
y = specificity.dHPC.aversive(PC.dHPC,1);
yy = specificity.dHPC.reward(PC.dHPC,1);
grp = [ones(length(y),1) ; ones(length(yy),1)*2];
ttest2(y , yy)

boxplot([y ; yy] , grp),ylim([0 2])


subplot(122)
y = specificity.vHPC.aversive(PC.vHPC,1);
yy = specificity.vHPC.reward(PC.vHPC,1);
grp = [ones(length(y),1) ; ones(length(yy),1)*2];
ttest2(y , yy)

boxplot([y ; yy] , grp),ylim([0 2])
