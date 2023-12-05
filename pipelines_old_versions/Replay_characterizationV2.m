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
sigma = 2;%round(15/(180/Xedges)); %defined for gauss kernel of 15cm
binSize = 0.001; % bin size for replay events detection
% Behavior
minimal_speed = 5; % minimal speed to detect quite periods
minimal_speed_time = 2; % minimal time to detect quite periods

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
        
        %% Firing maps calculation
        disp('dHPC Firing rate map calculation')
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
                if ~isnan(m)
                    % --- Aversive ---
                    spks_tmp = Restrict(spks , movement.aversive);
                    pos_tmp = Restrict(behavior.pos.aversive(:,1:2) , movement.aversive);
                    pos_tmp(:,2) = pos_tmp(:,2)-min(pos_tmp(:,2)); pos_tmp(:,2) = pos_tmp(:,2)/max(pos_tmp(:,2)); %normalization of position

                    %Firing curve construction
                    [curve , stats] = FiringCurve(pos_tmp , spks_tmp , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 2 , 'minPeak' , 0.2);
 
                    %Quantile for place-cell definition
                    q = SkaggsRandomFMT(spks_tmp, pos_tmp, sigma ,Xedges);
                    specificity_dHPC_A = [specificity_dHPC_A ; stats.specificity ,q , cond];
                    
                    % Store of place-field location and cluster of PCs
%                     if stats.specificity > q
                        [ff,f] = max(curve.rate);
                        fff = [1/Xedges:1/Xedges:1];
                        PC.dHPC.aversive=[PC.dHPC.aversive ; cluster , ff , fff(f)];
                        clear f ff fff
                        maps_dHPC_A = [maps_dHPC_A ; curve.rate];
                        peak_dHPC_A = [peak_dHPC_A ; stats.peak(1)] ;
                        field_dHPC_A = [field_dHPC_A ; stats.fieldX(1,:)];
%                     end
                    clear s q curve OccMap Nspikes spks_tmp pos_tmp m curve stats sp curve1
                    
                    % --- Reward ---
                    spks_tmp = Restrict(spks , movement.reward);
                    pos_tmp = Restrict(behavior.pos.reward(:,1:2) , movement.reward);
                    pos_tmp(:,2) = pos_tmp(:,2)-min(pos_tmp(:,2)); pos_tmp(:,2) = pos_tmp(:,2)/max(pos_tmp(:,2)); %normalization of position
 
                    %Firing curve construction
                    [curve , stats] = FiringCurve(pos_tmp , spks_tmp , 'smooth' , 2 , 'nBins' , Xedges , 'minSize' , 2 , 'minPeak' , 0.2);

                    %Quantile for place-cell definition
                    q = SkaggsRandomFMT(spks_tmp, pos_tmp, sigma, Xedges);
                    specificity_dHPC_R = [specificity_dHPC_R ; stats.specificity , q, cond];
                    
                    % Store of place-field location and cluster of PCs
%                     if stats.specificity > q
                        [ff,f] = max(curve.rate);
                        fff = [1/Xedges:1/Xedges:1];
                        PC.dHPC.reward=[PC.dHPC.reward ; cluster , ff , fff(f)];
                        clear f ff fff
                        maps_dHPC_R = [maps_dHPC_R ; curve.rate];
                        peak_dHPC_R = [peak_dHPC_R ; stats.peak(1)] ;
                        field_dHPC_R = [field_dHPC_R ; stats.fieldX(1,:)];
%                     end
                    clear s q curve OccMap Nspikes spks_tmp pos_tmp m curve stats sp curve1
                    cluster_dHPC = [cluster_dHPC ; cluster];
                end
            end
            clear celltype tmp b cluster
        end
        
        disp('vHPC Firing rate map calculation')
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
                    pos_tmp(:,2) = pos_tmp(:,2)-min(pos_tmp(:,2)); pos_tmp(:,2) = pos_tmp(:,2)/max(pos_tmp(:,2)); %normalization of position
                    
                    %Firing curve construction
                    [curve , stats] = FiringCurve(pos_tmp , spks_tmp , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4 , 'minPeak' , 0.2);
                    
                    %Quantile for place-cell definition
                    q = SkaggsRandomFMT(spks_tmp, pos_tmp, sigma, Xedges);
                    specificity_vHPC_A = [specificity_vHPC_A ; stats.specificity , q , cond];
                    
                    % Store of place-field location and cluster of PCs
%                     if stats.specificity > q
                        [ff,f] = max(curve.rate);
                        fff = [1/Xedges:1/Xedges:1];
                        PC.vHPC.aversive=[PC.vHPC.aversive ; cluster , ff , fff(f)];
                        clear f ff fff
                        maps_vHPC_A = [maps_vHPC_A ; curve.rate];
                        peak_vHPC_A = [peak_vHPC_A ; stats.peak(1)] ;
                        field_vHPC_A = [field_vHPC_A ; stats.fieldX(1,:)];
%                     end
                    clear s q curve OccMap Nspikes spks_tmp pos_tmp m curve stats sp curve1
                    
                    % --- Reward ---
                    spks_tmp = Restrict(spks , movement.reward);
                    pos_tmp = Restrict(behavior.pos.reward(:,1:2) , movement.reward);
                    pos_tmp(:,2) = pos_tmp(:,2)-min(pos_tmp(:,2)); pos_tmp(:,2) = pos_tmp(:,2)/max(pos_tmp(:,2)); %normalization of position
                    
                    %Firing curve construction
                    [curve , stats] = FiringCurve(pos_tmp , spks_tmp , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4 , 'minPeak' , 0.2);
                    
                    %Quantile for place-cell definition
                    q = SkaggsRandomFMT(spks_tmp, pos_tmp, sigma, Xedges);
                    specificity_vHPC_R = [specificity_vHPC_R ; stats.specificity , q , cond];
                    
                    % Store of place-field location and cluster of PCs
%                     if stats.specificity > q
                        [ff,f] = max(curve.rate);
                        fff = [1/Xedges:1/Xedges:1];
                        PC.vHPC.reward=[PC.vHPC.reward ; cluster , ff , fff(f)];
                        maps_vHPC_R = [maps_vHPC_R ; curve.rate];
                        peak_vHPC_R = [peak_vHPC_R ; stats.peak(1)] ;
                        field_vHPC_R = [field_vHPC_R ; stats.fieldX(1,:)];
                        clear f ff fff
%                     end
                    clear s q curve OccMap Nspikes spks_tmp pos_tmp m curve stats sp curve1
                    cluster_vHPC = [cluster_vHPC ; cluster];
                    %                     condition1 = or(isnan(specificity_vHPC_R(end,1)) , isnan(specificity_vHPC_A(end,1)));%store only those that contain spatial information
                    %                     condition2 = or(mean(maps_vHPC_R(end,:))<criteria_fr , mean(maps_vHPC_A(end,:))<criteria_fr); %store only those that fire at least 10 times
                    %                     if or(condition1 , condition2)
                    %                         specificity_vHPC_R(end,:) = [];
                    %                         peak_vHPC_R(end,:)        = [];
                    %                         field_vHPC_R(end,:)       = [];
                    %                         maps_vHPC_R(end,:)        = [];
                    %
                    %                         specificity_vHPC_A(end,:) = [];
                    %                         peak_vHPC_A(end,:)        = [];
                    %                         field_vHPC_A(end,:)       = [];
                    %                         maps_vHPC_A(end,:)        = [];
                    %
                    %                         cluster_vHPC(end,:) = [];
                    %
                    %                     end
                    %                     clear condition1 condition2
                end
            end
            clear celltype tmp b cluster
        end
        
%         c1 = and(size(PC.vHPC.aversive,1)>=3 , size(PC.vHPC.reward,1)>=3);
%         c2 = and(size(PC.dHPC.aversive,1)>=3 , size(PC.dHPC.reward,1)>=3);
        c1 = size(cluster_dHPC,1)>=3;
        c2 = size(cluster_vHPC,1)>=3;
        if and(c1 , c2)
            %% Replay detection
            disp('Replay detection in dHPC')
            freq = 1/0.001;
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
            [MUA.dHPC,Bins]=binspikes(sort(spks,'ascend'),freq,limits);clear spks
            replay.dHPC = FindReplay([Bins',MUA.dHPC],[1 3],[0.05 0.1 1],100);
            
            %filtering replay events by amount of PCs
            count = [];
            for i = 1 :  size(cluster_dHPC,1)
                cluster = cluster_dHPC(i,1);
                spks = spks_dHPC(spks_dHPC(:,1)==cluster,2);
                [status,interval,index] = InIntervals(spks,[replay.dHPC(:,1) replay.dHPC(:,3)]);
                interval = unique(interval);
                count = [count ; interval(interval~=0)];
                clear spks cluster status interval index
            end
            [gc,grps] = groupcounts(count);
            replay.dHPC = replay.dHPC(grps(gc>size(cluster_dHPC,1)*0.30),:);
            
            %Replay events in vHPC
            disp('Replay detection in vHPC')
            spks = [];
            for ii=1:size(group_vHPC,1)
                cluster = group_vHPC(ii,1);
                celltype = logical(cellulartype(cellulartype(:,1) == cluster,2));
                if celltype
                    spks = [spks;spks_vHPC(spks_vHPC(:,1)==cluster,2)];
                end
                clear celltype
            end
            [MUA.vHPC,Bins]=binspikes(sort(spks,'ascend'),freq,limits);clear spks
            replay.vHPC = FindReplay([Bins',MUA.vHPC],[1 3],[0.05 0.1 1],100);
            
            %filtering replay events by amount of PCs
            count = [];
            for i = 1 :  size(cluster_vHPC,1)
                cluster = cluster_vHPC(i,1);
                spks = spks_vHPC(spks_vHPC(:,1)==cluster,2);
                [status,interval,index] = InIntervals(spks,[replay.vHPC(:,1) , replay.vHPC(:,3)]);
                interval = unique(interval);
                count = [count ; interval(interval~=0)];
                clear spks cluster status interval index
            end
            [gc,grps] = groupcounts(count);
            replay.vHPC = replay.vHPC(grps(gc>size(cluster_vHPC,1)*0.30),:);
            
            %% Restriction to sleep replays
            ss = [REM.all ; NREM.all];
            [~ , s] = sort(ss(:,1));
            replay.vHPC = Restrict(replay.vHPC , ss(s,:));
            replay.dHPC = Restrict(replay.dHPC , ss(s,:));
            clear s ss
            
            %% Bayesian Decoding across conditions
            disp('Bayesian decoding across structures and conditions')
            % --- Aversive ---
            tmp.vHPC = replay.vHPC;
            tmp.dHPC = replay.dHPC;
            % for ventral hippocampus
            %         rZ.vHPC = [];
            r.vHPC.all = [];
            for i = 1 : length(tmp.vHPC)
                % Decoding using dHPC SUs
                start = tmp.vHPC(i,1);
                stop = tmp.vHPC(i,3);
                bin = ((stop-start)/2) + start;
                center = bin;
                %             bin = [start : 0.02 : stop];
                %             tmp1 = [];            tmp2 = [];
                %             for ii = 1 : length(bin)-1
                %             nSpks = count_spks(spks_vHPC, PC.vHPC.aversive(:,1), bin(ii), bin(ii+1));
                %             probability = bayesian_replay(maps_vHPC_A, nSpks, bin(ii), bin(ii+1));
                %
                %             tmp1 = [tmp1 , probability'];
                %             clear nSpks probability o p
                %
                %             nSpks = count_spks(spks_vHPC, PC.vHPC.reward(:,1), bin(ii), bin(ii+1));
                %             probability = bayesian_replay(maps_vHPC_R, nSpks, bin(ii), bin(ii+1));
                %
                %             tmp2 = [tmp2 , probability'];
                %             clear nSpks probability o p
                %             end
                %
                %             p = [1:Xedges];
                %             bin = bin(1:end-1);
                %
                %             % Aversive
                %             m.pos = sum(sum(tmp1.* p')) / sum(sum(tmp1));
                %             m.bin = sum(sum(tmp1.*bin)) / sum(sum(tmp1));
                %
                %             cova.pos_bin = sum(sum(tmp1'.*(p-m.pos).*(bin'-m.bin))) / sum(sum(tmp1));
                %             cova.pos_pos = sum(sum(tmp1'.*(p-m.pos).*(p-m.pos))) / sum(sum(tmp1));
                %             cova.bin_bin = sum(sum(tmp1.*(bin-m.bin).*(bin-m.bin))) / sum(sum(tmp1));
                %
                %             rA = cova.pos_bin / sqrt(cova.pos_pos * (cova.bin_bin));
                %             clear m cova
                %
                %             % Reward
                %             m.pos = sum(sum(tmp2.* p')) / sum(sum(tmp2));
                %             m.bin = sum(sum(tmp2.*bin)) / sum(sum(tmp2));
                %
                %             cova.pos_bin = sum(sum(tmp2'.*(p-m.pos).*(bin'-m.bin))) / sum(sum(tmp2));
                %             cova.pos_pos = sum(sum(tmp2'.*(p-m.pos).*(p-m.pos))) / sum(sum(tmp2));
                %             cova.bin_bin = sum(sum(tmp2.*(bin-m.bin).*(bin-m.bin))) / sum(sum(tmp2));
                %
                %             rR = cova.pos_bin / sqrt(cova.pos_pos * (cova.bin_bin));
                %             clear m cova
                %
                %             rZ.vHPC = [rZ.vHPC ; center , (abs(rA)-abs(rR))/(abs(rA)+abs(rR))];
                nSpks = count_spks(spks_vHPC, PC.vHPC.aversive(:,1), start , stop);
                probabilityA = bayesian_replay(maps_vHPC_A, nSpks, start , stop);
                nSpks = count_spks(spks_vHPC, PC.vHPC.reward(:,1), start , stop);
                probabilityR = bayesian_replay(maps_vHPC_R, nSpks, start , stop);
                r.vHPC.all = [r.vHPC.all ; center , (max(probabilityA)- max(probabilityR)) / (max(probabilityA)+max(probabilityR))];
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
                %                         realReplay.vHPC.aversive.backward = [realReplay.vHPC.aversive.backward ; start stop];
                %                         figure,plot(tmp1(:,3) , tmp1(:,2),'*')
                %                     end
                %                 else %for replay
                %                     p = sum(shuffle>c(1,2))/1000;
                %
                %                     if p<0.05
                %                         realReplay.vHPC.aversive.forward = [realReplay.vHPC.aversive.forward ; start stop];
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
                clear probability nSpks start stop p c bin m cova p probabilityA probabilityR
                clear bin center
            end
            disp('Finished ventral Replays events')
            %         figure,stem(r.vHPC.all(:,1),r.vHPC.all(:,2))
            %         hold on
            %         xline(segments.Var1(3)/1000,'--')
            %         xline(segments.Var1(4)/1000,'--')
            %         xline(segments.Var1(7)/1000,'--')
            %         xline(segments.Var1(8)/1000,'--')
            %
            r.vHPC.sws = Restrict(r.vHPC.all,NREM.all);
            %         figure,stem(r.vHPC.sws(:,1),r.vHPC.sws(:,2))
            %         hold on
            %         xline(segments.Var1(3)/1000,'--')
            %         xline(segments.Var1(4)/1000,'--')
            %         xline(segments.Var1(7)/1000,'--')
            %         xline(segments.Var1(8)/1000,'--')
            %
            r.vHPC.rem = Restrict(r.vHPC.all,REM.all);
            %         figure,stem(r.vHPC.sws(:,1),r.vHPC.sws(:,2))
            %         hold on
            %         xline(segments.Var1(3)/1000,'--')
            %         xline(segments.Var1(4)/1000,'--')
            %         xline(segments.Var1(7)/1000,'--')
            %         xline(segments.Var1(8)/1000,'--')
            %
            %
            r.vHPC.burst = Restrict(r.vHPC.all,[ripple_bursts.all(:,1) ripple_bursts.all(:,3)]);
            %         figure,stem(r.vHPC.burst(:,1),r.vHPC.burst(:,2))
            %         hold on
            %         xline(segments.Var1(3)/1000,'--')
            %         xline(segments.Var1(4)/1000,'--')
            %         xline(segments.Var1(7)/1000,'--')
            %         xline(segments.Var1(8)/1000,'--')
            %
            r.vHPC.ripple = Restrict(r.vHPC.all,[ripplesV(:,1) ripplesV(:,3)]);
            %         figure,stem(r.vHPC.ripple(:,1),r.vHPC.ripple(:,2))
            %         hold on
            %         xline(segments.Var1(3)/1000,'--')
            %         xline(segments.Var1(4)/1000,'--')
            %         xline(segments.Var1(7)/1000,'--')
            %         xline(segments.Var1(8)/1000,'--')
            %
            r.vHPC.coordinated = Restrict(r.vHPC.all,[coordinatedV(:,2)-0.1 coordinatedV(:,2)+0.1]);
            %         figure,stem(r.vHPC.coordinated(:,1),r.vHPC.coordinated(:,2))
            %         hold on
            %         xline(segments.Var1(3)/1000,'--')
            %         xline(segments.Var1(4)/1000,'--')
            %         xline(segments.Var1(7)/1000,'--')
            %         xline(segments.Var1(8)/1000,'--')
            %
            % for dorsal hippocampus
            %         rZ.dHPC = [];
            r.dHPC.all = [];
            for i = 1 : length(tmp.dHPC)
                % Decoding using dHPC SUs
                start = tmp.dHPC(i,1);
                stop = tmp.dHPC(i,3);
                bin = ((stop-start)/2) + start;
                center = bin;
                %             bin = [start : 0.02 : stop];
                %             tmp1 = [];            tmp2 = [];
                %             for ii = 1 : length(bin)-1
                %             nSpks = count_spks(spks_dHPC, PC.dHPC.aversive(:,1), bin(ii), bin(ii+1));
                %             probability = bayesian_replay(maps_dHPC_A, nSpks, bin(ii), bin(ii+1));
                %
                %             tmp1 = [tmp1 , probability'];
                %             clear nSpks probability o p
                %
                %             nSpks = count_spks(spks_dHPC, PC.dHPC.reward(:,1), bin(ii), bin(ii+1));
                %             probability = bayesian_replay(maps_dHPC_R, nSpks, bin(ii), bin(ii+1));
                %
                %             tmp2 = [tmp2 , probability'];
                %             clear nSpks probability o p
                %             end
                %
                %             p = [1:Xedges];
                %             bin = bin(1:end-1);
                %
                %             % Aversive
                %             m.pos = sum(sum(tmp1.* p')) / sum(sum(tmp1));
                %             m.bin = sum(sum(tmp1.*bin)) / sum(sum(tmp1));
                %
                %             cova.pos_bin = sum(sum(tmp1'.*(p-m.pos).*(bin'-m.bin))) / sum(sum(tmp1));
                %             cova.pos_pos = sum(sum(tmp1'.*(p-m.pos).*(p-m.pos))) / sum(sum(tmp1));
                %             cova.bin_bin = sum(sum(tmp1.*(bin-m.bin).*(bin-m.bin))) / sum(sum(tmp1));
                %
                %             rA = cova.pos_bin / sqrt(cova.pos_pos * (cova.bin_bin));
                %             clear m cova
                %
                %             % Reward
                %             m.pos = sum(sum(tmp2.* p')) / sum(sum(tmp2));
                %             m.bin = sum(sum(tmp2.*bin)) / sum(sum(tmp2));
                %
                %             cova.pos_bin = sum(sum(tmp2'.*(p-m.pos).*(bin'-m.bin))) / sum(sum(tmp2));
                %             cova.pos_pos = sum(sum(tmp2'.*(p-m.pos).*(p-m.pos))) / sum(sum(tmp2));
                %             cova.bin_bin = sum(sum(tmp2.*(bin-m.bin).*(bin-m.bin))) / sum(sum(tmp2));
                %
                %             rR = cova.pos_bin / sqrt(cova.pos_pos * (cova.bin_bin));
                %             clear m cova
                %
                %             rZ.dHPC = [rZ.dHPC ; center , (abs(rA)-abs(rR))/(abs(rA)+abs(rR))];
                nSpks = count_spks(spks_dHPC, PC.dHPC.aversive(:,1), start , stop);
                probabilityA = bayesian_replay(maps_dHPC_A, nSpks, start , stop);
                nSpks = count_spks(spks_dHPC, PC.dHPC.reward(:,1), start , stop);
                probabilityR = bayesian_replay(maps_dHPC_R, nSpks, start , stop);
                r.dHPC.all = [r.dHPC.all ; center , (max(probabilityA)- max(probabilityR)) / (max(probabilityA)+max(probabilityR))];
                clear probability nSpks start stop p c bin m cova p rR rA probabilityA probabilityR
                clear bin center
            end
            disp('Finished dorsal Replays events in aversive sleep')
            %         figure,stem(r.dHPC.all(:,1),r.dHPC.all(:,2))
            %         hold on
            %         xline(segments.Var1(3)/1000,'--')
            %         xline(segments.Var1(4)/1000,'--')
            %         xline(segments.Var1(7)/1000,'--')
            %         xline(segments.Var1(8)/1000,'--')
            %
            r.dHPC.sws = Restrict(r.dHPC.all,NREM.all);
            %         figure,stem(r.dHPC.sws(:,1),r.dHPC.sws(:,2))
            %         hold on
            %         xline(segments.Var1(3)/1000,'--')
            %         xline(segments.Var1(4)/1000,'--')
            %         xline(segments.Var1(7)/1000,'--')
            %         xline(segments.Var1(8)/1000,'--')
            %
            r.dHPC.rem = Restrict(r.dHPC.all,REM.all);
            %         figure,stem(r.dHPC.sws(:,1),r.dHPC.sws(:,2))
            %         hold on
            %         xline(segments.Var1(3)/1000,'--')
            %         xline(segments.Var1(4)/1000,'--')
            %         xline(segments.Var1(7)/1000,'--')
            %         xline(segments.Var1(8)/1000,'--')
            %
            r.dHPC.burst = Restrict(r.dHPC.all,[ripple_bursts.all(:,1) ripple_bursts.all(:,3)]);
            %         figure,stem(r.dHPC.burst(:,1),r.dHPC.burst(:,2))
            %         hold on
            %         xline(segments.Var1(3)/1000,'--')
            %         xline(segments.Var1(4)/1000,'--')
            %         xline(segments.Var1(7)/1000,'--')
            %         xline(segments.Var1(8)/1000,'--')
            %
            r.dHPC.ripple = Restrict(r.dHPC.all,[ripplesD(:,1) ripplesD(:,3)]);
            %         figure,stem(r.dHPC.ripple(:,1),r.dHPC.ripple(:,2))
            %         hold on
            %         xline(segments.Var1(3)/1000,'--')
            %         xline(segments.Var1(4)/1000,'--')
            %         xline(segments.Var1(7)/1000,'--')
            %         xline(segments.Var1(8)/1000,'--')
            %
            r.dHPC.coordinated = Restrict(r.dHPC.all,[coordinated(:,2)-0.1 coordinated(:,2)+0.1]);
            %         figure,stem(r.dHPC.coordinated(:,1),r.dHPC.coordinated(:,2))
            %         hold on
            %         xline(segments.Var1(3)/1000,'--')
            %         xline(segments.Var1(4)/1000,'--')
            %         xline(segments.Var1(7)/1000,'--')
            %         xline(segments.Var1(8)/1000,'--')
            %% Saving data
            disp('Saving outputs')
            if not(exist('maps','var'))
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
                %             WeightedCorrelation.dHPC.aversive = rZ.dHPC;
                %             WeightedCorrelation.vHPC.aversive = rZ.vHPC;
                % Emotional Ccndition Index
                %REM dHPC
                bias.dHPC.rem.baseline = Restrict(r.dHPC.all,REM.baseline);
                bias.dHPC.rem.reward = Restrict(r.dHPC.all,REM.reward);
                bias.dHPC.rem.aversive = Restrict(r.dHPC.all,REM.aversive);
                %REM vHPC
                bias.vHPC.rem.baseline = Restrict(r.vHPC.all,REM.baseline);
                bias.vHPC.rem.reward = Restrict(r.vHPC.all,REM.reward);
                bias.vHPC.rem.aversive = Restrict(r.vHPC.all,REM.aversive);
                %NREM dHPC
                bias.dHPC.sws.baseline = Restrict(r.dHPC.all,NREM.baseline);
                bias.dHPC.sws.reward = Restrict(r.dHPC.all,NREM.reward);
                bias.dHPC.sws.aversive = Restrict(r.dHPC.all,NREM.aversive);
                %NREM vHPC
                bias.vHPC.sws.baseline = Restrict(r.vHPC.all,NREM.baseline);
                bias.vHPC.sws.reward = Restrict(r.vHPC.all,NREM.reward);
                bias.vHPC.sws.aversive = Restrict(r.vHPC.all,NREM.aversive);
                %Ripple Bursts dHPC
                bias.dHPC.burst.baseline = Restrict(r.dHPC.burst,NREM.baseline);
                bias.dHPC.burst.reward = Restrict(r.dHPC.burst,NREM.reward);
                bias.dHPC.burst.aversive = Restrict(r.dHPC.burst,NREM.aversive);
                %Ripple Bursts vHPC
                bias.vHPC.burst.baseline = Restrict(r.vHPC.burst,NREM.baseline);
                bias.vHPC.burst.reward = Restrict(r.vHPC.burst,NREM.reward);
                bias.vHPC.burst.aversive = Restrict(r.vHPC.burst,NREM.aversive);
                %Cooridnated ripples dHPC
                bias.dHPC.coordinated.baseline = Restrict(r.dHPC.coordinated,NREM.baseline);
                bias.dHPC.coordinated.reward = Restrict(r.dHPC.coordinated,NREM.reward);
                bias.dHPC.coordinated.aversive = Restrict(r.dHPC.coordinated,NREM.aversive);
                %Cooridnated ripples vHPC
                bias.vHPC.coordinated.baseline = Restrict(r.vHPC.coordinated,NREM.baseline);
                bias.vHPC.coordinated.reward = Restrict(r.vHPC.coordinated,NREM.reward);
                bias.vHPC.coordinated.aversive = Restrict(r.vHPC.coordinated,NREM.aversive);
                %Ripples dHPC
                bias.dHPC.ripple.baseline = Restrict(r.dHPC.ripple,NREM.baseline);
                bias.dHPC.ripple.reward = Restrict(r.dHPC.ripple,NREM.reward);
                bias.dHPC.ripple.aversive = Restrict(r.dHPC.ripple,NREM.aversive);
                %Ripples vHPC
                bias.vHPC.ripple.baseline = Restrict(r.vHPC.ripple,NREM.baseline);
                bias.vHPC.ripple.reward = Restrict(r.vHPC.ripple,NREM.reward);
                bias.vHPC.ripple.aversive = Restrict(r.vHPC.ripple,NREM.aversive);
                
                %Saving means
                %REM dHPC
                m = Restrict(r.dHPC.all,REM.baseline);
                bias.dHPC.mean.rem.baseline = mean(m(:,2)); clear m
                m = Restrict(r.dHPC.all,REM.reward);
                bias.dHPC.mean.rem.reward = mean(m(:,2)); clear m
                m = Restrict(r.dHPC.all,REM.aversive);
                bias.dHPC.mean.rem.aversive = mean(m(:,2)); clear m
                
                %REM vHPC
                m = Restrict(r.vHPC.all,REM.baseline);
                bias.vHPC.mean.rem.baseline = mean(m(:,2)); clear m
                m = Restrict(r.vHPC.all,REM.reward);
                bias.vHPC.mean.rem.reward = mean(m(:,2)); clear m
                m = Restrict(r.vHPC.all,REM.aversive);
                bias.vHPC.mean.rem.aversive = mean(m(:,2)); clear m
                
                %NREM dHPC
                m = Restrict(r.dHPC.all,NREM.baseline);
                bias.dHPC.mean.sws.baseline = mean(m(:,2)); clear m
                m = Restrict(r.dHPC.all,NREM.reward);
                bias.dHPC.mean.sws.reward = mean(m(:,2)); clear m
                m = Restrict(r.dHPC.all,NREM.aversive);
                bias.dHPC.mean.sws.aversive = mean(m(:,2)); clear m
                
                %NREM vHPC
                m = Restrict(r.vHPC.all,NREM.baseline);
                bias.vHPC.mean.sws.baseline = mean(m(:,2)); clear m
                m = Restrict(r.vHPC.all,NREM.reward);
                bias.vHPC.mean.sws.reward = mean(m(:,2)); clear m
                m = Restrict(r.vHPC.all,NREM.aversive);
                bias.vHPC.mean.sws.aversive = mean(m(:,2)); clear m
                
                %Ripple Bursts dHPC
                m = Restrict(r.dHPC.burst,NREM.baseline);
                bias.dHPC.mean.burst.baseline = mean(m(:,2)); clear m
                m = Restrict(r.dHPC.burst,NREM.reward);
                bias.dHPC.mean.burst.reward = mean(m(:,2)); clear m
                m = Restrict(r.dHPC.burst,NREM.aversive);
                bias.dHPC.mean.burst.aversive = mean(m(:,2)); clear m
                
                %Ripple Bursts vHPC
                m = Restrict(r.vHPC.burst,NREM.baseline);
                bias.vHPC.mean.burst.baseline = mean(m(:,2)); clear m
                m = Restrict(r.vHPC.burst,NREM.reward);
                bias.vHPC.mean.burst.reward = mean(m(:,2)); clear m
                m = Restrict(r.vHPC.burst,NREM.aversive);
                bias.vHPC.mean.burst.aversive = mean(m(:,2)); clear m
                
                %Cooridnated ripples dHPC
                m = Restrict(r.dHPC.coordinated,NREM.baseline);
                bias.dHPC.mean.coordinated.baseline = mean(m(:,2)); clear m
                m = Restrict(r.dHPC.coordinated,NREM.reward);
                bias.dHPC.mean.coordinated.reward = mean(m(:,2)); clear m
                m = Restrict(r.dHPC.coordinated,NREM.aversive);
                bias.dHPC.mean.coordinated.aversive = mean(m(:,2)); clear m
                
                %Cooridnated ripples vHPC
                m = Restrict(r.vHPC.coordinated,NREM.baseline);
                bias.vHPC.mean.coordinated.baseline = mean(m(:,2)); clear m
                m = Restrict(r.vHPC.coordinated,NREM.reward);
                bias.vHPC.mean.coordinated.reward = mean(m(:,2)); clear m
                m = Restrict(r.vHPC.coordinated,NREM.aversive);
                bias.vHPC.mean.coordinated.aversive = mean(m(:,2)); clear m
                
                %Ripples dHPC
                m = Restrict(r.dHPC.ripple,NREM.baseline);
                bias.dHPC.mean.ripple.baseline = mean(m(:,2)); clear m
                m = Restrict(r.dHPC.ripple,NREM.reward);
                bias.dHPC.mean.ripple.reward = mean(m(:,2)); clear m
                m = Restrict(r.dHPC.ripple,NREM.aversive);
                bias.dHPC.mean.ripple.aversive = mean(m(:,2)); clear m
                
                %Ripples vHPC
                m = Restrict(r.vHPC.ripple,NREM.baseline);
                bias.vHPC.mean.ripple.baseline = mean(m(:,2)); clear m
                m = Restrict(r.vHPC.ripple,NREM.reward);
                bias.vHPC.mean.ripple.reward = mean(m(:,2)); clear m
                m = Restrict(r.vHPC.ripple,NREM.aversive);   
                bias.vHPC.mean.ripple.aversive =  mean(m(:,2)); clear m           
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
                %             WeightedCorrelation.dHPC = [WeightedCorrelation.dHPC ; rZ.dHPC];
                %             WeightedCorrelation.vHPC = [WeightedCorrelation.vHPC ; rZ.vHPC];
                % Emotional Ccndition Index
                %REM dHPC
                bias.dHPC.rem.baseline = [bias.dHPC.rem.baseline ; Restrict(r.dHPC.all,REM.baseline)];
                bias.dHPC.rem.reward = [bias.dHPC.rem.reward ; Restrict(r.dHPC.all,REM.reward)];
                bias.dHPC.rem.aversive = [bias.dHPC.rem.aversive ; Restrict(r.dHPC.all,REM.aversive)];
                %REM vHPC
                bias.vHPC.rem.baseline = [bias.vHPC.rem.baseline ; Restrict(r.vHPC.all,REM.baseline)];
                bias.vHPC.rem.reward = [bias.vHPC.rem.reward ; Restrict(r.vHPC.all,REM.reward)];
                bias.vHPC.rem.aversive = [bias.vHPC.rem.aversive ; Restrict(r.vHPC.all,REM.aversive)];
                %NREM dHPC
                bias.dHPC.sws.baseline = [bias.dHPC.sws.baseline ; Restrict(r.dHPC.all,NREM.baseline)];
                bias.dHPC.sws.reward = [bias.dHPC.sws.reward ; Restrict(r.dHPC.all,NREM.reward)];
                bias.dHPC.sws.aversive = [bias.dHPC.sws.aversive ; Restrict(r.dHPC.all,NREM.aversive)];
                %NREM vHPC
                bias.vHPC.sws.baseline = [bias.vHPC.sws.baseline ; Restrict(r.vHPC.all,NREM.baseline)];
                bias.vHPC.sws.reward = [bias.vHPC.sws.reward ; Restrict(r.vHPC.all,NREM.reward)];
                bias.vHPC.sws.aversive = [bias.vHPC.sws.aversive ; Restrict(r.vHPC.all,NREM.aversive)];
                %Ripple Bursts dHPC
                if and(and(not(isempty(Restrict(r.dHPC.burst,NREM.baseline))) , not(isempty(Restrict(r.dHPC.burst,NREM.reward)))) ,  not(isempty(Restrict(r.dHPC.burst,NREM.aversive))))
                    bias.dHPC.burst.baseline = [bias.dHPC.burst.baseline ; Restrict(r.dHPC.burst,NREM.baseline)];
                    bias.dHPC.burst.reward = [bias.dHPC.burst.reward ; Restrict(r.dHPC.burst,NREM.reward)];
                    bias.dHPC.burst.aversive = [bias.dHPC.burst.aversive ; Restrict(r.dHPC.burst,NREM.aversive)];
                end
                %Ripple Bursts vHPC
                if and(and(not(isempty(Restrict(r.vHPC.burst,NREM.baseline))) , not(isempty(Restrict(r.vHPC.burst,NREM.reward)))) ,  not(isempty(Restrict(r.vHPC.burst,NREM.aversive))))
                    bias.vHPC.burst.baseline = [bias.vHPC.burst.baseline ; Restrict(r.vHPC.burst,NREM.baseline)];
                    bias.vHPC.burst.reward = [bias.vHPC.burst.reward ; Restrict(r.vHPC.burst,NREM.reward)];
                    bias.vHPC.burst.aversive = [bias.vHPC.burst.aversive ; Restrict(r.vHPC.burst,NREM.aversive)];
                end
                %Cooridnated ripples dHPC
                if and(and(not(isempty(Restrict(r.dHPC.coordinated,NREM.baseline))) , not(isempty(Restrict(r.dHPC.coordinated,NREM.reward)))) ,  not(isempty(Restrict(r.dHPC.coordinated,NREM.aversive))))
                    bias.dHPC.coordinated.baseline = [bias.dHPC.coordinated.baseline ; Restrict(r.dHPC.coordinated,NREM.baseline)];
                    bias.dHPC.coordinated.reward = [bias.dHPC.coordinated.reward ; Restrict(r.dHPC.coordinated,NREM.reward)];
                    bias.dHPC.coordinated.aversive = [bias.dHPC.coordinated.aversive ; Restrict(r.dHPC.coordinated,NREM.aversive)];
                end
                %Cooridnated ripples vHPC
                if and(and(not(isempty(Restrict(r.vHPC.coordinated,NREM.baseline))) , not(isempty(Restrict(r.vHPC.coordinated,NREM.reward)))) ,  not(isempty(Restrict(r.vHPC.coordinated,NREM.aversive))))
                    bias.vHPC.coordinated.baseline = [bias.vHPC.coordinated.baseline ; Restrict(r.vHPC.coordinated,NREM.baseline)];
                    bias.vHPC.coordinated.reward = [bias.vHPC.coordinated.reward ; Restrict(r.vHPC.coordinated,NREM.reward)];
                    bias.vHPC.coordinated.aversive = [bias.vHPC.coordinated.aversive ; Restrict(r.vHPC.coordinated,NREM.aversive)];
                end
                %Ripples dHPC
                if and(and(not(isempty(Restrict(r.dHPC.ripple,NREM.baseline))) , not(isempty(Restrict(r.dHPC.ripple,NREM.reward)))) ,  not(isempty(Restrict(r.dHPC.ripple,NREM.aversive))))
                    bias.dHPC.ripple.baseline = [bias.dHPC.ripple.baseline ; Restrict(r.dHPC.ripple,NREM.baseline)];
                    bias.dHPC.ripple.reward = [bias.dHPC.ripple.reward ; Restrict(r.dHPC.ripple,NREM.reward)];
                    bias.dHPC.ripple.aversive = [bias.dHPC.ripple.aversive ; Restrict(r.dHPC.ripple,NREM.aversive)];
                end
                %Ripples vHPC
                if and(and(not(isempty(Restrict(r.vHPC.coordinated,NREM.baseline))) , not(isempty(Restrict(r.vHPC.coordinated,NREM.reward)))) ,  not(isempty(Restrict(r.vHPC.coordinated,NREM.aversive))))
                    bias.vHPC.ripple.baseline = [bias.vHPC.ripple.baseline ; Restrict(r.vHPC.ripple,NREM.baseline)];
                    bias.vHPC.ripple.reward = [bias.vHPC.ripple.reward ; Restrict(r.vHPC.ripple,NREM.reward)];
                    bias.vHPC.ripple.aversive = [bias.vHPC.ripple.aversive ; Restrict(r.vHPC.ripple,NREM.aversive)];
                end
                
                %Saving means
                %REM dHPC
                m = Restrict(r.dHPC.all,REM.baseline);
                bias.dHPC.mean.rem.baseline = [bias.dHPC.mean.rem.baseline ; mean(m(:,2))]; clear m
                m = Restrict(r.dHPC.all,REM.reward);
                bias.dHPC.mean.rem.reward = [bias.dHPC.mean.rem.reward ; mean(m(:,2))]; clear m
                m = Restrict(r.dHPC.all,REM.aversive);
                bias.dHPC.mean.rem.aversive = [bias.dHPC.mean.rem.aversive ; mean(m(:,2))]; clear m
                
                %REM vHPC
                m = Restrict(r.vHPC.all,REM.baseline);
                bias.vHPC.mean.rem.baseline = [bias.vHPC.mean.rem.baseline ; mean(m(:,2))]; clear m
                m = Restrict(r.vHPC.all,REM.reward);
                bias.vHPC.mean.rem.reward = [bias.vHPC.mean.rem.reward ; mean(m(:,2))]; clear m
                m = Restrict(r.vHPC.all,REM.aversive);
                bias.vHPC.mean.rem.aversive = [bias.vHPC.mean.rem.aversive ; mean(m(:,2))]; clear m
                
                %NREM dHPC
                m = Restrict(r.dHPC.all,NREM.baseline);
                bias.dHPC.mean.sws.baseline = [bias.dHPC.mean.sws.baseline ; mean(m(:,2))]; clear m
                m = Restrict(r.dHPC.all,NREM.reward);
                bias.dHPC.mean.sws.reward = [bias.dHPC.mean.sws.reward ; mean(m(:,2))]; clear m
                m = Restrict(r.dHPC.all,NREM.aversive);
                bias.dHPC.mean.sws.aversive = [bias.dHPC.mean.sws.aversive ; mean(m(:,2))]; clear m
                
                %NREM vHPC
                m = Restrict(r.vHPC.all,NREM.baseline);
                bias.vHPC.mean.sws.baseline = [bias.vHPC.mean.sws.baseline ; mean(m(:,2))]; clear m
                m = Restrict(r.vHPC.all,NREM.reward);
                bias.vHPC.mean.sws.reward = [bias.vHPC.mean.sws.reward ; mean(m(:,2))]; clear m
                m = Restrict(r.vHPC.all,NREM.aversive);
                bias.vHPC.mean.sws.aversive = [bias.vHPC.mean.sws.aversive ; mean(m(:,2))]; clear m
                
                %Ripple Bursts dHPC
                if and(and(not(isempty(Restrict(r.dHPC.burst,NREM.baseline))) , not(isempty(Restrict(r.dHPC.burst,NREM.reward)))) ,  not(isempty(Restrict(r.dHPC.burst,NREM.aversive))))
                    m = Restrict(r.dHPC.burst,NREM.baseline);
                    bias.dHPC.mean.burst.baseline = [bias.dHPC.mean.burst.baseline ; mean(m(:,2))]; clear m
                    m = Restrict(r.dHPC.burst,NREM.reward);
                    bias.dHPC.mean.burst.reward = [bias.dHPC.mean.burst.reward ; mean(m(:,2))]; clear m
                    m = Restrict(r.dHPC.burst,NREM.aversive);
                    bias.dHPC.mean.burst.aversive = [bias.dHPC.mean.burst.aversive ; mean(m(:,2))]; clear m
                end
                
                %Ripple Bursts vHPC
                if and(and(not(isempty(Restrict(r.vHPC.burst,NREM.baseline))) , not(isempty(Restrict(r.vHPC.burst,NREM.reward)))) ,  not(isempty(Restrict(r.vHPC.burst,NREM.aversive))))
                    m = Restrict(r.vHPC.burst,NREM.baseline);
                    bias.vHPC.mean.burst.baseline = [bias.vHPC.mean.burst.baseline ; mean(m(:,2))]; clear m
                    m = Restrict(r.vHPC.burst,NREM.reward);
                    bias.vHPC.mean.burst.reward = [bias.vHPC.mean.burst.reward ; mean(m(:,2))]; clear m
                    m = Restrict(r.vHPC.burst,NREM.aversive);
                    bias.vHPC.mean.burst.aversive = [bias.vHPC.mean.burst.aversive ; mean(m(:,2))]; clear m
                end
                
                %Cooridnated ripples dHPC
                if and(and(not(isempty(Restrict(r.dHPC.coordinated,NREM.baseline))) , not(isempty(Restrict(r.dHPC.coordinated,NREM.reward)))) ,  not(isempty(Restrict(r.dHPC.coordinated,NREM.aversive))))
                    m = Restrict(r.dHPC.coordinated,NREM.baseline);
                    bias.dHPC.mean.coordinated.baseline = [bias.dHPC.mean.coordinated.baseline ; mean(m(:,2))]; clear m
                    m = Restrict(r.dHPC.coordinated,NREM.reward);
                    bias.dHPC.mean.coordinated.reward = [bias.dHPC.mean.coordinated.reward ; mean(m(:,2))]; clear m
                    m = Restrict(r.dHPC.coordinated,NREM.aversive);
                    bias.dHPC.mean.coordinated.aversive = [bias.dHPC.mean.coordinated.aversive ; mean(m(:,2))]; clear m
                end
                
                %Cooridnated ripples vHPC
                if and(and(not(isempty(Restrict(r.vHPC.coordinated,NREM.baseline))) , not(isempty(Restrict(r.vHPC.coordinated,NREM.reward)))) ,  not(isempty(Restrict(r.vHPC.coordinated,NREM.aversive))))
                    m = Restrict(r.vHPC.coordinated,NREM.baseline);
                    bias.vHPC.mean.coordinated.baseline = [bias.vHPC.mean.coordinated.baseline ; mean(m(:,2))]; clear m
                    m = Restrict(r.vHPC.coordinated,NREM.reward);
                    bias.vHPC.mean.coordinated.reward = [bias.vHPC.mean.coordinated.reward ; mean(m(:,2))]; clear m
                    m = Restrict(r.vHPC.coordinated,NREM.aversive);
                    bias.vHPC.mean.coordinated.aversive = [bias.vHPC.mean.coordinated.aversive ; mean(m(:,2))]; clear m
                end
                
                %Ripples dHPC
                if and(and(not(isempty(Restrict(r.dHPC.ripple,NREM.baseline))) , not(isempty(Restrict(r.dHPC.ripple,NREM.reward)))) ,  not(isempty(Restrict(r.dHPC.ripple,NREM.aversive))))
                    m = Restrict(r.dHPC.ripple,NREM.baseline);
                    bias.dHPC.mean.ripple.baseline = [ bias.dHPC.mean.ripple.baseline ; mean(m(:,2))]; clear m
                    m = Restrict(r.dHPC.ripple,NREM.reward);
                    bias.dHPC.mean.ripple.reward = [bias.dHPC.mean.ripple.reward ; mean(m(:,2))]; clear m
                    m = Restrict(r.dHPC.ripple,NREM.aversive);
                    bias.dHPC.mean.ripple.aversive = [bias.dHPC.mean.ripple.aversive ; mean(m(:,2))]; clear m
                end
                
                %Ripples vHPC
                if and(and(not(isempty(Restrict(r.vHPC.ripple,NREM.baseline))) , not(isempty(Restrict(r.vHPC.ripple,NREM.reward)))) ,  not(isempty(Restrict(r.vHPC.ripple,NREM.aversive))))
                    m = Restrict(r.vHPC.ripple,NREM.baseline);
                    bias.vHPC.mean.ripple.baseline = [bias.vHPC.mean.ripple.baseline ; mean(m(:,2))]; clear m
                    m = Restrict(r.vHPC.ripple,NREM.reward);
                    bias.vHPC.mean.ripple.reward = [bias.vHPC.mean.ripple.reward ; mean(m(:,2))]; clear m
                    m = Restrict(r.vHPC.ripple,NREM.aversive);
                    bias.vHPC.mean.ripple.aversive =  [bias.vHPC.mean.ripple.aversive ; mean(m(:,2))]; clear m
                end
            end
        else
            disp('*** NOT ENOUGH PCs TO DO THE ANALYSIS ***')
        end
        clear specificity_dHPC_A specificity_dHPC_R specificity_vHPC_A specificity_vHPC_R
        clear field_dHPC_A field_dHPC_R field_vHPC_A field_vHPC_R
        clear peak_dHPC_A peak_dHPC_R peak_vHPC_A peak_vHPC_R
        clear maps_dHPC_A maps_dHPC_R maps_vHPC_A maps_vHPC_R
        clear aversiveTS aversiveTS_run rewardTS rewardTS_run baselineTS
        clear behavior camara Cell_type_classification cellulartype cluster
        clear group_dHPC group_vHPC K Kinfo leftvalve rightvalve movement R
        clear Rewards_filt Shocks_filt spks_dHPC spks_vHPC r rZ
        clear cluster cluster_dHPC cluster_vHPC coordinated coordinatedV coordinatedV_refined
        clear camaraA count dX dX_int dY dY_int gc grps i ii MUA segments tmp WAKE REM NREM
        clear ripple_bursts ripplesD ripplesV x w position_shocks posx posy ejeX ejeY PC replay Replay limits
        disp(['-- Finished folder # ---' , num2str(t) , ' from rat #' , num2str(tt)])
        disp('  ')
    end
    disp(['-------- Finished rat#' , num2str(tt) , ' --------'])
    disp(' ')
end


%% Plot section
grps = [ones(length(bias.dHPC.sws.baseline),1) ; ones(length(bias.dHPC.sws.reward),1)*2 ; ones(length(bias.dHPC.sws.aversive),1)*3];
x = [bias.dHPC.sws.baseline(:,2) ; bias.dHPC.sws.reward(:,2) ; bias.dHPC.sws.aversive(:,2)];
figure,subplot(121),boxplot(x,grps),ylim([-1 1])
title('dHPC')
grps = [ones(length(bias.vHPC.sws.baseline),1) ; ones(length(bias.vHPC.sws.reward),1)*2 ; ones(length(bias.vHPC.sws.aversive),1)*3];
x = [bias.vHPC.sws.baseline(:,2) ; bias.vHPC.sws.reward(:,2) ; bias.vHPC.sws.aversive(:,2)];
subplot(122),boxplot(x,grps),ylim([-1 1])
title('vHPC')
sgtitle('During SWS')

grps = [ones(length(bias.dHPC.ripple.baseline),1) ; ones(length(bias.dHPC.ripple.reward),1)*2 ; ones(length(bias.dHPC.ripple.aversive),1)*3];
x = [bias.dHPC.ripple.baseline(:,2) ; bias.dHPC.ripple.reward(:,2) ; bias.dHPC.ripple.aversive(:,2)];
figure,subplot(121),boxplot(x,grps),ylim([-1 1])
title('dHPC')
grps = [ones(length(bias.vHPC.ripple.baseline),1) ; ones(length(bias.vHPC.ripple.reward),1)*2 ; ones(length(bias.vHPC.ripple.aversive),1)*3];
x = [bias.vHPC.ripple.baseline(:,2) ; bias.vHPC.ripple.reward(:,2) ; bias.vHPC.ripple.aversive(:,2)];
subplot(122),boxplot(x,grps),ylim([-1 1])
title('vHPC')
sgtitle('During Ripples')

grps = [ones(length(bias.dHPC.coordinated.baseline),1) ; ones(length(bias.dHPC.coordinated.reward),1)*2 ; ones(length(bias.dHPC.coordinated.aversive),1)*3];
x = [bias.dHPC.coordinated.baseline(:,2) ; bias.dHPC.coordinated.reward(:,2) ; bias.dHPC.coordinated.aversive(:,2)];
figure,subplot(121),boxplot(x,grps),ylim([-1 1])
title('dHPC')
grps = [ones(length(bias.vHPC.coordinated.baseline),1) ; ones(length(bias.vHPC.coordinated.reward),1)*2 ; ones(length(bias.vHPC.coordinated.aversive),1)*3];
x = [bias.vHPC.coordinated.baseline(:,2) ; bias.vHPC.coordinated.reward(:,2) ; bias.vHPC.coordinated.aversive(:,2)];
subplot(122),boxplot(x,grps),ylim([-1 1])
title('vHPC')
sgtitle('During Coordinated Ripples')

grps = [ones(length(bias.dHPC.burst.baseline),1) ; ones(length(bias.dHPC.burst.reward),1)*2 ; ones(length(bias.dHPC.burst.aversive),1)*3];
x = [bias.dHPC.burst.baseline(:,2) ; bias.dHPC.burst.reward(:,2) ; bias.dHPC.burst.aversive(:,2)];
figure,subplot(121),boxplot(x,grps),ylim([-1 1])
title('dHPC')
grps = [ones(length(bias.vHPC.burst.baseline),1) ; ones(length(bias.vHPC.burst.reward),1)*2 ; ones(length(bias.vHPC.burst.aversive),1)*3];
x = [bias.vHPC.burst.baseline(:,2) ; bias.vHPC.burst.reward(:,2) ; bias.vHPC.burst.aversive(:,2)];
subplot(122),boxplot(x,grps),ylim([-1 1])
title('vHPC')
sgtitle('During Coordinated Ripple bursts')