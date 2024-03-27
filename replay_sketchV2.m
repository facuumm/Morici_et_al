clear
clc
close all

%% Parameters
path = {'E:\Rat126\Ephys\in_Pyr';'E:\Rat103\usable';'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat132\recordings\in_pyr';'E:\Rat165\in_pyr\'};%List of folders from the path

% for SU
criteria_fr = 0; %criteria to include or not a SU into the analysis
criteria_n = 3; % minimal number of neurons from each structure
criteria_type = 0; %criteria for celltype (0:pyr, 1:int, 2:all)

binSize = 0.005; % bin size for replay events detection

% Behavior
minimal_speed = 7; % minimal speed to detect active periods
minimal_speed_time = 2; % minimal time to detect active periods

Replay.candidate.dHPC.baseline=[]; Replay.candidate.dHPC.reward=[]; Replay.candidate.dHPC.aversive=[];
Replay.candidate.vHPC.baseline=[]; Replay.candidate.vHPC.reward=[]; Replay.candidate.vHPC.aversive=[];

Replay.selected.dHPC.baseline=[]; Replay.selected.dHPC.reward=[]; Replay.selected.dHPC.aversive=[];
Replay.selected.vHPC.baseline=[]; Replay.selected.vHPC.reward=[]; Replay.selected.vHPC.aversive=[];

Replay.Shocks.vHPC.coordinated = []; Replay.Shocks.vHPC.uncoordinated = [];
Replay.Shocks.dHPC.coordinated = []; Replay.Shocks.dHPC.uncoordinated = [];

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
        
        % Definition of the place where the shocks were admninistrated
        position_shocks = [];
        for i = 1 : length(Shocks_filt)
            [~ , ii] = min(abs(behavior.pos.aversive(:,1)-Shocks_filt(i)));
            position_shocks = [position_shocks ; behavior.pos.aversive(ii,2)];
            clear ii
        end
        clear i
        template = linspace(min(behavior.pos.aversive(:,2)) , max(behavior.pos.aversive(:,2)) , 60); % define bins to estimate what position is more reactivated
        [~,position_shocks] = histc(position_shocks,template);
        
        % Merging events close in space (<5 spatial bins)
        position_shocks = sort(position_shocks);
        D = find(diff(position_shocks)<=5);
        tmp = [];
        tmp1 = [];
        for i = 1 : size(position_shocks,1)
            if any(D == i) || any(D +1 == i)
                tmp = [tmp ; position_shocks(i)];
            else
                tmp1 = [tmp1 ; position_shocks(i)];
            end
        end
        
        if not(isempty(tmp))
            tmp3 = [];
            suma = tmp(1);
            c = 1;
            for i = 2 : size(tmp,1)
                if tmp(i) - tmp(i-1) <=5
                    suma = suma+tmp(i);
                    c = c+1;
                    if i == size(tmp,1)
                        tmp3 = [tmp3 ; suma/c];
                        suma = tmp(i);
                        c = 1;
                    end
                else
                    tmp3 = [tmp3 ; suma/c];
                    suma = tmp(i);
                    c = 1;
                end
            end
            
            position_shocks = sort([tmp1 ; tmp3]); clear tmp tmp1 tmp3 D i
        else
            clear tmp tmp1 tmp3 D i
        end
%         position_shocks = position_shocks-min(behavior.pos.aversive(:,2));
%         position_shocks = position_shocks./max(behavior.pos.aversive(:,2));

        
        %% load sleep states
        disp('Uploading sleep scoring')
        x = dir([cd,'\*-states.mat']);    states = load([cd,'\',x.name]);    states = states.states;
        REM.all = ToIntervals(states==5);    NREM.all = ToIntervals(states==3);    WAKE.all = ToIntervals(states==1);
        %         REM.all(REM.all(:,2)-REM.all(:,1)<60,:) = [];
        clear x states
        
        %         NREM.all(NREM.all(:,2)-NREM.all(:,1)<60*time_criteria,:)=[];
        NREM.baseline = Restrict(NREM.all,baselineTS./1000);
        NREM.aversive = Restrict(NREM.all,aversiveTS./1000);
        NREM.reward = Restrict(NREM.all,rewardTS./1000);
        
        REM.baseline = Restrict(REM.all,baselineTS./1000);
        REM.aversive = Restrict(REM.all,aversiveTS./1000);
        REM.reward = Restrict(REM.all,rewardTS./1000);
        
        WAKE.baseline = Restrict(WAKE.all,baselineTS./1000);
        WAKE.aversive = Restrict(WAKE.all,aversiveTS./1000);
        WAKE.reward = Restrict(WAKE.all,rewardTS./1000);
        
        
        %% Load ripples
        if exist('ripplesD_customized2.csv')
            ripplesD = table2array(readtable('ripplesD_customized2.csv'));
            RD = true;
        else
            RD = false;
        end
        
        if exist('ripplesV_customized2.csv')
            ripplesV = table2array(readtable('ripplesV_customized2.csv'));
            RV = true;
        else
            RV = false;
        end
        
            
        if isfile('coordinated_ripple_bursts.mat')
            load('coordinated_ripple_bursts.mat')
        % Detection of coordinated ripples
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
            [C,IA,IC] = unique(coordinatedV(:,1));
            coordinatedV  = coordinatedV(IA,:); clear C IA IC
            % Detection of uncoordinated ripples
            uncoordinated = ripplesD(~ismember(ripplesD(:,1),coordinated(:,1)),:);
            uncoordinatedV = ripplesV(~ismember(ripplesV(:,1),unique(coordinatedV(:,1))),:);            
            
            
            % Spikes
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
            
            
            %% Firing maps related to each event
            if and(isfile('dHPC_pc.mat') , isfile('vHPC_pc.mat'))
                disp('Loading maps from dHPC and vHPC')
                load('dHPC_pc.mat')              
                clusters.pc.dHPC = [];
                maps.aversive.dHPC = [];
                maps.reward.dHPC = [];
                for i = 1 : size(dHPC,2)
                    clusters.pc.dHPC = [clusters.pc.dHPC ; dHPC{i}.id];
                    maps.aversive.dHPC = [maps.aversive.dHPC ; dHPC{i}.frMap_ave];
                    maps.reward.dHPC = [maps.reward.dHPC ; dHPC{i}.frMap_rew];
                end
                
                load('vHPC_pc.mat')
                clusters.pc.vHPC = [];
                maps.aversive.vHPC = [];
                maps.reward.vHPC = [];
                for i = 1 : size(vHPC,2)
                    clusters.pc.vHPC = [clusters.pc.vHPC ; vHPC{i}.id];
                    maps.aversive.vHPC = [maps.aversive.vHPC ; vHPC{i}.frMap_ave];
                    maps.reward.vHPC = [maps.reward.vHPC ; vHPC{i}.frMap_rew];
                end
                
                %% Detection of putative replay events
                if and(size(dHPC,2)>criteria_n , size(vHPC,2)>criteria_n)
                    % --- dHPC replay detection ---
                    disp('Detetion of dorsal putarive-replay events')
                    limits = [0 segments.Var1(end)/1000];
                    events = [];
                    
                    [MUA.dHPC , bins] = spike_train_construction_counts(spks_dHPC, clusters.pc.dHPC, binSize, limits);
                    [replay.dHPC] = FindReplay([bins' , MUA.dHPC], [0.5 3] , [0.04 0.05 0.8] , 10 , []);

                    
                    % --- vHPC replay detection ---
                    disp('Detetion of ventral putarive-replay events')
                    [MUA.vHPC , bins] = spike_train_construction_counts(spks_vHPC, clusters.pc.vHPC, binSize, limits);
                    [replay.vHPC] = FindReplay([bins' , MUA.vHPC], [0.5 3] , [0.04 0.05 0.8] , 10 , []);

%                     
                    %% Decoding
                    % vHPC aversive
                    disp('Decoding using vHPC aversive in coordinated')
%                     limits = [0 : 0.02 : segments.Var1(end)/1000];
                    R = replay.vHPC;
                    r = Restrict(ripplesV(:,2),NREM.aversive);
                    R = Restrict(R,[r-0.1 r+0.1]);
%                     R = Restrict(replay.vHPC,R);
                    Prob.vHPC = [];
                    for i = 1 : size(R,1)
%                         limits = [R(i)-0.05 : 0.02 : R(i)+0.05];
%                         tmp = [];
%                         for ii = 1 : size(limits,2)-1
%                             nSpks = count_spks(spks_vHPC, clusters.pc.vHPC, limits(ii), limits(ii+1));
%                             probability = bayesian_replay(maps.aversive.vHPC, nSpks, limits(ii), limits(ii+1));
%                             p = probability';
%                             p = p;
%                             tmp = [tmp , p]; clear p probability
%                         end
                        nSpks = count_spks(spks_vHPC, clusters.pc.vHPC, R(i,1), R(i,3));
                        [R H P] = bayesian_replay(spks_vHPC, clusters.pc.vHPC , maps.aversive.vHPC, [R(i,1), R(i,3)] , 0.02 , 15);
                        
                        
                        probability = bayesian_decoder(maps.aversive.vHPC, nSpks, R(i,1), R(i,3));
%                         probability = probability - min(probability);
%                         probability = probability./max(probability);
                        [ii iii] = findpeaks(probability,[1 : 1 :60],'MinPeakHeight',nanmean(probability),'NPeaks',1);
%                         figure,imagesc(tmp), colorbar
                        Prob.vHPC = [Prob.vHPC ; ii' iii']; clear nSpks probability ii iii
                    end
                    
                    C = [];
                    for i = 1 : size(position_shocks,1)
                        C = [C ; abs(position_shocks(i) - Prob.vHPC(:,2))];
                    end
                    Replay.Shocks.vHPC.coordinated = [Replay.Shocks.vHPC.coordinated ; (C)]; clear C

                    % Uncooridnated
                    disp('Decoding using vHPC aversive in uncooridnated')
                    R = replay.vHPC;
                    r = Restrict(uncoordinatedV(:,2),NREM.aversive);
                    R = Restrict(R,[r-0.1 r+0.1]);
                    Prob.vHPC = [];
                    for i = 1 : size(R,1)
                        nSpks = count_spks(spks_vHPC, clusters.pc.vHPC, R(i,1), R(i,3));
                        probability = bayesian_decoder(maps.aversive.vHPC, nSpks, R(i,1), R(i,3));
                        [ii iii] = findpeaks(probability,[1 : 1 :60],'MinPeakHeight',nanmean(probability),'NPeaks',1);
                        Prob.vHPC = [Prob.vHPC ; ii' iii']; clear nSpks probability ii iii
                    end
                    
                    C = [];
                    for i = 1 : size(position_shocks,1)
                        C = [C ; abs(position_shocks(i) - Prob.vHPC(:,2))];
                    end
                    Replay.Shocks.vHPC.uncoordinated = [Replay.Shocks.vHPC.uncoordinated ; (C)]; clear C
                    
                    
                    
                    % dHPC aversive
                    disp('Decoding using dHPC aversive in coordinated')
%                     limits = [0 : 0.02 : segments.Var1(end)/1000];
                    R = replay.dHPC;
                    r = Restrict(coordinatedV(:,2),NREM.aversive);
                    R = Restrict(R,[r-0.1 r+0.1]);
%                     R = Restrict(replay.vHPC,R);
                    Prob.dHPC = [];
                    for i = 1 : size(R,1)
%                         limits = [R(i)-0.05 : 0.02 : R(i)+0.05];
%                         tmp = [];
%                         for ii = 1 : size(limits,2)-1
%                             nSpks = count_spks(spks_vHPC, clusters.pc.vHPC, limits(ii), limits(ii+1));
%                             probability = bayesian_replay(maps.aversive.vHPC, nSpks, limits(ii), limits(ii+1));
%                             p = probability';
%                             p = p;
%                             tmp = [tmp , p]; clear p probability
%                         end
                        nSpks = count_spks(spks_dHPC, clusters.pc.dHPC, R(i,1), R(i,3));
                        probability = bayesian_decoder(maps.aversive.dHPC, nSpks, R(i,1), R(i,3));
%                         probability = probability - min(probability);
%                         probability = probability./max(probability);
                        [ii iii] = findpeaks(probability,[1 : 1 :60],'MinPeakHeight',nanmean(probability),'NPeaks',1);
%                         figure,imagesc(tmp), colorbar
                        Prob.dHPC = [Prob.dHPC ; ii' iii']; clear nSpks probability ii iii
                    end
                    
                    C = [];
                    for i = 1 : size(position_shocks,1)
                        C = [C ; abs(position_shocks(i) - Prob.dHPC(:,2))];
                    end
                    Replay.Shocks.dHPC.coordinated = [Replay.Shocks.dHPC.coordinated ; (C)]; clear C

                    % Uncooridnated
                    disp('Decoding using dHPC aversive in uncooridnated')
                    R = replay.dHPC;
                    r = Restrict(uncoordinated(:,2),NREM.aversive);
                    R = Restrict(R,[r-0.1 r+0.1]);
                    Prob.vHPC = [];
                    for i = 1 : size(R,1)
                        nSpks = count_spks(spks_dHPC, clusters.pc.dHPC, R(i,1), R(i,3));
                        probability = bayesian_decoder(maps.aversive.dHPC, nSpks, R(i,1), R(i,3));
                        [ii iii] = findpeaks(probability,[1 : 1 :60],'MinPeakHeight',nanmean(probability),'NPeaks',1);
                        Prob.vHPC = [Prob.vHPC ; ii' iii']; clear nSpks probability ii iii
                    end
                    
                    C = [];
                    for i = 1 : size(position_shocks,1)
                        C = [C ; abs(position_shocks(i) - Prob.vHPC(:,2))];
                    end
                    Replay.Shocks.dHPC.uncoordinated = [Replay.Shocks.dHPC.uncoordinated ; (C)]; clear C                    
                    
                    
%                     % dHPC aversive
% %                     limits = [0 : 0.02 : segments.Var1(end)/1000];
%                     R = replay.dHPC(:,3);
%                     r = Restrict(coordinated(:,2),NREM.aversive);
%                     R = Restrict(R,[r-0.1 r+0.1]);
% %                     R = Restrict(replay.vHPC,R);
%                     Prob.dHPC = [];
%                     for i = 1 : size(R,1)
%                         limits = [R(i)-0.05 : 0.02 : R(i)+0.05];
%                         tmp = [];
%                         for ii = 1 : size(limits,2)-1
%                             nSpks = count_spks(spks_dHPC, clusters.pc.dHPC, limits(ii), limits(ii+1));
%                             probability = bayesian_replay(maps.aversive.dHPC, nSpks, limits(ii), limits(ii+1));
%                             p = probability';
%                             tmp = [tmp , p]; clear p probability
%                         end
%                         tmp = tmp - min(min(tmp));
%                         tmp = tmp./max(max(tmp));
%                         figure,imagesc(tmp),colorbar
% %                         Prob.vHPC = [Prob.vHPC , probability']; clear nSpks probability
%                         i
%                     end
%                     
%                     
%                   % vHPC Reward
% %                     limits = [0 : 0.02 : segments.Var1(end)/1000];
%                     R = ((replay.vHPC(:,2)-replay.vHPC(:,1))/2)+replay.vHPC(:,1);
%                     R = Restrict(R,NREM.reward);
% %                     R = Restrict(replay.vHPC,R);
%                     Prob.vHPC = [];
%                     for i = 1 : size(R,1)
%                         limits = [R(i)-0.2 : 0.02 : R(i)+0.2];
%                         tmp = [];
%                         for ii = 1 : size(limits,2)-1
%                             nSpks = count_spks(spks_vHPC, clusters.pc.vHPC, limits(ii), limits(ii+1));
%                             probability = bayesian_replay(maps.reward.vHPC, nSpks, limits(ii), limits(ii+1));
%                             p = probability';
%                             p = p;
%                             tmp = [tmp , p]; clear p probability
%                         end
%                         figure,imagesc(tmp), colorbar
% %                         Prob.vHPC = [Prob.vHPC , probability']; clear nSpks probability
%                         i
%                     end
%                     
%                     % dHPC reward
% %                     limits = [0 : 0.02 : segments.Var1(end)/1000];
%                     R = ((replay.dHPC(:,2)-replay.dHPC(:,1))/2)+replay.dHPC(:,1);
%                     R = Restrict(R,NREM.reward);
%                     Prob.dHPC = [];
%                     for i = 1 : size(R,1)
%                         limits = [R(i)-0.2 : 0.02 : R(i)+0.2];
%                         tmp = [];
%                         for ii = 1 : size(limits,2)-1
%                             nSpks = count_spks(spks_dHPC, clusters.pc.dHPC, limits(ii), limits(ii+1));
%                             probability = bayesian_replay(maps.reward.dHPC, nSpks, limits(ii), limits(ii+1));
%                             p = probability';
%                             tmp = [tmp , p]; clear p probability
%                         end
%                         figure,imagesc(tmp),colorbar
% %                         Prob.vHPC = [Prob.vHPC , probability']; clear nSpks probability
%                         i
%                     end
                    
                    
                end
                clear dHPC vHPC
            end
            
        end
        disp(' ')
        clear specificity_dHPC_A specificity_dHPC_R specificity_vHPC_A specificity_vHPC_R
        clear field_dHPC_A field_dHPC_R field_vHPC_A field_vHPC_R
        clear peak_dHPC_A peak_dHPC_R peak_vHPC_A peak_vHPC_R
        clear maps_dHPC_A maps_dHPC_R maps_vHPC_A maps_vHPC_R
        clear aversiveTS aversiveTS_run rewardTS rewardTS_run baselineTS
        clear behavior camara Cell_type_classification cellulartype cluster
        clear group_dHPC group_vHPC K Kinfo leftvalve rightvalve movement R
        clear Rewards_filt Shocks_filt spks_dHPC spks_vHPC
        clear NREM REM WAKE coordinated coordinatedV coordinatedV_refined uncoordinated uncoordinatedV
        clear posx posy numberD numberV position_shocks segments
        clear spks template limits camaraA bins  replay ripplesD ripplesV MUA
    end
    disp(['-- Finishing analysis from rat #',num2str(tt) , ' --'])
    disp('  ')
end


%% Plot
figure,
subplot(221),histogram(Replay.Shocks.vHPC.coordinated,20,'Normalization','probability')
subplot(222),histogram(Replay.Shocks.vHPC.uncoordinated,20,'Normalization','probability')
subplot(223),histogram(Replay.Shocks.dHPC.coordinated,20,'Normalization','probability')
subplot(224),histogram(Replay.Shocks.dHPC.uncoordinated,20,'Normalization','probability')


