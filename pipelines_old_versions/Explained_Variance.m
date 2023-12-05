clear
clc
close all

%% Parameters
path = {'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat103\usable'};%List of folders from the path

%Sleep
time_criteria = 3600; %time criteria to define the maximal time of sleep to include

% for SU
criteria_fr = 0; %criteria to include or not a SU into the analysis
criteria_n = 6; % minimal number of neurons from each structure
criteria_type = 0; %criteria for celltype (0:pyr, 1:int, 2:all)
pval = 0.01; % p value to define if SU are ripple modulated
ss = 2; %smooth level of CCG
n_SU_V = 0;
n_SU_D = 0;
FR_B_V = []; FR_R_V = [];  FR_A_V = []; % Firing Rate during NREM
FR_B_D = []; FR_R_D = [];  FR_A_D = []; % Firing Rate during NREM
poisson_dHPC_split = []; poisson_vHPC_split = []; %poisson results split by conditions
Xedges = 60; % number of spatial bins for firing curve construction

%For EV and REV
binSize = 0.05;
EV_A = []; REV_A = [];
EV_R = []; REV_R = [];

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
        % Load digitalin:mat
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
            behavior.speed.reward = LinearVelocity(pos);
            behavior.quiet.reward = QuietPeriods( behavior.speed.reward , minimal_speed , minimal_speed_time);
            clear pos camaraR posx posy
            clear ejeX ejeY dX dY
            
            load('laps2.mat','posx','posy');
            [camaraA,~] = find((camara(:,1)-aversiveTS_run(1)/1000)>0,1,'first'); %TimeStamp of the begining of aversive
            pos = [camara(camaraA : camaraA+length(posx)-1),posx,posy];
            
            %interpolation of dropped frames
            ejeX = pos(~isnan(pos(:,2)),1); dX = pos(~isnan(pos(:,2)),2); dX_int = interp1(ejeX , dX , pos(:,1));
            ejeY = pos(~isnan(pos(:,3)),1); dY = pos(~isnan(pos(:,3)),3); dY_int = interp1(ejeY , dY , pos(:,1));
            pos(:,2) =dX_int; pos(:,3) =dY_int;%saving corrected pos
            
            behavior.pos.aversive = [pos];
            behavior.speed.aversive = LinearVelocity(pos);
            behavior.quiet.aversive = QuietPeriods(behavior.speed.aversive , minimal_speed , minimal_speed_time);
            clear pos camaraA posx posy
            clear ejeX ejeY dX dY
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
            behavior.speed.reward = LinearVelocity(pos,15);
            behavior.quiet.reward = QuietPeriods( behavior.speed.reward , minimal_speed , minimal_speed_time);
            clear pos camaraR posx posy
            clear ejeX ejeY dX dY
            
            load('laps1.mat','posx','posy');
            [camaraA,~] = find((camara(:,1)-aversiveTS_run(1)/1000)>0,1,'first'); %TimeStamp of the begining of aversive
            pos = [camara(camaraA : camaraA+length(posx)-1),posx,posy];
            
            %interpolation of dropped frames
            ejeX = pos(~isnan(pos(:,2)),1); dX = pos(~isnan(pos(:,2)),2); dX_int = interp1(ejeX , dX , pos(:,1));
            ejeY = pos(~isnan(pos(:,3)),1); dY = pos(~isnan(pos(:,3)),3); dY_int = interp1(ejeY , dY , pos(:,1));
            pos(:,2) =dX_int; pos(:,3) =dY_int; %saving corrected pos            
            
            behavior.pos.aversive = [pos];
            behavior.speed.aversive = LinearVelocity(pos,15);
            behavior.quiet.aversive = QuietPeriods(behavior.speed.aversive , minimal_speed , minimal_speed_time);
            clear pos camaraA posx posy
            clear ejeX ejeY dX dY
        end
        
        % Generation of no-movements periods
        start = behavior.speed.reward(1,1);   stop = behavior.speed.reward(end,1);
        movement.reward = InvertIntervals(behavior.quiet.reward , start , stop); %keep only those higher than criteria
%         [movement.reward , ~] =  ExcludeIntervals(movement.reward ,[Rewards_filt Rewards_filt+4]);
        clear tmp start stop
        
        start = behavior.speed.aversive(1,1);   stop = behavior.speed.aversive(end,1);
        movement.aversive = InvertIntervals(behavior.quiet.aversive , start , stop);%keep only those higher than criteria
%         [movement.aversive , ~] =  ExcludeIntervals(movement.aversive ,[Shocks_filt-0.5 Shocks_filt+2]);
        clear tmp start stop
        
        %% Sleep
        
        x = dir([cd,'\*-states.mat']);    states = load([cd,'\',x.name]);    states = states.states;
        
        REM = ToIntervals(states==5);    NREM = ToIntervals(states==3);    WAKE = ToIntervals(states==1);
        clear x states
        
        NREM_B = NREM(NREM(:,1)<baselineTS(1,2)/1000,:);
        NREM_R = NREM(and(NREM(:,1)>rewardTS(1,1)/1000 , NREM(:,1)<rewardTS(1,2)/1000),:);
        NREM_A = NREM(and(NREM(:,1)>aversiveTS(1,1)/1000 ,NREM(:,1)<aversiveTS(1,2)/1000),:);

        REM_B = REM(REM(:,1)<baselineTS(1,2)/1000,:);
        REM_R = REM(and(REM(:,1)>rewardTS(1,1)/1000 , REM(:,1)<rewardTS(1,2)/1000),:);
        REM_A = REM(and(REM(:,1)>aversiveTS(1,1)/1000 ,REM(:,1)<aversiveTS(1,2)/1000),:);
        
        ripplesD = table2array(readtable('ripplesD_customized2.csv'));
        ripplesV = table2array(readtable('ripplesV_customized2.csv'));
        
        % Coordination
        coordinated = [];
        coordinatedV = [];
        coordinatedV_refined = [];
        centroide = [];        centroide = [];

        for i = 1:length(ripplesD)
            r = ripplesD(i,:);
            tmp = sum(and(ripplesV(:,2)>= r(1,2)-0.1, ripplesV(:,2)<= r(1,2)+0.1));
            if tmp>0
                z = ripplesV(and(ripplesV(:,2)>= r(1,2)-0.1, ripplesV(:,2)<= r(1,2)+0.1),:);
                coordinatedV = [coordinatedV ; z];
                [p,indice] = min(abs(r(2)-z(:,2)));
                coordinatedV_refined = [coordinatedV_refined ; z(indice,:)];
                coordinated = [coordinated ; r];
                
                z = z(indice,:);
                
                if z(2) > r(2)
                    c = ((z(2) - r(2))/2)+r(2);
                    centroide = [centroide ; r(1) c z(3)];
                    clear c
                else
                    c = ((r(2) - z(2))/2)+z(2);
                    centroide = [centroide ; z(1) c r(3)];
                    clear c
                end
                                
                clear tmp2 tmp1 p indice z
            end
            clear r
        end
        clear x tmp i
        
        coordinatedB = coordinated(coordinated(:,2)<baselineTS(1,2)/1000,:);
        coordinatedR = coordinated(and(coordinated(:,2)>rewardTS(1,1)/1000 , coordinated(:,2)<rewardTS(1,2)/1000),:);
        coordinatedA = coordinated(and(coordinated(:,2)>aversiveTS(1,1)/1000 , coordinated(:,2)<aversiveTS(1,2)/1000),:);
      
        coordinatedB_V = coordinatedV_refined(coordinatedV_refined(:,2)<baselineTS(1,2)/1000,:);
        coordinatedR_V = coordinatedV_refined(and(coordinatedV_refined(:,2)>rewardTS(1,1)/1000 , coordinatedV_refined(:,2)<rewardTS(1,2)/1000),:);
        coordinatedA_V = coordinatedV_refined(and(coordinatedV_refined(:,2)>aversiveTS(1,1)/1000 , coordinatedV_refined(:,2)<aversiveTS(1,2)/1000),:);
              
        
        centroideB = centroide(centroide(:,2)<baselineTS(1,2)/1000,:);
%         centroideB = Restrict(centroideB,[NREM_B(1,1) NREM_B(1,1)+600]);
        centroideR = centroide(and(centroide(:,2)>rewardTS(1,1)/1000 , centroide(:,2)<rewardTS(1,2)/1000),:);
%         centroideR = Restrict(centroideR,[NREM_R(1,1) NREM_R(1,1)+600]);
        centroideA = centroide(and(centroide(:,2)>aversiveTS(1,1)/1000 , centroide(:,2)<aversiveTS(1,2)/1000),:);
%         centroideA = Restrict(centroideA,[NREM_A(1,1) NREM_A(1,1)+600]);
        
%         ripplesD_B = ripplesD(ripplesD(:,2)<baselineTS(1,2)/1000,:);
%         ripplesD_R = ripplesD(and(ripplesD(:,2)>rewardTS(1,1)/1000 , ripplesD(:,2)<rewardTS(1,2)/1000),:);
%         ripplesD_A = ripplesD(and(ripplesD(:,2)>aversiveTS(1,1)/1000 , ripplesD(:,2)<aversiveTS(1,2)/1000),:);
% 
%         ripplesV_B = ripplesV(ripplesV(:,2)<baselineTS(1,2)/1000,:);
%         ripplesV_R = ripplesV(and(ripplesV(:,2)>rewardTS(1,1)/1000 , ripplesV(:,2)<rewardTS(1,2)/1000),:);
%         ripplesV_A = ripplesV(and(ripplesV(:,2)>aversiveTS(1,1)/1000 , ripplesV(:,2)<aversiveTS(1,2)/1000),:);
%         
%         
%         % Load coordinated ripple bursts
%         load('coordinated_ripple_bursts.mat')
%         Coordinated_Bursts_B = coordinated_ripple_bursts(coordinated_ripple_bursts(:,2)<baselineTS(1,2)/1000,:);
%         Coordinated_Bursts_R = coordinated_ripple_bursts(and(coordinated_ripple_bursts(:,2)>rewardTS(1,1)/1000 , coordinated_ripple_bursts(:,2)<rewardTS(1,2)/1000),:);
%         Coordinated_Bursts_A = coordinated_ripple_bursts(and(coordinated_ripple_bursts(:,2)>aversiveTS(1,1)/1000 , coordinated_ripple_bursts(:,2)<aversiveTS(1,2)/1000),:);
%         
%         Coordinated_Bursts_B = [Coordinated_Bursts_B(:,2)-0.2 , Coordinated_Bursts_B(:,2) , Coordinated_Bursts_B(:,2)+0.2];
%         Coordinated_Bursts_R = [Coordinated_Bursts_R(:,2)-0.2 , Coordinated_Bursts_R(:,2) , Coordinated_Bursts_R(:,2)+0.2];
%         Coordinated_Bursts_A = [Coordinated_Bursts_A(:,2)-0.2 , Coordinated_Bursts_A(:,2) , Coordinated_Bursts_A(:,2)+0.2];
%         Coordinated_Bursts_B = Coordinated_Bursts_B(Coordinated_Bursts_B(:,2) < Coordinated_Bursts_B(1,1)+1200 , :);
%         Coordinated_Bursts_R = Coordinated_Bursts_R(Coordinated_Bursts_R(:,2) < Coordinated_Bursts_R(1,1)+1200 , :);
%         Coordinated_Bursts_A = Coordinated_Bursts_A(Coordinated_Bursts_A(:,2) < Coordinated_Bursts_A(1,1)+1200 , :);

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

        %% Firing maps calculation
        maps_dHPC_A = [];        maps_dHPC_R = [];
        PC.dHPC.aversive=[]; PC.dHPC.reward=[];
        for ii=1:length(group_dHPC)
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
                    
                    if not(isempty(spks_tmp))
                        %Firing curve construction
                        [curve , stats] = FiringCurve(pos_tmp , spks_tmp , 'smooth' , 2 , 'nBins' , Xedges , 'minSize' , 2 , 'minPeak' , 0.2);
                        maps_dHPC_A = [maps_dHPC_A ; curve.rate];
                        
%                         q = SkaggsRandom(spks_tmp, pos_tmp, [] , [] , Xedges);
                        
                        % Store of place-field location and cluster of PCs
                        if stats.specificity > 0.25
                            [ff,f] = max(curve.rate);
                            fff = [1/Xedges:1/Xedges:1];
                            PC.dHPC.aversive=[PC.dHPC.aversive ; cluster , ff , fff(f)];
                            clear f ff fff
                        end
                    end
                    clear s q curve OccMap Nspikes spks_tmp pos_tmp m curve stats sp curve1 q
                    clear curve1 OccMap Nspikes x

                    % --- Reward ---
                    spks_tmp = Restrict(spks , movement.reward);
                    pos_tmp = Restrict(behavior.pos.reward(:,1:2) , movement.reward);
                    %                     pos_tmp(or(pos_tmp(:,2)<30 , pos_tmp(:,2)>170),2) = nan; %Eliminate platform periods
                    pos_tmp(:,2) = pos_tmp(:,2)-min(pos_tmp(:,2)); pos_tmp(:,2) = pos_tmp(:,2)/max(pos_tmp(:,2)); %normalization of position
                    if not(isempty(spks_tmp))
                        [curve , stats] = FiringCurve(pos_tmp , spks_tmp , 'smooth' , 2 , 'nBins' , Xedges , 'minSize' , 2 , 'minPeak' , 0.2);
                        maps_dHPC_R = [maps_dHPC_R ; curve.rate];
                        
%                         q = SkaggsRandom(spks_tmp, pos_tmp, [] , [] , Xedges);
                        
                        % Store of place-field location and cluster of PCs
                        if stats.specificity > 0.25
                            [ff,f] = max(curve.rate);
                            fff = [1/Xedges:1/Xedges:1];
                            PC.dHPC.reward=[PC.dHPC.reward ; cluster , ff , fff(f)];
                            clear f ff fff
                        end
                    end
                    clear s q curve OccMap Nspikes spks_tmp pos_tmp m curve stats sp curve1 q 
                    clear curve1 OccMap Nspikes x
                end
            end
            clear celltype tmp b cluster
        end
        
        
        maps_vHPC_A = [];        maps_vHPC_R = [];
        PC.vHPC.aversive=[];     PC.vHPC.reward=[];
        for ii=1:length(group_vHPC)
            cluster = group_vHPC(ii,1);
            celltype = logical(cellulartype(cellulartype(:,1) == cluster,2));
            if celltype
                spks = spks_vHPC(spks_vHPC(:,1)==cluster,2);
                b = [movement.aversive ; movement.reward];
                tmp = Restrict(spks,b);
                m = length(tmp) / sum(b(:,2) - b(:,1));
                if ~isnan(m)
                    % --- Aversive ---
                    spks_tmp = Restrict(spks , movement.aversive);
                    pos_tmp = Restrict(behavior.pos.aversive(:,1:2) , movement.aversive);
                    pos_tmp(or(pos_tmp(:,2)<30 , pos_tmp(:,2)>150),2) = nan; %Eliminate platform periods
                    pos_tmp(:,2) = pos_tmp(:,2)-min(pos_tmp(:,2)); pos_tmp(:,2) = pos_tmp(:,2)/max(pos_tmp(:,2)); %normalization of position
                    
                    if not(isempty(spks_tmp))
                        %Firing curve construction
                        [curve , stats] = FiringCurve(pos_tmp , spks_tmp , 'smooth' , 2 , 'nBins' , Xedges , 'minSize' , 4 , 'minPeak' , 0.2);
                        maps_vHPC_A = [maps_vHPC_A ; curve.rate];
                        
%                         q = SkaggsRandom(spks_tmp, pos_tmp, [] , [], Xedges);
                        
                        % Store of place-field location and cluster of PCs
                        if stats.specificity > 0.25
                            [ff,f] = max(curve.rate);
                            fff = [1/Xedges:1/Xedges:1];
                            PC.vHPC.aversive=[PC.vHPC.aversive ; cluster , ff , fff(f)];
                            clear f ff fff
                        end
                    end
                    clear s q curve OccMap Nspikes spks_tmp pos_tmp m curve stats sp curve1 q
                    
                    % --- Reward ---
                    spks_tmp = Restrict(spks , movement.reward);
                    pos_tmp = Restrict(behavior.pos.reward(:,1:2) , movement.reward);
                    %                     pos_tmp(or(pos_tmp(:,2)<30 , pos_tmp(:,2)>150),2) = nan; %Eliminate platform periods
                    pos_tmp(:,2) = pos_tmp(:,2)-min(pos_tmp(:,2)); pos_tmp(:,2) = pos_tmp(:,2)/max(pos_tmp(:,2)); %normalization of position
                    
                    if not(isempty(spks_tmp))
                        %Firing curve construction
                        [curve , stats] = FiringCurve(pos_tmp , spks_tmp , 'smooth' , 2 , 'nBins' , Xedges , 'minSize' , 4 , 'minPeak' , 0.2);
                        maps_vHPC_R = [maps_vHPC_R ; curve.rate];
                        
%                         q = SkaggsRandom(spks_tmp, pos_tmp, [] , [] , Xedges);
                        
                        % Store of place-field location and cluster of PCs
                        if stats.specificity > 0.25
                            [ff,f] = max(curve.rate);
                            fff = [1/Xedges:1/Xedges:1];
                            PC.vHPC.reward=[PC.vHPC.reward ; cluster , ff , fff(f)];
                            clear f ff fff
                        end
                    end
                    clear s q curve OccMap Nspikes spks_tmp pos_tmp m curve stats sp curve1 q
                end
            end
            clear celltype tmp b cluster
        end
        
        %% Constructing Spiketrains
        freq = 1/binSize;
        limits = [0 segments.Var1(end)/1000];
        spiketrains_dHPC = [];
        spiketrains_vHPC = [];
        clusters.dHPC = [];        clusters.vHPC = [];
        
        for ii=1:length(group_dHPC)
            cluster = group_dHPC(ii,1);
            celltype = logical(cellulartype(cellulartype(:,1) == cluster,2));
            if celltype
                spks = spks_dHPC(spks_dHPC(:,1)==cluster,2);
                [tmp,bins]=binspikes(spks,freq,limits);
                spiketrains_dHPC = [spiketrains_dHPC , tmp./binSize];
                clusters.dHPC = [clusters.dHPC;cluster];
                clear spks tmp m1 m2 m3 m4
            end
            clear celltype cluster
        end
        
        for ii=1:length(group_vHPC)
            cluster = group_vHPC(ii,1);
            celltype = logical(cellulartype(cellulartype(:,1) == cluster,2));
            if celltype
                spks = spks_vHPC(spks_vHPC(:,1)==cluster,2);
                [tmp,bins]=binspikes(spks,freq,limits);
                
                spiketrains_vHPC = [spiketrains_vHPC , tmp./binSize];
                clusters.vHPC = [clusters.vHPC;cluster];
            end
            clear spks tmp m1 m2 m3 m4
            clear celltype cluster
        end
        clear  freq limits
        
        if and(size(spiketrains_vHPC,2) >= criteria_n,size(spiketrains_dHPC,2) >= criteria_n)
%           if and(and(length(PC.vHPC.aversive)>criteria_n , length(PC.vHPC.reward)>criteria_n),and(length(PC.dHPC.aversive)>criteria_n , length(PC.dHPC.reward)>criteria_n))
            %Restricting bins inside each condition
            %coordinated bursts
            %             is.coordinatedBurst.baseline = InIntervals(bins,[Coordinated_Bursts_B(:,1),Coordinated_Bursts_B(:,3)]);
            %             is.coordinatedBurst.reward = InIntervals(bins,[Coordinated_Bursts_R(:,1),Coordinated_Bursts_R(:,3)]);
            %             is.coordinatedBurst.aversive = InIntervals(bins,[Coordinated_Bursts_A(:,1),Coordinated_Bursts_A(:,3)]);
            %
            
            is.baseline.coordinated.dHPC = InIntervals(bins,[coordinatedB(:,2)-0.1 coordinatedB(:,2)+0.1]);
            is.aversive.coordinated.dHPC = InIntervals(bins,[coordinatedA(:,2)-0.1 coordinatedA(:,2)+0.1]);
            is.reward.coordinated.dHPC = InIntervals(bins,[coordinatedR(:,2)-0.1 coordinatedR(:,2)+0.1]);
            is.baseline.coordinated.vHPC = InIntervals(bins,[coordinatedB_V(:,2)-0.1 coordinatedB_V(:,2)+0.1]);
            is.aversive.coordinated.vHPC = InIntervals(bins,[coordinatedA_V(:,2)-0.1 coordinatedA_V(:,2)+0.1]);
            is.reward.coordinated.vHPC = InIntervals(bins,[coordinatedR_V(:,2)-0.1 coordinatedR_V(:,2)+0.1]);
            
            is.baseline.centroide = InIntervals(bins,[centroideB(:,2)-0.2 centroideB(:,2)+0.2]);
            is.aversive.centroide = InIntervals(bins,[centroideA(:,2)-0.2 centroideA(:,2)+0.2]);
            is.reward.centroide = InIntervals(bins,[centroideR(:,2)-0.2 centroideR(:,2)+0.2]);

            
            is.baseline.sws = InIntervals(bins,NREM_B);
            is.aversive.sws = InIntervals(bins,NREM_A);
            is.reward.sws = InIntervals(bins,NREM_R);
            
            is.baseline.rem = InIntervals(bins,REM_B);
            is.aversive.rem = InIntervals(bins,REM_A);
            is.reward.rem = InIntervals(bins,REM_R);
            
            is.aversive.run = InIntervals(bins,aversiveTS_run./1000);
            is.reward.run = InIntervals(bins,rewardTS_run./1000);
            
            % selection of place cells for each condition
%             is.aversive.PC.dHPC = ismember(clusters.dHPC,PC.dHPC.aversive(:,1));
%             is.aversive.PC.vHPC = ismember(clusters.vHPC,PC.vHPC.aversive(:,1));
%             is.reward.PC.dHPC = ismember(clusters.dHPC,PC.dHPC.reward(:,1));
%             is.reward.PC.vHPC = ismember(clusters.vHPC,PC.vHPC.reward(:,1));

            %% Explained variance calculation
            % Aversive
            if aversiveTS_run(1) < rewardTS_run(1)
                % Aversive
                m = mean(spiketrains_dHPC(InIntervals(bins,[baselineTS(1)/1000 aversiveTS(2)/1000]),:));
                s = std(spiketrains_dHPC(InIntervals(bins,[baselineTS(1)/1000 aversiveTS(2)/1000]),:));
                m1 = mean(spiketrains_vHPC(InIntervals(bins,[baselineTS(1)/1000 aversiveTS(2)/1000]),:));
                s1 = std(spiketrains_vHPC(InIntervals(bins,[baselineTS(1)/1000 aversiveTS(2)/1000]),:));
                
                x = [(spiketrains_dHPC(is.baseline.centroide,:)-m)./s];
                y = [(spiketrains_vHPC(is.baseline.centroide,:)-m1)./s1];
                [S1 , PvalMatrixB_NREM]=corr(x,y); clear x y
%                 S1 = S1 - diag(diag(S1));
                x = [(spiketrains_dHPC(is.aversive.centroide,:)-m)./s];
                y = [(spiketrains_vHPC(is.aversive.centroide,:)-m1)./s1];
                [S2 , PvalMatrixA_NREM]=corr(x,y); clear x y
%                 S2 = S2 - diag(diag(S2));
                x = [(spiketrains_dHPC(is.aversive.run,:)-m)./s];
                y = [(spiketrains_vHPC(is.aversive.run,:)-m1)./s1];
                [T , PvalMatrixR_RUN]=corr(x , y); clear c y m m1 s s1
%                 T = T - diag(diag(T));
                
%                 figure,
%                 m = max(max([S1;S2;T]));
%                 mm = min(min([S1;S2;T]));
%                 subplot(131),imagesc(S1),clim([mm m])
%                 subplot(132),imagesc(T),clim([mm m])
%                 subplot(133),imagesc(S2),clim([mm m])
%                 sgtitle('Aversive'),hold on
%                 clear m mm
                
                Sx = corrcoef(S2,T,'rows','complete');     Sx = Sx(1,2);
                Sy = corrcoef(S1,T,'rows','complete');     Sy = Sy(1,2);
                Sz = corrcoef(S1,S2,'rows','complete');     Sz = Sz(1,2);
                
                EV = (Sx-Sy*Sz);
                EV = (EV/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
                EV_A = [EV_A ; EV];
                
                REV = (Sy-Sx*Sz);
                REV = (REV/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
                REV_A = [REV_A ; REV];
                clear Sx Sy Sz REV x y S1 S2 EV T
                
                % Reward
                m = mean(spiketrains_dHPC(InIntervals(bins,[aversiveTS(1)/1000 rewardTS(2)/1000]),:));
                s = std(spiketrains_dHPC(InIntervals(bins,[aversiveTS(1)/1000 rewardTS(2)/1000]),:));
                m1 = mean(spiketrains_vHPC(InIntervals(bins,[aversiveTS(1)/1000 rewardTS(2)/1000]),:));
                s1 = std(spiketrains_vHPC(InIntervals(bins,[aversiveTS(1)/1000 rewardTS(2)/1000]),:));
                
                x = [(spiketrains_dHPC(is.aversive.centroide,:)-m)./s];
                y = [(spiketrains_vHPC(is.aversive.centroide,:)-m1)./s1];
                [S1 , PvalMatrixB_NREM]=corr(x,y); clear x y
%                 S1 = S1 - diag(diag(S1));

                x = [(spiketrains_dHPC(is.reward.centroide,:)-m)./s];
                y = [(spiketrains_vHPC(is.reward.centroide,:)-m1)./s1];
                [S2 , PvalMatrixA_NREM]=corr(x,y); clear x y
%                 S2 = S2 - diag(diag(S2));
                
                x = [(spiketrains_dHPC(is.reward.run,:)-m)./s];
                y = [(spiketrains_vHPC(is.reward.run,:)-m1)./s1];
                [T , PvalMatrixR_RUN]=corr(x , y); clear c y m m1 s s1
%                 T = T - diag(diag(T));

%                 figure,
%                 m = max(max([S1;S2;T]));
%                 mm = min(min([S1;S2;T]));
%                 subplot(131),imagesc(S1),clim([mm m])
%                 subplot(132),imagesc(T),clim([mm m])
%                 subplot(133),imagesc(S2),clim([mm m])
%                 sgtitle('Reward'),hold on
%                 clear m mm
                
                Sx = corrcoef(S2,T,'rows','complete');     Sx = Sx(1,2);
                Sy = corrcoef(S1,T,'rows','complete');     Sy = Sy(1,2);
                Sz = corrcoef(S1,S2,'rows','complete');     Sz = Sz(1,2);
                
                EV = (Sx-Sy*Sz);
                EV = (EV/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
                EV_R = [EV_R ; EV];
                
                REV = (Sy-Sx*Sz);
                REV = (REV/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
                REV_R = [REV_R ; REV];
                clear Sx Sy Sz REV x y S1 S2 EV T
                
            else
                % Reward
                m = mean(spiketrains_dHPC(InIntervals(bins,[baselineTS(1)/1000 rewardTS(2)/1000]),:));
                s = std(spiketrains_dHPC(InIntervals(bins,[baselineTS(1)/1000 rewardTS(2)/1000]),:));
                m1 = mean(spiketrains_vHPC(InIntervals(bins,[baselineTS(1)/1000 rewardTS(2)/1000]),:));
                s1 = std(spiketrains_vHPC(InIntervals(bins,[baselineTS(1)/1000 rewardTS(2)/1000]),:));
                
                x = [(spiketrains_dHPC(is.baseline.centroide,:)-m)./s];
                y = [(spiketrains_vHPC(is.baseline.centroide,:)-m1)./s1];
                [S1 , PvalMatrixB_NREM]=corr(x,y); clear x y
%                 S1 = S1 - diag(diag(S1));

                x = [(spiketrains_dHPC(is.reward.centroide,:)-m)./s];
                y = [(spiketrains_vHPC(is.reward.centroide,:)-m1)./s1];
                [S2 , PvalMatrixA_NREM]=corr(x,y); clear x y
%                 S2 = S2 - diag(diag(S2));
                
                x = [(spiketrains_dHPC(is.reward.run,:)-m)./s];
                y = [(spiketrains_vHPC(is.reward.run,:)-m1)./s1];
                [T , PvalMatrixR_RUN]=corr(x , y); clear c y m m1 s s1
%                 T = T - diag(diag(T));
                
%                 figure,
%                 m = max(max([S1;S2;T]));
%                 mm = min(min([S1;S2;T]));
%                 subplot(131),imagesc(S1),clim([mm m])
%                 subplot(132),imagesc(T),clim([mm m])
%                 subplot(133),imagesc(S2),clim([mm m])
%                 sgtitle('Reward'),hold on
%                 clear m mm
                
                Sx = corrcoef(S2,T,'rows','complete');     Sx = Sx(1,2);
                Sy = corrcoef(S1,T,'rows','complete');     Sy = Sy(1,2);
                Sz = corrcoef(S1,S2,'rows','complete');     Sz = Sz(1,2);
                
                EV = (Sx-Sy*Sz);
                EV = (EV/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
                EV_R = [EV_R ; EV];
                
                
                REV = (Sy-Sx*Sz);
                REV = (REV/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
                REV_R = [REV_R ; REV];
                clear Sx Sy Sz REV x y S1 S2 EV T
                
                % Aversive
                m = mean(spiketrains_dHPC(InIntervals(bins,[rewardTS(1)/1000 aversiveTS(2)/1000]),:));
                s = std(spiketrains_dHPC(InIntervals(bins,[rewardTS(1)/1000 aversiveTS(2)/1000]),:));
                m1 = mean(spiketrains_vHPC(InIntervals(bins,[rewardTS(1)/1000 aversiveTS(2)/1000]),:));
                s1 = std(spiketrains_vHPC(InIntervals(bins,[rewardTS(1)/1000 aversiveTS(2)/1000]),:));
                
                x = [(spiketrains_dHPC(is.reward.centroide,:)-m)./s];
                y = [(spiketrains_vHPC(is.reward.centroide,:)-m1)./s1];
                [S1 , PvalMatrixB_NREM]=corr(x,y); clear x y
%                 S1 = S1 - diag(diag(S1));

                x = [(spiketrains_dHPC(is.aversive.centroide,:)-m)./s];
                y = [(spiketrains_vHPC(is.aversive.centroide,:)-m1)./s1];
                [S2 , PvalMatrixA_NREM]=corr(x,y); clear x y
%                 S2 = S2 - diag(diag(S2));
                
                x = [(spiketrains_dHPC(is.aversive.run,:)-m)./s];
                y = [(spiketrains_vHPC(is.aversive.run,:)-m1)./s1];
                [T , PvalMatrixR_RUN]=corr(x , y); clear c y m m1 s s1
%                 T = T - diag(diag(T));

%                 figure,
%                 m = max(max([S1;S2;T]));
%                 mm = min(min([S1;S2;T]));
%                 subplot(131),imagesc(S1),clim([mm m])
%                 subplot(132),imagesc(T),clim([mm m])
%                 subplot(133),imagesc(S2),clim([mm m])
%                 sgtitle('Aversive'),hold on
%                 clear m mm
                
                
                Sx = corrcoef(S2,T,'rows','complete');     Sx = Sx(1,2);
                Sy = corrcoef(S1,T,'rows','complete');     Sy = Sy(1,2);
                Sz = corrcoef(S1,S2,'rows','complete');     Sz = Sz(1,2);
                
                EV = (Sx-Sy*Sz);
                EV = (EV/sqrt((1-Sy^2)*(1-Sz^2)))^2;
                EV_A = [EV_A ; EV];
                
                REV = (Sy-Sx*Sz);
                REV = (REV/sqrt((1-Sy^2)*(1-Sz^2)))^2;
                REV_A = [REV_A ; REV];
                clear Sx Sy Sz REV x y S1 S2 EV T 
                
            end
            
        end
        
        clear spiketrains_dHPC spiketrains_vHPC is
        clear Sx Sy Sz T S1 S2 CorrMatrixR_RUN CorrMatrixA_RUN
        clear coordinated coordinatedA coordinatedB coordinatedR coordinatedV coordinatedV_refined
        clear ii K Kinfo leftvalve rightvalve movement NREM NREM_A NREM_B NREM_R
        clear REM REM_A REM_B REM_R PvalMatrixA_NREM PvalMatrixB_NREM PvalMatrixR_RUN
        clear Shocks_filt Rewards_filt rewardTS rewardTS_run aversiveTS aversiveTS_run baselineTS
        clear camara Cell_type_classification cellulartype behavior
        clear group_dHPC group_vHPC maps_dHPC_A maps_dHPC_R maps_vHPC_A maps_vHPC_R
        clear ripplesD ripplesV spks_dHPC spks_tmp spks_vHPC WAKE clusters
        clear dX_int dY_int
        t
    end
    % end
    tt
    
end

figure,
subplot(121)
boxplot([REV_A*100 , EV_A*100] , {'REV' , 'EV'}),ylim([0 10]),ylabel('Explained Variance (%)'),title('Aversive')
subplot(122)
boxplot([REV_R*100 , EV_R*100], {'REV' , 'EV'}),ylim([0 10]),ylabel('Explained Variance (%)'),title('Reward')

M = [EV_A ; REV_A];
MM = [EV_R ; REV_R];
MMM = [M MM];
[~,~,stats] = anova2(MMM,20);
c1 = multcompare(stats);
c2 = multcompare(stats,"Estimate","row");
