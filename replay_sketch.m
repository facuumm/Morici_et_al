clear
clc
close all

%% Parameters
path = {'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat103\usable'};%List of folders from the path

% for SU
criteria_fr = 0.01; %criteria to include or not a SU into the analysis
criteria_n = 4; % minimal number of neurons from each structure
criteria_type = 0; %criteria for celltype (0:pyr, 1:int, 2:all)
criteria_distance = 600; %criteria to define if animal explored
n_SU_V = [];
n_SU_D = [];
Xedges = 36; %number of bins for RateMap construction
SD = 2; %SD for gauss kernel for placemaps

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
            behavior.speed.reward = LinearVelocity(pos,15);
            behavior.quiet.reward = QuietPeriods( behavior.speed.reward , 5 , 2);
            clear pos camaraR posx posy
            
            load('laps2.mat','posx','posy');
            [camaraA,~] = find((camara(:,1)-aversiveTS_run(1)/1000)>0,1,'first'); %TimeStamp of the begining of aversive
            pos = [camara(camaraA : camaraA+length(posx)-1),posx,posy];
            
            %interpolation of dropped frames
            ejeX = pos(~isnan(pos(:,2)),1); dX = pos(~isnan(pos(:,2)),2); dX_int = interp1(ejeX , dX , pos(:,1));
            ejeY = pos(~isnan(pos(:,3)),1); dY = pos(~isnan(pos(:,3)),3); dY_int = interp1(ejeY , dY , pos(:,1));
            pos(:,2) =dX_int; pos(:,3) =dY_int;%saving corrected pos
            
            behavior.pos.aversive = [pos];
            behavior.speed.aversive = LinearVelocity(pos,15);
            behavior.quiet.aversive = QuietPeriods(behavior.speed.aversive , 5 , 2);
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
            behavior.speed.reward = LinearVelocity(pos,15);
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
            behavior.speed.aversive = LinearVelocity(pos,15);
            behavior.quiet.aversive = QuietPeriods(behavior.speed.aversive , minimal_speed , minimal_speed_time);
            clear pos camaraA posx posy
        end
        
        % Generation of no-movements periods
        start = behavior.speed.reward(1,1);   stop = behavior.speed.reward(end,1);
        tmp = InvertIntervals(behavior.quiet.reward , start , stop); %keep only those higher than criteria
        [movement.reward , ~] =  ExcludeIntervals(tmp ,[Rewards_filt Rewards_filt+4]);
        clear tmp start stop
        
        start = behavior.speed.aversive(1,1);   stop = behavior.speed.aversive(end,1);
        tmp = InvertIntervals(behavior.quiet.aversive , start , stop);%keep only those higher than criteria
        [movement.aversive , ~] =  ExcludeIntervals(tmp ,[Shocks_filt-0.5 Shocks_filt+2]);
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
        
        
        if and(length(group_dHPC)>criteria_n , length(group_vHPC)>criteria_n)
            %% Firing maps related to each event
            % Detection of space close to the event
            pos_normalized = [behavior.pos.aversive(:,1),behavior.pos.aversive(:,2)];
            pos_normalized(:,2) =  pos_normalized(:,2) - min( pos_normalized(:,2));
            pos_normalized(:,2) =  pos_normalized(:,2) ./ max( pos_normalized(:,2));
            x = and(pos_normalized(:,2)<0.9 , pos_normalized(:,2)>0.1);
            y = ToIntervals(pos_normalized(:,1),x); y = y(y(:,2)-y(:,1)>1,:);
            laps.aversive.TS = y; clear x y
            laps.aversive.AB = [];             laps.aversive.BA = [];
            for i = 1:length(Shocks_filt)
                origin = Shocks_filt(i);
                [~,ind] = min(abs(pos_normalized(:,1)-origin));
                if origin < pos_normalized(end,1)
                    lap = laps.aversive.TS(laps.aversive.TS(:,2) < origin, :);
                    lap = lap(end,:);
                    lapAB =  Restrict(pos_normalized , lap);
                    lapAB(:,2) = abs(lapAB(:,2) - pos_normalized(ind,2));
                    speedAB = Restrict(behavior.speed.aversive , lap);
                    lapAB = lapAB(speedAB(:,2)>minimal_speed,:);
                    laps.aversive.AB = [laps.aversive.AB ; lapAB]; clear lapAB speedAB lap
                    
                    lap = laps.aversive.TS(laps.aversive.TS(:,1) > origin, :);
                    lap = lap(1,:);
                    lapBA =  Restrict(pos_normalized , lap);
                    lapBA(:,2) = abs(lapBA(:,2) - pos_normalized(ind,2));
                    speedBA = Restrict(behavior.speed.aversive , lap);
                    lapBA = lapBA(speedBA(:,2)>minimal_speed,:);
                    laps.aversive.BA = [laps.aversive.BA ; lapBA]; clear lapBA speedBA lap                    

                end
                %             plot(x(:,1) , x(:,2),'*')
                clear lap ind x
            end
            clear pos_normalized i
            
            pos_normalized = [behavior.pos.reward(:,1),behavior.pos.reward(:,2)];
            pos_normalized(:,2) =  pos_normalized(:,2) - min( pos_normalized(:,2));
            pos_normalized(:,2) =  pos_normalized(:,2) ./ max( pos_normalized(:,2));
            x = and(pos_normalized(:,2)<0.9 , pos_normalized(:,2)>0.1);
            y = ToIntervals(pos_normalized(:,1),x); y = y(y(:,2)-y(:,1)>1,:);
            laps.reward.TS = y; clear x y
            laps.reward.AB = [];
            Rewards_filt = sort(Rewards_filt);
            for i = 2:length(Rewards_filt)
                start = Rewards_filt(i-1,2);
                stop = Rewards_filt(i,1);
                if stop-start>1
                    if and(start>pos_normalized(1,1) , stop<pos_normalized(end,1))
                        [~,ind] = min(abs(pos_normalized(:,1) - start));
                        speed = Restrict(behavior.speed.reward,[start stop]);
                        speed = QuietPeriods(speed,minimal_speed,minimal_speed_time);
                        [x , ~] = SubtractIntervals([start,stop],speed);
                        lap =  Restrict(pos_normalized , x);
                        lap(:,2) = abs((lap(:,2) - pos_normalized(ind,2)));
                        laps.reward.AB = [laps.reward.AB ; lap];
                        clear lap ind x
                    end
                end
            end
            clear pos_normalized i
            if and(length(laps.reward.AB)>criteria_distance , and(length(laps.aversive.AB)>criteria_distance,length(laps.aversive.BA)>criteria_distance))
                %% Firing maps calculation
                PHIST.dHPC.aversive = [];        PHIST.dHPC.reward = [];
                clusterdHPC.aversive = [];      clusterdHPC.reward = [];
                for ii=1:length(group_dHPC)
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
                for ii=1:length(group_vHPC)
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
                
                %% Calculation of maps normalized to the event
                PHIST.dHPC.normalizado.aversive.AB = [];                PHIST.dHPC.normalizado.aversive.BA = [];
                PHIST.dHPC.normalizado.reward = [];
                clusterdHPC.normalizado.aversive = [];      clusterdHPC.normalizado.reward = [];
                for ii=1:length(group_dHPC)
                    cluster = group_dHPC(ii,1);
                    celltype = logical(cellulartype(cellulartype(:,1) == cluster,2));
                    if celltype
                        spks = spks_dHPC(spks_dHPC(:,1)==cluster,2);
                        b = [Shocks_filt(:,1)-2 Shocks_filt(:,1)+2];
                        cond1 = Restrict(spks,b);
                        
                        b = [Rewards_filt(:,1)-2 Rewards_filt(:,1)+2];
                        cond2 = Restrict(spks,b);
                        
                        if length(cond1)>5
                            % --- Aversive ---
                            spks_tmp = Restrict(spks ,[min(laps.aversive.AB(:,1)) max(laps.aversive.AB(:,1))]);
                            [curve , stats] = FiringCurve(laps.aversive.AB , spks_tmp , 'smooth' , SD , 'nBins' , Xedges , 'minSize' , 2 , 'minPeak' , 0.2);
                            
                            PHIST.dHPC.normalizado.aversive.AB = [PHIST.dHPC.normalizado.aversive.AB ; curve.rate];
                            clear curve stats
                                
                            spks_tmp = Restrict(spks ,[min(laps.aversive.BA(:,1)) max(laps.aversive.BA(:,1))]);
                            [curve , stats] = FiringCurve(laps.aversive.BA , spks_tmp , 'smooth' , SD , 'nBins' , Xedges , 'minSize' , 2 , 'minPeak' , 0.2);
                            
                            PHIST.dHPC.normalizado.aversive.BA = [PHIST.dHPC.normalizado.aversive.BA ; curve.rate];
                            clear curve stats
                            
                            clusterdHPC.normalizado.aversive = [clusterdHPC.normalizado.aversive ; cluster];
                        end
                        
                        
                        if length(cond2)>5
                            % --- Reward ---
                            spks_tmp = Restrict(spks ,[min(laps.reward.AB(:,1)) max(laps.reward.AB(:,1))]);
                            [curve , stats] = FiringCurve(laps.reward.AB , spks_tmp , 'smooth' , SD , 'nBins' , Xedges , 'minSize' , 2 , 'minPeak' , 0.2);
                            
                            PHIST.dHPC.normalizado.reward = [PHIST.dHPC.normalizado.reward ; curve.rate];
                            clear curve stats
                            clusterdHPC.normalizado.reward = [clusterdHPC.normalizado.reward ; cluster];
                        end
                    end
                    
                    clear celltype tmp b cluster
                end
                
                PHIST.vHPC.normalizado.aversive.AB = [];                PHIST.vHPC.normalizado.aversive.BA = [];
                PHIST.vHPC.normalizado.reward = [];
                clustervHPC.normalizado.aversive = [];      clustervHPC.normalizado.reward = [];
                for ii=1:length(group_vHPC)
                    cluster = group_vHPC(ii,1);
                    celltype = logical(cellulartype(cellulartype(:,1) == cluster,2));
                    if celltype
                        spks = spks_vHPC(spks_vHPC(:,1)==cluster,2);
                        b = [Shocks_filt(:,1)-2 Shocks_filt(:,1)+2];
                        cond1 = Restrict(spks,b);
                        
                        b = [Rewards_filt(:,1)-2 Rewards_filt(:,1)+2];
                        cond2 = Restrict(spks,b);
                        
                        if length(cond1)>5
                            % --- Aversive ---
                            spks_tmp = Restrict(spks ,[min(laps.aversive.AB(:,1)) max(laps.aversive.AB(:,1))]);
                            [curve , stats] = FiringCurve(laps.aversive.AB , spks_tmp , 'smooth' , SD , 'nBins' , Xedges , 'minSize' , 2 , 'minPeak' , 0.2);
                            
                            PHIST.vHPC.normalizado.aversive.AB = [PHIST.vHPC.normalizado.aversive.AB ; curve.rate];
                            clear curve stats
                                
                            spks_tmp = Restrict(spks ,[min(laps.aversive.BA(:,1)) max(laps.aversive.BA(:,1))]);
                            [curve , stats] = FiringCurve(laps.aversive.BA , spks_tmp , 'smooth' , SD , 'nBins' , Xedges , 'minSize' , 2 , 'minPeak' , 0.2);
                            
                            PHIST.vHPC.normalizado.aversive.BA = [PHIST.vHPC.normalizado.aversive.BA ; curve.rate];
                            clear curve stats    
                            
                            clustervHPC.normalizado.aversive = [clustervHPC.normalizado.aversive ; cluster];
                        end
                        if length(cond2)>5
                            % --- Reward ---
                            spks_tmp = Restrict(spks ,[min(laps.reward.AB(:,1)) max(laps.reward.AB(:,1))]);
                            [curve , stats] = FiringCurve(laps.reward.AB , spks_tmp , 'smooth' , SD , 'nBins' , Xedges , 'minSize' , 2 , 'minPeak' , 0.2);
                            
                            PHIST.vHPC.normalizado.reward = [PHIST.vHPC.normalizado.reward ; curve.rate];
                            clear curve stats
                            clustervHPC.normalizado.reward = [clustervHPC.normalizado.reward ; cluster];
                        end
                    end
                    
                    clear celltype tmp b cluster
                end
                
                
                %% Replay detection
                freq = 1/binSize;
                limits = [0 segments.Var1(end)/1000];
                spks = [];
                %Replay events in dHPC
                for ii=1:length(group_dHPC)
                    cluster = group_dHPC(ii,1);
                    celltype = logical(cellulartype(cellulartype(:,1) == cluster,2));
                    if celltype
                        spks = [spks;spks_dHPC(spks_dHPC(:,1)==cluster,2)];
                    end
                    clear celltype
                end
                [MUA.dHPC,bins]=binspikes(sort(spks,'ascend'),freq,limits);clear spks
                MUA.dHPC = Smooth(MUA.dHPC , 30 ,'kernel','gaussian');
                replay.dHPC = MUA.dHPC>mean(MUA.dHPC)+std(MUA.dHPC)*2;
                replay.dHPC = ToIntervals(bins',replay.dHPC);
                replay.dHPC = replay.dHPC(replay.dHPC(:,2)-replay.dHPC(:,1)>0.1,:);
                replay.dHPC = replay.dHPC(replay.dHPC(:,2)-replay.dHPC(:,1)<0.8,:);
                [replay.dHPC] = merge_events(replay.dHPC, 0.04);
                
                
                %filtering replay events by amount of PCs
                count = [];
                cluster_dHPC = unique([clusterdHPC.normalizado.aversive ; clusterdHPC.normalizado.reward]);
                for i = 1 :  length(cluster_dHPC)
                    cluster = cluster_dHPC(i,1);
                    spks = spks_dHPC(spks_dHPC(:,1)==cluster,2);
                    [status,interval,index] = InIntervals(spks,replay.dHPC);
                    interval = unique(interval);
                    count = [count ; interval(interval~=0)];
                    clear spks cluster status interval index
                end
                [gc,grps] = groupcounts(count);
                replay.dHPC = replay.dHPC(grps(gc>length(cluster_dHPC)*0.30),:);
                
                %Replay events in vHPC
                spks = [];
                for ii=1:length(group_vHPC)
                    cluster = group_vHPC(ii,1);
                    celltype = logical(cellulartype(cellulartype(:,1) == cluster,2));
                    if celltype
                        spks = [spks;spks_vHPC(spks_vHPC(:,1)==cluster,2)];
                    end
                    clear celltype
                end
                [MUA.vHPC,bins]=binspikes(sort(spks,'ascend'),freq,limits);clear spks
                MUA.vHPC = Smooth(MUA.vHPC , 30 ,'kernel','gaussian');
                replay.vHPC = MUA.vHPC>mean(MUA.vHPC)+std(MUA.vHPC)*2;
                replay.vHPC = ToIntervals(bins',replay.vHPC);
                replay.vHPC = replay.vHPC(replay.vHPC(:,2)-replay.vHPC(:,1)>0.1,:);
                replay.vHPC = replay.vHPC(replay.vHPC(:,2)-replay.vHPC(:,1)<0.8,:);
                [replay.vHPC] = merge_events(replay.vHPC, 0.04);
                
                %filtering replay events by amount of PCs
                count = [];
                cluster_vHPC = unique([clustervHPC.normalizado.aversive; clustervHPC.normalizado.reward]);
                for i = 1 :  length(group_vHPC)
                    cluster = group_vHPC(i,1);
                    spks = spks_vHPC(spks_vHPC(:,1)==cluster,2);
                    [status,interval,index] = InIntervals(spks,replay.vHPC);
                    interval = unique(interval);
                    count = [count ; interval(interval~=0)];
                    clear spks cluster status interval index
                end
                [gc,grps] = groupcounts(count);
                replay.vHPC = replay.vHPC(grps(gc>length(group_vHPC)*0.30),:);
                
                %% Save Replay percentage
                % outside ripples
                % --- Baseline ---
                x = length(Restrict(replay.dHPC ,baselineTS./1000));
                y = baselineTS(2)/1000 - baselineTS(1)/1000;
                Replay.candidate.dHPC.baseline=[Replay.candidate.dHPC.baseline ; x/y]; clear x y
                % --- Reward ---
                x = length(Restrict(replay.dHPC,rewardTS./1000));
                y = rewardTS(2)/1000 - rewardTS(1)/1000;
                Replay.candidate.dHPC.reward=[Replay.candidate.dHPC.reward ; x/y]; clear x y
                % --- Aversive ---
                x = length(Restrict(replay.dHPC ,aversiveTS./1000));
                y = aversiveTS(2)/1000 - aversiveTS(1)/1000;
                Replay.candidate.dHPC.aversive=[Replay.candidate.dHPC.aversive ; x/y]; clear x y
                
                % --- Baseline ---
                x = length(Restrict(replay.vHPC , baselineTS./1000));
                y = baselineTS(2)/1000 - baselineTS(1)/1000;
                Replay.candidate.vHPC.baseline=[Replay.candidate.vHPC.baseline ; x/y]; clear x y
                % --- Reward ---
                x = length(Restrict(replay.vHPC , rewardTS./1000));
                y = rewardTS(2)/1000 - rewardTS(1)/1000;
                Replay.candidate.vHPC.reward=[Replay.candidate.vHPC.reward ; x/y]; clear x y
                % --- Aversive ---
                x = length(Restrict(replay.vHPC , aversiveTS./1000));
                y = aversiveTS(2)/1000 - aversiveTS(1)/1000;
                Replay.candidate.vHPC.aversive=[Replay.candidate.vHPC.aversive ; x/y]; clear x y
                
                % within ripples
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
                Replay.selected.dHPC.aversive=[Replay.selected.dHPC.aversive ; x/y]; clear x y
                
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
                Replay.selected.vHPC.aversive=[Replay.selected.vHPC.aversive ; x/y]; clear x y                
                
                %% Bayesian Decoding across conditions
                % Position
                decoded_positions.dHPC.aversive = [];
                decoded_positions.vHPC.aversive = [];
                decoded_positions.dHPC.reward = [];
                decoded_positions.vHPC.reward = [];
                
                % --- Aversive ---
                tmp.vHPC = Restrict(replay.vHPC ,aversiveTS./1000);
                tmp.dHPC = Restrict(replay.dHPC ,aversiveTS./1000);
                
                % for ventral hippocampus
                realReplay.vHPC.aversive.forward = [];
                realReplay.vHPC.aversive.backward = [];
                for i = 1 : length(tmp.vHPC)
                    % Decoding using dHPC SUs
                    start = tmp.vHPC(i,1);
                    stop = tmp.vHPC(i,2);
                    bin = [start : 0.01 :stop];
                    tmp1 = [];
                    for ii = 2 : length(bin)-1
                        nSpks = count_spks(spks_vHPC, clustervHPC.normalizado.aversive, bin(ii-1), bin(ii+1));
                        probability = bayesian_replay(PHIST.vHPC.normalizado.aversive.AB, nSpks, bin(ii-1), bin(ii+1));
                        [o,p] = findpeaks(probability,'WidthReference','halfheight','MinPeakProminence',mean(probability),'SortStr','descend');
                        if not(isempty(p))
                            tmp1 = [tmp1 ; o(1) p(1) ii];
                        end
                        clear nSpks probability o p
                    end
                    
                    
                    if length(tmp1) > 3
                        c = corrcoef(tmp1(:,3) , tmp1(:,2));
                        
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
                            clear probability nSpks start stop p c
                        end
                    end
                end
                disp('--- Done Ventral Aversive ---')
                
                % for dorsal hippocampus
                realReplay.dHPC.aversive.forward = [];
                realReplay.dHPC.aversive.backward = [];
                for i = 1 : length(tmp.dHPC)
                    % Decoding using dHPC SUs
                    start = tmp.dHPC(i,1);
                    stop = tmp.dHPC(i,2);
%                     bin = ((stop-start)/2) + start;
                    bin = [start : 0.01 :stop];
                    %                     if length(bin)>=4
                    tmp1 = [];
                    for ii = 2 : length(bin)-1
                        nSpks = count_spks(spks_dHPC, clusterdHPC.normalizado.aversive, bin(ii-1), bin(ii+1));
                        probability = bayesian_replay(PHIST.dHPC.normalizado.aversive, nSpks, bin(ii-1), bin(ii+1));
                        [o,p] = findpeaks(probability,'WidthReference','halfheight','MinPeakProminence',mean(probability),'SortStr','descend');
                        %                         [o p] = max(probability);
                        if not(isempty(p))
                            tmp1 = [tmp1 ; o(1) p(1) ii];
%                         plot(probability),hold on
                        end
                        clear nSpks probability o p
                    end
                    
                    
                    if length(tmp1) > 3
                        c = corrcoef(tmp1(:,3) , tmp1(:,2));
                        
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
                                    realReplay.vHPC.aversive.forward = [realReplay.dHPC.aversive.forward ; start stop];
                                    plot(tmp1(:,3) , tmp1(:,2),'*')
                                end
                                
                            end
                            clear probability nSpks start stop p c
                        end
                    end
                end
                disp('--- Done Dorsal Aversive ---')
                
                % --- Reward ---
                tmp.vHPC = Restrict(Restrict(replay.vHPC , [ripplesV(:,2)-0.1 ripplesV(:,2)+0.1]),rewardTS./1000);
                tmp.dHPC = Restrict(Restrict(replay.dHPC , [ripplesD(:,2)-0.1 ripplesD(:,2)+0.1]),rewardTS./1000);
                
                % for ventral hippocampus
                realReplay.vHPC.reward.forward = [];
                realReplay.vHPC.reward.backward = [];
                for i = 1 : length(tmp.vHPC)
                    % Decoding using dHPC SUs
                    start = tmp.vHPC(i,1);
                    stop = tmp.vHPC(i,2);
%                     bin = ((stop-start)/2) + start;
                    bin = [start : 0.01 :stop];
                    %                     if length(bin)>=4
                    tmp1 = [];
                    for ii = 2 : length(bin)-1
                        nSpks = count_spks(spks_vHPC, clustervHPC.normalizado.reward, bin(ii-1), bin(ii+1));
                        probability = bayesian_replay(PHIST.vHPC.normalizado.reward, nSpks, bin(ii-1), bin(ii+1));
                        [o,p] = findpeaks(probability,'WidthReference','halfheight','MinPeakProminence',mean(probability),'SortStr','descend');
                        %                         [o p] = max(probability);
                        if not(isempty(p))
                            tmp1 = [tmp1 ; o(1) p(1) ii];
%                         plot(probability),hold on
                        end
                        clear nSpks probability o p
                    end
                    
                    
                    if length(tmp1) > 3
                        c = corrcoef(tmp1(:,3) , tmp1(:,2));
                        
                        shuffle = [];
                        for ii = 1:1000
                            p = corrcoef(tmp1(randperm(length(tmp1)),3) , tmp1(:,2));
                            shuffle = [shuffle ; p(1,2)];
                        end
                        
                        if not(isnan(c(1,2))) %if correlation exist
                            if sign(c(1,2))<0 %for reverse replay
                                p = sum(shuffle<c(1,2))/1000;
                                if p<0.05
                                    realReplay.vHPC.aversive.backward = [realReplay.vHPC.reward.backward ; start stop];
                                    figure,plot(tmp1(:,3) , tmp1(:,2),'*')
                                end
                            else %for replay
                                p = sum(shuffle>c(1,2))/1000;
                                
                                if p<0.05
                                    realReplay.vHPC.aversive.forward = [realReplay.vHPC.reward.forward ; start stop];
                                    plot(tmp1(:,3) , tmp1(:,2),'*')
                                end
                                
                            end
                            clear probability nSpks start stop p c
                        end
                    end
                end
                disp('--- Done Ventral Reward ---')

                % for dorsal hippocampus
                realReplay.dHPC.reward.forward = [];
                realReplay.dHPC.reward.backward = [];
                for i = 1 : length(tmp.dHPC)
                    % Decoding using dHPC SUs
                    start = tmp.dHPC(i,1);
                    stop = tmp.dHPC(i,2);
                    bin = [start : 0.01 :stop];
                    tmp1 = [];
                    for ii = 2 : length(bin)-1
                        nSpks = count_spks(spks_dHPC, clusterdHPC.normalizado.reward, bin(ii-1), bin(ii+1));
                        probability = bayesian_replay(PHIST.dHPC.normalizado.reward, nSpks, bin(ii-1), bin(ii+1));
                        [o,p] = findpeaks(probability,'WidthReference','halfheight','MinPeakProminence',mean(probability),'SortStr','descend');
                        %                         [o p] = max(probability);
                        if not(isempty(p))
                            tmp1 = [tmp1 ; o(1) p(1) ii];
%                         plot(probability),hold on
                        end
                        clear nSpks probability o p
                    end
                    
                    
                    if length(tmp1) > 3
                        c = corrcoef(tmp1(:,3) , tmp1(:,2));
                        
                        shuffle = [];
                        for ii = 1:1000
                            p = corrcoef(tmp1(randperm(length(tmp1)),3) , tmp1(:,2));
                            shuffle = [shuffle ; p(1,2)];
                        end
                        
                        if not(isnan(c(1,2))) %if correlation exist
                            if sign(c(1,2))<0 %for reverse replay
                                p = sum(shuffle<c(1,2))/1000;
                                if p<0.05
                                    realReplay.dHPC.aversive.backward = [realReplay.dHPC.reward.backward ; start stop];
                                    figure,plot(tmp1(:,3) , tmp1(:,2),'*')
                                end
                            else %for replay
                                p = sum(shuffle>c(1,2))/1000;
                                
                                if p<0.05
                                    realReplay.vHPC.aversive.forward = [realReplay.dHPC.reward.forward ; start stop];
                                    plot(tmp1(:,3) , tmp1(:,2),'*')
                                end
                            end
                            clear probability nSpks start stop p c
                        end
                    end
                end
                disp('--- Done Dorsal Reward ---')
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



subplot(121),boxplot([Replay.selected.dHPC.baseline Replay.selected.dHPC.reward Replay.selected.dHPC.aversive])
subplot(122),boxplot([Replay.selected.vHPC.baseline Replay.selected.vHPC.reward Replay.selected.vHPC.aversive])


figure
subplot(121),boxplot([Replay.candidate.dHPC.baseline Replay.candidate.dHPC.reward Replay.candidate.dHPC.aversive])
subplot(122),boxplot([Replay.candidate.vHPC.baseline Replay.candidate.vHPC.reward Replay.candidate.vHPC.aversive])