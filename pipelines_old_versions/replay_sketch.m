clear
clc
close all

%% Parameters
path = {'E:\Rat126\Ephys\in_Pyr';'E:\Rat103\usable';'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat132\recordings\in_pyr';'E:\Rat165\in_pyr\'};%List of folders from the path

% for SU
criteria_fr = 0; %criteria to include or not a SU into the analysis
criteria_n = 3; % minimal number of neurons from each structure
criteria_type = 0; %criteria for celltype (0:pyr, 1:int, 2:all)

binSize = 0.001; % bin size for replay events detection

% Behavior
minimal_speed = 7; % minimal speed to detect active periods
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
        
%         % Definition of the place where the shocks were admninistrated
%         position_shocks = [];
%         for i = 1 : length(Shocks_filt)
%             [~ , ii] = min(abs(behavior.pos.aversive(:,1)-Shocks_filt(i)));
%             position_shocks = [position_shocks ; behavior.pos.aversive(ii,2)];
%             clear ii
%         end
%         clear i
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
            
            if and(length(group_dHPC)>criteria_n , length(group_vHPC)>criteria_n)
                %% Firing maps related to each event
                if and(isfile('dHPC_pc.mat') , isfile('vHPC_pc.mat'))
                    load('dHPC_pc.mat')
                    load('vHPC_pc.mat')
                    
                    
                end
            end
            
            
        end
        
            

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