clear
clc
close all

%% Parameters
path = {'E:\Rat103\usable';'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat132\recordings\in_pyr';'E:\Rat165\in_pyr\'};%List of folders from the path

%Sleep
time_criteria = 1; % minimal time to include a NREM epoch (in min)

% What par of the code I want to run
S = logical(1);   % Reactivation Strength Calculation
MUAselection = logical(0); % to select ripples by their MUA
W = 'E'; % to select what kind of ripples I want to check
% E= all coordinated ripples, DV dRipple-vRipple, VD vRipple-dRipple
% CB = cooridnated bursts
% N= NREM, R= REM
TA =  logical(0); % Trigger Reactivation Strength
TPks = logical(0); %trigger CCG assemblies peaks to events
REC = logical(0); % Assemblie Recruitment during cooridnated events
C = logical(0);   % Check if coordinated ripples are ocurring in a manner
SRC = logical(0); % If I want to calculate the tuning curve for shocks

% for SU
criteria_fr = 0; %criteria to include or not a SU into the analysis
criteria_n = [3 3]; % minimal number of neurons from each structure [vHPC dHPC]
criteria_type = 0; %criteria for celltype (0:pyr, 1:int, 2:all)
binSize = 0.025; %for qssemblie detection qnd qxctivity strength
binSizeT = 1; %for speed strength correlations
n_SU_V = 0;
n_SU_D = 0;

win = 300; % time window for bin construction

% Behavior
minimal_speed = 7; % minimal speed to detect quite periods
minimal_speed_time = 2; % minimal time to detect quite periods


counts.both = [];
counts.dHPC = [];
counts.vHPC = [];

reactivation.aversive.dvHPC = [];
reactivation.reward.dvHPC = [];
reactivation.aversive.dHPC = [];
reactivation.reward.dHPC = [];
reactivation.aversive.vHPC = [];
reactivation.reward.vHPC = [];

normalization = true; % to define if normalization over Reactivation Strength is applied or not
th = 5; % threshold for peak detection

for_scatter.aversive = [];
for_scatter.reward = [];
Rsquared.aversive = [];
Rsquared.reward = [];


gain.both.reward.pre = [];     gain.both.reward.post = [];
gain.both.aversive.pre = [];   gain.both.aversive.post = [];

percentages = [];

weigths.up = [];
weigths.down = [];

curveShock = [];
ResponsiveS = [];

Number_of_assemblies.aversive = [];
Number_of_assemblies.reward = [];

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
        
        %% load coordinated ripple bursts
        load('coordinated_ripple_bursts.mat')
        coordinated_ripple_bursts = [coordinated_ripple_bursts(:,1)  coordinated_ripple_bursts(:,3)];
        %         [coordinated_ripple_bursts] = merge_events(coordinated_ripple_bursts, 0.05);
        
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
        cooridnated_eventDV = [];
        cooridnated_eventVD = [];
        for i = 1:length(ripplesD)
            r = ripplesD(i,:);
            tmp = sum(and(ripplesV(:,2)>= r(1,2)-0.1, ripplesV(:,2)<= r(1,2)+0.1));
            if tmp>0
                z = ripplesV(and(ripplesV(:,2)>= r(1,2)-0.1, ripplesV(:,2)<= r(1,2)+0.1),:);
                coordinatedV = [coordinatedV ; z];
                [p,indice] = min(abs(r(2)-z(:,2)));
                coordinatedV_refined = [coordinatedV_refined ; z(indice,:)];
                coordinated = [coordinated ; r];
                
                cooridnated_event = [cooridnated_event ; min([r(1) , z(indice,1)]) , min(min(z(indice,2),r(2)))+abs(z(indice,2)-r(2))/2 , max([r(3) , z(indice,3)])];
                
                if r(2)<z(indice,2)
                    cooridnated_eventDV = [cooridnated_eventDV ; min([r(1) , z(indice,1)]) , min(min(z(indice,2),r(2)))+abs(z(indice,2)-r(2))/2 , max([r(3) , z(indice,3)])];
                else
                    cooridnated_eventVD = [cooridnated_eventVD ; min([r(1) , z(indice,1)]) , min(min(z(indice,2),r(2)))+abs(z(indice,2)-r(2))/2 , max([r(3) , z(indice,3)])];
                end
                
                clear tmp2 tmp1 p indice z
            end
            clear r
        end
        clear x tmp i
        
        % Store events time stamps
        % dRipples
        ripples.dHPC.baseline = Restrict(ripplesD , NREM.baseline);
        ripples.dHPC.reward = Restrict(ripplesD , NREM.reward);
        ripples.dHPC.aversive = Restrict(ripplesD , NREM.aversive);
        % vRipples
        ripples.vHPC.baseline = Restrict(ripplesV , NREM.baseline);
        ripples.vHPC.reward = Restrict(ripplesV , NREM.reward);
        ripples.vHPC.aversive = Restrict(ripplesV , NREM.aversive);
        % coordinated dRipples
        ripples.dHPC.coordinated.baseline = Restrict(coordinated , NREM.baseline);
        ripples.dHPC.coordinated.reward = Restrict(coordinated , NREM.reward);
        ripples.dHPC.coordinated.aversive = Restrict(coordinated , NREM.aversive);
        % coordinated vRipples
        ripples.vHPC.coordinated.baseline = Restrict(coordinatedV_refined , NREM.baseline);
        ripples.vHPC.coordinated.reward = Restrict(coordinatedV_refined , NREM.reward);
        ripples.vHPC.coordinated.aversive = Restrict(coordinatedV_refined , NREM.aversive);
        %coordinated event
        cooridnated_event((cooridnated_event(:,3)-cooridnated_event(:,1)<0.04),:) = [];
        ripple_event.baseline = Restrict(cooridnated_event,baselineTS./1000);
        ripple_event.reward = Restrict(cooridnated_event,rewardTS./1000);
        ripple_event.aversive = Restrict(cooridnated_event,aversiveTS./1000);
        ripple_event.all = cooridnated_event;
        % coordinated event when dRipple was first
        ripple_event.DV.baseline = Restrict(cooridnated_eventDV,baselineTS./1000);
        ripple_event.DV.reward = Restrict(cooridnated_eventDV,rewardTS./1000);
        ripple_event.DV.aversive = Restrict(cooridnated_eventDV,aversiveTS./1000);
        ripple_event.DV.all = cooridnated_eventDV;
        % coordinated event when vRipple was first
        ripple_event.VD.baseline = Restrict(cooridnated_eventVD,baselineTS./1000);
        ripple_event.VD.reward = Restrict(cooridnated_eventVD,rewardTS./1000);
        ripple_event.VD.aversive = Restrict(cooridnated_eventVD,aversiveTS./1000);
        ripple_event.VD.all = cooridnated_eventVD;
        
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
        
        %% Shock responsive cells
        if SRC
            [curve , binsS , responsive] = SU_responsivness(spks,clusters.all,Shocks_filt,aversiveTS_run./1000,[0 1],6,0.1,1,'gain',2);
            curveShock = [curveShock , curve];
            i = [ones(size(clusters.dHPC,1),1) ; ones(size(clusters.vHPC,1),1)*2];
            ResponsiveS = [ResponsiveS ; i , responsive'];
            clear curve responsive i
        end
        
        criteria_n = [3 3];
        
        %% Assemblies detection
        if and(numberD >= criteria_n(1),numberV >= criteria_n(2))
            if MUAselection
                disp('Selection of coordinated events with high MUA')
                freq = 1/0.001;
                limits = [0 segments.Var1(end)/1000];
                tmp = [];
                for ii=1:length(cellulartype) % selection of SU for spks train
                    if logical(cellulartype(ii,2))
                        tmp = [tmp ; spks(spks(:,1)==cellulartype(ii,1),2)];
                    end
                end
                [MUA,b]=binspikes(sort(tmp,'ascend'),freq,limits); % spiketrain
                %             MUA = Smooth(MUA , 30 ,'kernel','gaussian'); % smooth
                MUA = smoothdata(MUA , 'gaussian' , 30); % smooth
                MUA = MUA>mean(MUA)+std(MUA)*2; % detection of high MUA
                MUA = ToIntervals(b',MUA); % definition of periods
                MUA = merge_events(MUA, 0.04); % merge periods with short IEI
                
                if strcmp(W,'E')
                    a = ripple_event.all(:,2); % center of the event
                elseif strcmp(W,'DV')
                    a = ripple_event.DV.all(:,2); % center of the event
                elseif strcmp(W,'VD')
                    a = ripple_event.VD.all(:,2); % center of the event
                end
                
                a = Restrict(a , MUA); % restriction of events to MUA periods
                ripple_event.filtered.all = Restrict(a,NREM.all); % store all events
                ripple_event.filtered.baseline = Restrict(a,NREM.baseline); % store baseline events
                ripple_event.filtered.reward = Restrict(a,NREM.reward); % store reward events
                ripple_event.filtered.aversive = Restrict(a,NREM.aversive); % store aversive events
                clear a b  spks
                
            else
                a = ripplesD(:,1:3); % center of the event
                ripple_event.filtered.all = Restrict(a,NREM.all); % store all events
                ripple_event.filtered.baseline = Restrict(a,NREM.baseline); % store baseline events
                ripple_event.filtered.reward = Restrict(a,NREM.reward); % store reward events
                ripple_event.filtered.aversive = Restrict(a,NREM.aversive); % store aversive events
                clear a b  spks
            end
            
            %%
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
            
            
            percentages = [percentages  ; sum(cond.dHPC.aversive) , sum(cond.dHPC.reward) , sum(sum(AR.dHPC)) , sum(cond.vHPC.aversive) , sum(cond.vHPC.reward) , sum(sum(AR.vHPC)) , sum(cond.both.aversive) , sum(cond.both.reward) , sum(sum(AR.both))];
            clear A R r p
            
            % to save the clusters I used for further analysis
            %             save([cd,'\clusters_included_in_assemblies.mat'],'clusters')
            
            
            %% SpikeTrains construction
            limits = [0 segments.Var1(end)/1000];
            events = [];
            [Spikes , bins , Clusters] = spike_train_construction([spks_dHPC;spks_vHPC], clusters.all, cellulartype, binSize, limits, events, true, true);
            clear limits events
            
            %% Assemblies activation in the entier recording
            
            if S
                % NREM sleep
                if strcmp(W,'N')
                    is.sws.baseline = InIntervals(bins,NREM.baseline);
                    is.sws.reward = InIntervals(bins,NREM.reward);
                    is.sws.aversive = InIntervals(bins,NREM.aversive);
                elseif strcmp(W,'R')
                    is.sws.baseline = InIntervals(bins,REM.baseline);
                    is.sws.reward = InIntervals(bins,REM.reward);
                    is.sws.aversive = InIntervals(bins,REM.aversive);
                elseif or(or(strcmp(W,'E'),strcmp(W,'DV')),strcmp(W,'VD'))
                    % coordinated event
                    is.sws.baseline = InIntervals(bins,[ripple_event.filtered.baseline(:,1)-0.2 ripple_event.filtered.baseline(:,1)+0.2]);
                    is.sws.reward = InIntervals(bins,[ripple_event.filtered.reward(:,1)-0.2 ripple_event.filtered.reward(:,1)+0.2]);
                    is.sws.aversive = InIntervals(bins,[ripple_event.filtered.aversive(:,1)-0.2 ripple_event.filtered.aversive(:,1)+0.2]);
                elseif strcmp(W,'CB')
                    % coordinated event
                    is.sws.baseline = InIntervals(bins,ripple_bursts.baseline);
                    is.sws.reward = InIntervals(bins,ripple_bursts.reward);
                    is.sws.aversive = InIntervals(bins,ripple_bursts.aversive);
                end
                
                is.sws.runaversive = InIntervals(bins,movement.aversive);
                is.sws.runreward = InIntervals(bins,movement.reward);
                
                is.aversive = InIntervals(bins,aversiveTS_run./1000);
                is.reward = InIntervals(bins,rewardTS_run./1000);
                
                %% Reactivation Strenght
                
                %                 if aversiveTS_run(1) < rewardTS_run(1)
                if sum(cond.both.aversive)>=1
                    [R] = reactivation_strength(patterns.all.aversive , cond.both.aversive , [bins' , Spikes] , is.sws , th , 'A' , config , normalization);
                    reactivation.aversive.dvHPC = [reactivation.aversive.dvHPC ; R];
                    RBA = R; clear R
                end
                
                if sum(cond.dHPC.aversive)>=1
                    [R] = reactivation_strength(patterns.all.aversive , cond.dHPC.aversive , [bins' , Spikes] , is.sws , th , 'A' , config , normalization);
                    reactivation.aversive.dHPC = [reactivation.aversive.dHPC ; R];
                    RDA = R; clear R
                end
                
                if sum(cond.vHPC.aversive)>=1
                    [R] = reactivation_strength(patterns.all.aversive , cond.vHPC.aversive , [bins' , Spikes] , is.sws , th , 'A' , config , normalization);
                    reactivation.aversive.vHPC = [reactivation.aversive.vHPC ; R];
                    RVA = R; clear R
                end
                %                 end
                
                
                %                 if aversiveTS_run(1) > rewardTS_run(1)
                if sum(cond.both.reward)>=1
                    [R] = reactivation_strength(patterns.all.reward , cond.both.reward , [bins' , Spikes] , is.sws , th , 'R' , config , normalization);
                    reactivation.reward.dvHPC = [reactivation.reward.dvHPC ; R];
                    RBR = R; clear R
                end
                
                if sum(cond.dHPC.reward)>=1
                    [R] = reactivation_strength(patterns.all.reward , cond.dHPC.reward , [bins' , Spikes] , is.sws , th , 'R' , config , normalization);
                    reactivation.reward.dHPC = [reactivation.reward.dHPC ; R];
                    RDR = R; clear R
                end
                
                if sum(cond.vHPC.reward)>=1
                    [R] = reactivation_strength(patterns.all.reward , cond.vHPC.reward , [bins' , Spikes] , is.sws , th , 'R' , config , normalization);
                    reactivation.reward.vHPC = [reactivation.reward.vHPC ; R];
                    RVR = R; clear R
                end
                %                 end
            end
            
            %% Coordinated events triggered Assemblies activity
            if TPks
                disp('Triggered assemblies activity sourrounding cooridnated events')
                % Aversive
                if sum(RBA(:,end))>=1
                    condd = logical(RBA(:,end));
                    if sum(cond.both.aversive) >= 1
                        if aversiveTS_run(1) < rewardTS_run(1)
                            [R.baseline , ttttt] = triggered_CCG(patterns.all.aversive , cond.both.aversive , [bins' , Spikes], 2 , ripple_event.filtered.baseline , 5 , NREM.baseline);
                            [R.aversive , ttttt] = triggered_CCG(patterns.all.aversive , cond.both.aversive , [bins' , Spikes], 2 , ripple_event.filtered.aversive , 5 , NREM.aversive);
                            
                            if isempty(gain.both.aversive.pre)
                                gain.both.aversive.pre = [gain.both.aversive.pre , R.baseline(:,condd)];
                                gain.both.aversive.post = [gain.both.aversive.post , R.aversive(:,condd)];
                            else
                                gain.both.aversive.pre = [gain.both.aversive.pre , R.baseline(:,condd)];
                                gain.both.aversive.post = [gain.both.aversive.post , R.aversive(:,condd)];
                            end
                            
                            clear R
                        else
                            [R.reward , ttttt] = triggered_CCG(patterns.all.aversive , cond.both.aversive , [bins' , Spikes], 2 , ripple_event.filtered.reward , 5 , NREM.reward);
                            [R.aversive , ttttt] = triggered_CCG(patterns.all.aversive , cond.both.aversive , [bins' , Spikes], 2 , ripple_event.filtered.aversive , 5 , NREM.aversive);
                            
                            if isempty(gain.both.aversive.pre)
                                gain.both.aversive.pre = [gain.both.aversive.pre , R.reward(:,condd)];
                                gain.both.aversive.post = [gain.both.aversive.post , R.aversive(:,condd)];
                            else
                                gain.both.aversive.pre = [gain.both.aversive.pre , R.reward(:,condd)];
                                gain.both.aversive.post = [gain.both.aversive.post , R.aversive(:,condd)];
                            end
                            clear R
                        end
                    end
                    clear condd
                end
                
                % Reward
                if RBR(:,end)>=1
                    condd = logical(RBR(:,end));
                    if sum(cond.both.reward) >= 1
                        if aversiveTS_run(1) > rewardTS_run(1)
                            [R.baseline , ttttt] = triggered_CCG(patterns.all.reward , cond.both.reward , [bins' , Spikes], 2 , ripple_event.filtered.baseline , 5 , NREM.baseline);
                            [R.reward , ttttt] = triggered_CCG(patterns.all.reward , cond.both.reward , [bins' , Spikes], 2 , ripple_event.filtered.reward , 5 , NREM.reward);
                            
                            if isempty(gain.both.reward.pre)
                                gain.both.reward.pre = [gain.both.reward.pre , R.baseline(:,condd)];
                                gain.both.reward.post = [gain.both.reward.post , R.reward(:,condd)];
                            else
                                gain.both.reward.pre = [gain.both.reward.pre , R.baseline(:,condd)];
                                gain.both.reward.post = [gain.both.reward.post , R.reward(:,condd)];
                            end
                            clear R
                        else
                            [R.aversive , ttttt] = triggered_CCG(patterns.all.reward , cond.both.reward , [bins' , Spikes], 2 , ripple_event.filtered.aversive , 5 , NREM.aversive);
                            [R.reward , ttttt] = triggered_CCG(patterns.all.reward , cond.both.reward , [bins' , Spikes], 2 , ripple_event.filtered.reward , 5 , NREM.reward);
                            
                            if isempty(gain.both.reward.pre)
                                gain.both.reward.pre = [gain.both.reward.pre , R.aversive(:,condd)];
                                gain.both.reward.post = [gain.both.reward.post , R.reward(:,condd)];
                            else
                                gain.both.reward.pre = [gain.both.reward.pre , R.aversive(:,condd)];
                                gain.both.reward.post = [gain.both.reward.post , R.reward(:,condd)];
                            end
                            clear R
                        end
                    end
                    clear cond
                end
            end
            
            if TA
                disp('Triggered assemblies activity sourrounding cooridnated events')
                % Aversive
                if RBA(:,end)>=1
                    condd = logical(RBA(:,end));
                    if sum(cond.both.aversive) >= 1
                        if aversiveTS_run(1) < rewardTS_run(1)
                            [R.baseline] = triggered_activity(patterns.all.aversive , cond.both.aversive , [bins' , Spikes], 2 , ripple_event.filtered.baseline , 0 , 1);
                            [R.aversive] = triggered_activity(patterns.all.aversive , cond.both.aversive , [bins' , Spikes], 2 , ripple_event.filtered.aversive , 0 , 1);
                            
                            if isempty(gain.both.aversive.pre)
                                gain.both.aversive.pre = [gain.both.aversive.pre ; R.baseline(condd,:)];
                                gain.both.aversive.post = [gain.both.aversive.post ; R.aversive(condd,:)];
                            else
                                gain.both.aversive.pre = [gain.both.aversive.pre ; R.baseline(condd,1:size(gain.both.aversive.pre,2))];
                                gain.both.aversive.post = [gain.both.aversive.post ; R.aversive(condd,1:size(gain.both.aversive.pre,2))];
                            end
                            
                            clear R
                        else
                            [R.reward] = triggered_activity(patterns.all.aversive , cond.both.aversive , [bins' , Spikes], 2 , ripple_event.filtered.reward , 0 , 1);
                            [R.aversive] = triggered_activity(patterns.all.aversive , cond.both.aversive , [bins' , Spikes], 2 , ripple_event.filtered.aversive , 0 , 1);
                            
                            if isempty(gain.both.aversive.pre)
                                gain.both.aversive.pre = [gain.both.aversive.pre ; R.reward(condd,:)];
                                gain.both.aversive.post = [gain.both.aversive.post ; R.aversive(condd,:)];
                            else
                                gain.both.aversive.pre = [gain.both.aversive.pre ; R.reward(condd,1:size(gain.both.aversive.pre,2))];
                                gain.both.aversive.post = [gain.both.aversive.post ; R.aversive(condd,1:size(gain.both.aversive.pre,2))];
                            end
                            clear R
                        end
                    end
                    clear condd
                end
                
                % Reward
                if RBR(:,end)>=1
                    condd = logical(RBR(:,end));
                    if sum(cond.both.reward) >= 1
                        if aversiveTS_run(1) > rewardTS_run(1)
                            [R.baseline] = triggered_activity(patterns.all.reward , cond.both.reward , [bins' , Spikes], 2 , ripple_event.filtered.baseline , 0 , 1);
                            [R.reward] = triggered_activity(patterns.all.reward , cond.both.reward , [bins' , Spikes], 2 , ripple_event.filtered.reward , 0 , 1);
                            
                            if isempty(gain.both.reward.pre)
                                gain.both.reward.pre = [gain.both.reward.pre ; R.baseline(condd,:)];
                                gain.both.reward.post = [gain.both.reward.post ; R.reward(condd,:)];
                            else
                                gain.both.reward.pre = [gain.both.reward.pre ; R.baseline(condd,1:size(gain.both.reward.pre,2))];
                                gain.both.reward.post = [gain.both.reward.post ; R.reward(condd,1:size(gain.both.reward.pre,2))];
                            end
                            clear R
                        else
                            [R.aversive] = triggered_activity(patterns.all.reward , cond.both.reward , [bins' , Spikes], 2 , ripple_event.filtered.aversive , 0 , 1);
                            [R.reward] = triggered_activity(patterns.all.reward , cond.both.reward , [bins' , Spikes], 2 , ripple_event.filtered.reward , 0 , 1);
                            
                            if isempty(gain.both.reward.pre)
                                gain.both.reward.pre = [gain.both.reward.pre ; R.aversive(condd,:)];
                                gain.both.reward.post = [gain.both.reward.post ; R.reward(condd,:)];
                            else
                                gain.both.reward.pre = [gain.both.reward.pre ; R.aversive(condd,1:size(gain.both.reward.pre,2))];
                                gain.both.reward.post = [gain.both.reward.post ; R.reward(condd,1:size(gain.both.reward.pre,2))];
                            end
                            clear R
                        end
                    end
                    clear cond
                end
            end
            
            
            %% Recruitment of Assemblies during coordinated Events
            if REC
                disp('Recruitment of Assemblies during coordinated Events')
                % Aversive
                if sum(cond.both.aversive) >= 1
                    if aversiveTS_run(1) < rewardTS_run(1)
                        [R.baseline] = assembly_recruitment(patterns.all.aversive , cond.both.aversive , [bins' , Spikes], 0.2 , ripple_event.filtered.baseline , th , 1);
                        [R.aversive] = assembly_recruitment(patterns.all.aversive , cond.both.aversive , [bins' , Spikes], 0.2 , ripple_event.filtered.aversive , th , 1);
                        
                        gain.both.aversive.pre = [gain.both.aversive.pre ; R.baseline];
                        gain.both.aversive.post = [gain.both.aversive.post ; R.aversive];
                        
                        clear R
                    else
                        [R.reward] = assembly_recruitment(patterns.all.aversive , cond.both.aversive , [bins' , Spikes], 0.2 , ripple_event.filtered.reward , th , 1);
                        [R.aversive] = assembly_recruitment(patterns.all.aversive , cond.both.aversive , [bins' , Spikes], 0.2 , ripple_event.filtered.aversive , th , 1);
                        
                        gain.both.aversive.pre = [gain.both.aversive.pre ; R.reward];
                        gain.both.aversive.post = [gain.both.aversive.post ; R.aversive];
                        clear R
                    end
                end
                
                % Reward
                if sum(cond.both.reward) >= 1
                    if aversiveTS_run(1) > rewardTS_run(1)
                        [R.baseline] = assembly_recruitment(patterns.all.reward , cond.both.reward , [bins' , Spikes], 0.2 , ripple_event.filtered.baseline , th , 1);
                        [R.reward] = assembly_recruitment(patterns.all.reward , cond.both.reward , [bins' , Spikes], 0.2 , ripple_event.filtered.reward , th , 1);
                        
                        gain.both.reward.pre = [gain.both.reward.pre ; R.baseline];
                        gain.both.reward.post = [gain.both.reward.post ; R.reward];
                        
                        clear R
                    else
                        [R.aversive] = assembly_recruitment(patterns.all.reward , cond.both.reward , [bins' , Spikes], 0.2 , ripple_event.filtered.aversive , th , 1);
                        [R.reward] = assembly_recruitment(patterns.all.reward , cond.both.reward , [bins' , Spikes], 0.2 , ripple_event.filtered.reward , th , 1);
                        
                        gain.both.aversive.pre = [gain.both.reward.pre ; R.aversive];
                        gain.both.aversive.post = [gain.both.reward.post ; R.reward];
                        clear R
                    end
                end
            end
            
            %% Temporality of the cooridnated events
            if C
                disp('Cross-Correlograms of coordinated ripple events')
                count = 0;
                x = Restrict(ripples.dHPC.coordinated.baseline(:,2),MUA);
                y = Restrict(ripples.vHPC.coordinated.baseline(:,2),MUA);
                if and(length(x) > 10 , length(y) > 10)
                    [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                    [ccg.b,tttt] = CCG(s,ids,'binSize',0.005,'duration',0.4,'smooth',1,'mode','ccv'); %ccg calculation
                    %                 plot(tttt,ccg(:,1,2)./sum(ccg(:,1,2)),'k'),hold on
                    count = count+1;
                end
                clear x y s ids groups
                
                x = Restrict(ripples.dHPC.coordinated.reward(:,2),MUA);
                y = Restrict(ripples.vHPC.coordinated.reward(:,2),MUA);
                if and(length(x) > 10 , length(y) > 10)
                    [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                    [ccg.r,tttt] = CCG(s,ids,'binSize',0.005,'duration',0.4,'smooth',1,'mode','ccv'); %ccg calculation
                    %                 plot(tttt,ccg(:,1,2)./sum(ccg(:,1,2)),'b')
                    count = count+1;
                end
                clear x y s ids groups
                
                
                x = Restrict(ripples.dHPC.coordinated.aversive(:,2),MUA);
                y = Restrict(ripples.vHPC.coordinated.aversive(:,2),MUA);
                if and(length(x) > 10 , length(y) > 10)
                    
                    [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                    [ccg.a,tttt] = CCG(s,ids,'binSize',0.005,'duration',0.4,'smooth',1,'mode','ccv'); %ccg calculation
                    %                 plot(tttt,ccg(:,1,2)./sum(ccg(:,1,2)),'r')
                    count = count+1;
                end
                clear x y s ids groups
                
                if count ==3
                    if not(exist('crosscorrelograms','var'))
                        crosscorrelograms.baseline = ccg.b(:,1,2);
                        crosscorrelograms.aversive = ccg.a(:,1,2);
                        crosscorrelograms.reward = ccg.r(:,1,2);
                        clear ccg
                    else
                        crosscorrelograms.baseline = [crosscorrelograms.baseline , ccg.b(:,1,2)];
                        crosscorrelograms.aversive = [crosscorrelograms.aversive , ccg.a(:,1,2)];
                        crosscorrelograms.reward = [crosscorrelograms.reward , ccg.r(:,1,2)];
                        clear ccg
                    end
                end
                clear count
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


%% cleaning data
% joint
x = or(isnan(reactivation.reward.dvHPC(:,2)) , isnan(reactivation.reward.dvHPC(:,3)));
reactivation.reward.dvHPC(x,:) = [];

x = or(isnan(reactivation.aversive.dvHPC(:,2)) , isnan(reactivation.aversive.dvHPC(:,3)));
reactivation.aversive.dvHPC(x,:) = [];

% dHPC
x = or(isnan(reactivation.reward.dHPC(:,2)) , isnan(reactivation.reward.dHPC(:,3)));
reactivation.reward.dHPC(x,:) = [];

x = or(isnan(reactivation.aversive.dHPC(:,2)) , isnan(reactivation.aversive.dHPC(:,3)));
reactivation.aversive.dHPC(x,:) = [];

% vHPC
x = or(isnan(reactivation.reward.vHPC(:,2)) , isnan(reactivation.reward.vHPC(:,3)));
reactivation.reward.vHPC(x,:) = [];

x = or(isnan(reactivation.aversive.vHPC(:,2)) , isnan(reactivation.aversive.vHPC(:,3)));
reactivation.aversive.vHPC(x,:) = [];

%% Plot Strenght Reactivation
% for joint assemblies

figure
x = logical(reactivation.reward.dvHPC(:,6));
y = logical(reactivation.aversive.dvHPC(:,6));

reactivation.reward.dvHPC(isnan(reactivation.reward.dvHPC(:,1)),:) = [];
reactivation.aversive.dvHPC(isnan(reactivation.aversive.dvHPC(:,1)),:) = [];
x = reactivation.reward.dvHPC(:,1);
y = reactivation.aversive.dvHPC(:,1);

kstest(x)
kstest(y)
[h, p] = ttest2(x,y)  
[h, p] = ttest(y)
[h, p] = ttest(x)


xx = [1 2];
yy = [mean(x) mean(y)];
err = [(std(x)/sqrt(length(x))) (std(y)/sqrt(length(y)))];

subplot(131),
bar(xx,yy),hold on
er = errorbar(xx,yy,err);ylim([-0.07 0.07])
er.Color = [0 0 0];                            
er.LineStyle = 'none';
hold off

% for dHPC assemblies
x = logical(reactivation.reward.dHPC(:,6));
y = logical(reactivation.aversive.dHPC(:,6));
reactivation.reward.dHPC(isnan(reactivation.reward.dHPC(:,1)),:) = [];
reactivation.aversive.dHPC(isnan(reactivation.aversive.dHPC(:,1)),:) = [];
x = reactivation.reward.dHPC(:,1);
y = reactivation.aversive.dHPC(:,1);
kstest(x)
kstest(y)
[h, p] = ttest2(x,y)  
[h, p] = ttest(y)
[h, p] = ttest(x)

xx = [1 2];
yy = [mean(x) mean(y)];
err = [(std(x)/sqrt(length(x))) (std(y)/sqrt(length(y)))];

subplot(132),
bar(xx,yy),hold on
er = errorbar(xx,yy,err);ylim([-0.07 0.07])
er.Color = [0 0 0];                            
er.LineStyle = 'none';
hold off

% for vHPC assemblies
x = logical(reactivation.reward.vHPC(:,6));
y = logical(reactivation.aversive.vHPC(:,6));
reactivation.reward.vHPC(isnan(reactivation.reward.vHPC(:,1)),:) = [];
reactivation.aversive.vHPC(isnan(reactivation.aversive.vHPC(:,1)),:) = [];
x = reactivation.reward.vHPC(:,1);
y = reactivation.aversive.vHPC(:,1);
kstest(x)
kstest(y)
[h, p] = ttest2(x,y)  
[h, p] = ttest(y)
[h, p] = ttest(x)

xx = [1 2];
yy = [mean(x) mean(y)];
err = [(std(x)/sqrt(length(x))) (std(y)/sqrt(length(y)))];

subplot(133),
bar(xx,yy),hold on
er = errorbar(xx,yy,err);ylim([-0.07 0.07])
er.Color = [0 0 0];                            
er.LineStyle = 'none';
hold off

%%
percentage.dvHPC.aversive = sum(reactivation.aversive.dvHPC(:,6))/length(reactivation.aversive.dvHPC(:,4));
percentage.dvHPC.reward = sum(reactivation.reward.dvHPC(:,6))/length(reactivation.reward.dvHPC(:,4));

percentage.dHPC.aversive = sum(reactivation.aversive.dHPC(:,6))/length(reactivation.aversive.dHPC(:,4));
percentage.dHPC.reward = sum(reactivation.reward.dHPC(:,6))/length(reactivation.reward.dHPC(:,4));

percentage.vHPC.aversive = sum(reactivation.aversive.vHPC(:,6))/length(reactivation.aversive.vHPC(:,4));
percentage.vHPC.reward = sum(reactivation.reward.vHPC(:,6))/length(reactivation.reward.vHPC(:,4));


%% boxplot
figure
reactivation.reward.dvHPC(isnan(reactivation.reward.dvHPC(:,1)),:) = [];
reactivation.aversive.dvHPC(isnan(reactivation.aversive.dvHPC(:,1)),:) = [];
x = reactivation.reward.dvHPC(:,1);
y = reactivation.aversive.dvHPC(:,1);

[h, p] = ttest2(x,y)  
[h, p] = ttest(y)
[h, p] = ttest(x)
kstest(x)
kstest(y)


subplot(131),
scatter([ones(length(x),1) ; ones(length(y),1)*2],[x;y],'filled'),xlim([0 3])%,ylim([-0.5 0.5])
hold on
scatter([1 2],[mean(x);mean(y)],'filled'),xlim([0 3])%,ylim([-0.5 0.5])
yline(0,'--')
% boxplot([x;y],[ones(length(x),1) ; ones(length(y),1)*2]),xlim([0 3]),ylim([-0.5 0.5])

% for dHPC assemblies
reactivation.reward.dHPC(isnan(reactivation.reward.dHPC(:,1)),:) = [];
reactivation.aversive.dHPC(isnan(reactivation.aversive.dHPC(:,1)),:) = [];
x = reactivation.reward.dHPC(:,1);
y = reactivation.aversive.dHPC(:,1);
[h, p] = ttest2(x,y)             
[h, p] = ttest(x)
[h, p] = ttest(y)

subplot(132),
scatter([ones(length(x),1) ; ones(length(y),1)*2],[x;y],'filled'),xlim([0 3]),ylim([-10 10])
hold on
scatter([1 2],[mean(x);mean(y)],'filled'),xlim([0 3]),ylim([-0.5 0.5])
% boxplot([x;y],[ones(length(x),1) ; ones(length(y),1)*2]),xlim([0 3]),ylim([-0.5 0.5])
hold on
yline(0,'--')

% for vHPC assemblies
reactivation.reward.vHPC(isnan(reactivation.reward.vHPC(:,1)),:) = [];
reactivation.aversive.vHPC(isnan(reactivation.aversive.vHPC(:,1)),:) = [];
x = reactivation.reward.vHPC(:,1);
y = reactivation.aversive.vHPC(:,1);
[h, p] = ttest2(x,y)             
[h, p] = ttest(x)
[h, p] = ttest(y)

subplot(133),
scatter([ones(length(x),1) ; ones(length(y),1)*2],[x;y],'filled'),xlim([0 3]),ylim([-10 10])
hold on
scatter([1 2],[mean(x);mean(y)],'filled'),xlim([0 3]),ylim([-0.5 0.5])
yline(0,'--')
% boxplot([x;y],[ones(length(x),1) ; ones(length(y),1)*2]),xlim([0 3]),ylim([-0.5 0.5])


%% violin plot
reactivation.reward.dvHPC(isnan(reactivation.reward.dvHPC(:,1)),:) = [];
reactivation.aversive.dvHPC(isnan(reactivation.aversive.dvHPC(:,1)),:) = [];
x = reactivation.reward.dvHPC(:,1);
y = reactivation.aversive.dvHPC(:,1);

m = max(size(x,1),size(y,1));
data_box = nan(m,2);
data_box(1:size(x,1),1) = x;
data_box(1:size(y,1),2) = y;

figure
subplot(131)
label = {'Reward', 'Aversive'}; 
[h,L,MX,MED,bw] = violin(data_box, 'xlabel', label,'facecolor',[0 0 0;1 0 0]);
ylabel('Reactivation Strength','FontSize',14);
title('Joint');


reactivation.reward.dHPC(isnan(reactivation.reward.dHPC(:,1)),:) = [];
reactivation.aversive.dHPC(isnan(reactivation.aversive.dHPC(:,1)),:) = [];
x = reactivation.reward.dHPC(:,1);
y = reactivation.aversive.dHPC(:,1);

m = max(size(x,1),size(y,1));
data_box = nan(m,2);
data_box(1:size(x,1),1) = x;
data_box(1:size(y,1),2) = y;

subplot(132)
label = {'Reward', 'Aversive'}; 
[h,L,MX,MED,bw] = violin(data_box, 'xlabel', label,'facecolor',[0 0 0;1 0 0]);
ylabel('Reactivation Strength','FontSize',14);
title('dHPC');


reactivation.reward.vHPC(isnan(reactivation.reward.vHPC(:,1)),:) = [];
reactivation.aversive.vHPC(isnan(reactivation.aversive.vHPC(:,1)),:) = [];
x = reactivation.reward.vHPC(:,1);
y = reactivation.aversive.vHPC(:,1);

m = max(size(x,1),size(y,1));
data_box = nan(m,2);
data_box(1:size(x,1),1) = x;
data_box(1:size(y,1),2) = y;

subplot(133)
label = {'Reward', 'Aversive'}; 
[h,L,MX,MED,bw] = violin(data_box, 'xlabel', label,'facecolor',[0 0 0;1 0 0]);
ylabel('Reactivation Strength','FontSize',14);
title('vHPC');

%%

figure
scatter(reactivation.reward.dvHPC(:,2) , reactivation.reward.dvHPC(:,1),'filled','MarkerFaceAlpha',0.5),hold on
scatter(reactivation.aversive.dvHPC(:,2) , reactivation.aversive.dvHPC(:,1),'filled','MarkerFaceAlpha',0.5,'MarkerFaceColor','r'),hold on

xx = or(isnan(reactivation.aversive.dvHPC(:,2)) , isnan(reactivation.aversive.dvHPC(:,1)));
x = reactivation.aversive.dvHPC(not(xx),2);
y = reactivation.aversive.dvHPC(not(xx),1);

xx = or(isnan(reactivation.reward.dvHPC(:,2)) , isnan(reactivation.reward.dvHPC(:,1)));
x = reactivation.reward.dvHPC(not(xx),2);
y = reactivation.reward.dvHPC(not(xx),1);
corr(x,y)


%% Plotting the number of peaks
x = reactivation.reward.dvHPC(:,2)./reactivation.reward.dvHPC(:,3);
y = reactivation.aversive.dvHPC(:,2)./reactivation.aversive.dvHPC(:,3);
[h, p] = ttest2(x,y)             

xx = [1 2];
yy = [nanmean(x) nanmean(y)];
err = [(nanstd(x)/sqrt(sum(not(isnan(x))))) (nanstd(y)/sqrt(sum(not(isnan(y)))))];

subplot(131),
bar(xx,yy),hold on
er = errorbar(xx,yy,err);ylim([0 3])
er.Color = [0 0 0];                            
er.LineStyle = 'none';
hold off
ylabel('#PeaksPost / #PeaksPre')

% for dHPC assemblies
x = reactivation.reward.dHPC(:,2)./reactivation.reward.dHPC(:,3);
y = reactivation.aversive.dHPC(:,2)./reactivation.aversive.dHPC(:,3);
[h, p] = ttest2(x,y)             

xx = [1 2];
yy = [nanmean(x) nanmean(y)];
err = [(nanstd(x)/sqrt(sum(not(isnan(x))))) (nanstd(y)/sqrt(sum(not(isnan(y)))))];

subplot(132),
bar(xx,yy),hold on
er = errorbar(xx,yy,err);ylim([0 3])
er.Color = [0 0 0];                            
er.LineStyle = 'none';
hold off

% for vHPC assemblies
x = reactivation.reward.vHPC(:,2)./reactivation.reward.vHPC(:,3);
y = reactivation.aversive.vHPC(:,2)./reactivation.aversive.vHPC(:,3);
[h, p] = ttest2(x,y)             

xx = [1 2];
yy = [nanmean(x) nanmean(y)];
err = [(nanstd(x)/sqrt(sum(not(isnan(x))))) (nanstd(y)/sqrt(sum(not(isnan(y)))))];

subplot(133),
bar(xx,yy),hold on
er = errorbar(xx,yy,err);ylim([0 3])
er.Color = [0 0 0];                            
er.LineStyle = 'none';
hold off

%% Plotting the mean of the activity strengtht
x = reactivation.reward.dvHPC(:,2)-reactivation.reward.dvHPC(:,3);
y = reactivation.aversive.dvHPC(:,2)-reactivation.aversive.dvHPC(:,3);
[h, p] = ttest2(x,y)             

xx = [1 2];
yy = [nanmean(x) nanmean(y)];
err = [(nanstd(x)/sqrt(sum(not(isnan(x))))) (nanstd(y)/sqrt(sum(not(isnan(y)))))];

subplot(131),
bar(xx,yy),hold on
er = errorbar(xx,yy,err);ylim([-1 1])
er.Color = [0 0 0];                            
er.LineStyle = 'none';
hold off
ylabel('#PeaksPost / #PeaksPre')

% for dHPC assemblies
x = reactivation.reward.dHPC(:,4)-reactivation.reward.dHPC(:,5);
y = reactivation.aversive.dHPC(:,4)-reactivation.aversive.dHPC(:,5);
[h, p] = ttest2(x,y)             

xx = [1 2];
yy = [mean(x) mean(y)];
err = [(std(x)/sqrt(length(x))) (std(y)/sqrt(length(y)))];

subplot(132),
bar(xx,yy),hold on
er = errorbar(xx,yy,err);ylim([-1 1])
er.Color = [0 0 0];                            
er.LineStyle = 'none';
hold off

% for vHPC assemblies
x = reactivation.reward.vHPC(:,4)-reactivation.reward.vHPC(:,5);
y = reactivation.aversive.vHPC(:,4)-reactivation.aversive.vHPC(:,5);
[h, p] = ttest2(x,y)             

xx = [1 2];
yy = [mean(x) mean(y)];
err = [(std(x)/sqrt(length(x))) (std(y)/sqrt(length(y)))];

subplot(133),
bar(xx,yy),hold on
er = errorbar(xx,yy,err);ylim([-1 1])
er.Color = [0 0 0];                            
er.LineStyle = 'none';
hold off

%% 
[i ii] = max(gain.both.aversive.post);
[i , iii] = sort(ii,'descend');

figure,
subplot(131),imagesc(ttttt,[1:1:size(gain.both.aversive.pre,2)],gain.both.aversive.pre(:,iii)')%,xlim([-2 2]),caxis([-2 2]),colormap 'jet'
subplot(132),imagesc(ttttt,[1:1:size(gain.both.aversive.post,2)],gain.both.aversive.post(:,iii)')%,xlim([-2 2]),caxis([-2 2]),colormap 'jet'
clear i ii iii

subplot(133),plot(ttttt,mean(gain.both.aversive.pre,2,'omitnan'),'k','LineWidth',2),hold on
sem = std(gain.both.aversive.pre,0,2,'omitnan')/sqrt(size(gain.both.aversive.pre,2));
ciplot(mean(gain.both.aversive.pre,2,'omitnan')-sem , mean(gain.both.aversive.pre,2,'omitnan')+sem,ttttt,'k'),alpha 0.5
plot(ttttt,mean(gain.both.aversive.post,2,'omitnan'),'r','LineWidth',1),hold on
sem = std(gain.both.aversive.post,0,2,'omitnan')/sqrt(size(gain.both.aversive.post,2));
ciplot(mean(gain.both.aversive.post,2,'omitnan')-sem , mean(gain.both.aversive.post,2,'omitnan')+sem,ttttt,'r'),alpha 0.5,xlim([-0.5 0.5])


% 
[i ii] = max(gain.both.reward.post);
[i , iii] = sort(ii,'descend');
figure,
subplot(131),imagesc(ttttt,[1:1:size(gain.both.reward.post,2)],gain.both.reward.pre(:,iii)')%,xlim([-0.2 0.2]),caxis([-3 3]),colormap 'jet'
subplot(132),imagesc(ttttt,[1:1:size(gain.both.reward.post,2)],gain.both.reward.post(:,iii)')%,xlim([-0.2 0.2]),caxis([-3 3]),colormap 'jet'
clear i ii iii



subplot(133),plot(ttttt,mean(gain.both.reward.pre,2,'omitnan'),'k','LineWidth',2),hold on
sem = std(gain.both.reward.pre,0,2,'omitnan')/sqrt(size(gain.both.reward.pre,2));
ciplot(mean(gain.both.reward.pre,2,'omitnan')-sem , mean(gain.both.reward.pre,2,'omitnan')+sem,ttttt,'k'),alpha 0.5
plot(ttttt,mean(gain.both.reward.post,2,'omitnan'),'r','LineWidth',2),hold on
sem = std(gain.both.reward.post,0,2,'omitnan')/sqrt(size(gain.both.reward.post,2));
ciplot(mean(gain.both.reward.post,2,'omitnan')-sem , mean(gain.both.reward.post,2,'omitnan')+sem,ttttt,'r'),alpha 0.5,%ylim([-0.3 0.6]),xlim([-0.2 0.2])

%% PHIST sourrounding the shock
d = ResponsiveS(:,1)==1;
v = ResponsiveS(:,1)==2;

d = curveShock(:,d);
v = curveShock(:,v);

[h , p] = min(abs(binsS - 0));
[h , pp] = min(abs(binsS - 1));

[i ii] = sort(mean(d(p:pp,:),1),'descend');
subplot(121),imagesc(binsS , [1:1:size(d,2)] , d(:,ii)'),clim([0 6]),colormap 'jet',hold on
xline(0,'--'),xline(1,'--')
title('dHPC Pyr')

[i ii] = sort(mean(v(p:pp,:),1),'descend');
subplot(122),imagesc(binsS , [1:1:size(v,2)] , v(:,ii)'),clim([0 6]),colormap 'jet',hold on
xline(0,'--'),xline(1,'--')
title('vHPC Pyr')

r = [];
[h pp] = min(abs(binsS -0));
[h ppp] = min(abs(binsS -1));

for i = 1 : size(curveShock,2)
    if ResponsiveS(i,1) == 1
        if mean(curveShock(pp : ppp , i)) >= 2
            r = [r ; 1 1];
        else
            r = [r ; 1 0];
        end
    else
        if mean(curveShock(pp : ppp , i)) >= 2
            r = [r ; 2 1];
        else
            r = [r ; 2 0];
        end
    end
end

d = and(r(:,1)==1 , r(:,2)==1);
d = curveShock(:,d);
m = mean(d,2,'omitnan');
s = std(d,0,2,'omitnan');
s = s./sqrt(sum(ResponsiveS(:,1)==1));
subplot(121), plot(binsS,m','k','LineWidth',2),hold on
ciplot(m-s , m+s , binsS,'k')
alpha 0.5

d = and(r(:,1)==1 , r(:,2)==0);
d = curveShock(:,d);
m = mean(d,2,'omitnan');
s = std(d,0,2,'omitnan');
s = s./sqrt(sum(ResponsiveS(:,1)==1));
plot(binsS,m','k','LineWidth',2),hold on
ciplot(m-s , m+s , binsS,'k')
alpha 0.5
ylim([0 30])



d = and(r(:,1)==2 , r(:,2)==1);
d = curveShock(:,d);
m = mean(d,2,'omitnan');
s = std(d,0,2,'omitnan');
s = s./sqrt(sum(ResponsiveS(:,1)==1));
subplot(122), plot(binsS,m','k','LineWidth',2),hold on
ciplot(m-s , m+s , binsS,'k')
alpha 0.5


d = and(r(:,1)==2 , r(:,2)==0);
d = curveShock(:,d);
m = mean(d,2,'omitnan');
s = std(d,0,2,'omitnan');
s = s./sqrt(sum(ResponsiveS(:,1)==1));
plot(binsS,m','k','LineWidth',2),hold on
ciplot(m-s , m+s , binsS,'k')
alpha 0.5
xlim([-3 3])
ylim([0 30])
