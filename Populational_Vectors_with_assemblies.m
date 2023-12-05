clear
clc
close all

%% Parameters
path = {'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat103\usable';'E:\Rat132\recordings\in_pyr'};%List of folders from the path

%Sleep
time_criteria = 1; % minimal time to include a NREM epoch (in min)

% for SU
criteria_fr = 0.1; %criteria to include or not a SU into the analysis
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

FiringRate.dHPC.baseline = [];
FiringRate.dHPC.reward = [];
FiringRate.dHPC.aversive = [];
FiringRate.vHPC.baseline = [];
FiringRate.vHPC.reward = [];
FiringRate.vHPC.aversive = [];

% to save data
output.aversive.first.aversive = [];
output.aversive.second.aversive = [];
output.aversive.first.reward = [];
output.aversive.second.reward = [];


counts.both = [];
counts.dHPC = [];
counts.vHPC = [];

reactivation.aversive.dvHPC = [];
reactivation.reward.dvHPC = [];
reactivation.aversive.dHPC = [];
reactivation.reward.dHPC = [];
reactivation.aversive.vHPC = [];
reactivation.reward.vHPC = [];

th = 10;
for_scatter.aversive = [];
for_scatter.reward = [];
Rsquared.aversive = [];
Rsquared.reward = [];

Act_aversive = [];
Act_reward = [];

gain.both.reward.pre = [];
gain.both.reward.post = [];
gain.both.aversive.pre = [];
gain.both.aversive.post = [];

percentages = [];

store.RA = [];
store.AR = [];

Cross.aversive.pre = []; Cross.aversive.post = [];
Cross.reward.pre = []; Cross.reward.post = [];

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
        REM.all(REM.all(:,2)-REM.all(:,1)<60,:) = [];
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
        
        cooridnated_event((cooridnated_event(:,2)-cooridnated_event(:,1)<0.04),:) = [];
        ripple_event.baseline = Restrict(cooridnated_event,baselineTS./1000);
        ripple_event.reward = Restrict(cooridnated_event,rewardTS./1000);
        ripple_event.aversive = Restrict(cooridnated_event,aversiveTS./1000);
        ripple_event.all = cooridnated_event;
                
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
                numberD = numberD+1;
                clusters.all = [clusters.all ; cluster];
                clusters.dHPC = [clusters.dHPC ; cluster];
            end
            clear tmp cluster Cellulartype fr1 fr2 fr3 fr4 fr5
        end
        
        numberV = 0;
        clusters.vHPC = [];
        for ii=1:size(group_vHPC,1)
            cluster = group_vHPC(ii,1);
            Cellulartype = logical(cellulartype(cellulartype(:,1) == cluster,2));
            if Cellulartype
                numberV = numberV+1;
                clusters.all = [clusters.all ; cluster];
                clusters.vHPC = [clusters.vHPC ; cluster];
            end
            clear tmp cluster Cellulartype fr1 fr2 fr3 fr4 fr5
        end
        clear freq limits
        clear camara shock rightvalve leftvalve
        clear ejeX ejeY dX dY dX_int dY_int
        
        %% Assemblies detection
        if and(numberD >= criteria_n(1),numberV >= criteria_n(2))
            disp('Lets go for the SUs')
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
            MUA = Smooth(MUA , 30 ,'kernel','gaussian'); % smooth
            MUA = MUA>mean(MUA)+std(MUA)*3; % detection of high MUA
            MUA = ToIntervals(b',MUA); % definition of periods
            MUA = merge_events(MUA, 0.04); % merge periods with short IEI
            
            a = ((ripple_event.all(:,2) - ripple_event.all(:,1))/2) + ripple_event.all(:,1); % center of the event
            a = Restrict(a , MUA); % restriction of events to MUA periods
            ripple_event.filtered.all = Restrict(a,NREM.all); % store all events
            ripple_event.filtered.baseline = Restrict(a,NREM.baseline); % store baseline events
            ripple_event.filtered.reward = Restrict(a,NREM.reward); % store reward events
            ripple_event.filtered.aversive = Restrict(a,NREM.aversive); % store aversive events
            clear a MUA b  spks
            
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
            load('dorsalventral_assemblies_aversive.mat')
            
            %             disp('Detection of assemblies using Aversive template')
            %             limits = aversiveTS_run./1000;
            %             events = [];
            %             events = movement.aversive;
            %             [SpksTrains.all.aversive , Bins.aversive , Cluster.all.aversive] = spike_train_construction([spks_dHPC;spks_vHPC], clusters.all, cellulartype, binSize, limits, events, false,false);
            %             [Th , pat] = assembly_patternsJFM([SpksTrains.all.aversive'],opts);
            %             save([cd,'\dorsalventral_assemblies_aversive.mat'],'Th' , 'pat')
            
            Thresholded.aversive.all = Th;
            patterns.all.aversive = pat;
            clear cond Th pat
            
            % Detection of members
            cond1 =  sum(Thresholded.aversive.all(1:size(clusters.dHPC,1),:),1)>0; %checking of dHPC SU
            cond2 =  sum(Thresholded.aversive.all(size(clusters.dHPC,1)+1:end,:),1)>0; %checking of vHPC SU
            cond.dHPC.aversive = and(cond1 , not(cond2));
            cond.vHPC.aversive = and(cond2 , not(cond1));
            cond.both.aversive = and(cond1 , cond2); clear cond1 cond2
            
            % --- Reward ---
            disp('Loading Reward template')
            load('dorsalventral_assemblies_reward.mat')
            
            %             disp('Detection of assemblies using Rewarded template')
            %             limits = rewardTS_run./1000;
            %             events = [];
            %             events = movement.reward;
            %             [SpksTrains.all.reward , Bins.reward , Cluster.all.reward] = spike_train_construction([spks_dHPC;spks_vHPC], clusters.all, cellulartype, binSize, limits, events, false,false);
            %             [Th , pat] = assembly_patternsJFM([SpksTrains.all.reward'],opts);
            %             save([cd,'\dorsalventral_assemblies_reward.mat'],'Th' , 'pat')
            
            Thresholded.reward.all = Th;
            patterns.all.reward = pat;
            clear Th pat
            
            % Detection of members using
            cond1 =  sum(Thresholded.reward.all(1:size(clusters.dHPC,1),:),1)>0; %checking of dHPC SU
            cond2 =  sum(Thresholded.reward.all(size(clusters.dHPC,1)+1:end,:),1)>0; %checking of vHPC SU
            cond.dHPC.reward = and(cond1 , not(cond2));
            cond.vHPC.reward = and(cond2 , not(cond1));
            cond.both.reward = and(cond1 , cond2); clear cond1 cond2
            
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
            
            %%
            a = sum(sum(and(cond.dHPC.aversive .* cond.both.reward' , p.AR)));
            b = sum(sum(and(cond.vHPC.aversive .* cond.both.reward' , p.AR)));
            c = sum(sum(and(cond.both.aversive .* cond.dHPC.reward' , p.AR)));
            d = sum(sum(and(cond.both.aversive .* cond.vHPC.reward' , p.AR)));
            store.AR = [store.AR ; a , b , c , d];
            clear a b c d
            
            a = sum(sum(and(cond.dHPC.reward' .* cond.both.aversive , p.AR)));
            b = sum(sum(and(cond.vHPC.reward' .* cond.both.aversive , p.AR)));
            c = sum(sum(and(cond.both.reward' .* cond.dHPC.aversive , p.AR)));
            d = sum(sum(and(cond.both.reward' .* cond.vHPC.aversive , p.AR)));
            store.RA = [store.RA ; a , b , c , d];
            clear a b c d
            
            %% SpikeTrains construction
            limits = [0 segments.Var1(end)/1000];
            events = [];
            [Spikes , bins , Clusters] = spike_train_construction([spks_dHPC;spks_vHPC], clusters.all, cellulartype, binSize, limits, events, true, true);
            clear limits events
            
            %% Assemblies activation in the entier recording
%             is.sws.baseline = InIntervals(bins,NREM.baseline);
%             is.sws.reward = InIntervals(bins,NREM.reward);
%             is.sws.aversive = InIntervals(bins,NREM.aversive);
            
            is.sws.baseline = InIntervals(bins,[ripple_event.filtered.baseline(:,1)-0.2 ripple_event.filtered.baseline(:,1)+0.2]);
            is.sws.reward = InIntervals(bins,[ripple_event.filtered.reward(:,1)-0.2 ripple_event.filtered.reward(:,1)+0.2]);
            is.sws.aversive = InIntervals(bins,[ripple_event.filtered.aversive(:,1)-0.2 ripple_event.filtered.aversive(:,1)+0.2]);
            
            is.burst.baseline = InIntervals(bins,ripple_bursts.baseline);
            is.burst.reward = InIntervals(bins,ripple_bursts.reward);
            is.burst.aversive = InIntervals(bins,ripple_bursts.aversive);
            
            is.rem.baseline = InIntervals(bins,REM.baseline);
            is.rem.reward = InIntervals(bins,REM.reward);
            is.rem.aversive = InIntervals(bins,REM.aversive);
            
            is.coordinated.baseline = InIntervals(bins,ripple_event.baseline);
            is.coordinated.reward = InIntervals(bins,ripple_event.reward);
            is.coordinated.aversive = InIntervals(bins,ripple_event.aversive);
            
            is.run.aversive = InIntervals(bins,movement.aversive);
            is.run.reward = InIntervals(bins,movement.reward);
            
            is.aversive = InIntervals(bins,aversiveTS_run./1000);
            is.reward = InIntervals(bins,rewardTS_run./1000);
            
            %% Reactivation Strenght
            % --- for joint assemblies ---
            if sum(cond.both.aversive)>=1
                a = assembly_activity(patterns.all.aversive(:,cond.both.aversive) , Spikes');
                aa = a;
                a = zscore(a,1,2);
                for_reactivation.aversive = [];
                AssemblyTime.joint.aversive.baseline = cell(size(a,1),1);
                AssemblyTime.joint.aversive.reward = cell(size(a,1),1);
                AssemblyTime.joint.aversive.aversive = cell(size(a,1),1);
                for i = 1:size(a,1)
                    [pks.baseline,loc.baseline] = findpeaks(a(i,is.sws.baseline),bins(is.sws.baseline),'MinPeakHeight',th);
                    [pks.reward,loc.reward] = findpeaks(a(i,is.sws.reward),bins(is.sws.reward),'MinPeakHeight',th);
                    [pks.aversive,loc.aversive] = findpeaks(a(i,is.sws.aversive),bins(is.sws.aversive),'MinPeakHeight',th);
                    
                    %                     [pks.baseline,loc.baseline] = findpeaks(a(i,is.burst.baseline),bins(is.burst.baseline),'MinPeakHeight',th);
                    %                     [pks.reward,loc.reward] = findpeaks(a(i,is.burst.reward),bins(is.burst.reward),'MinPeakHeight',th);
                    %                     [pks.aversive,loc.aversive] = findpeaks(a(i,is.burst.aversive),bins(is.burst.aversive),'MinPeakHeight',th);
                    
                    [pks.run,loc.run] = findpeaks(a(i,is.run.aversive),bins(is.run.aversive),'MinPeakHeight',th);
                    
                    [pks.a,loc.a] = findpeaks(a(i,is.run.aversive),bins(is.run.aversive),'MinPeakHeight',th);
                    [pks.r,loc.r] = findpeaks(a(i,is.run.reward),bins(is.run.reward),'MinPeakHeight',th);
                    
                    Rsquared.aversive = [Rsquared.aversive ; mean(pks.a,'omitnan')  mean(pks.r,'omitnan') length(Shocks_filt) length(Rewards_filt)];
                    
                    if aversiveTS_run(1) < rewardTS_run(1)
                        tmp = [sum(aa(i,is.sws.aversive)) , sum(aa(i,is.sws.baseline))];
                        reactivation.aversive.dvHPC = [reactivation.aversive.dvHPC ; (mean(pks.aversive) - mean(pks.baseline)) length(pks.aversive) length(pks.baseline) , tmp]; clear tmp
                        for_scatter.aversive = [for_scatter.aversive ; mean(pks.run) (mean(pks.aversive) - mean(pks.baseline)) length(Shocks_filt)];
                        for_reactivation.aversive = [for_reactivation.aversive ; mean(pks.baseline)];
                    else
                        tmp = [sum(aa(i,is.sws.aversive)) , sum(aa(i,is.sws.reward))];
                        reactivation.aversive.dvHPC = [reactivation.aversive.dvHPC ; (mean(pks.aversive) - mean(pks.reward)) length(pks.aversive) length(pks.reward) , tmp]; clear tmp
                        for_scatter.aversive = [for_scatter.aversive ; mean(pks.run) (mean(pks.aversive) - mean(pks.reward)) length(Shocks_filt)];
                        for_reactivation.aversive = [for_reactivation.aversive ; mean(pks.reward)];
                    end
                        AssemblyTime.joint.aversive.baseline{i} = loc.baseline;
                        AssemblyTime.joint.aversive.reward{i} = loc.reward;
                        AssemblyTime.joint.aversive.aversive{i} = loc.aversive;
                        clear pks loc
                end
                clear a MUA aa
            end
            
            if sum(cond.both.reward)>=1
                a = assembly_activity(patterns.all.reward(:,cond.both.reward) , Spikes');
                aa = a;
                a = zscore(a,1,2);
                for_reactivation.reward = [];
                AssemblyTime.joint.reward.baseline = cell(size(a,1),1);
                AssemblyTime.joint.reward.reward = cell(size(a,1),1);
                AssemblyTime.joint.reward.aversive = cell(size(a,1),1);
                for i = 1:size(a,1)
                    [pks.baseline,loc.baseline] = findpeaks(a(i,is.sws.baseline),bins(is.sws.baseline),'MinPeakHeight',th);
                    [pks.reward,loc.reward] = findpeaks(a(i,is.sws.reward),bins(is.sws.reward),'MinPeakHeight',th);
                    [pks.aversive,loc.aversive] = findpeaks(a(i,is.sws.aversive),bins(is.sws.aversive),'MinPeakHeight',th);
                    
                    [pks.run,loc.run] = findpeaks(a(i,is.run.reward),bins(is.run.reward),'MinPeakHeight',th);
                    
                    [pks.a,loc.a] = findpeaks(a(i,is.run.aversive),bins(is.run.aversive),'MinPeakHeight',th);
                    [pks.r,loc.r] = findpeaks(a(i,is.run.reward),bins(is.run.reward),'MinPeakHeight',th);
                    
                    Rsquared.reward = [Rsquared.reward ; mean(pks.r,'omitnan')  mean(pks.a,'omitnan') length(Shocks_filt) length(Rewards_filt)];
                    
                    if aversiveTS_run(1) > rewardTS_run(1)
                        %                         disp(['Calculation of Reactivation Strenght for Reward assembly #',num2str(i)])
                        tmp = [sum(aa(i,is.sws.reward)) , sum(aa(i,is.sws.baseline))];
                        reactivation.reward.dvHPC = [reactivation.reward.dvHPC ; (mean(pks.reward) - mean(pks.baseline)) length(pks.reward) length(pks.baseline) , tmp]; clear tmp
                        for_scatter.reward = [for_scatter.reward ; mean(pks.run) (mean(pks.reward) - mean(pks.baseline)) length(Rewards_filt)];
                        for_reactivation.reward = [for_reactivation.reward ; mean(pks.baseline)];                
                    else
                        %                         disp(['Calculation of Reactivation Strenght for Reward assembly #',num2str(i)])
                        tmp = [sum(aa(i,is.sws.reward)) , sum(aa(i,is.sws.aversive))];
                        reactivation.reward.dvHPC = [reactivation.reward.dvHPC ; ((mean(pks.reward) - mean(pks.aversive))) length(pks.reward) length(pks.aversive) , tmp]; clear tmp
                        for_scatter.reward = [for_scatter.reward ; mean(pks.run) ((mean(pks.reward) - mean(pks.aversive))) length(Rewards_filt)];
                    end
                        AssemblyTime.joint.reward.baseline{i} = loc.baseline;
                        AssemblyTime.joint.reward.reward{i} = loc.reward;
                        AssemblyTime.joint.reward.aversive{i} = loc.aversive;
                    clear pks loc
                end
                clear a MUA aa
            end
            
            % --- for dHPC assemblies ---
            % Aversive
            if sum(cond.dHPC.aversive)>=1
                a = assembly_activity(patterns.all.aversive(:,cond.dHPC.aversive) , Spikes');
                aa = a;
                a = zscore(a,1,2);
                AssemblyTime.dHPC.aversive.baseline = cell(size(a,1),1);
                AssemblyTime.dHPC.aversive.reward = cell(size(a,1),1);
                AssemblyTime.dHPC.aversive.aversive = cell(size(a,1),1);
                for i = 1:size(a,1)
                    [pks.baseline,loc.baseline] = findpeaks(a(i,is.sws.baseline),bins(is.sws.baseline),'MinPeakHeight',th);
                    [pks.reward,loc.reward] = findpeaks(a(i,is.sws.reward),bins(is.sws.reward),'MinPeakHeight',th);
                    [pks.aversive,loc.aversive] = findpeaks(a(i,is.sws.aversive),bins(is.sws.aversive),'MinPeakHeight',th);
                    
                    if aversiveTS_run(1) < rewardTS_run(1)
                        tmp = [sum(aa(i,is.sws.aversive)) , sum(aa(i,is.sws.baseline))];
                        reactivation.aversive.dHPC = [reactivation.aversive.dHPC ; (mean(pks.aversive) - mean(pks.baseline)) length(pks.aversive) length(pks.baseline) , tmp]; clear tmp
                    else
                        tmp = [sum(aa(i,is.sws.aversive)) , sum(aa(i,is.sws.reward))];
                        reactivation.aversive.dHPC = [reactivation.aversive.dHPC ; (mean(pks.aversive) - mean(pks.reward)) length(pks.aversive) length(pks.reward) , tmp]; clear tmp
                    end
                    AssemblyTime.dHPC.aversive.baseline{i} = loc.baseline;
                    AssemblyTime.dHPC.aversive.reward{i} = loc.reward;
                    AssemblyTime.dHPC.aversive.aversive{i} = loc.aversive;
                    clear pks loc
                end
                clear a MUA aa
            end
            
            % Reward
            if sum(cond.dHPC.reward)>=1
                a = assembly_activity(patterns.all.reward(:,cond.dHPC.reward) , Spikes');
                aa = a;
                a = zscore(a,1,2);
                AssemblyTime.dHPC.reward.baseline = cell(size(a,1),1);
                AssemblyTime.dHPC.reward.reward = cell(size(a,1),1);
                AssemblyTime.dHPC.reward.aversive = cell(size(a,1),1);
                for i = 1:size(a,1)
                    [pks.baseline,loc.baseline] = findpeaks(a(i,is.sws.baseline),bins(is.sws.baseline),'MinPeakHeight',th);
                    [pks.reward,loc.reward] = findpeaks(a(i,is.sws.reward),bins(is.sws.reward),'MinPeakHeight',th);
                    [pks.aversive,loc.aversive] = findpeaks(a(i,is.sws.aversive),bins(is.sws.aversive),'MinPeakHeight',th);
                    
                    if aversiveTS_run(1) > rewardTS_run(1)
                        tmp = [sum(aa(i,is.sws.reward)) , sum(aa(i,is.sws.baseline))];
                        reactivation.reward.dHPC = [reactivation.reward.dHPC ; (mean(pks.reward) - mean(pks.baseline)) length(pks.reward) length(pks.baseline) , tmp]; clear tmp
                    else
                        tmp = [sum(aa(i,is.sws.reward)) , sum(aa(i,is.sws.aversive))];
                        reactivation.reward.dHPC = [reactivation.reward.dHPC ; ((mean(pks.reward) - mean(pks.aversive))) length(pks.reward) length(pks.aversive) , tmp]; clear tmp
                    end
                    AssemblyTime.dHPC.reward.baseline{i} = loc.baseline;
                    AssemblyTime.dHPC.reward.aversive{i} = loc.aversive;
                    AssemblyTime.dHPC.reward.reward{i} = loc.reward;
                        clear pks loc
                end
                clear a MUA aa
            end
            
            % --- for vHPC assemblies ---
            % Aversive
            if sum(cond.vHPC.aversive)>=1
                a = assembly_activity(patterns.all.aversive(:,cond.vHPC.aversive) , Spikes');
                aa = a;
                a = zscore(a,1,2);
                AssemblyTime.vHPC.aversive.baseline = cell(size(a,1),1);
                AssemblyTime.vHPC.aversive.reward = cell(size(a,1),1);
                AssemblyTime.vHPC.aversive.aversive = cell(size(a,1),1);
                for i = 1:size(a,1)
                    [pks.baseline,loc.baseline] = findpeaks(a(i,is.sws.baseline),bins(is.sws.baseline),'MinPeakHeight',th);
                    [pks.reward,loc.reward] = findpeaks(a(i,is.sws.reward),bins(is.sws.reward),'MinPeakHeight',th);
                    [pks.aversive,loc.aversive] = findpeaks(a(i,is.sws.aversive),bins(is.sws.aversive),'MinPeakHeight',th);
                    
                    if aversiveTS_run(1) < rewardTS_run(1)
                        tmp = [sum(aa(i,is.sws.aversive)) , sum(aa(i,is.sws.baseline))];
                        reactivation.aversive.vHPC = [reactivation.aversive.vHPC ; (mean(pks.aversive) - mean(pks.baseline)) length(pks.aversive) length(pks.baseline) , tmp]; clear tmp
                    else
                        tmp = [sum(aa(i,is.sws.aversive)) , sum(aa(i,is.sws.reward))];
                        reactivation.aversive.vHPC = [reactivation.aversive.vHPC ; (mean(pks.aversive) - mean(pks.reward)) length(pks.aversive) length(pks.reward) , tmp]; clear tmp
                    end
                    AssemblyTime.vHPC.aversive.baseline{i} = loc.baseline;
                    AssemblyTime.vHPC.aversive.reward{i} = loc.reward;
                    AssemblyTime.vHPC.aversive.aversive{i} = loc.aversive;
                    clear pks loc
                end
                clear a MUA aa
            end
            
            % Reward
            if sum(cond.vHPC.reward)>=1
                a = assembly_activity(patterns.all.reward(:,cond.vHPC.reward) , Spikes');
                aa = a;
                a = zscore(a,1,2);
                AssemblyTime.vHPC.reward.baseline = cell(size(a,1),1);
                AssemblyTime.vHPC.reward.reward = cell(size(a,1),1);
                AssemblyTime.vHPC.reward.aversive = cell(size(a,1),1);
                for i = 1:size(a,1)
                    [pks.baseline,loc.baseline] = findpeaks(a(i,is.sws.baseline),bins(is.sws.baseline),'MinPeakHeight',th);
                    [pks.reward,loc.reward] = findpeaks(a(i,is.sws.reward),bins(is.sws.reward),'MinPeakHeight',th);
                    [pks.aversive,loc.aversive] = findpeaks(a(i,is.sws.aversive),bins(is.sws.aversive),'MinPeakHeight',th);
                    
                    if aversiveTS_run(1) > rewardTS_run(1)
                        tmp = [sum(aa(i,is.sws.reward)) , sum(aa(i,is.sws.baseline))];
                        reactivation.reward.vHPC = [reactivation.reward.vHPC ; (mean(pks.reward) - mean(pks.baseline)) length(pks.reward) length(pks.baseline) , tmp]; clear tmp
                    else
                        tmp = [sum(aa(i,is.sws.reward)) , sum(aa(i,is.sws.aversive))];
                        reactivation.reward.vHPC = [reactivation.reward.vHPC ; ((mean(pks.reward) - mean(pks.aversive))) length(pks.reward) length(pks.aversive) , tmp]; clear tmp
                    end
                    AssemblyTime.vHPC.reward.baseline{i} = loc.baseline;
                    AssemblyTime.vHPC.reward.aversive{i} = loc.aversive;
                    AssemblyTime.vHPC.reward.reward{i} = loc.reward;
                    clear pks loc
                end
                clear a MUA aa
            end
            
            %% CCG between assemblies
            
%             tmp.b = [];
%             tmp.r = [];
%             tmp.a = [];
%             for i = 1 : size(AssemblyTime.joint.reward.baseline,1)
%                 x.b = AssemblyTime.joint.reward.baseline{i}';
%                 x.r = AssemblyTime.joint.reward.reward{i}';
%                 x.a = AssemblyTime.joint.reward.aversive{i}';
%                 for ii = 1 : size(AssemblyTime.joint.aversive.baseline,1)
%                     y.b = AssemblyTime.joint.aversive.baseline{ii}';
%                     y.r = AssemblyTime.joint.aversive.reward{ii}';
%                     y.a = AssemblyTime.joint.aversive.aversive{ii}';
%                     
%                     [s,ids,groups] = CCGParameters(x.b,ones(length(x.b),1),y.b,ones(length(y.b),1)*2);
%                     [ccg,ttttt] = CCG(s,ids,'binSize',0.01,'duration',1,'smooth',2,'mode','ccg'); %ccg calculation
%                     tmp.b = [tmp.b , ccg(:,1,2)./sum(ccg(:,1,2))];
%                     clear ccg s ids groups
%                     
%                     [s,ids,groups] = CCGParameters(x.r,ones(length(x.r),1),y.r,ones(length(y.r),1)*2);
%                     [ccg,ttttt] = CCG(s,ids,'binSize',0.01,'duration',1,'smooth',2,'mode','ccg'); %ccg calculation
%                     tmp.r = [tmp.r , ccg(:,1,2)./sum(ccg(:,1,2))];
%                     clear ccg s ids groups
%                     
%                     [s,ids,groups] = CCGParameters(x.a,ones(length(x.a),1),y.a,ones(length(y.a),1)*2);
%                     [ccg,ttttt] = CCG(s,ids,'binSize',0.01,'duration',1,'smooth',2,'mode','ccg'); %ccg calculation
%                     tmp.a = [tmp.a , ccg(:,1,2)./sum(ccg(:,1,2))];          
%                     clear ccg s ids groups y
%                 end
%                 clear x
%             end
%             
%             subplot(131)
%             plot(ttttt,tmp.b)
%             subplot(132)
%             plot(ttttt,tmp.r)
%             subplot(133)
%             plot(ttttt,tmp.a)
            
            %% Definition of SU id that are members of my assemblies
            disp('Members definition using 1/sqrt(N)thresold')
            % --- Aversive Condition ---
            if sum(cond.both.aversive)>=1
                % When they have SUs from both regions
                members.both.aversive = cell(1,sum(cond.both.aversive));
                nonmembers.both.aversive = cell(1,sum(cond.both.aversive));
                ii = 0;
                for i = 1:size(Thresholded.aversive.all,2)
                    if cond.both.aversive(i)
                        ii = ii+1;
                        members.both.aversive{ii} = clusters.all(Thresholded.aversive.all(:,i));
                        nonmembers.both.aversive{ii} = clusters.all(not(Thresholded.aversive.all(:,i)));
                    end
                end
                % When they have only dorsal SUs
                members.dHPC.aversive = cell(1,sum(cond.dHPC.aversive));
                ii = 0;
                for i = 1:size(Thresholded.aversive.all,2)
                    if cond.dHPC.aversive(i)
                        ii = ii+1;
                        members.dHPC.aversive{ii} = clusters.all(Thresholded.aversive.all(:,i));
                    end
                end
                % When they have only ventral SUs
                members.vHPC.aversive = cell(1,sum(cond.vHPC.aversive));
                ii = 0;
                for i = 1:size(Thresholded.aversive.all,2)
                    if cond.vHPC.aversive(i)
                        ii = ii+1;
                        members.vHPC.aversive{ii} = clusters.all(Thresholded.aversive.all(:,i));
                    end
                end
            end
            
            % --- Rewarded Condition ---
            if sum(cond.both.reward)>=1
                % When they have SUs from both regions
                members.both.reward = cell(1,sum(cond.both.reward));
                nonmembers.both.reward = cell(1,sum(cond.both.reward));
                ii = 0;
                for i = 1:size(Thresholded.reward.all,2)
                    if cond.both.reward(i)
                        ii = ii+1;
                        members.both.reward{ii} = clusters.all(Thresholded.reward.all(:,i));
                        nonmembers.both.reward{ii} = clusters.all(not(Thresholded.reward.all(:,i)));
                    end
                end
                % When they have only dorsal SUs
                members.dHPC.reward = cell(1,sum(cond.dHPC.reward));
                ii = 0;
                for i = 1:size(Thresholded.reward.all,2)
                    if cond.dHPC.reward(i)
                        ii = ii+1;
                        members.dHPC.reward{ii} = clusters.all(Thresholded.reward.all(:,i));
                    end
                end
                % When they have only ventral SUs
                members.vHPC.reward = cell(1,sum(cond.vHPC.reward));
                ii = 0;
                for i = 1:size(Thresholded.reward.all,2)
                    if cond.vHPC.reward(i)
                        ii = ii+1;
                        members.vHPC.reward{ii} = clusters.all(Thresholded.reward.all(:,i));
                    end
                end
            end
            
            % Save members ID for further analysis
            save([cd,'\assemblies_members.mat'],'members')
            
     
            %                 %% CCG
            %                 %Only members
            %                 % --- Reward ---
            %                 if sum(cond.both.reward)>=1
            %                     cross.reward= [];
            %                     for i = 1:size(members.both.reward,2)
            %                         x = members.both.reward{i}(ismember(members.both.reward{i},group_dHPC(:,1)));
            %                         y = members.both.reward{i}(ismember(members.both.reward{i},group_vHPC(:,1)));
            %                         for ii = 1 : size(x,1)
            %                             xx.reward = Restrict(spks_dHPC(ismember(spks_dHPC(:,1),x(ii)),2),movement.reward);
            %                             for iii = 1 : size(y,1)
            %                                 yy.reward =  Restrict(spks_vHPC(ismember(spks_vHPC(:,1),y(iii)),2),movement.reward);
            %                                 if and(length(yy.reward )>10 , length(xx.reward )>10)
            %                                             [s,ids,groups] = CCGParameters(yy.reward,ones(length(yy.reward),1),xx.reward,ones(length(xx.reward),1)*2);
            %                                             [ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',0.6,'smooth',2,'mode','ccg'); %ccg calculation
            %                                             cross.reward = [cross.reward , zscore(ccg(:,1,2)./sum(ccg(:,1,2)))];
            %                                             clear ccg s ids groups
            %                                 end
            %                                 clear yy
            %                             end
            %                             clear xx
            %                         end
            %                         clear x y
            %                     end
            %                 end
            % %                 imagesc(tttt,size(cross.reward,2),cross.reward')
            %
            %
            %                 % --- Aversive ---
            %                 if sum(cond.both.aversive)>=1
            %                     cross.aversive= [];
            %                     for i = 1:size(members.both.aversive,2)
            %                         x = members.both.aversive{i}(ismember(members.both.aversive{i},group_dHPC(:,1)));
            %                         y = members.both.aversive{i}(ismember(members.both.aversive{i},group_vHPC(:,1)));
            %                         for ii = 1 : size(x,1)
            %                             xx.aversive = Restrict(spks_dHPC(ismember(spks_dHPC(:,1),x(ii)),2),movement.aversive);
            %                             for iii = 1 : size(y,1)
            %                                 yy.aversive =  Restrict(spks_vHPC(ismember(spks_vHPC(:,1),y(iii)),2),movement.aversive);
            %                                 if and(length(yy.aversive )>10 , length(xx.aversive )>10)
            %                                             [s,ids,groups] = CCGParameters(yy.aversive,ones(length(yy.aversive),1),xx.aversive,ones(length(xx.aversive),1)*2);
            %                                             [ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',0.6,'smooth',2,'mode','ccg'); %ccg calculation
            %                                             cross.aversive = [cross.aversive , zscore(ccg(:,1,2)./sum(ccg(:,1,2)))];
            %                                             clear ccg s ids groups
            %                                 end
            %                                 clear yy
            %                             end
            %                             clear xx
            %                         end
            %                         clear x y
            %                     end
            %                 end
            % %                imagesc(tttt,size(cross.aversive,2),cross.aversive')
            %
            %
            %                 % Only non-members
            %                 % --- Aversive ---
            %                 if sum(cond.both.aversive)>=1
            %                     crossN.aversive = [];
            %                     for i = 1:size(members.both.aversive,2)
            %                         x = nonmembers.both.aversive{i}(ismember(nonmembers.both.aversive{i},group_dHPC(:,1)));
            %                         y = nonmembers.both.aversive{i}(ismember(nonmembers.both.aversive{i},group_vHPC(:,1)));
            %                         for ii = 1 : size(x,1)
            %                             xx.aversive = Restrict(spks_dHPC(ismember(spks_dHPC(:,1),x(ii)),2),movement.aversive);
            %                             for iii = 1 : size(y,1)
            %                                 yy.aversive =  Restrict(spks_vHPC(ismember(spks_vHPC(:,1),y(iii)),2),movement.aversive);
            %                                 if and(length(yy.aversive )>10 , length(xx.aversive )>10)
            %                                     [s,ids,groups] = CCGParameters(yy.aversive,ones(length(yy.aversive),1),xx.aversive,ones(length(xx.aversive),1)*2);
            %                                     [ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',0.6,'smooth',2,'mode','ccg'); %ccg calculation
            %                                     crossN.aversive = [crossN.aversive , zscore(ccg(:,1,2)./sum(ccg(:,1,2)))];
            %                                     clear ccg s ids groups
            %                                 end
            %                                 clear yy
            %                             end
            %                             clear xx
            %                         end
            %                         clear x y
            %                     end
            %                 end
            %
            %                 % --- Reward ---
            %                 if sum(cond.both.reward)>=1
            %                     crossN.reward = [];
            %                     for i = 1:size(members.both.reward,2)
            %                         x = nonmembers.both.reward{i}(ismember(nonmembers.both.reward{i},group_dHPC(:,1)));
            %                         y = nonmembers.both.reward{i}(ismember(nonmembers.both.reward{i},group_vHPC(:,1)));
            %                         for ii = 1 : size(x,1)
            %                             xx.reward = Restrict(spks_dHPC(ismember(spks_dHPC(:,1),x(ii)),2),movement.reward);
            %                             for iii = 1 : size(y,1)
            %                                 yy.reward =  Restrict(spks_vHPC(ismember(spks_vHPC(:,1),y(iii)),2),movement.reward);
            %                                 if and(length(yy.reward )>10 , length(xx.reward )>10)
            %                                     [s,ids,groups] = CCGParameters(yy.reward,ones(length(yy.reward),1),xx.reward,ones(length(xx.reward),1)*2);
            %                                     [ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',0.6,'smooth',2,'mode','ccg'); %ccg calculation
            %                                     crossN.reward = [crossN.reward , zscore(ccg(:,1,2)./sum(ccg(:,1,2)))];
            %                                     clear ccg s ids groups
            %                                 end
            %                                 clear yy
            %                             end
            %                             clear xx
            %                         end
            %                         clear x y
            %                     end
            %                 end
            %
            %
            %
            %                 if not(exist('CrossCorrelograms'))
            %                     if sum(cond.both.reward)>=1
            %                         CrossCorrelograms.reward = cross.reward;
            %                         CrossCorrelograms.non.reward = crossN.reward;
            %                     end
            %
            %                     if sum(cond.both.aversive)>=1
            %                         CrossCorrelograms.aversive = cross.aversive;
            %                         CrossCorrelograms.non.aversive = crossN.aversive;
            %                     end
            %
            %                 else
            %                     if sum(cond.both.reward)>=1
            %                         CrossCorrelograms.reward = [CrossCorrelograms.reward , cross.reward];
            %                         CrossCorrelograms.non.reward = [CrossCorrelograms.non.reward , crossN.reward];
            %                     end
            %
            %                     if sum(cond.both.aversive)>=1
            %                         CrossCorrelograms.aversive = [CrossCorrelograms.aversive , cross.aversive];
            %                         CrossCorrelograms.non.aversive = [CrossCorrelograms.non.aversive , crossN.aversive];
            %                     end
            %                 end
            %
            %% Assembly pattern reactivation sourrounding ripple cooridnated events
            if sum(cond.both.aversive)>=1
                a = assembly_activity(patterns.all.aversive(:,cond.both.aversive) , Spikes');
                a = zscore(a,1,2);
                x.pre = [];
                x.post = [];
                for i = 1 : size(a,1)
                    [pks,loc] = findpeaks(a(i,:),bins,'MinPeakHeight',th);
                    
                    xx = Restrict(cooridnated_event(:,1)+(cooridnated_event(:,2) - cooridnated_event(:,1))/2 , NREM.aversive);
                    yy = Restrict(loc , NREM.aversive);
                    [s,ids,groups] = CCGParameters(xx,ones(length(xx),1),yy,ones(length(yy),1)*2);
                    [ccg,tttt] = CCG(s,ids,'binSize',0.025,'duration',1,'smooth',0,'mode','ccg'); %ccg calculation
                    
                    x.post = [x.post , ccg(:,1,2)./sum(ccg(:,1,2))];
                    clear yy xx ccg s ids groups
                    if aversiveTS_run(1)<rewardTS_run(1)
                        xx = Restrict(cooridnated_event(:,1)+(cooridnated_event(:,2) - cooridnated_event(:,1))/2 , NREM.baseline);
                        yy = Restrict(loc , NREM.baseline);
                        [s,ids,groups] = CCGParameters(xx,ones(length(xx),1),yy,ones(length(yy),1)*2);
                        [ccg,tttt] = CCG(s,ids,'binSize',0.025,'duration',1,'smooth',0,'mode','ccg'); %ccg calculation
                        x.pre = [x.pre , ccg(:,1,2)./sum(ccg(:,1,2))];
                        clear yy xx ccg s ids groups
                    else
                        xx = Restrict(cooridnated_event(:,1)+(cooridnated_event(:,2) - cooridnated_event(:,1))/2 , NREM.reward);
                        yy = Restrict(loc , NREM.reward);
                        [s,ids,groups] = CCGParameters(xx,ones(length(xx),1),yy,ones(length(yy),1)*2);
                        [ccg,tttt] = CCG(s,ids,'binSize',0.025,'duration',1,'smooth',0,'mode','ccg'); %ccg calculation
                        x.pre = [x.pre , ccg(:,1,2)./sum(ccg(:,1,2))];
                        clear yy xx ccg s ids groups
                    end
                    clear loc pks
                end
                
                % save CCG
                    Cross.aversive.pre = [Cross.aversive.pre , x.pre];
                    Cross.aversive.post = [Cross.aversive.post , x.post];
                    clear x a             
            end
            
            
            
            if sum(cond.both.reward)>=1
                a = assembly_activity(patterns.all.reward(:,cond.both.reward) , Spikes');
                a = zscore(a,1,2);
                x.pre = [];
                x.post = [];
                for i = 1 : size(a,1)
                    [pks,loc] = findpeaks(a(i,:),bins,'MinPeakHeight',th);
                    
                    xx = Restrict(cooridnated_event(:,1)+(cooridnated_event(:,2) - cooridnated_event(:,1))/2 , NREM.reward);
                    yy = Restrict(loc , NREM.reward);
                    [s,ids,groups] = CCGParameters(xx,ones(length(xx),1),yy,ones(length(yy),1)*2);
                    [ccg,tttt] = CCG(s,ids,'binSize',0.025,'duration',1,'smooth',0,'mode','ccg'); %ccg calculation
                    
                    x.post = [x.post , ccg(:,1,2)./sum(ccg(:,1,2))];
                    clear yy xx ccg s ids groups
                    
                    if aversiveTS_run(1) > rewardTS_run(1)
                        xx = Restrict(cooridnated_event(:,1)+(cooridnated_event(:,2) - cooridnated_event(:,1))/2 , NREM.baseline);
                        yy = Restrict(loc , NREM.baseline);
                        [s,ids,groups] = CCGParameters(xx,ones(length(xx),1),yy,ones(length(yy),1)*2);
                        [ccg,tttt] = CCG(s,ids,'binSize',0.025,'duration',1,'smooth',0,'mode','ccg'); %ccg calculation
                        x.pre = [x.pre , ccg(:,1,2)./sum(ccg(:,1,2))];
                        clear yy xx ccg s ids groups
                    else
                        xx = Restrict(cooridnated_event(:,1)+(cooridnated_event(:,2) - cooridnated_event(:,1))/2 , NREM.aversive);
                        yy = Restrict(loc , NREM.aversive);
                        [s,ids,groups] = CCGParameters(xx,ones(length(xx),1),yy,ones(length(yy),1)*2);
                        [ccg,tttt] = CCG(s,ids,'binSize',0.025,'duration',1,'smooth',0,'mode','ccg'); %ccg calculation
                        x.pre = [x.pre , ccg(:,1,2)./sum(ccg(:,1,2))];
                        clear yy xx ccg s ids groups
                    end
                    clear loc pks
                end
                
                % save CCG
                    Cross.reward.pre = [Cross.reward.pre , x.pre];
                    Cross.reward.post = [Cross.reward.post , x.post];
                    clear x a
            end
            
            
        end
        clear A aversiveTS aversiveTS_run baselineTS rewardTS rewardTS_run
        clear behavior bins Cell_type_classification cellulartype cond
        clear is K Kinfo group_dHPC group_vHPC matrixC matrixD matrixV
        clear NREM REM WAKE segmentation tmp
        clear spiketrains_dHPC spiketrains_vHPC opts
        clear patterns Thresholded i  ii numberD numberV movement cross crossN
        clear Spikes bins Clusters Shocks_filt Rewards_filt
    end
end

%% Para el facu del futuro, mejorar esto.
%  Hacerlo mas automatico
Rsquared.reward(or(isnan(Rsquared.reward(:,1)),isnan(Rsquared.reward(:,2))),:) = [];
Rsquared.aversive(or(isnan(Rsquared.aversive(:,1)),isnan(Rsquared.aversive(:,2))),:) = [];


grps = [ones(length(Rsquared.reward),1) ; ones(length(Rsquared.reward),1)*2 ; ones(length(Rsquared.aversive),1)*3 ; ones(length(Rsquared.aversive),1)*4];
y = [Rsquared.reward(:,1) ; Rsquared.reward(:,2) ; Rsquared.aversive(:,2) ; Rsquared.aversive(:,1)];
boxplot(y,grps),ylim([0 40])

reps = [[ones(length(Rsquared.reward),1) ; ones(length(Rsquared.aversive),1)*2] , [ones(length(Rsquared.reward),1) ; ones(length(Rsquared.aversive),1)*2]];
g1 = [ones(length(Rsquared.reward),1) ; ones(length(Rsquared.reward),1) ; ones(length(Rsquared.aversive),1)*2 ; ones(length(Rsquared.aversive),1)*2];
g2 = [ones(length(Rsquared.reward),1) ; ones(length(Rsquared.reward),1)*2 ; ones(length(Rsquared.aversive),1) ; ones(length(Rsquared.aversive),1)*2];
y = [Rsquared.reward(:,1) ; Rsquared.reward(:,2) ; Rsquared.aversive(:,2) ; Rsquared.aversive(:,1)];

[~,~,stats] = anovan(y,{g1 g2},"Model","interaction", "Varnames",["g1","g2"]);
[results,~,~,gnames] = multcompare(stats,'CType','hsd','Dimension',[1 2]);
[h ,p] = anovan(y,reps)

x = Rsquared.reward(:,1)./Rsquared.reward(:,2);
x(isnan(x)) = [];
y = Rsquared.aversive(:,1) ./ Rsquared.aversive(:,2);
y(isnan(y)) = [];
xx = [1 2];
yy = [mean(x,'omitnan') mean(y,'omitnan')];
err = [(std(x)/sqrt(length(x))) (std(y)/sqrt(length(y)))];
figure
bar(xx,yy),hold on
er = errorbar(xx,yy,err);ylim([-1.5 2.5])

grps = [ones(length(x),1) ; ones(length(y),1)*2];
boxplot([x;y] , grps),ylim([0 2.5])


%% Plot Strenght Reactivation
% for joint assemblies
reactivation.reward.dvHPC(isnan(reactivation.reward.dvHPC(:,1)),:) = [];
reactivation.aversive.dvHPC(isnan(reactivation.aversive.dvHPC(:,1)),:) = [];
x = reactivation.reward.dvHPC(:,1);
y = reactivation.aversive.dvHPC(:,1);
[h, p] = ttest2(x,y)             

xx = [1 2];
yy = [mean(x) mean(y)];
err = [(std(x)/sqrt(length(x))) (std(y)/sqrt(length(y)))];

subplot(131),
bar(xx,yy),hold on
er = errorbar(xx,yy,err);ylim([-10 10])
er.Color = [0 0 0];                            
er.LineStyle = 'none';
hold off

% for dHPC assemblies
reactivation.reward.dHPC(isnan(reactivation.reward.dHPC(:,1)),:) = [];
reactivation.aversive.dHPC(isnan(reactivation.aversive.dHPC(:,1)),:) = [];
x = reactivation.reward.dHPC(:,1);
y = reactivation.aversive.dHPC(:,1);
[h, p] = ttest2(x,y)             

xx = [1 2];
yy = [mean(x) mean(y)];
err = [(std(x)/sqrt(length(x))) (std(y)/sqrt(length(y)))];

subplot(132),
bar(xx,yy),hold on
er = errorbar(xx,yy,err);ylim([-10 10])
er.Color = [0 0 0];                            
er.LineStyle = 'none';
hold off

% for vHPC assemblies
reactivation.reward.vHPC(isnan(reactivation.reward.vHPC(:,1)),:) = [];
reactivation.aversive.vHPC(isnan(reactivation.aversive.vHPC(:,1)),:) = [];
x = reactivation.reward.vHPC(:,1);
y = reactivation.aversive.vHPC(:,1);
[h, p] = ttest2(x,y)             

xx = [1 2];
yy = [mean(x) mean(y)];
err = [(std(x)/sqrt(length(x))) (std(y)/sqrt(length(y)))];

subplot(133),
bar(xx,yy),hold on
er = errorbar(xx,yy,err);ylim([-10 10])
er.Color = [0 0 0];                            
er.LineStyle = 'none';
hold off

%% Plotting the number of peaks
x = reactivation.reward.dvHPC(:,2)./reactivation.reward.dvHPC(:,3);
y = reactivation.aversive.dvHPC(:,2)./reactivation.aversive.dvHPC(:,3);
[h, p] = ttest2(x,y)             

xx = [1 2];
yy = [mean(x) mean(y)];
err = [(std(x)/sqrt(length(x))) (std(y)/sqrt(length(y)))];

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
yy = [mean(x) mean(y)];
err = [(std(x)/sqrt(length(x))) (std(y)/sqrt(length(y)))];

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
yy = [mean(x) mean(y)];
err = [(std(x)/sqrt(length(x))) (std(y)/sqrt(length(y)))];

subplot(133),
bar(xx,yy),hold on
er = errorbar(xx,yy,err);ylim([0 3])
er.Color = [0 0 0];                            
er.LineStyle = 'none';
hold off

%% Plotting the mean of the activity strengtht
x = reactivation.reward.dvHPC(:,4)-reactivation.reward.dvHPC(:,5);
y = reactivation.aversive.dvHPC(:,4)-reactivation.aversive.dvHPC(:,5);
[h, p] = ttest2(x,y)             

xx = [1 2];
yy = [mean(x) mean(y)];
err = [(std(x)/sqrt(length(x))) (std(y)/sqrt(length(y)))];

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


figure,
x = for_scatter.aversive;
x(x(:,1)>60,:)=[];
scatter(x(:,1),x(:,2),'filled','r')%,ylim([-12 12]),xlim([8 52]),hold on
figure,plot(fitlm(x(:,1),x(:,2)))%,ylim([-12 12]),xlim([8 52])


figure
y = for_scatter.reward;
y(y(:,1)>60,:)=[];
scatter(y(:,1),y(:,2),'filled','b')%,ylim([-12 12]),xlim([8 52]),
figure,plot(fitlm(y(:,1),y(:,2)))%,ylim([-12 12]),xlim([8 52])


%% For venn graph

num = sum(percentages);
            


            if or(sum(cond.both.reward)>=1 , sum(cond.both.aversive)>=1)

                
                %% Characterization of joint assemblies
                if and(sum(cond.both.aversive)>=1 , sum(cond.both.reward)>=1)
                    characterization.both = [];
                    for i = 1 : size(cond.dHPC.aversive,2)
                        if cond.both.aversive(i)
                            for ii = 1 : size(cond.dHPC.reward,2)
                                if cond.both.reward(ii)
                                    if sum(and(Thresholded.aversive.all(:,i) , Thresholded.reward.all(:,ii)))>0
                                        %                                                                 break
                                        iii = not(and(Thresholded.aversive.all(:,i) , Thresholded.reward.all(:,ii)));
                                        c = corrcoef(patterns.all.aversive(:,i) , patterns.all.reward(:,ii));
                                        cc = dot(patterns.all.aversive(:,i) , patterns.all.reward(:,ii)) / (norm(patterns.all.aversive(:,i)) * norm(patterns.all.reward(:,ii)));
                                        ccc = corrcoef(patterns.all.aversive(iii,i) , patterns.all.reward(iii,ii));
                                        characterization.both = [characterization.both ; i ii c(1,2) cc ccc(1,2) sum(abs(patterns.all.aversive(:,i))) sum(abs(patterns.all.reward(:,ii)))];
                                        clear c cc ccc iii
                                    end
                                end
                            end
                            %                     break
                        end
                    end
                    
                    characterization.dHPC = [];
                    for i = 1 : size(cond.dHPC.aversive,2)
                        if cond.dHPC.aversive(i)
                            for ii = 1 : size(cond.dHPC.reward,2)
                                if cond.dHPC.reward(ii)
                                    if sum(and(Thresholded.aversive.all(:,i) , Thresholded.reward.all(:,ii)))>0
                                        %                                                                 break
                                        iii = not(and(Thresholded.aversive.all(:,i) , Thresholded.reward.all(:,ii)));
                                        c = corrcoef(patterns.all.aversive(:,i) , patterns.all.reward(:,ii));
                                        cc = dot(patterns.all.aversive(:,i) , patterns.all.reward(:,ii)) / (norm(patterns.all.aversive(:,i)) * norm(patterns.all.reward(:,ii)));
                                        ccc = corrcoef(patterns.all.aversive(iii,i) , patterns.all.reward(iii,ii));
                                        characterization.dHPC = [characterization.dHPC ; i ii c(1,2) cc ccc(1,2) sum(abs(patterns.all.aversive(:,i))) sum(abs(patterns.all.reward(:,ii)))];
                                        clear c cc ccc iii
                                    end
                                end
                            end
                            %                         break
                        end
                    end
                    
                    characterization.vHPC = [];
                    for i = 1 : size(cond.vHPC.aversive,2)
                        if cond.vHPC.aversive(i)
                            for ii = 1 : size(cond.vHPC.reward,2)
                                if cond.vHPC.reward(ii)
                                    if sum(and(Thresholded.aversive.all(:,i) , Thresholded.reward.all(:,ii)))>0
                                        %                                                                 break
                                        iii = not(and(Thresholded.aversive.all(:,i) , Thresholded.reward.all(:,ii)));
                                        c = corrcoef(patterns.all.aversive(:,i) , patterns.all.reward(:,ii));
                                        cc = dot(patterns.all.aversive(:,i) , patterns.all.reward(:,ii)) / (norm(patterns.all.aversive(:,i)) * norm(patterns.all.reward(:,ii)));
                                        ccc = corrcoef(patterns.all.aversive(iii,i) , patterns.all.reward(iii,ii));
                                        characterization.vHPC = [characterization.vHPC ; i ii c(1,2) cc ccc(1,2) sum(abs(patterns.all.aversive(:,i))) sum(abs(patterns.all.reward(:,ii)))];
                                        clear c cc ccc iii
                                    end
                                end
                            end
                            %                                             break
                        end
                    end
                end
                


[i ii] = max(CrossCorrelograms.reward);
[i , iii] = sort(ii,'descend');
figure,
subplot(131),imagesc(tttt,[1:1:size(CrossCorrelograms.reward,2)],CrossCorrelograms.reward(:,iii)'),caxis([-3 3]),xlim([-0.3 0.3]),colormap 'jet'
clear i ii iii

CrossCorrelograms.non.reward(:,isnan(sum(CrossCorrelograms.non.reward,1))) = [];
o = randperm(size(CrossCorrelograms.non.reward,2));
x = CrossCorrelograms.non.reward(:,o(1:size(CrossCorrelograms.reward,2)));
[i ii] = max(x);
[i , iii] = sort(ii,'descend');
subplot(132),imagesc(tttt,[1:1:size(CrossCorrelograms.reward,2)],x(:,iii)'),caxis([-3 3]),xlim([-0.3 0.3]),colormap 'jet'


subplot(133),plot(tttt,mean(CrossCorrelograms.reward,2,'omitnan'),'b','LineWidth',2),hold on
sem = std(CrossCorrelograms.reward,0,2,'omitnan')/sqrt(size(CrossCorrelograms.reward,2));
ciplot(mean(CrossCorrelograms.reward,2,'omitnan')-sem , mean(CrossCorrelograms.reward,2,'omitnan')+sem,tttt,'b'),alpha 0.5
plot(tttt,mean(x,2,'omitnan'),'k','LineWidth',2),hold on
sem = std(x,0,2,'omitnan')/sqrt(size(x,2));
ciplot(mean(x,2,'omitnan')-sem , mean(x,2,'omitnan')+sem,tttt,'k'),alpha 0.5,ylim([-0.3 0.6]),xlim([-0.3 0.3])

[i , ii] = min(abs(tttt-(-0.025)));
[i , iii] = min(abs(tttt-(-0.025)));

xx = mean(CrossCorrelograms.reward(ii:iii,:),1);
xxx = mean(CrossCorrelograms.non.reward(ii:iii,o),1);

[h , p] = ranksum(xx,xxx)



[i ii] = max(CrossCorrelograms.aversive);
[i , iii] = sort(ii,'descend');
figure,
subplot(131),imagesc(tttt,[1:1:size(CrossCorrelograms.aversive,2)],CrossCorrelograms.aversive(:,iii)'),caxis([-3 3]),xlim([-0.3 0.3]),colormap 'jet'
clear i ii iii

CrossCorrelograms.non.aversive(:,isnan(sum(CrossCorrelograms.non.aversive,1))) = [];
o = randperm(size(CrossCorrelograms.non.aversive,2));
x = CrossCorrelograms.non.aversive(:,o(1:size(CrossCorrelograms.aversive,2)));
[i ii] = max(x);
[i , iii] = sort(ii,'descend');
subplot(132),imagesc(tttt,[1:1:size(CrossCorrelograms.aversive,2)],x(:,iii)'),caxis([-3 3]),xlim([-0.3 0.3]),colormap 'jet'


subplot(133),plot(tttt,mean(CrossCorrelograms.aversive,2,'omitnan'),'r','LineWidth',2),hold on
sem = std(CrossCorrelograms.aversive,0,2,'omitnan')/sqrt(size(CrossCorrelograms.aversive,2));
ciplot(mean(CrossCorrelograms.aversive,2,'omitnan')-sem , mean(CrossCorrelograms.aversive,2,'omitnan')+sem,tttt,'r'),alpha 0.5
plot(tttt,mean(x,2,'omitnan'),'k','LineWidth',2),hold on
sem = std(x,0,2,'omitnan')/sqrt(size(x,2));
ciplot(mean(x,2,'omitnan')-sem , mean(x,2,'omitnan')+sem,tttt,'k'),alpha 0.5,ylim([-0.3 0.6]),xlim([-0.3 0.3])
xline(-0.025,'--'),xline(0.025,'--'),


[i , ii] = min(abs(tttt-(-0.025)));
[i , iii] = min(abs(tttt-(-0.025)));

xx = mean(CrossCorrelograms.aversive(ii:iii,:),1);
xxx = mean(CrossCorrelograms.non.aversive(ii:iii,o),1);

[h , p] = ranksum(xx,xxx)

% figure
% plot(tttt,mean(CrossCorrelograms.reward.baseline,2,'omitnan'),'k','LineWidth',2),hold on
% plot(tttt,mean(CrossCorrelograms.reward.reward,2,'omitnan'),'b','LineWidth',2),hold on
% plot(tttt,mean(CrossCorrelograms.reward.aversive,2,'omitnan'),'r','LineWidth',2),hold on


% figure
% plot(tttt,mean(CrossCorrelograms.aversive.baseline,2,'omitnan'),'k','LineWidth',2),hold on
% plot(tttt,mean(CrossCorrelograms.aversive.reward,2,'omitnan'),'b','LineWidth',2),hold on
% plot(tttt,mean(CrossCorrelograms.aversive.aversive,2,'omitnan'),'r','LineWidth',2),hold on



scatter3(counts.both(:,5),counts.both(:,6),counts.both(:,3)),xlabel('Aversive'),ylabel('Reward'),zlabel('R')

scatter(counts.both(:,3),counts.both(:,5)) , xlabel('With'), ylabel('Without'),xline(0,'--'),yline(0,'--')
scatter(counts.dHPC(:,3),counts.dHPC(:,5))
scatter(counts.vHPC(:,3),counts.vHPC(:,5))

            
            %% Binning of sessions in 100 setps sleep sessions
            segmentation = [0 : win : segments.Var1(end)/1000]; clear dt
            tmp = [];
            for i = 2 : size(segmentation,2)-1
                tmp = [tmp , InIntervals(bins,[segmentation(i-1) segmentation(i+1)])];
            end
            segmentation = logical(tmp);
            clear tmp i
            
            correlations.cross.aversive = cell(1,size(segmentation,2));
            correlations.cross.reward = cell(1,size(segmentation,2));
            correlations.cross.dHPC.reward = cell(1,size(segmentation,2));
            correlations.cross.vHPC.reward = cell(1,size(segmentation,2));
            correlations.cross.dHPC.aversive = cell(1,size(segmentation,2));
            correlations.cross.vHPC.aversive = cell(1,size(segmentation,2));            
            for i = 1 : size(segmentation,2)
                x = ActivityAssembles.dHPC.aversive(:,segmentation(:,i));
                y = ActivityAssembles.vHPC.aversive(:,segmentation(:,i));
                correlations.cross.aversive{i} = corr(x',y');
                clear x y
                
                x = ActivityAssembles.dHPC.reward(:,segmentation(:,i));
                y = ActivityAssembles.vHPC.reward(:,segmentation(:,i));
                correlations.cross.reward{i} = corr(x',y');
                clear x y
                
                x = ActivityAssembles.dHPC.reward(:,segmentation(:,i));
                y = ActivityAssembles.dHPC.reward(:,segmentation(:,i));
                correlations.cross.dHPC.reward{i} = corr(x',y');
                clear x y  
                
                x = ActivityAssembles.vHPC.reward(:,segmentation(:,i));
                y = ActivityAssembles.vHPC.reward(:,segmentation(:,i));
                correlations.cross.vHPC.reward{i} = corr(x',y');
                clear x y                 
                
                x = ActivityAssembles.dHPC.aversive(:,segmentation(:,i));
                y = ActivityAssembles.dHPC.aversive(:,segmentation(:,i));
                correlations.cross.dHPC.aversive{i} = corr(x',y');
                clear x y  
                
                x = ActivityAssembles.vHPC.aversive(:,segmentation(:,i));
                y = ActivityAssembles.vHPC.aversive(:,segmentation(:,i));
                correlations.cross.vHPC.aversive{i} = corr(x',y');
                clear x y                                 
            end
            clear i
            
            matrixC.aversive = [];
            matrixC.reward = [];
            matrixC.cross.dHPC.reward = [];
            matrixC.cross.vHPC.reward = [];
            matrixC.cross.dHPC.aversive = [];
            matrixC.cross.vHPC.aversive = [];
            for i = 1 : size(segmentation,2)
                tmpA = [];
                tmpR = [];
                tmp1 = []; tmp2 = []; tmp3 = []; tmp4 = [];
                for ii = 1 : size(segmentation,2)
                    x = correlations.cross.aversive{i};
                    y = correlations.cross.aversive{ii};
                    c = corrcoef(x,y,'rows','complete');
                    tmpA = [tmpA c(1,2)];
                    clear x y c
                    
                    x = correlations.cross.dHPC.reward{i};
                    y = correlations.cross.dHPC.reward{ii};
                    c = corrcoef(x,y,'rows','complete');
                    tmp1 = [tmp1 c(1,2)];
                    clear x y c  
                    
                    x = correlations.cross.vHPC.reward{i};
                    y = correlations.cross.vHPC.reward{ii};
                    c = corrcoef(x,y,'rows','complete');
                    tmp2 = [tmp2 c(1,2)];
                    clear x y c 
                    
                    x = correlations.cross.dHPC.aversive{i};
                    y = correlations.cross.dHPC.aversive{ii};
                    c = corrcoef(x,y,'rows','complete');
                    tmp3 = [tmp3 c(1,2)];
                    clear x y c  
                    
                    x = correlations.cross.vHPC.aversive{i};
                    y = correlations.cross.vHPC.aversive{ii};
                    c = corrcoef(x,y,'rows','complete');
                    tmp4 = [tmp4 c(1,2)];
                    clear x y c                     
                end
                
                matrixC.aversive = [matrixC.aversive ; tmpA];
                matrixC.reward = [matrixC.reward ; tmpR];
                matrixC.cross.dHPC.reward = [matrixC.cross.dHPC.reward ; tmp1];
                matrixC.cross.vHPC.reward = [matrixC.cross.vHPC.reward ; tmp2];
                matrixC.cross.dHPC.aversive = [matrixC.cross.dHPC.aversive ; tmp3];
                matrixC.cross.vHPC.aversive = [matrixC.cross.vHPC.aversive ; tmp4];
                
                clear tmpA tmpR tmp1 tmp2 tmp3 tmp4
            end
            
            figure,
            subplot(121)
            matrixC.cross.vHPC(eye(size(matrixC.aversive))==1) = nan;
            imagesc([300 : 300 : segments.Var1(end)/1000],[0 : 300 : segments.Var1(end)/1000],matrixC.cross.vHPC)
            hold on
            xline(aversiveTS_run(1)/1000,'r')
            xline(aversiveTS_run(2)/1000,'r')
            xline(rewardTS_run(1)/1000,'b')
            xline(rewardTS_run(2)/1000,'b')
            yline(aversiveTS_run(1)/1000,'r')
            yline(aversiveTS_run(2)/1000,'r')
            yline(rewardTS_run(1)/1000,'b')
            yline(rewardTS_run(2)/1000,'b')
            title('Cross A')
            
            
            subplot(122)
            matrixC.reward(eye(size(matrixC.reward))==1) = nan;
            imagesc([300 : 300 : segments.Var1(end)/1000],[0 : 300 : segments.Var1(end)/1000],matrixC.reward)
            hold on
            xline(aversiveTS_run(1)/1000,'r')
            xline(aversiveTS_run(2)/1000,'r')
            xline(rewardTS_run(1)/1000,'b')
            xline(rewardTS_run(2)/1000,'b')
            yline(aversiveTS_run(1)/1000,'r')
            yline(aversiveTS_run(2)/1000,'r')
            yline(rewardTS_run(1)/1000,'b')
            yline(rewardTS_run(2)/1000,'b')
            title('Cross R')            
            
            segmentation = [win : win : segments.Var1(end)/1000-win]; 
            i = InIntervals(segmentation,aversiveTS_run./1000);
            figure,
            i = matrixC.aversive(i,:);
            mean_cross_aversive = mean(i);
            plot(segmentation,mean(i),'r')
            hold on
            
            i = InIntervals(segmentation,rewardTS_run./1000);
            i = matrixV(i,:);
            mean_cross_reward = mean(i);
            plot(segmentation,mean(i),'b')
            hold on
            xline(aversiveTS_run(1)/1000,'r')
            xline(aversiveTS_run(2)/1000,'r')
            xline(rewardTS_run(1)/1000,'b')
            xline(rewardTS_run(2)/1000,'b')
            title('Cross')            
            
            
            % Normalization of the graphs
            if aversiveTS_run(2) > rewardTS_run(2)
                index = [baselineTS ; rewardTS_run ; rewardTS ; aversiveTS_run ; aversiveTS] ./ 1000;
            else
                index = [baselineTS ; aversiveTS_run ; aversiveTS ; rewardTS_run ; rewardTS] ./ 1000;
            end
            
            tmp = (index(:,2) - index(:,1));
            tmp(1) = tmp(1)/40;
            tmp(2) = tmp(2)/4;
            tmp(3) = tmp(3)/40;
            tmp(4) = tmp(4)/4;
            tmp(5) = tmp(5)/40;
            index = [index(:,1) , tmp , index(:,2)]; clear tmp
            
            time = [];
            for i = 1 : size(index,1)
                time = [time , index(i,1) : index(i,2) : index(i,3)-index(i,2)];
            end
            
            tmp = [];
            for i = 1 : size(time,2)-1
                ii = InIntervals(segmentation,[time(i) time(i+1)]);
                iii = mean(mean_cross_aversive(ii));
                iiii =  mean(mean_cross_reward(ii));
                ii = [iii ; iiii];
                tmp = [tmp , ii]; clear ii iii iiii
            end
            clear mean_cross_aversive mean_cross_reward i time index
            
            if aversiveTS_run(1)<rewardTS_run(1)
                output.aversive.first.aversive = [output.aversive.first.aversive ; tmp(1,:)];
                output.aversive.first.reward = [output.aversive.first.reward ; tmp(2,:)];
            else
                output.aversive.second.aversive = [output.aversive.second.aversive ; tmp(1,:)];
                output.aversive.second.reward = [output.aversive.second.reward ; tmp(2,:)];
            end

%             subplot(132)
%             matrixD(eye(size(matrixD))==1) = nan;
%             imagesc([300 : 300 : segments.Var1(end)/1000],[0 : 300 : segments.Var1(end)/1000],matrixD)
%             hold on
%             xline(aversiveTS_run(1)/1000,'r')
%             xline(aversiveTS_run(2)/1000,'r')
%             xline(rewardTS_run(1)/1000,'b')
%             xline(rewardTS_run(2)/1000,'b')
%             yline(aversiveTS_run(1)/1000,'r')
%             yline(aversiveTS_run(2)/1000,'r')
%             yline(rewardTS_run(1)/1000,'b')
%             yline(rewardTS_run(2)/1000,'b')
%             title('dHPC')          
%             
%             subplot(133)
%             matrixV(eye(size(matrixV))==1) = nan;
%             imagesc([300 : 300 : segments.Var1(end)/1000],[0 : 300 : segments.Var1(end)/1000],matrixV)
%             hold on
%             xline(aversiveTS_run(1)/1000,'r')
%             xline(aversiveTS_run(2)/1000,'r')
%             xline(rewardTS_run(1)/1000,'b')
%             xline(rewardTS_run(2)/1000,'b')
%             yline(aversiveTS_run(1)/1000,'r')
%             yline(aversiveTS_run(2)/1000,'r')
%             yline(rewardTS_run(1)/1000,'b')
%             yline(rewardTS_run(2)/1000,'b')
%             title('vHPC')             


        clear A aversiveTS aversiveTS_run baselineTS rewardTS rewardTS_run
        clear behavior bins Cell_type_classification cellulartype cond
        clear is K Kinfo group_dHPC group_vHPC matrixC matrixD matrixV
        clear NREM REM WAKE segmentation tmp
        clear spiketrains_dHPC spiketrains_vHPC
        end
    end
end

figure,
subplot(221),plot(mean(output.aversive.first.aversive,'omitnan'),'r')
subplot(222),plot(mean(output.aversive.first.reward,'omitnan'),'b')
subplot(223),plot(mean(output.aversive.second.aversive,'omitnan'),'r')
subplot(224),plot(mean(output.aversive.second.reward,'omitnan'),'b')
