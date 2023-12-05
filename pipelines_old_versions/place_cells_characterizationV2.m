clear
clc
close all
%% Parameters
path = {'Z:\Facundo\Rat127\Ephys\pyr';'Z:\Facundo\Rat128\Ephys\in_pyr\ready';'Z:\Facundo\Rat103\usable';'Z:\Facundo\Rat132\recordings\in_pyr'};%List of folders from the path
% path = {'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat103\usable';'E:\Rat132\recordings\in_pyr'};%List of folders from the path
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
minimal_speed = 7; % minimal speed to detect quite periods
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
        maps_dHPC_SA = [];        maps_dHPC_SR = [];
        peak_dHPC_A = [] ;       peak_dHPC_R = [];
        field_dHPC_A = [];       field_dHPC_R = [];
        PC.dHPC.aversive=[]; PC.dHPC.reward=[];
        cluster_dHPC = [];
        for ii=1:size(group_dHPC,1)
            cluster = group_dHPC(ii,1);
            celltype = logical(cellulartype(cellulartype(:,1) == cluster,2));
            if celltype % check if pyr
                spks = spks_dHPC(spks_dHPC(:,1)==cluster,2);
                b = [movement.aversive ; movement.reward];
                tmp = Restrict(spks,b);
                m = length(tmp) / sum(b(:,2) - b(:,1)); % check if the SU fires
                
                if ~isnan(m)
                    % --- Aversive ---
                    spks_tmp = Restrict(spks , movement.aversive);
                    pos_tmp = Restrict(behavior.pos.aversive(:,1:2) , movement.aversive); % Restrict to movement periods
                    pos_tmp(:,2) = pos_tmp(:,2)-min(pos_tmp(:,2)); pos_tmp(:,2) = pos_tmp(:,2)/max(pos_tmp(:,2)); %normalization of position
                    
                    %Firing curve construction
                    [curveA , statsA] = FiringCurve(pos_tmp , spks_tmp , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 2 , 'minPeak' , 0.2);
                    [curveSA , statsSA] = FiringCurve(pos_tmp , ShuffleSpks(spks_tmp) , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 2 , 'minPeak' , 0.2);
                    %Quantile for place-cell definition
                    qA = SkaggsRandomFMT(spks_tmp, pos_tmp, sigma ,Xedges);
                    
                    % Store of place-field location and cluster of PCs
                    clear s q curve OccMap Nspikes spks_tmp pos_tmp m curve stats sp curve1
                    
                    % --- Reward ---
                    spks_tmp = Restrict(spks , movement.reward); % Restrict to movement periods
                    pos_tmp = Restrict(behavior.pos.reward(:,1:2) , movement.reward);
                    pos_tmp(:,2) = pos_tmp(:,2)-min(pos_tmp(:,2)); pos_tmp(:,2) = pos_tmp(:,2)/max(pos_tmp(:,2)); %normalization of position
                    
                    %Firing curve construction
                    [curveR , statsR] = FiringCurve(pos_tmp , spks_tmp , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 2 , 'minPeak' , 0.2);
                    [curveSR , statsSR] = FiringCurve(pos_tmp , ShuffleSpks(spks_tmp) , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 2 , 'minPeak' , 0.2);
                    %Quantile for place-cell definition
                    qR = SkaggsRandomFMT(spks_tmp, pos_tmp, sigma, Xedges);
                    
                    % Store of place-field location and cluster of PCs
                    
                    if or(statsA.specificity > qA , statsR.specificity > qR)
                        
                        [ff,f] = max(curveA.rate);
                        fff = [1/Xedges:1/Xedges:1];
                        PC.dHPC.aversive=[PC.dHPC.aversive ; cluster , ff , fff(f)];
                        clear f ff fff
                        maps_dHPC_A = [maps_dHPC_A ; curveA.rate];
                        maps_dHPC_SA = [maps_dHPC_SA ; curveSA.rate];
                        peak_dHPC_A = [peak_dHPC_A ; statsA.peak(1)] ;
                        field_dHPC_A = [field_dHPC_A ; statsA.fieldX(1,:)];
                        
                        [ff,f] = max(curveR.rate);
                        fff = [1/Xedges:1/Xedges:1];
                        PC.dHPC.reward=[PC.dHPC.reward ; cluster , ff , fff(f)];
                        clear f ff fff
                        maps_dHPC_R = [maps_dHPC_R ; curveR.rate];
                        maps_dHPC_SR = [maps_dHPC_SR ; curveSR.rate];
                        peak_dHPC_R = [peak_dHPC_R ; statsR.peak(1)] ;
                        field_dHPC_R = [field_dHPC_R ; statsR.fieldX(1,:)];
                        cluster_dHPC = [cluster_dHPC ; cluster , cond];
                    end
                    clear s q curve OccMap Nspikes spks_tmp pos_tmp m curve stats sp curve1
                    clear curveA curveR qA qR curveSA curveSR statsSA statsSR
                end
            end
            clear celltype tmp b cluster
        end
        
        disp('vHPC Firing rate map calculation')
        maps_vHPC_A = [];        maps_vHPC_R = [];
        maps_vHPC_SA = [];        maps_vHPC_SR = [];
        peak_vHPC_A = [] ;       peak_vHPC_R = [];
        field_vHPC_A = [];       field_vHPC_R = [];
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
                if ~isnan(m)
                    % --- Aversive ---
                    spks_tmp = Restrict(spks , movement.aversive);
                    pos_tmp = Restrict(behavior.pos.aversive(:,1:2) , movement.aversive);
                    pos_tmp(:,2) = pos_tmp(:,2)-min(pos_tmp(:,2)); pos_tmp(:,2) = pos_tmp(:,2)/max(pos_tmp(:,2)); %normalization of position
                    
                    %Firing curve construction
                    [curveA , statsA] = FiringCurve(pos_tmp , spks_tmp , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 2 , 'minPeak' , 0.2);
                    [curveSA , statsSA] = FiringCurve(pos_tmp , ShuffleSpks(spks_tmp) , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 2 , 'minPeak' , 0.2);
                    %Quantile for place-cell definition
                    qA = SkaggsRandomFMT(spks_tmp, pos_tmp, sigma ,Xedges);
                    
                    % Store of place-field location and cluster of PCs
                    clear s q curve OccMap Nspikes spks_tmp pos_tmp m curve stats sp curve1
                    
                    % --- Reward ---
                    spks_tmp = Restrict(spks , movement.reward);
                    pos_tmp = Restrict(behavior.pos.reward(:,1:2) , movement.reward);
                    pos_tmp(:,2) = pos_tmp(:,2)-min(pos_tmp(:,2)); pos_tmp(:,2) = pos_tmp(:,2)/max(pos_tmp(:,2)); %normalization of position
                    
                    %Firing curve construction
                    [curveR , statsR] = FiringCurve(pos_tmp , spks_tmp , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 2 , 'minPeak' , 0.2);
                    [curveSR , statsSR] = FiringCurve(pos_tmp , ShuffleSpks(spks_tmp) , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 2 , 'minPeak' , 0.2);
                    %Quantile for place-cell definition
                    qR = SkaggsRandomFMT(spks_tmp, pos_tmp, sigma, Xedges);
                    
                    % Store of place-field location and cluster of PCs
                    if or(statsA.specificity > qA , statsR.specificity > qR)
                        [ff,f] = max(curveA.rate);
                        fff = [1/Xedges:1/Xedges:1];
                        PC.vHPC.aversive=[PC.vHPC.aversive ; cluster , ff , fff(f)];
                        clear f ff fff
                        maps_vHPC_A = [maps_vHPC_A ; curveA.rate];
                        maps_vHPC_SA = [maps_vHPC_SA ; curveSA.rate];
                        peak_vHPC_A = [peak_vHPC_A ; statsA.peak(1)] ;
                        field_vHPC_A = [field_vHPC_A ; statsA.fieldX(1,:)];
                        
                        [ff,f] = max(curveR.rate);
                        fff = [1/Xedges:1/Xedges:1];
                        PC.vHPC.reward=[PC.vHPC.reward ; cluster , ff , fff(f)];
                        clear f ff fff
                        maps_vHPC_R = [maps_vHPC_R ; curveR.rate];
                        maps_vHPC_SR = [maps_vHPC_SR ; curveSR.rate];
                        peak_vHPC_R = [peak_vHPC_R ; statsR.peak(1)] ;
                        field_vHPC_R = [field_vHPC_R ; statsR.fieldX(1,:)];
                        cluster_vHPC = [cluster_vHPC ; cluster , cond];
                    end
                    clear s q curve OccMap Nspikes spks_tmp pos_tmp m curve stats sp curve1
                    clear curveA curveR qA qR curveSA curveSR statsSA statsSR
                end
            end
            clear celltype tmp b cluster
        end
        

            %% Saving data
            disp('Saving outputs')
            if not(exist('maps','var'))
                maps.dHPC.aversive = maps_dHPC_A;
                maps.vHPC.aversive = maps_vHPC_A;
                maps.dHPC.reward   = maps_dHPC_R;
                maps.vHPC.reward   = maps_vHPC_R;
                maps.dHPC.aversiveS = maps_dHPC_SA;
                maps.vHPC.aversiveS = maps_vHPC_SA;
                maps.dHPC.rewardS   = maps_dHPC_SR;
                maps.vHPC.rewardS   = maps_vHPC_SR;
                peak.dHPC.aversive = peak_dHPC_A;
                peak.vHPC.aversive = peak_vHPC_A;
                peak.dHPC.reward   = peak_dHPC_R;
                peak.vHPC.reward   = peak_vHPC_R;
                field.dHPC.aversive = field_dHPC_A;
                field.vHPC.aversive = field_vHPC_A;
                field.dHPC.reward   = field_dHPC_R;
                field.vHPC.reward   = field_vHPC_R;
                clus.dHPC = cluster_dHPC;
                clus.vHPC = cluster_vHPC;
            else
                maps.dHPC.aversive = [maps.dHPC.aversive ; maps_dHPC_A];
                maps.vHPC.aversive = [maps.vHPC.aversive ; maps_vHPC_A];
                maps.dHPC.reward   = [maps.dHPC.reward ; maps_dHPC_R];
                maps.vHPC.reward   = [maps.vHPC.reward ; maps_vHPC_R];
                maps.dHPC.aversiveS = [maps.dHPC.aversiveS ; maps_dHPC_SA];
                maps.vHPC.aversiveS = [maps.vHPC.aversiveS ; maps_vHPC_SA];
                maps.dHPC.rewardS   = [maps.dHPC.rewardS ; maps_dHPC_SR];
                maps.vHPC.rewardS   = [maps.vHPC.rewardS ; maps_vHPC_SR];                
                peak.dHPC.aversive = [peak.dHPC.aversive ; peak_dHPC_A];
                peak.vHPC.aversive = [peak.vHPC.aversive ; peak_vHPC_A];
                peak.dHPC.reward   = [peak.dHPC.reward ; peak_dHPC_R];
                peak.vHPC.reward   = [peak.vHPC.reward ; peak_vHPC_R];
                field.dHPC.aversive = [field.dHPC.aversive ; field_dHPC_A];
                field.vHPC.aversive = [field.vHPC.aversive ; field_vHPC_A];
                field.dHPC.reward   = [field.dHPC.reward ; field_dHPC_R];
                field.vHPC.reward   = [field.vHPC.reward ; field_vHPC_R];
                clus.dHPC = [clus.dHPC ; cluster_dHPC];
                clus.vHPC = [clus.vHPC ; cluster_vHPC];                
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
        disp(['-- Finished folder #' , num2str(t) , ' from rat #' , num2str(tt) , ' ---'])
        disp('  ')
    end
    disp(['-------- Finished rat#' , num2str(tt) , ' --------'])
    disp(' ')
end


%% Aversive  --> Reward
criteria = clus.dHPC(:,2) == 1;
tmpA = [];
MA = [];
tmpR = [];
MR = [];
SpatialCorr.dHPC = [];
for i = 1 : size(maps.dHPC.aversive,1)
    if logical(criteria(i))
        ii = maps.dHPC.aversive(i,:) - min(maps.dHPC.aversive(i,:));
        ii = ii ./ max(ii);
        [trash , iii] = max(ii);
        MA = [MA ; iii];
        tmpA = [tmpA ; ii]; clear ii trash iii
        
        
        ii = maps.dHPC.reward(i,:) - min(maps.dHPC.reward(i,:));
        ii = ii ./ max(ii);
        [trash , iii] = max(ii);
        MR = [MR ; iii];
        tmpR = [tmpR ; ii]; clear ii trash iii
        c = corrcoef(maps.dHPC.aversive(i,:) , maps.dHPC.reward(i,:));
        ca = corrcoef(maps.dHPC.aversive(i,:) , maps.dHPC.rewardS(i,:));
        cr = corrcoef(maps.dHPC.reward(i,:) , maps.dHPC.aversiveS(i,:));
        SpatialCorr.dHPC = [SpatialCorr.dHPC ; c(1,2)/(ca(1,2) + cr(1,2))]; clear c ca cr
    end
end
clear i
[i ii] = sort(MA);

figure,
subplot(121),imagesc(tmpA(ii,:))
subplot(122),imagesc(tmpR(ii,:))


criteria = clus.vHPC(:,2) == 1;
tmpA = [];
MA = [];
tmpR = [];
MR = [];
SpatialCorr.vHPC = [];
for i = 1 : size(maps.vHPC.aversive,1)
    if logical(criteria(i))
        ii = maps.vHPC.aversive(i,:) - min(maps.vHPC.aversive(i,:));
        ii = ii ./ max(ii);
        [trash , iii] = max(ii);
        MA = [MA ; iii];
        tmpA = [tmpA ; ii]; clear ii trash iii
        
        
        ii = maps.vHPC.reward(i,:) - min(maps.vHPC.reward(i,:));
        ii = ii ./ max(ii);
        [trash , iii] = max(ii);
        MR = [MR ; iii];
        tmpR = [tmpR ; ii]; clear ii trash iii
        
        c = corrcoef(maps.vHPC.aversive(i,:) , maps.vHPC.reward(i,:));
        ca = corrcoef(maps.vHPC.aversive(i,:) , maps.vHPC.rewardS(i,:));
        cr = corrcoef(maps.vHPC.reward(i,:) , maps.vHPC.aversiveS(i,:));
        SpatialCorr.vHPC = [SpatialCorr.dHPC ; c(1,2)/(ca(1,2) + cr(1,2))]; clear c ca cr
    end
end
clear i

[i ii] = sort(MA);

figure,
subplot(121),imagesc(tmpA(ii,:))
subplot(122),imagesc(tmpR(ii,:))

figure,
x = ones(size(SpatialCorr.dHPC));
y = ones(size(SpatialCorr.vHPC))*2;
boxplot([SpatialCorr.dHPC ; SpatialCorr.vHPC],[x;y]),ylim([-3 3])

%% Reward --> Aversive
criteria = clus.dHPC(:,2) == 2;
tmpA = [];
MA = [];
tmpR = [];
MR = [];
SpatialCorr.dHPC = [];
for i = 1 : size(maps.dHPC.aversive,1)
    if logical(criteria(i))
        ii = maps.dHPC.aversive(i,:) - min(maps.dHPC.aversive(i,:));
        ii = ii ./ max(ii);
        [trash , iii] = max(ii);
        MA = [MA ; iii];
        tmpA = [tmpA ; ii]; clear ii trash iii
        
        
        ii = maps.dHPC.reward(i,:) - min(maps.dHPC.reward(i,:));
        ii = ii ./ max(ii);
        [trash , iii] = max(ii);
        MR = [MR ; iii];
        tmpR = [tmpR ; ii]; clear ii trash iii
        c = corrcoef(maps.dHPC.aversive(i,:) , maps.dHPC.reward(i,:));
        SpatialCorr.dHPC = [SpatialCorr.dHPC ; c(1,2)]; clear c
    end
end
clear i
[i ii] = sort(MA);

figure,
subplot(121),imagesc(tmpA(ii,:))
subplot(122),imagesc(tmpR(ii,:))


criteria = clus.vHPC(:,2) == 2;
tmpA = [];
MA = [];
tmpR = [];
MR = [];
SpatialCorr.vHPC = [];
for i = 1 : size(maps.vHPC.aversive,1)
    if logical(criteria(i))
        ii = maps.vHPC.aversive(i,:) - min(maps.vHPC.aversive(i,:));
        ii = ii ./ max(ii);
        [trash , iii] = max(ii);
        MA = [MA ; iii];
        tmpA = [tmpA ; ii]; clear ii trash iii
        
        
        ii = maps.vHPC.reward(i,:) - min(maps.vHPC.reward(i,:));
        ii = ii ./ max(ii);
        [trash , iii] = max(ii);
        MR = [MR ; iii];
        tmpR = [tmpR ; ii]; clear ii trash iii
        
        c = corrcoef(maps.vHPC.aversive(i,:) , maps.vHPC.reward(i,:));
        SpatialCorr.vHPC = [SpatialCorr.vHPC ; c(1,2)]; clear c
    end
end
clear i

[i ii] = sort(MA);

figure,
subplot(121),imagesc(tmpA(ii,:))
subplot(122),imagesc(tmpR(ii,:))

figure,
x = ones(size(SpatialCorr.dHPC));
y = ones(size(SpatialCorr.vHPC))*2;
boxplot([SpatialCorr.dHPC ; SpatialCorr.vHPC],[x;y])
[h p] = ttest2(SpatialCorr.dHPC , SpatialCorr.vHPC);

%% All
tmpA = [];
MA = [];
tmpR = [];
MR = [];
SpatialCorr.dHPC = [];
for i = 1 : size(maps.dHPC.aversive,1)
    ii = maps.dHPC.aversive(i,:) - min(maps.dHPC.aversive(i,:));
    ii = ii ./ max(ii);
    [trash , iii] = max(ii);
    MA = [MA ; iii];
    tmpA = [tmpA ; ii]; clear ii trash iii
    
    
    ii = maps.dHPC.reward(i,:) - min(maps.dHPC.reward(i,:));
    ii = ii ./ max(ii);
    [trash , iii] = max(ii);
    MR = [MR ; iii];
    tmpR = [tmpR ; ii]; clear ii trash iii  
    c = corrcoef(maps.dHPC.aversive(i,:) , maps.dHPC.reward(i,:));
    SpatialCorr.dHPC = [SpatialCorr.dHPC ; c(1,2)]; clear c
end
clear i
[i ii] = sort(MA);

figure,
subplot(121),imagesc(tmpA(ii,:))
subplot(122),imagesc(tmpR(ii,:))

tmpA = [];
MA = [];
tmpR = [];
MR = [];
SpatialCorr.vHPC = [];
for i = 1 : size(maps.vHPC.aversive,1)
    ii = maps.vHPC.aversive(i,:) - min(maps.vHPC.aversive(i,:));
    ii = ii ./ max(ii);
    [trash , iii] = max(ii);
    MA = [MA ; iii];
    tmpA = [tmpA ; ii]; clear ii trash iii
    
    
    ii = maps.vHPC.reward(i,:) - min(maps.vHPC.reward(i,:));
    ii = ii ./ max(ii);
    [trash , iii] = max(ii);
    MR = [MR ; iii];
    tmpR = [tmpR ; ii]; clear ii trash iii   
    
    c = corrcoef(maps.vHPC.aversive(i,:) , maps.vHPC.reward(i,:));
    SpatialCorr.vHPC = [SpatialCorr.vHPC ; c(1,2)]; clear c
end
clear i

[i ii] = sort(MA);

figure,
subplot(121),imagesc(tmpA(ii,:))
subplot(122),imagesc(tmpR(ii,:))

figure,
x = ones(size(SpatialCorr.dHPC));
y = ones(size(SpatialCorr.vHPC))*2;
boxplot([SpatialCorr.dHPC ; SpatialCorr.vHPC],[x;y])
