clear
clc
close all

%% Parameters
path = {'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat103\usable';'E:\Rat132\recordings\in_pyr'};%List of folders from the path

%Sleep
time_criteria = 600; %time criteria to define the maximal time of sleep to include

% for SU
criteria_fr = 0.1; %criteria to include or not a SU into the analysis
criteria_n = 6; % minimal number of neurons from each structure
criteria_type = 0; %criteria for celltype (0:pyr, 1:int, 2:all)
binSize = 0.05;

% Behavior
minimal_speed = 7; % minimal speed to detect quite periods
minimal_speed_time = 2; % minimal time to detect quite periods

n_SU_V = 0;
n_SU_D = 0;

EV.aversive.dHPC = [];    EV.reward.dHPC = [];       
EV.aversive.vHPC = [];    EV.reward.vHPC = [];
EV.aversive.dHPCi = [];    EV.reward.dHPCi = [];       
EV.aversive.vHPCi = [];    EV.reward.vHPCi = [];
EV.aversive.dvHPC = [];   EV.reward.dvHPC = [];
EV.aversive.dvHPCi = [];   EV.reward.dvHPCi = [];
FiringRate.dHPC.baseline = [];
FiringRate.dHPC.reward = [];
FiringRate.dHPC.aversive = [];
FiringRate.vHPC.baseline = [];
FiringRate.vHPC.reward = [];
FiringRate.vHPC.aversive = [];
%% Main lo op, to iterate across sessions
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
        
                %% Constructing Spiketrains
        freq = 1/binSize;
        limits = [0 segments.Var1(end)/1000];
        spiketrains_dHPC.pyr = [];        spiketrains_dHPC.int = [];
        spiketrains_vHPC.pyr = [];        spiketrains_vHPC.int = [];
        for ii=1:size(group_dHPC,1)
            cluster = group_dHPC(ii,1);
            Cellulartype = logical(Cell_type_classification(Cell_type_classification(:,1) == cluster,6));
            if Cellulartype
                spks = spks_dHPC(spks_dHPC(:,1)==cluster,2);
                [tmp,bins]=binspikes(spks,freq,limits);
                is.sleep.baseline = InIntervals(bins,baselineTS./1000);                 is.sleep.reward = InIntervals(bins,rewardTS./1000);
                is.sleep.aversive = InIntervals(bins,aversiveTS./1000);                 is.run.reward = InIntervals(bins,rewardTS_run./1000);
                is.run.aversive = InIntervals(bins,aversiveTS_run./1000);
                fr1 = sum(tmp(is.sleep.baseline))/(baselineTS(2)/1000);
                fr2 = sum(tmp(is.sleep.reward))/((rewardTS(2)-rewardTS(1))/1000);
                fr3 = sum(tmp(is.sleep.aversive))/((aversiveTS(2)-aversiveTS(1))/1000);
                fr4 = sum(tmp(is.run.reward))/((rewardTS_run(2)-rewardTS_run(1))/1000);
                fr5 = sum(tmp(is.run.aversive))/((aversiveTS_run(2)-aversiveTS_run(1))/1000);
                
                FiringRate.dHPC.baseline = [FiringRate.dHPC.baseline ; fr1];
                FiringRate.dHPC.reward = [FiringRate.dHPC.reward ; fr2];
                FiringRate.dHPC.aversive = [FiringRate.dHPC.aversive ; fr3];
                
                
                if and(and(fr1>criteria_fr , fr2 > criteria_fr),fr3>criteria_fr)
                    spiketrains_dHPC.pyr = [spiketrains_dHPC.pyr , zscore(tmp)];
                end
            else
                spks = spks_dHPC(spks_dHPC(:,1)==cluster,2);
                [tmp,bins]=binspikes(spks,freq,limits);
                spiketrains_dHPC.int = [spiketrains_dHPC.int , zscore(tmp)];
            end
            clear spks tmp cluster Cellulartype fr1 fr2 fr3 fr4 fr5
        end
        
        for ii=1:size(group_vHPC,1)
            cluster = group_vHPC(ii,1);
            Cellulartype = logical(Cell_type_classification(Cell_type_classification(:,1) == cluster,6));
            if Cellulartype
                spks = spks_vHPC(spks_vHPC(:,1)==cluster,2);
                [tmp,bins]=binspikes(spks,freq,limits);
                is.sleep.baseline = InIntervals(bins,baselineTS./1000);                 is.sleep.reward = InIntervals(bins,rewardTS./1000);
                is.sleep.aversive = InIntervals(bins,aversiveTS./1000);                 is.run.reward = InIntervals(bins,rewardTS_run./1000);
                is.run.aversive = InIntervals(bins,aversiveTS_run./1000);
                fr1 = sum(tmp(is.sleep.baseline))/(baselineTS(2)/1000);
                fr2 = sum(tmp(is.sleep.reward))/((rewardTS(2)-rewardTS(1))/1000);
                fr3 = sum(tmp(is.sleep.aversive))/((aversiveTS(2)-aversiveTS(1))/1000);
                fr4 = sum(tmp(is.run.reward))/((rewardTS_run(2)-rewardTS_run(1))/1000);
                fr5 = sum(tmp(is.run.aversive))/((aversiveTS_run(2)-aversiveTS_run(1))/1000);
                
                FiringRate.vHPC.baseline = [FiringRate.vHPC.baseline ; fr1];
                FiringRate.vHPC.reward = [FiringRate.vHPC.reward ; fr2];
                FiringRate.vHPC.aversive = [FiringRate.vHPC.aversive ; fr3];                
                
                if and(and(fr1>criteria_fr , fr2 > criteria_fr),fr3>criteria_fr)
                    spiketrains_vHPC.pyr = [spiketrains_vHPC.pyr , zscore(tmp)];
                end
            else
                spks = spks_vHPC(spks_vHPC(:,1)==cluster,2);
                [tmp,bins]=binspikes(spks,freq,limits);
                spiketrains_vHPC.int = [spiketrains_vHPC.int , zscore(tmp)];
            end
            clear spks tmp cluster Cellulartype fr1 fr2 fr3 fr4 fr5
        end
        clear freq limits
        clear spks spks_dHPC spks_vHPC camara shock rightvalve leftvalve
        clear ejeX ejeY dX dY dX_int dY_int
        
        %% Explained variance calculation
        if and(size(spiketrains_vHPC.pyr,2) >= criteria_n,size(spiketrains_dHPC.pyr,2) >= criteria_n)
            disp('Lets go for the SUs')
            %Restricting bins inside each condition
            is.baseline.sws = InIntervals(bins,NREM.baseline); clear tmp
            is.reward.sws = InIntervals(bins,NREM.reward);  clear tmp
            is.aversive.sws = InIntervals(bins,NREM.aversive);  clear tmp
            
%             is.baseline.sws = InIntervals(bins,[ripple_bursts.baseline(:,1) ripple_bursts.baseline(:,3)]);
%             is.aversive.sws = InIntervals(bins,[ripple_bursts.aversive(:,1) ripple_bursts.aversive(:,3)]);
%             is.reward.sws = InIntervals(bins,[ripple_bursts.reward(:,1) ripple_bursts.reward(:,3)]);            
%             
%             is.baseline.sws = InIntervals(bins,ripple_event.baseline);
%             is.aversive.sws = InIntervals(bins,ripple_event.aversive);
%             is.reward.sws = InIntervals(bins,ripple_event.reward);  
            
%             is.baseline.sws = InIntervals(bins,[ripples.vHPC.baseline(:,1) ripples.vHPC.baseline(:,3)]);
%             is.aversive.sws = InIntervals(bins,[ripples.vHPC.aversive(:,1) ripples.vHPC.aversive(:,3)]);
%             is.reward.sws = InIntervals(bins,[ripples.vHPC.reward(:,1) ripples.vHPC.reward(:,3)]);  
            
            is.baseline.rem = InIntervals(bins,REM.baseline);
            is.aversive.rem = InIntervals(bins,REM.aversive);
            is.reward.rem = InIntervals(bins,REM.reward);
            
            is.aversive.run = InIntervals(bins,movement.aversive);
            is.reward.run = InIntervals(bins,movement.reward);
            
            if and(and(~isempty(NREM.aversive),~isempty(NREM.reward)),~isempty(NREM.baseline))
                if aversiveTS_run(1)>rewardTS_run(1)
%                     % dHPC
%                     % Reward
%                     x = [spiketrains_dHPC.pyr(is.baseline.sws,:)];
%                     y = [spiketrains_dHPC.pyr(is.baseline.sws,:)];
%                     S1 = x'*y/size(x,1);
% %                     [S1 , p]=corr(x,y);
% %                     S1 = S1 - diag(diag(S1));
%                     clear x y p
%                     
%                     x = [spiketrains_dHPC.pyr(is.reward.run,:)];
%                     y = [spiketrains_dHPC.pyr(is.reward.run,:)];
%                     S2= x'*y/size(x,1);
% %                     [S2 , p]=corr(x,y);
% %                     S2 = S2 - diag(diag(S2));
%                     clear x y p
%                     
%                     x = [spiketrains_dHPC.pyr(is.reward.sws,:)];
%                     y = [spiketrains_dHPC.pyr(is.reward.sws,:)];
%                     S3= x'*y/size(x,1);
% %                     [S3 , p]=corr(x,y);
% %                     S3 = S3 - diag(diag(S3));
%                     clear x y p
%                     
%                     % EV and REV
%                     Sx = corrcoef(S2,S3,'rows','complete');     Sx = Sx(1,2);
%                     Sy = corrcoef(S2,S1,'rows','complete');     Sy = Sy(1,2);
%                     Sz = corrcoef(S3,S1,'rows','complete');     Sz = Sz(1,2);
%                     
%                     ev = (Sx-Sy*Sz);
%                     ev = (ev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
%                     
%                     rev = (Sy-Sx*Sz);
%                     rev = (rev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
%                     
%                     EV.reward.dHPC = [EV.reward.dHPC ; rev , ev];
%                     clear Sx Sy Sz ev rev
%                     
%                     % Aversive
%                     x = [spiketrains_dHPC.pyr(is.aversive.run,:)];
%                     y = [spiketrains_dHPC.pyr(is.aversive.run,:)];
%                     S4= x'*y/size(x,1);
% %                     [S4 , p]=corr(x,y);
% %                     S4 = S4 - diag(diag(S4));
%                     clear x y p
%                     
%                     x = [spiketrains_dHPC.pyr(is.aversive.sws,:)];
%                     y = [spiketrains_dHPC.pyr(is.aversive.sws,:)];
%                     S5= x'*y/size(x,1);
% %                     [S5 , p]=corr(x,y);
% %                     S5 = S5 - diag(diag(S5));
%                     clear x y p
%                     
%                     % EV and REV
%                     Sx = corrcoef(S4,S5,'rows','complete');     Sx = Sx(1,2);
%                     Sy = corrcoef(S4,S3,'rows','complete');     Sy = Sy(1,2);
%                     Sz = corrcoef(S5,S3,'rows','complete');     Sz = Sz(1,2);
%                     
%                     ev = (Sx-Sy*Sz);
%                     ev = (ev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
%                     
%                     rev = (Sy-Sx*Sz);
%                     rev = (rev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
%                     
%                     EV.aversive.dHPC = [EV.aversive.dHPC ; rev , ev];
%                     clear Sx Sy Sz ev rev
%                     
%                     % Reward in Aversive
%                     Sx = corrcoef(S2,S5,'rows','complete');     Sx = Sx(1,2);
%                     Sy = corrcoef(S2,S1,'rows','complete');     Sy = Sy(1,2);
%                     Sz = corrcoef(S5,S1,'rows','complete');     Sz = Sz(1,2);
%                     
%                     ev = (Sx-Sy*Sz);
%                     ev = (ev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
%                     
%                     rev = (Sy-Sx*Sz);
%                     rev = (rev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
%                     
%                     EV.reward.dHPCi = [EV.reward.dHPCi ; rev , ev];
%                     clear Sx Sy Sz ev rev S1 S2 S3 S4 S5                    
%                     
%                     % vHPC
%                     % Reward
%                     x = [spiketrains_vHPC.pyr(is.baseline.sws,:)];
%                     y = [spiketrains_vHPC.pyr(is.baseline.sws,:)];
%                     S1= x'*y/size(x,1);
% %                     [S1 , p]=corr(x,y);
% %                     S1 = S1 - diag(diag(S1));
%                     clear x y p
%                     
%                     x = [spiketrains_vHPC.pyr(is.reward.run,:)];
%                     y = [spiketrains_vHPC.pyr(is.reward.run,:)];
%                     S2= x'*y/size(x,1);
% %                     [S2 , p]=corr(x,y);
% %                     S2 = S2 - diag(diag(S2));
%                     clear x y p
%                     
%                     x = [spiketrains_vHPC.pyr(is.reward.sws,:)];
%                     y = [spiketrains_vHPC.pyr(is.reward.sws,:)];
%                     S3= x'*y/size(x,1);
% %                     [S3 , p]=corr(x,y);
% %                     S3 = S3 - diag(diag(S3));
%                     clear x y p
%                     
%                     % EV and REV
%                     Sx = corrcoef(S2,S3,'rows','complete');     Sx = Sx(1,2);
%                     Sy = corrcoef(S2,S1,'rows','complete');     Sy = Sy(1,2);
%                     Sz = corrcoef(S3,S1,'rows','complete');     Sz = Sz(1,2);
%                     
%                     ev = (Sx-Sy*Sz);
%                     ev = (ev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
%                     
%                     rev = (Sy-Sx*Sz);
%                     rev = (rev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
%                     
%                     EV.reward.vHPC = [EV.reward.vHPC ; rev , ev];
%                     clear Sx Sy Sz ev rev
%                                         
%                     % Aversive
%                     x = [spiketrains_vHPC.pyr(is.aversive.run,:)];
%                     y = [spiketrains_vHPC.pyr(is.aversive.run,:)];
%                     S4= x'*y/size(x,1);
% %                     [S4 , p]=corr(x,y);
% %                     S4 = S4 - diag(diag(S4));
%                     clear x y p
%                     
%                     x = [spiketrains_vHPC.pyr(is.aversive.sws,:)];
%                     y = [spiketrains_vHPC.pyr(is.aversive.sws,:)];
%                     S5= x'*y/size(x,1);
% %                     [S5 , p]=corr(x,y);
% %                     S5 = S5 - diag(diag(S5));
%                     clear x y p
%                     
%                     % EV and REV
%                     Sx = corrcoef(S4,S5,'rows','complete');     Sx = Sx(1,2);
%                     Sy = corrcoef(S4,S3,'rows','complete');     Sy = Sy(1,2);
%                     Sz = corrcoef(S5,S3,'rows','complete');     Sz = Sz(1,2);
%                     
%                     ev = (Sx-Sy*Sz);
%                     ev = (ev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
%                     
%                     rev = (Sy-Sx*Sz);
%                     rev = (rev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
%                     
%                     EV.aversive.vHPC = [EV.aversive.vHPC ; rev , ev];
%                     clear Sx Sy Sz ev rev
% 
%                     % Reward in Aversive
%                     Sx = corrcoef(S2,S5,'rows','complete');     Sx = Sx(1,2);
%                     Sy = corrcoef(S2,S1,'rows','complete');     Sy = Sy(1,2);
%                     Sz = corrcoef(S5,S1,'rows','complete');     Sz = Sz(1,2);
%                     
%                     ev = (Sx-Sy*Sz);
%                     ev = (ev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
%                     
%                     rev = (Sy-Sx*Sz);
%                     rev = (rev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
%                     
%                     EV.reward.vHPCi = [EV.reward.vHPCi ; rev , ev];
%                     clear Sx Sy Sz ev rev S1 S2 S3 S4 S5                    
                    
                    % dHPC and vHPC
                    % Reward
                    x = [spiketrains_dHPC.pyr(is.baseline.sws,:)];
                    y = [spiketrains_vHPC.pyr(is.baseline.sws,:)];                  
                    S1= x'*y./size(x,1);
%                     [S1 , p]=corr(x,y);
                    clear x y p
                    
                    x = [spiketrains_dHPC.pyr(is.reward.run,:)];
                    y = [spiketrains_vHPC.pyr(is.reward.run,:)];
                    S2= x'*y/size(x,1);
%                     [S2 , p]=corr(x,y);
                    clear x y p
                    
                    x = [spiketrains_dHPC.pyr(is.reward.sws,:)];
                    y = [spiketrains_vHPC.pyr(is.reward.sws,:)];
                    S3= x'*y/size(x,1);
%                     [S3 , p]=corr(x,y);
                    clear x y p
                    
                    % EV and REV
                    Sx = corrcoef(S2,S3,'rows','complete');     Sx = Sx(1,2);
                    Sy = corrcoef(S2,S1,'rows','complete');     Sy = Sy(1,2);
                    Sz = corrcoef(S3,S1,'rows','complete');     Sz = Sz(1,2);
                    
                    ev = (Sx-Sy*Sz);
                    ev = (ev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
                    
                    rev = (Sy-Sx*Sz);
                    rev = (rev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
                    
                    EV.reward.dvHPC = [EV.reward.dvHPC ; rev , ev];
                    clear Sx Sy Sz ev rev
                                        
                    % Aversive
                    x = [spiketrains_dHPC.pyr(is.aversive.run,:)];
                    y = [spiketrains_vHPC.pyr(is.aversive.run,:)];
                    S4= x'*y/size(x,1);
%                     [S4 , p]=corr(x,y);
                    clear x y p
                    
                    x = [spiketrains_dHPC.pyr(is.aversive.sws,:)];
                    y = [spiketrains_vHPC.pyr(is.aversive.sws,:)];
                    S5= x'*y/size(x,1);
%                     [S5 , p]=corr(x,y);
                    clear x y p
                    
                    % EV and REV
                    Sx = corrcoef(S4,S5,'rows','complete');     Sx = Sx(1,2);
                    Sy = corrcoef(S4,S3,'rows','complete');     Sy = Sy(1,2);
                    Sz = corrcoef(S5,S3,'rows','complete');     Sz = Sz(1,2);
                    
                    ev = (Sx-Sy*Sz);
                    ev = (ev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
                    
                    rev = (Sy-Sx*Sz);
                    rev = (rev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
                    
%                     EV.aversive.dvHPC = [EV.aversive.dvHPC ; rev , ev];
                    clear Sx Sy Sz ev rev

                    % Reward in Aversive
                    Sx = corrcoef(S2,S5,'rows','complete');     Sx = Sx(1,2);
                    Sy = corrcoef(S2,S1,'rows','complete');     Sy = Sy(1,2);
                    Sz = corrcoef(S5,S1,'rows','complete');     Sz = Sz(1,2);
                    
                    ev = (Sx-Sy*Sz);
                    ev = (ev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
                    
                    rev = (Sy-Sx*Sz);
                    rev = (rev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
                    
                    EV.reward.dvHPCi = [EV.reward.dvHPCi ; rev , ev];
                    clear Sx Sy Sz ev rev S1 S2 S3 S4 S5                     
                    
                else
%                     % dHPC
%                     % Aversive
%                     x = [spiketrains_dHPC.pyr(is.baseline.sws,:)];
%                     y = [spiketrains_dHPC.pyr(is.baseline.sws,:)];
%                     S1= x'*y/size(x,1);
% %                     [S1 , p]=corr(x,y);
% %                     S1 = S1 - diag(diag(S1));
%                     clear x y p
%                     
%                     x = [spiketrains_dHPC.pyr(is.aversive.run,:)];
%                     y = [spiketrains_dHPC.pyr(is.aversive.run,:)];
%                     S2= x'*y/size(x,1);
% %                     [S2 , p]=corr(x,y);
% %                     S2 = S2 - diag(diag(S2));
%                     clear x y p
%                     
%                     x = [spiketrains_dHPC.pyr(is.aversive.sws,:)];
%                     y = [spiketrains_dHPC.pyr(is.aversive.sws,:)];
%                     S3= x'*y/size(x,1);
% %                     [S3 , p]=corr(x,y);
% %                     S3 = S3 - diag(diag(S3));
%                     clear x y p
%                     
%                     % EV and REV
%                     Sx = corrcoef(S2,S3,'rows','complete');     Sx = Sx(1,2);
%                     Sy = corrcoef(S2,S1,'rows','complete');     Sy = Sy(1,2);
%                     Sz = corrcoef(S3,S1,'rows','complete');     Sz = Sz(1,2);
%                     
%                     ev = (Sx-Sy*Sz);
%                     ev = (ev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
%                     
%                     rev = (Sy-Sx*Sz);
%                     rev = (rev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
%                     
%                     EV.aversive.dHPC = [EV.aversive.dHPC ; rev , ev];
%                     clear Sx Sy Sz ev rev
%                     
%                     % Reward
%                     x = [spiketrains_dHPC.pyr(is.reward.run,:)];
%                     y = [spiketrains_dHPC.pyr(is.reward.run,:)];
%                     S4= x'*y/size(x,1);
% %                     [S4 , p]=corr(x,y);
% %                     S4 = S4 - diag(diag(S4));
%                     clear x y p
%                     
%                     x = [spiketrains_dHPC.pyr(is.reward.sws,:)];
%                     y = [spiketrains_dHPC.pyr(is.reward.sws,:)];
%                     S5= x'*y/size(x,1);
% %                     [S5 , p]=corr(x,y);
% %                     S5 = S5 - diag(diag(S5));
%                     clear x y p
%                     
%                     % EV and REV
%                     Sx = corrcoef(S4,S5,'rows','complete');     Sx = Sx(1,2);
%                     Sy = corrcoef(S4,S3,'rows','complete');     Sy = Sy(1,2);
%                     Sz = corrcoef(S5,S3,'rows','complete');     Sz = Sz(1,2);
%                     
%                     ev = (Sx-Sy*Sz);
%                     ev = (ev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
%                     
%                     rev = (Sy-Sx*Sz);
%                     rev = (rev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
%                     
%                     EV.reward.dHPC = [EV.reward.dHPC ; rev , ev];
%                     clear Sx Sy Sz ev rev
% 
%                     % Aversive in Reward
%                     Sx = corrcoef(S2,S5,'rows','complete');     Sx = Sx(1,2);
%                     Sy = corrcoef(S2,S1,'rows','complete');     Sy = Sy(1,2);
%                     Sz = corrcoef(S5,S1,'rows','complete');     Sz = Sz(1,2);
%                     
%                     ev = (Sx-Sy*Sz);
%                     ev = (ev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
%                     
%                     rev = (Sy-Sx*Sz);
%                     rev = (rev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
%                     
%                     EV.aversive.dHPCi = [EV.aversive.dHPCi ; rev , ev];                    
%                     clear Sx Sy Sz ev rev S1 S2 S3 S4 S5
%                     
%                     % vHPC
%                     % Aversive
%                     x = [spiketrains_vHPC.pyr(is.baseline.sws,:)];
%                     y = [spiketrains_vHPC.pyr(is.baseline.sws,:)];
%                     S1= x'*y/size(x,1);
% %                     [S1 , p]=corr(x,y);
% %                     S1 = S1 - diag(diag(S1));
%                     clear x y p
%                     
%                     x = [spiketrains_vHPC.pyr(is.aversive.run,:)];
%                     y = [spiketrains_vHPC.pyr(is.aversive.run,:)];
%                     S2= x'*y/size(x,1);
% %                     [S2 , p]=corr(x,y);
% %                     S2 = S2 - diag(diag(S2));
%                     clear x y p
%                     
%                     x = [spiketrains_vHPC.pyr(is.aversive.sws,:)];
%                     y = [spiketrains_vHPC.pyr(is.aversive.sws,:)];
%                     S3= x'*y/size(x,1);
% %                     [S3 , p]=corr(x,y);
% %                     S3 = S3 - diag(diag(S3));
%                     clear x y p
%                     
%                     % EV and REV
%                     Sx = corrcoef(S2,S3,'rows','complete');     Sx = Sx(1,2);
%                     Sy = corrcoef(S2,S1,'rows','complete');     Sy = Sy(1,2);
%                     Sz = corrcoef(S3,S1,'rows','complete');     Sz = Sz(1,2);
%                     
%                     ev = (Sx-Sy*Sz);
%                     ev = (ev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
%                     
%                     rev = (Sy-Sx*Sz);
%                     rev = (rev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
%                     
%                     EV.aversive.vHPC = [EV.aversive.vHPC ; rev , ev];
%                     clear Sx Sy Sz ev rev
%                     
%                     % Reward
%                     x = [spiketrains_vHPC.pyr(is.reward.run,:)];
%                     y = [spiketrains_vHPC.pyr(is.reward.run,:)];
%                     S4= x'*y/size(x,1);
% %                     [S4 , p]=corr(x,y);
% %                     S4 = S4 - diag(diag(S4));
%                     clear x y p
%                     
%                     x = [spiketrains_vHPC.pyr(is.reward.sws,:)];
%                     y = [spiketrains_vHPC.pyr(is.reward.sws,:)];
%                     S5= x'*y/size(x,1);
% %                     [S5 , p]=corr(x,y);
% %                     S5 = S5 - diag(diag(S5));
%                     clear x y p
%                     
%                     % EV and REV
%                     Sx = corrcoef(S4,S5,'rows','complete');     Sx = Sx(1,2);
%                     Sy = corrcoef(S4,S3,'rows','complete');     Sy = Sy(1,2);
%                     Sz = corrcoef(S5,S3,'rows','complete');     Sz = Sz(1,2);
%                     
%                     ev = (Sx-Sy*Sz);
%                     ev = (ev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
%                     
%                     rev = (Sy-Sx*Sz);
%                     rev = (rev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
%                     
%                     EV.reward.vHPC = [EV.reward.vHPC ; rev , ev];
%                     clear Sx Sy Sz ev rev
% 
%                     % Aversive in Reward
%                     Sx = corrcoef(S2,S5,'rows','complete');     Sx = Sx(1,2);
%                     Sy = corrcoef(S2,S1,'rows','complete');     Sy = Sy(1,2);
%                     Sz = corrcoef(S5,S1,'rows','complete');     Sz = Sz(1,2);
%                     
%                     ev = (Sx-Sy*Sz);
%                     ev = (ev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
%                     
%                     rev = (Sy-Sx*Sz);
%                     rev = (rev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
%                     
%                     EV.aversive.vHPCi = [EV.aversive.vHPCi ; rev , ev];
%                     clear Sx Sy Sz ev rev S1 S2 S3 S4 S5

                    % dHPC and vHPC
                    % Aversive
                    x = [spiketrains_dHPC.pyr(is.baseline.sws,:)];
                    y = [spiketrains_vHPC.pyr(is.baseline.sws,:)];
                    S1= x'*y/size(x,1);
%                     [S1 , p]=corr(x,y);
                    clear x y p
                    
                    x = [spiketrains_dHPC.pyr(is.aversive.run,:)];
                    y = [spiketrains_vHPC.pyr(is.aversive.run,:)];
                    S2= x'*y/size(x,1);
%                     [S2 , p]=corr(x,y);
                    clear x y p
                    
                    x = [spiketrains_dHPC.pyr(is.aversive.sws,:)];
                    y = [spiketrains_vHPC.pyr(is.aversive.sws,:)];
                    S3= x'*y/size(x,1);
%                     [S3 , p]=corr(x,y);
                    clear x y p
                    
                    % EV and REV
                    Sx = corrcoef(S2,S3,'rows','complete');     Sx = Sx(1,2);
                    Sy = corrcoef(S2,S1,'rows','complete');     Sy = Sy(1,2);
                    Sz = corrcoef(S3,S1,'rows','complete');     Sz = Sz(1,2);
                    
                    ev = (Sx-Sy*Sz);
                    ev = (ev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
                    
                    rev = (Sy-Sx*Sz);
                    rev = (rev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
                    
                    EV.aversive.dvHPC = [EV.aversive.dvHPC ; rev , ev];
                    clear Sx Sy Sz ev rev
                    
                    % Reward
                    x = [spiketrains_dHPC.pyr(is.reward.run,:)];
                    y = [spiketrains_vHPC.pyr(is.reward.run,:)];
                    S4= x'*y/size(x,1);
%                     [S4 , p]=corr(x,y);
                    clear x y p
                    
                    x = [spiketrains_dHPC.pyr(is.reward.sws,:)];
                    y = [spiketrains_vHPC.pyr(is.reward.sws,:)];
                    S5= x'*y/size(x,1);
%                     [S5 , p]=corr(x,y);
                    clear x y p
                    
                    % EV and REV
                    Sx = corrcoef(S4,S5,'rows','complete');     Sx = Sx(1,2);
                    Sy = corrcoef(S4,S3,'rows','complete');     Sy = Sy(1,2);
                    Sz = corrcoef(S5,S3,'rows','complete');     Sz = Sz(1,2);
                    
                    ev = (Sx-Sy*Sz);
                    ev = (ev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
                    
                    rev = (Sy-Sx*Sz);
                    rev = (rev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
                    
%                     EV.reward.dvHPC = [EV.reward.dvHPC ; rev , ev];
                    clear Sx Sy Sz ev rev

                    % Aversive in Reward
                    Sx = corrcoef(S2,S5,'rows','complete');     Sx = Sx(1,2);
                    Sy = corrcoef(S2,S1,'rows','complete');     Sy = Sy(1,2);
                    Sz = corrcoef(S5,S1,'rows','complete');     Sz = Sz(1,2);
                    
                    ev = (Sx-Sy*Sz);
                    ev = (ev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
                    
                    rev = (Sy-Sx*Sz);
                    rev = (rev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
                    
                    EV.aversive.dvHPCi = [EV.aversive.dvHPCi ; rev , ev];
                    
                    clear Sx Sy Sz ev rev S1 S2 S3 S4 S5                    
                end
            end
            clear spiketrains_dHPC_int spiketrains_dHPC_pyr spiketrains_vHPC_int spiketrains_vHPC_pyr
        end
        disp(['-- Analysis of folder #' , num2str(t) , ' from rat #',num2str(tt) , ' finished --'])
        disp(' ')
    end
    
end

% figure
% subplot(3,4,1),boxplot([EV.reward.dHPC(:,1) , EV.reward.dHPC(:,2)]*100) , ylim([0 100]), [h p] = ranksum(EV.reward.dHPC(:,1) , EV.reward.dHPC(:,2))
% subplot(3,4,2),boxplot([EV.aversive.dHPC(:,1) , EV.aversive.dHPC(:,2)]*100) , ylim([0 100]), [h p] = ranksum(EV.aversive.dHPC(:,1) , EV.aversive.dHPC(:,2))
% subplot(3,4,3),boxplot([EV.reward.dHPCi(:,1) , EV.reward.dHPCi(:,2)]*100) , ylim([0 100]), [h p] = ranksum(EV.reward.dHPCi(:,1) , EV.reward.dHPCi(:,2))
% subplot(3,4,4),boxplot([EV.aversive.dHPCi(:,1) , EV.aversive.dHPCi(:,2)]*100) , ylim([0 100]), [h p] = ranksum(EV.aversive.dHPCi(:,1) , EV.aversive.dHPCi(:,2))

% subplot(3,4,5),boxplot([EV.reward.vHPC(:,1) , EV.reward.vHPC(:,2)]*100) , ylim([0 100]), [h p] =ranksum(EV.reward.vHPC(:,1) , EV.reward.vHPC(:,2))
% subplot(3,4,6),boxplot([EV.aversive.vHPC(:,1) , EV.aversive.vHPC(:,2)]*100) , ylim([0 100]), [h p] =ranksum(EV.aversive.vHPC(:,1) , EV.aversive.vHPC(:,2))
% subplot(3,4,7),boxplot([EV.reward.vHPCi(:,1) , EV.reward.vHPCi(:,2)]*100) , ylim([0 100]), [h p] =ranksum(EV.reward.vHPCi(:,1) , EV.reward.vHPCi(:,2))
% subplot(3,4,8),boxplot([EV.aversive.vHPCi(:,1) , EV.aversive.vHPCi(:,2)]*100) , ylim([0 100]), [h p] =ranksum(EV.aversive.vHPCi(:,1) , EV.aversive.vHPCi(:,2))

% subplot(3,4,9),boxplot([EV.reward.dvHPC(:,1) , EV.reward.dvHPC(:,2)]*100) , ylim([0 10]), [h p] =ranksum(EV.reward.dvHPC(:,1) , EV.reward.dvHPC(:,2))
% subplot(3,4,10),boxplot([EV.aversive.dvHPC(:,1) , EV.aversive.dvHPC(:,2)]*100) , ylim([0 10]), [h p] =ranksum(EV.aversive.dvHPC(:,1) , EV.aversive.dvHPC(:,2))
% subplot(3,4,11),boxplot([EV.reward.dvHPCi(:,1) , EV.reward.dvHPCi(:,2)]*100) , ylim([0 10]), [h p] =ranksum(EV.reward.dvHPCi(:,1) , EV.reward.dvHPCi(:,2))
% subplot(3,4,12),boxplot([EV.aversive.dvHPCi(:,1) , EV.aversive.dvHPCi(:,2)]*100) , ylim([0 10]), [h p] =ranksum(EV.aversive.dvHPCi(:,1) , EV.aversive.dvHPCi(:,2))


figure
subplot(1,4,1),boxplot([EV.reward.dvHPC(:,1) , EV.reward.dvHPC(:,2)]*100) , ylim([0 6]), [h p] =ranksum(EV.reward.dvHPC(:,1) , EV.reward.dvHPC(:,2))
subplot(1,4,2),boxplot([EV.aversive.dvHPC(:,1) , EV.aversive.dvHPC(:,2)]*100) , ylim([0 6]), [h p] =ranksum(EV.aversive.dvHPC(:,1) , EV.aversive.dvHPC(:,2))
subplot(1,4,3),boxplot([EV.reward.dvHPCi(:,1) , EV.reward.dvHPCi(:,2)]*100) , ylim([0 6]), [h p] =ranksum(EV.reward.dvHPCi(:,1) , EV.reward.dvHPCi(:,2))
subplot(1,4,4),boxplot([EV.aversive.dvHPCi(:,1) , EV.aversive.dvHPCi(:,2)]*100) , ylim([0 6]), [h p] =ranksum(EV.aversive.dvHPCi(:,1) , EV.aversive.dvHPCi(:,2))
