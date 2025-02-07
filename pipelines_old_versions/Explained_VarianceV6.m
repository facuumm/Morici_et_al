clear
clc
% close all

%% Parameters
path = {'E:\Rat103\usable';'E:\Rat126\Ephys\in_Pyr';'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat132\recordings\in_pyr';'E:\Rat165\in_pyr\'};%List of folders from the path

%Sleep
time_criteria = 1; % minimal time to include a NREM epoch (in min)

% for SU
criteria_fr = 0; %criteria to include or not a SU into the analysis
criteria_n = [6 6]; % minimal number of neurons from each structure [vHPC dHPC]
criteria_type = 0; %criteria for celltype (0:pyr, 1:int, 2:all)
binSize = 0.05;
n_SU_V = 0;
n_SU_D = 0;

% Behavior
minimal_speed = 7; % minimal speed to detect quite periods
minimal_speed_time = 2; % minimal time to detect quite periods

EV.aversive.dvHPCi = [];   EV.reward.dvHPCi = [];
EV.aversive.dvHPC = [];   EV.reward.dvHPC = [];

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
%         movement.reward(movement.reward(:,2) - movement.reward(:,1) <1,:)=[]; %eliminate 1sec segments
        clear tmp start stop
        % Aversive
        start = behavior.speed.aversive(1,1);   stop = behavior.speed.aversive(end,1);
        movement.aversive = InvertIntervals(behavior.quiet.aversive , start , stop);%keep only those higher than criteria
%         movement.aversive(movement.aversive(:,2) - movement.aversive(:,1) <1,:)=[];
        clear tmp start stop
        
        %% load sleep states
        disp('Uploading sleep scoring')
        x = dir([cd,'\*-states.mat']);    states = load([cd,'\',x.name]);    states = states.states;
        REM.all = ToIntervals(states==5);    NREM.all = ToIntervals(states==3);    WAKE.all = ToIntervals(states==1);
        clear x states
        
%         NREM.all(NREM.all(:,2)-NREM.all(:,1)<60*time_criteria,:)=[];
        NREM.baseline = Restrict(NREM.all,baselineTS./1000);
        NREM.aversive = Restrict(NREM.all,aversiveTS./1000);
        NREM.reward = Restrict(NREM.all,rewardTS./1000);    
        
        REM.baseline = Restrict(REM.all,baselineTS./1000);
        REM.aversive = Restrict(REM.all,aversiveTS./1000);
        REM.reward = Restrict(REM.all,rewardTS./1000);
                
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
            cellulartype = [K(:,1) , K(:,4)];
        elseif criteria_type == 1 % int
            cellulartype = [K(:,1) , not(K(:,4))];
        elseif criteria_type == 2 % all
            cellulartype = [K(:,1) , ones(size(K,1),1)];
        end
        % Definition of the emotional transition
        if aversiveTS_run(1)<rewardTS_run(1)
            config = 1; % 1 if the order was Aversive -> Reward
        else
            config = 2;% 1 if the order was Reward -> Aversive
        end
        
        %% Constructing Spiketrains
        freq = 1/binSize;
        limits = [0 segments.Var1(end)/1000];
        spiketrains_dHPC.pyr = [];        spiketrains_dHPC.int = [];
        spiketrains_vHPC.pyr = [];        spiketrains_vHPC.int = [];
        for ii=1:size(group_dHPC,1)
            cluster = group_dHPC(ii,1);
            Cellulartype = logical(Cell_type_classification(Cell_type_classification(:,1) == cluster,7));
            if Cellulartype
                spks = spks_dHPC(spks_dHPC(:,1)==cluster,2);
                [tmp,bins]=binspikes(spks,freq,limits);
%                 is.sleep.baseline = InIntervals(bins,baselineTS./1000);                 is.sleep.reward = InIntervals(bins,rewardTS./1000);
%                 is.sleep.aversive = InIntervals(bins,aversiveTS./1000);                 is.run.reward = InIntervals(bins,rewardTS_run./1000);
%                 is.run.aversive = InIntervals(bins,aversiveTS_run./1000);
%                 fr1 = sum(tmp(is.sleep.baseline))/(baselineTS(2)/1000);
%                 fr2 = sum(tmp(is.sleep.reward))/((rewardTS(2)-rewardTS(1))/1000);
%                 fr3 = sum(tmp(is.sleep.aversive))/((aversiveTS(2)-aversiveTS(1))/1000);
%                 fr4 = sum(tmp(is.run.reward))/((rewardTS_run(2)-rewardTS_run(1))/1000);
%                 fr5 = sum(tmp(is.run.aversive))/((aversiveTS_run(2)-aversiveTS_run(1))/1000);
%                 
%                 FiringRate.dHPC.baseline = [FiringRate.dHPC.baseline ; fr1];
%                 FiringRate.dHPC.reward = [FiringRate.dHPC.reward ; fr2];
%                 FiringRate.dHPC.aversive = [FiringRate.dHPC.aversive ; fr3];
%                 
%                 if and(and(fr1>criteria_fr , fr2 > criteria_fr),fr3>criteria_fr)
                    spiketrains_dHPC.pyr = [spiketrains_dHPC.pyr , (tmp)];
%                 end
            end
            clear spks tmp cluster Cellulartype fr1 fr2 fr3 fr4 fr5
        end
        
        for ii=1:size(group_vHPC,1)
            cluster = group_vHPC(ii,1);
            Cellulartype = logical(Cell_type_classification(Cell_type_classification(:,1) == cluster,7));
            if Cellulartype
                spks = spks_vHPC(spks_vHPC(:,1)==cluster,2);
                [tmp,bins]=binspikes(spks,freq,limits);
%                 is.sleep.baseline = InIntervals(bins,baselineTS./1000);                 is.sleep.reward = InIntervals(bins,rewardTS./1000);
%                 is.sleep.aversive = InIntervals(bins,aversiveTS./1000);                 is.run.reward = InIntervals(bins,rewardTS_run./1000);
%                 is.run.aversive = InIntervals(bins,aversiveTS_run./1000);
%                 fr1 = sum(tmp(is.sleep.baseline))/(baselineTS(2)/1000);
%                 fr2 = sum(tmp(is.sleep.reward))/((rewardTS(2)-rewardTS(1))/1000);
%                 fr3 = sum(tmp(is.sleep.aversive))/((aversiveTS(2)-aversiveTS(1))/1000);
%                 fr4 = sum(tmp(is.run.reward))/((rewardTS_run(2)-rewardTS_run(1))/1000);
%                 fr5 = sum(tmp(is.run.aversive))/((aversiveTS_run(2)-aversiveTS_run(1))/1000);
%                 
%                 FiringRate.vHPC.baseline = [FiringRate.vHPC.baseline ; fr1];
%                 FiringRate.vHPC.reward = [FiringRate.vHPC.reward ; fr2];
%                 FiringRate.vHPC.aversive = [FiringRate.vHPC.aversive ; fr3];                
%                 
%                 if and(and(fr1>criteria_fr , fr2 > criteria_fr),fr3>criteria_fr)
                    spiketrains_vHPC.pyr = [spiketrains_vHPC.pyr , (tmp)];
%                 end
            end
            clear spks tmp cluster Cellulartype fr1 fr2 fr3 fr4 fr5
        end
        clear freq limits
        clear spks spks_dHPC spks_vHPC camara shock rightvalve leftvalve
        clear ejeX ejeY dX dY dX_int dY_int
        
        %% Explained variance calculation
        if and(size(spiketrains_vHPC.pyr,2) >= criteria_n(1),size(spiketrains_dHPC.pyr,2) >= criteria_n(2))
            disp('Lets go for the SUs')
            %Restricting bins inside each condition
            is.baseline.sws = InIntervals(bins,NREM.baseline);
            is.reward.sws = InIntervals(bins,NREM.reward);
            is.aversive.sws = InIntervals(bins,NREM.aversive);
            
%             is.baseline.sws = InIntervals(bins,REM.baseline);
%             is.aversive.sws = InIntervals(bins,REM.aversive);
%             is.reward.sws = InIntervals(bins,REM.reward);
            
            is.aversive.run = InIntervals(bins,movement.aversive);
            is.reward.run = InIntervals(bins,movement.reward);
            
%             is.aversive.run = InIntervals(bins,aversiveTS_run./1000);
%             is.reward.run = InIntervals(bins,rewardTS_run./1000);

            %%
            if and(and(~isempty(NREM.aversive),~isempty(NREM.reward)),~isempty(NREM.baseline))
                if aversiveTS_run(1)>rewardTS_run(1)
                    %% Reward
                    % Correlation Matrix Calculation
                    x = [spiketrains_dHPC.pyr(is.baseline.sws,:)];
                    y = [spiketrains_vHPC.pyr(is.baseline.sws,:)];      
                    [S1 , p] = corr(x,y);
%                     S1 = S1 - diag(diag(S1));
                    clear x y p
                    
                    x = [spiketrains_dHPC.pyr(is.reward.run,:)];
                    y = [spiketrains_vHPC.pyr(is.reward.run,:)];
                    [S2 , p] = corr(x,y);
%                     S2 = S2 - diag(diag(S2));
                    clear x y p
                    
                    x = [spiketrains_dHPC.pyr(is.reward.sws,:)];
                    y = [spiketrains_vHPC.pyr(is.reward.sws,:)];
                    [S3 , p] = corr(x,y);
%                     S3 = S3 - diag(diag(S3));
                    clear x y p
                    
                    % EV and REV
                    Sx = corrcoef(S2,S3,'rows','complete');     Sx = Sx(1,2);
                    Sy = corrcoef(S2,S1,'rows','complete');     Sy = Sy(1,2);
                    Sz = corrcoef(S3,S1,'rows','complete');     Sz = Sz(1,2);
                    
                    ev = (Sx-Sy*Sz);
                    ev = (ev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
                    
                    rev = (Sy-Sx*Sz);
                    rev = (rev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
                    
                    EV.reward.dvHPC = [EV.reward.dvHPC ; rev*100 , ev*100 tt t config];
                    clear Sx Sy Sz ev rev
                                        
                    %% Aversive
                    x = [spiketrains_dHPC.pyr(is.reward.sws,:)];
                    y = [spiketrains_vHPC.pyr(is.reward.sws,:)];
                    [S3 , p] = corr(x,y);
%                     S3 = S3 - diag(diag(S3));
                    clear x y p
                    
                    % Correlation Matrix Calculation
                    x = [spiketrains_dHPC.pyr(is.aversive.run,:)];
                    y = [spiketrains_vHPC.pyr(is.aversive.run,:)];
                    [S4 , p] = corr(x,y);
%                     S4 = S4 - diag(diag(S4));
                    clear x y p
                    
                    x = [spiketrains_dHPC.pyr(is.aversive.sws,:)];
                    y = [spiketrains_vHPC.pyr(is.aversive.sws,:)];
                    [S5 , p] = corr(x,y);
%                     S5 = S5 - diag(diag(S5));
                    clear x y p
                    
                    % EV and REV
                    Sx = corrcoef(S4,S5,'rows','complete');     Sx = Sx(1,2);
                    Sy = corrcoef(S4,S3,'rows','complete');     Sy = Sy(1,2);
                    Sz = corrcoef(S5,S3,'rows','complete');     Sz = Sz(1,2);
                    
                    ev = (Sx-Sy*Sz);
                    ev = (ev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
                    
                    rev = (Sy-Sx*Sz);
                    rev = (rev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
                    
                    EV.aversive.dvHPC = [EV.aversive.dvHPC ; rev*100 , ev*100 tt t config];
                    clear Sx Sy Sz ev rev
                    
%                     EV.reward.aversive.dvHPC = [EV.reward.aversive.dvHPC ; rev1 ev1 rev2 ev2];
                    clear rev1 ev1 rev2 ev2
                    %% Reward in Aversive
                    Sx = corrcoef(S2,S5,'rows','complete');     Sx = Sx(1,2);
                    Sy = corrcoef(S2,S1,'rows','complete');     Sy = Sy(1,2);
                    Sz = corrcoef(S5,S1,'rows','complete');     Sz = Sz(1,2);
                    
                    ev = (Sx-Sy*Sz);
                    ev = (ev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
                    
                    rev = (Sy-Sx*Sz);
                    rev = (rev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
                    
                    EV.reward.dvHPCi = [EV.reward.dvHPCi ; rev*100 , ev*100 tt t config];
                    clear Sx Sy Sz ev rev S1 S2 S3 S4 S5                     
                    
                else
                    %% Aversive
                    % Correlation Matrix Calculation
                    x = [spiketrains_dHPC.pyr(is.baseline.sws,:)];
                    y = [spiketrains_vHPC.pyr(is.baseline.sws,:)];
                    [S1 , p] = corr(x,y);
%                     S1 = S1 - diag(diag(S1));
                    clear x y p
                    
                    x = [spiketrains_dHPC.pyr(is.aversive.run,:)];
                    y = [spiketrains_vHPC.pyr(is.aversive.run,:)];
                    [S2 , p] = corr(x,y);
%                     S2 = S2 - diag(diag(S2));
                    clear x y p
                    
                    x = [spiketrains_dHPC.pyr(is.aversive.sws,:)];
                    y = [spiketrains_vHPC.pyr(is.aversive.sws,:)];
                    [S3 , p] = corr(x,y);
%                     S3 = S3 - diag(diag(S3));
                    clear x y p
                    
                    % EV and REV
                    Sx = corrcoef(S2,S3,'rows','complete');     Sx = Sx(1,2);
                    Sy = corrcoef(S2,S1,'rows','complete');     Sy = Sy(1,2);
                    Sz = corrcoef(S3,S1,'rows','complete');     Sz = Sz(1,2);
                    
                    ev = (Sx-Sy*Sz);
                    ev = (ev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
                    
                    rev = (Sy-Sx*Sz);
                    rev = (rev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
                    
                    EV.aversive.dvHPC = [EV.aversive.dvHPC ; rev*100 , ev*100 tt t config];
                    clear Sx Sy Sz ev rev
                    
                    %% Reward
                    x = [spiketrains_dHPC.pyr(is.aversive.sws,:)];
                    y = [spiketrains_vHPC.pyr(is.aversive.sws,:)];
                    [S3 , p] = corr(x,y);
%                     S3 = S3 - diag(diag(S3));
                    clear x y p
                    
                    % Correlation Matrix Calculation      
                    x = [spiketrains_dHPC.pyr(is.reward.run,:)];
                    y = [spiketrains_vHPC.pyr(is.reward.run,:)];
                    [S4 , p] = corr(x,y);
%                     S4 = S4 - diag(diag(S4));
                    clear x y p
                    
                    x = [spiketrains_dHPC.pyr(is.reward.sws,:)];
                    y = [spiketrains_vHPC.pyr(is.reward.sws,:)];
                    [S5 , p] = corr(x,y);
%                     S5 = S5 - diag(diag(S5));
                    clear x y p
                    
                    % EV and REV
                    Sx = corrcoef(S4,S5,'rows','complete');     Sx = Sx(1,2);
                    Sy = corrcoef(S4,S3,'rows','complete');     Sy = Sy(1,2);
                    Sz = corrcoef(S5,S3,'rows','complete');     Sz = Sz(1,2);
                    
                    ev = (Sx-Sy*Sz);
                    ev = (ev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
                    
                    rev = (Sy-Sx*Sz);
                    rev = (rev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
                    
                    EV.reward.dvHPC = [EV.reward.dvHPC ; rev*100 , ev*100 tt t config];
                    clear Sx Sy Sz ev rev
                    
%                     EV.aversive.reward.dvHPC = [EV.aversive.reward.dvHPC ; rev1 ev1 rev2 ev2];
%                     clear rev1 ev1 rev2 ev2
                    
                    %% Aversive in Reward
                    Sx = corrcoef(S2,S5,'rows','complete');     Sx = Sx(1,2);
                    Sy = corrcoef(S2,S1,'rows','complete');     Sy = Sy(1,2);
                    Sz = corrcoef(S5,S1,'rows','complete');     Sz = Sz(1,2);
                    
                    ev = (Sx-Sy*Sz);
                    ev = (ev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
                    
                    rev = (Sy-Sx*Sz);
                    rev = (rev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
                    
                    EV.aversive.dvHPCi = [EV.aversive.dvHPCi ; rev*100 , ev*100 tt t config];
                    clear Sx Sy Sz ev rev S1 S2 S3 S4 S5                     

                end
            end
            
        end
        disp(['-- Analysis of folder #' , num2str(t) , ' from rat #',num2str(tt) , ' finished --'])
        disp(' ')
        clear spiketrains_dHPC_int spiketrains_dHPC_pyr spiketrains_vHPC_int spiketrains_vHPC_pyr
        clear aversiveTS aversiveTS_run rewardTS rewardTS_run behavior bins Cell_type_classification cellulartype
        clear group_dHPC group_vHPC is K Kinfo NREM REM WAKE ripplesD ripplesV segments
        clear spiketrains_dHPC spiketrains_vHPC ripple_bursts ripple_event baselineTS cond
        clear coordinated coordinatedV coordinatedV_refined movement
    end
end

% save([cd,'\Explained_Variance_REM.mat'] , 'EV')

figure
subplot(1,2,1),boxplot([EV.reward.dvHPC(:,1) , EV.reward.dvHPC(:,2)]) , ylim([0 4]), [h p] =ranksum(EV.reward.dvHPC(:,1) , EV.reward.dvHPC(:,2),'tail','left'),hold on
x = [[EV.reward.dvHPC(:,1) ; EV.reward.dvHPC(:,2)] , [ones(length(EV.reward.dvHPC),1) ; ones(length(EV.reward.dvHPC),1)*2]];
scatter(x(:,2),x(:,1),"filled",'jitter','on', 'jitterAmount',0.1)

subplot(1,2,2),boxplot([EV.aversive.dvHPC(:,1) , EV.aversive.dvHPC(:,2)]) , ylim([0 4]), [h p] =ranksum(EV.aversive.dvHPC(:,1) , EV.aversive.dvHPC(:,2),'tail','left'),hold on
x = [[EV.aversive.dvHPC(:,1) ; EV.aversive.dvHPC(:,2)] , [ones(length(EV.aversive.dvHPC),1) ; ones(length(EV.aversive.dvHPC),1)*2]];
scatter(x(:,2),x(:,1),"filled",'jitter','on', 'jitterAmount',0.1)

% subplot(1,4,3),boxplot([EV.reward.dvHPCi(:,1) , EV.reward.dvHPCi(:,2)]) , ylim([0 3]), [h p] =ranksum(EV.reward.dvHPCi(:,1) , EV.reward.dvHPCi(:,2),'tail','left')
% subplot(1,4,4),boxplot([EV.aversive.dvHPCi(:,1) , EV.aversive.dvHPCi(:,2)]) , ylim([0 3]), [h p] =ranksum(EV.aversive.dvHPCi(:,1) , EV.aversive.dvHPCi(:,2),'tail','left')

% figure
% subplot(1,2,1),scatter([ones(length(EV.reward.dvHPC(:,1)),1) ; ones(length(EV.reward.dvHPC(:,2)),1)*2],[EV.reward.dvHPC(:,1) ; EV.reward.dvHPC(:,2)],"filled",'jitter','on', 'jitterAmount',0.1),xlim([0 3]),hold on
% scatter([1 2] , [nanmean(EV.reward.dvHPC(:,1)) nanmean(EV.reward.dvHPC(:,2))],'filled')
% subplot(1,2,2),scatter([ones(length(EV.aversive.dvHPC(:,1)),1) ; ones(length(EV.aversive.dvHPC(:,2)),1)*2],[EV.aversive.dvHPC(:,1) ; EV.aversive.dvHPC(:,2)],"filled",'jitter','on', 'jitterAmount',0.1),xlim([0 3]),hold on
% scatter([1 2] , [nanmean(EV.aversive.dvHPC(:,1)) nanmean(EV.aversive.dvHPC(:,2))],'filled')


x = [EV.reward.dvHPC(:,1) ; EV.reward.dvHPC(:,2) ; EV.aversive.dvHPC(:,1) ; EV.aversive.dvHPC(:,2)];
y = [ones(length([EV.reward.dvHPC(:,1) ; EV.reward.dvHPC(:,2)]),1) ; ones(length([EV.aversive.dvHPC(:,1) ; EV.aversive.dvHPC(:,2)]),1)*2];
y = [y , [ones(length(EV.reward.dvHPC(:,1)),1) ; ones(length(EV.reward.dvHPC(:,2)),1)*2 ; ones(length(EV.aversive.dvHPC(:,1)),1) ; ones(length(EV.aversive.dvHPC(:,2)),1)*2]]

[~,~,stats] = anovan(x,y,'model','interaction','varnames',{'condition','tipo'})
[results,~,~,gnames] = multcompare(stats,"Dimension",[1 2]);

x = [EV.reward.dvHPC(:,1) ; EV.reward.dvHPC(:,2) ; EV.aversive.dvHPC(:,1) ; EV.aversive.dvHPC(:,2)];
y = {'REV';'REV';'REV';'REV';'REV';'REV';'REV';'REV';'REV';'REV';'REV';'REV';'REV';'REV';'EV';'EV';'EV';'EV';'EV';'EV';'EV';'EV';'EV';'EV';'EV';'EV';'EV';'EV';'REV';'REV';'REV';'REV';'REV';'REV';'REV';'REV';'EV';'EV';'EV';'EV';'EV';'EV';'EV';'EV'};
z = {'R';'R';'R';'R';'R';'R';'R';'R';'R';'R';'R';'R';'R';'R';'R';'R';'R';'R';'R';'R';'R';'R';'R';'R';'R';'R';'R';'R';'A';'A';'A';'A';'A';'A';'A';'A';'A';'A';'A';'A';'A';'A';'A';'A'};
[~,~,stats] = anovan(x,{y,z},'interaction')
[results,~,~,gnames] = multcompare(stats,"Dimension",[1 2]);


x = EV.aversive.reward.dvHPC .* 100;
y = EV.reward.aversive.dvHPC .* 100;

figure
subplot(121),boxplot(x),ylim([0 8]), ranksum(x(:,3),x(:,4),'tail','left')
subplot(122),boxplot(y),ylim([0 8]), ranksum(y(:,3),y(:,4),'tail','left')

