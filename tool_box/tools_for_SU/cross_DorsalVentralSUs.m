function [cross time] = cross_DorsalVentralSUs(path)
% This function a CCG between dorsal and ventral SUs. It does it during RUN
% periods, and if they are coordinated (>95th percentile of a surrogate)
% they store their CCG during PRE and POST sleep.
%
% INPUTS
% path: cell, inside each cell you should put the path of each session you
%       want to analyze
%
% OUTPUT
% cross: structure, it contains the CCG for aversive and reward
%                .pre.aversive
%                    .reward
%
%                .post.aversive
%                     .reward
%
%
%                .run.reward
%                .run.aversive.up
%                             .down
%
% time: column vector, it contains the time axis for the CCG plot.
%
% other functions: CCG from FMA toolbox
% Morci Juan Facundo 01/2024

% variables for behavrio
minimal_speed = 7;
minimal_speed_time = 2;

% variables for CCG
sm = 2;
d = 0.5;
b = 0.005;

% Initialization of structures that wil contain the outputs
cross.pre.aversive = []; cross.post.aversive = []; cross.run.aversive.down = [];  cross.run.aversive.up = [];
cross.pre.reward = []; cross.post.reward = []; cross.run.reward.down = [];  cross.run.reward.up = [];

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
        
        if aversiveTS_run(1) < rewardTS_run(1)
            %             PreA = Restrict(NREM.baseline, [aversiveTS_run(1,1)/1000 - 3600 aversiveTS_run(1,2)/1000]);
            %             PostA = Restrict(NREM.aversive, [aversiveTS_run(1,2)/1000 aversiveTS_run(1,2)/1000 + 3600]);
            %             PreR = Restrict(NREM.aversive, [rewardTS_run(1,1)/1000 - 3600 rewardTS_run(1,2)/1000]);
            %             PostR =Restrict(NREM.reward, [rewardTS_run(1,2)/1000 rewardTS_run(1,2)/1000 + 3600]);
            PreA = NREM.baseline;
            PostA = NREM.aversive;
            PreR = NREM.aversive;
            PostR = NREM.reward;
        else
            %             PreA = Restrict(NREM.reward, [aversiveTS_run(1,1)/1000 - 3600 aversiveTS_run(1,2)/1000]);
            %             PostA = Restrict(NREM.aversive, [aversiveTS_run(1,1)/1000 - 3600 aversiveTS_run(1,2)/1000]);
            %             PreR = Restrict(NREM.baseline, [aversiveTS_run(1,2)/1000 aversiveTS_run(1,2)/1000 + 3600]);
            %             PostR = Restrict(NREM.reward, [rewardTS_run(1,2)/1000 rewardTS_run(1,2)/1000 + 3600]);
            PreA = NREM.reward;
            PostA = NREM.aversive;
            PreR = NREM.baseline;
            PostR = NREM.reward;
        end
        
        
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
        cellulartype = [K(:,1) , K(:,4)];
        
        
        %% Counting the Number f SU
        numberD = 0;
        clusters.all = [];
        clusters.dHPC = [];
        clusters.int.dHPC = [];
        for ii=1:size(group_dHPC,1)
            cluster = group_dHPC(ii,1);
            Cellulartype = logical(cellulartype(cellulartype(:,1) == cluster,2));
            if Cellulartype
                a = length(Restrict(spks(spks(:,1)==cluster,2),aversiveTS_run./1000)) / ((aversiveTS_run(2)-aversiveTS_run(1))/1000);
                r = length(Restrict(spks(spks(:,1)==cluster,2),rewardTS_run./1000)) / ((rewardTS_run(2)-rewardTS_run(1))/1000);
                if or(a > 0 ,  r > 0)
                    numberD = numberD+1;
                    clusters.all = [clusters.all ; cluster];
                    clusters.dHPC = [clusters.dHPC ; cluster];
                end
            else
                clusters.int.dHPC = [clusters.int.dHPC ; cluster];
            end
            clear tmp cluster Cellulartype fr1 fr2 fr3 fr4 fr5 r a
        end
        
        numberV = 0;
        clusters.vHPC = [];
        clusters.int.vHPC = [];
        for ii=1:size(group_vHPC,1)
            cluster = group_vHPC(ii,1);
            Cellulartype = logical(cellulartype(cellulartype(:,1) == cluster,2));
            if Cellulartype
                a = length(Restrict(spks(spks(:,1)==cluster,2),aversiveTS_run./1000)) / ((aversiveTS_run(2)-aversiveTS_run(1))/1000);
                r = length(Restrict(spks(spks(:,1)==cluster,2),rewardTS_run./1000)) / ((rewardTS_run(2)-rewardTS_run(1))/1000);
                if or(a > 0 ,  r > 0)
                    numberV = numberV+1;
                    clusters.all = [clusters.all ; cluster];
                    clusters.vHPC = [clusters.vHPC ; cluster];
                end
            else
                clusters.int.vHPC = [clusters.int.vHPC ; cluster];
            end
            clear tmp cluster Cellulartype fr1 fr2 fr3 fr4 fr5 r a
        end
        clear freq limits
        clear camara shock rightvalve leftvalve
        clear ejeX ejeY dX dY dX_int dY_int
        
        if and(numberD >0 , numberV > 0)
            for i = 1 : size(clusters.dHPC,1)
                SpikesD = spks(spks(:,1)==clusters.dHPC(i),2);
                for ii = 1 : size(clusters.vHPC,1)
                    SpikesV = spks(spks(:,1)==clusters.vHPC(ii),2);
                    
                    % Aversive
                    % Run Aversive
                    Times1 = Restrict(SpikesD,movement.aversive);
                    Times2 = Restrict(SpikesV,movement.aversive);
                    if and(size(Times1,1)>5 , size(Times2,1)>5)
                        [s,ids,groups] = CCGParameters(Times1,ones(length(Times1),1),Times2,ones(length(Times2),1)*2);
                        [ccg,T] = CCG(s,ids,'binSize',b,'duration',d,'smooth',sm,'mode','ccg');
                        ccgR = ccg(:,1,2);
                        
                        surrogated = [];
                        [o beg] = min(abs(T - (-0.025)));
                        [o sto] = min(abs(T - (0.025)));
                        for iii = 1 : 100
                            Times1 = Restrict(SpikesD,movement.aversive);
                            Times2 = Restrict(ShuffleSpks(SpikesV),movement.aversive);
                            if and(size(Times1,1)>5 , size(Times2,1)>5)
                                [s,ids,groups] = CCGParameters(Times1,ones(length(Times1),1),Times2,ones(length(Times2),1)*2);
                                [ccg,T] = CCG(s,ids,'binSize',b,'duration',d,'smooth',sm,'mode','ccg');
                                ccg = (ccg(:,1,2));
                                surrogated = [surrogated ; mean(ccg(beg:sto))];
                            end
                        end
                        
                        if nanmean(ccgR(beg:sto)) > prctile(surrogated,99)
                            cross.run.aversive.up = [cross.run.aversive.up , ccgR];
                            % Pre
                            Times1 = Restrict(SpikesD,PreA);
                            Times2 = Restrict(SpikesV,PreA);
                            if and(size(Times1,1)>5 , size(Times2,1)>5)
                                [s,ids,groups] = CCGParameters(Times1,ones(length(Times1),1),Times2,ones(length(Times2),1)*2);
                                [ccg,T] = CCG(s,ids,'binSize',b,'duration',d,'smooth',sm,'mode','ccg');
                                %                         ccg = (ccg(:,1,2)./sum(ccg(:,1,2)));
                                ccg = (ccg(:,1,2));
                                cross.pre.aversive = [cross.pre.aversive , ccg]; clear ccg s ids groups Times1 Times2
                            end
                            
                            %Post
                            Times1 = Restrict(SpikesD,PostA);
                            Times2 = Restrict(SpikesV,PostA);
                            if and(size(Times1,1)>5 , size(Times2,1)>5)
                                [s,ids,groups] = CCGParameters(Times1,ones(length(Times1),1),Times2,ones(length(Times2),1)*2);
                                [ccg,T] = CCG(s,ids,'binSize',b,'duration',d,'smooth',sm,'mode','ccg');
                                %                         ccg = (ccg(:,1,2)./sum(ccg(:,1,2)));
                                ccg = (ccg(:,1,2));
                                cross.post.aversive = [cross.post.aversive , ccg]; clear ccg s ids groups Times1 Times2
                            end
                        else
                            cross.run.aversive.down = [cross.run.aversive.down , ccgR];
                        end
                        clear ccgR
                    end
                    
                    %Reward
                    % Run Reward
                    Times1 = Restrict(SpikesD,movement.reward);
                    Times2 = Restrict(SpikesV,movement.reward);
                    if and(size(Times1,1)>5 , size(Times2,1)>5)
                        [s,ids,groups] = CCGParameters(Times1,ones(length(Times1),1),Times2,ones(length(Times2),1)*2);
                        [ccg,T] = CCG(s,ids,'binSize',b,'duration',d,'smooth',sm,'mode','ccg');
                        %                         ccg = (ccg(:,1,2)./sum(ccg(:,1,2)));
                        ccgR = (ccg(:,1,2));
                        
                        surrogated = [];
                        [o beg] = min(abs(T - (-0.025)));
                        [o sto] = min(abs(T - (0.025)));
                        for iii = 1 : 100
                            Times1 = Restrict(SpikesD,movement.reward);
                            Times2 = Restrict(ShuffleSpks(SpikesV),movement.reward);
                            if and(size(Times1,1)>5 , size(Times2,1)>5)
                                [s,ids,groups] = CCGParameters(Times1,ones(length(Times1),1),Times2,ones(length(Times2),1)*2);
                                [ccg,T] = CCG(s,ids,'binSize',b,'duration',d,'smooth',sm,'mode','ccg');
                                ccg = (ccg(:,1,2));
                                surrogated = [surrogated ; mean(ccg(beg:sto))];
                            end
                        end
                        
                        if nanmean(ccgR(beg:sto)) > prctile(surrogated,99)
                            cross.run.reward.up = [cross.run.reward.up , ccgR];
                            % Pre
                            Times1 = Restrict(SpikesD,PreR);
                            Times2 = Restrict(SpikesV,PreR);
                            if and(size(Times1,1)>5 , size(Times2,1)>5)
                                [s,ids,groups] = CCGParameters(Times1,ones(length(Times1),1),Times2,ones(length(Times2),1)*2);
                                [ccg,T] = CCG(s,ids,'binSize',b,'duration',d,'smooth',sm,'mode','ccg');
                                %                         ccg = (ccg(:,1,2)./sum(ccg(:,1,2)));
                                ccg = (ccg(:,1,2));
                                cross.pre.reward = [cross.pre.reward , ccg]; clear ccg s ids groups Times1 Times2
                            end
                            
                            %Post
                            Times1 = Restrict(SpikesD,PostR);
                            Times2 = Restrict(SpikesV,PostR);
                            if and(size(Times1,1)>5 , size(Times2,1)>5)
                                [s,ids,groups] = CCGParameters(Times1,ones(length(Times1),1),Times2,ones(length(Times2),1)*2);
                                [ccg,T] = CCG(s,ids,'binSize',b,'duration',d,'smooth',sm,'mode','ccg');
                                %                         ccg = (ccg(:,1,2)./sum(ccg(:,1,2)));
                                ccg = (ccg(:,1,2));
                                cross.post.reward = [cross.post.reward , ccg]; clear ccg s ids groups Times1 Times2
                            end
                        else
                            cross.run.reward.down = [cross.run.reward.down , ccgR];
                        end
                    end
                    clear SpikesV ccgR
                end
                clear SpikesD
            end
        end
        
        clear aversiveTS aversiveTS_run baselineTS bins Cell_type_classification
        clear cellulartype clusters Events group_dHPC group_vHPC i ii iii K
        clear Kinfo NREM REM numberD numberV p parameter patterns post pre
        clear SpikesD SpikesV spks spks_dHPC spks_vHPC Thresholded WAKE
        clear Kinfo K ii iii ripples ripplesD ripplesV ripple_event Spikes
        clear x y loc1 loc2 pks1 pks2 cond Is rewardTS rewardTS_run aversiveTS aversiveTS_run
        clear segments c PreA PostA PreR PostR
    end
end
time = T;

 
 figure,
 subplot(131),
 plot(time,nanmean((cross.pre.aversive./sum(cross.pre.aversive)),2),'k'),hold on
 plot(time,nanmean((cross.post.aversive./sum(cross.post.aversive)),2),'r'),
 
 subplot(132),
 plot(time,nanmean(cross.run.aversive.down./sum(cross.run.aversive.down),2),'k'),hold on
 plot(time,nanmean(cross.run.aversive.up./sum(cross.run.aversive.up),2),'r'),
 
 subplot(133),
 plot(time,nanmean(cross.pre.reward./sum(cross.pre.reward),2),'k'),hold on
 plot(time,nanmean(cross.post.reward./sum(cross.post.reward),2),'b'),
 
end