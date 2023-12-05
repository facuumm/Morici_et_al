clear
clc
% close all

%% Parameters
path = {'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat103\usable';'E:\Rat132\recordings\in_pyr'};%List of folders from the path

%Sleep
time_criteria = 1; % minimal time to include a NREM epoch (in min)

% for SU
criteria_fr = 0.1; %criteria to include or not a SU into the analysis
criteria_n = [6 6]; % minimal number of neurons from each structure [vHPC dHPC]
criteria_type = 0; %criteria for celltype (0:pyr, 1:int, 2:all)
binSize = 0.1;
n_SU_V = 0;
n_SU_D = 0;

win = 120; % time window for bin construction

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
        clear x states
        
        NREM.all(NREM.all(:,2)-NREM.all(:,1)<60*time_criteria,:)=[];
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
                    spiketrains_dHPC.pyr = [spiketrains_dHPC.pyr , (tmp)];
                end
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
                    spiketrains_vHPC.pyr = [spiketrains_vHPC.pyr , (tmp)];
                end
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
            tmp = Restrict(NREM.baseline,[NREM.baseline(end,2)-1800 , NREM.baseline(end,2)]);
            is.baseline.sws2 = InIntervals(bins,tmp); clear tmp
            
            tmp = Restrict(NREM.reward,[NREM.reward(1,1) , NREM.reward(1,1)+1800]);
            is.reward.sws1 = InIntervals(bins,tmp);  clear tmp
            
            tmp = Restrict(NREM.reward,[NREM.reward(end,2)-1800 , NREM.reward(end,2)]);
            is.reward.sws2 = InIntervals(bins,tmp);  clear tmp
            
            tmp = Restrict(NREM.aversive,[NREM.aversive(1,1) , NREM.aversive(1,1)+1800]);
            is.aversive.sws1 = InIntervals(bins,tmp);  clear tmp
            
            tmp = Restrict(NREM.aversive,[NREM.aversive(end,2)-1800 , NREM.aversive(end,2)]);
            is.aversive.sws2 = InIntervals(bins,tmp);  clear tmp
            
            is.baseline.rem = InIntervals(bins,REM.baseline);
            is.aversive.rem = InIntervals(bins,REM.aversive);
            is.reward.rem = InIntervals(bins,REM.reward);
            
            is.aversive.run = InIntervals(bins,movement.aversive);
            is.reward.run = InIntervals(bins,movement.reward);
            
            is.aversive.quiet  = InIntervals(bins,behavior.quiet.aversive);
            is.reward.quiet  = InIntervals(bins,behavior.quiet.reward);
            
%             is.aversive.run = InIntervals(bins,aversiveTS_run./1000);
%             is.reward.run = InIntervals(bins,rewardTS_run./1000);

            % Eliminations of bins when the animal is not mooving
             spiketrains_dHPC.pyr(is.aversive.quiet,:) = 0;
             spiketrains_vHPC.pyr(is.aversive.quiet,:) = 0;
             
             spiketrains_dHPC.pyr(is.reward.quiet,:) = 0;
             spiketrains_vHPC.pyr(is.reward.quiet,:) = 0;             

            %% Binning of sessions in 100 setps sleep sessions
            segmentation = [0 : win : segments.Var1(end)/1000]; clear dt
            tmp = [];
            for i = 2 : size(segmentation,2)-1
                tmp = [tmp , InIntervals(bins,[segmentation(i-1) segmentation(i+1)])];
            end
            segmentation = logical(tmp);
            clear tmp i
            
            correlations.cross = cell(1,size(segmentation,2));
            correlations.dHPC = cell(1,size(segmentation,2));
            correlations.vHPC = cell(1,size(segmentation,2));
            for i = 1 : size(segmentation,2)
                x = spiketrains_dHPC.pyr(segmentation(:,i),:);
                y = spiketrains_vHPC.pyr(segmentation(:,i),:);
                correlations.cross{i} = corr(x,y);
                clear x y
                
                x = spiketrains_dHPC.pyr(segmentation(:,i),:);
                y = spiketrains_dHPC.pyr(segmentation(:,i),:);
                correlations.dHPC{i} = corr(x,y);
                clear x y                
                
                x = spiketrains_vHPC.pyr(segmentation(:,i),:);
                y = spiketrains_vHPC.pyr(segmentation(:,i),:);
                correlations.vHPC{i} = corr(x,y);
                clear x y                          
                
            end
            clear i
            
            matrixC = [];
            matrixD = [];
            matrixV = [];
            for i = 1 : size(segmentation,2)
                tmpC = [];
                tmpD = [];
                tmpV = [];
                for ii = 1 : size(segmentation,2)
                    x = correlations.cross{i};
                    y = correlations.cross{ii};
                    c = corrcoef(x,y,'rows','complete');
                    tmpC = [tmpC c(1,2)];
                    clear x y c
                    
                    x = correlations.dHPC{i};
                    y = correlations.dHPC{ii};
                    c = corrcoef(x,y,'rows','complete');
                    tmpD = [tmpD c(1,2)];
                    clear x y c
                    
                    x = correlations.vHPC{i};
                    y = correlations.vHPC{ii};
                    c = corrcoef(x,y,'rows','complete');
                    tmpV = [tmpV c(1,2)];
                    clear x y c
                
                end
                matrixC = [matrixC ; tmpC];
                matrixD = [matrixD ; tmpD];
                matrixV = [matrixV ; tmpV];
                
                clear tmpC tmpD tmpV
            end
            
%             figure,
%             subplot(131)
%             matrixC(eye(size(matrixC))==1) = nan;
%             imagesc([300 : 300 : segments.Var1(end)/1000],[0 : 300 : segments.Var1(end)/1000],matrixC)
%             hold on
%             xline(aversiveTS_run(1)/1000,'r')
%             xline(aversiveTS_run(2)/1000,'r')
%             xline(rewardTS_run(1)/1000,'b')
%             xline(rewardTS_run(2)/1000,'b')
%             yline(aversiveTS_run(1)/1000,'r')
%             yline(aversiveTS_run(2)/1000,'r')
%             yline(rewardTS_run(1)/1000,'b')
%             yline(rewardTS_run(2)/1000,'b')
%             title('Cross')
            
            segmentation = [win : win : segments.Var1(end)/1000-win]; 
            i = InIntervals(segmentation,aversiveTS_run./1000);
            figure,
            i = matrixV(i,:);
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
