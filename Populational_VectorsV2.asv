clear
clc
% close all

%% Parameters
path = {'E:\Rat103\usable';'E:\Rat126\Ephys\in_Pyr';'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat132\recordings\in_pyr';'E:\Rat165\in_pyr\'};%List of folders from the path

% for SU
criteria_fr = 0; %criteria to include or not a SU into the analysis
criteria_n = [3 3]; % minimal number of neurons from each structure [vHPC dHPC]
criteria_type = 0; %criteria for celltype (0:pyr, 1:int, 2:all)
win = 120; %for assemblie detection qnd qxctivity strength
binSize = 0.025;

% Behavior
minimal_speed = 7; % minimal speed to detect quite periods
minimal_speed_time = 2; % minimal time to detect quite periods

reactivation.aversive.dvHPC = [];
reactivation.reward.dvHPC = [];
reactivation.aversive.dHPC = [];
reactivation.reward.dHPC = [];
reactivation.aversive.vHPC = [];
reactivation.reward.vHPC = [];

normalization = false; % to define if normalization over Reactivation Strength is applied or not
th = 5; % threshold for peak detection

gain.both.reward.pre = [];     gain.both.reward.post = [];
gain.both.aversive.pre = [];   gain.both.aversive.post = [];


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
        
%         %% Awake
%         disp('Uploading digital imputs')
%         % Load digitalin.mat
%         load('digitalin.mat')
%         
%         % Behavioral calculations
%         disp('Uploading DLC outputs')
%         camara = ((camara(:,2)-camara(:,1))/2)+camara(:,1);
%         % periods of movment during eacj condition
%         if rewardTS_run(1) < aversiveTS_run(1)
%             load('laps1.mat','posx','posy');
%             [camaraR,~] = find((camara(:,1)-rewardTS_run(1)/1000)>0,1,'first'); %TimeStamp of the begining of aversive
%             pos = [camara(camaraR : camaraR+length(posx)-1),posx,posy];
%             %interpolation of dropped frames
%             ejeX = pos(~isnan(pos(:,2)),1); dX = pos(~isnan(pos(:,2)),2); dX_int = interp1(ejeX , dX , pos(:,1));
%             ejeY = pos(~isnan(pos(:,3)),1); dY = pos(~isnan(pos(:,3)),3); dY_int = interp1(ejeY , dY , pos(:,1));
%             pos(:,2) =dX_int; pos(:,3) =dY_int;%saving corrected pos
%             behavior.pos.reward = [pos];
%             behavior.speed.reward = LinearVelocity(pos,0);
%             behavior.speed.reward(:,2) = smoothdata(behavior.speed.reward(:,2),'gaussian',1,'SamplePoints',behavior.speed.reward(:,1));
%             behavior.quiet.reward = QuietPeriods( behavior.speed.reward , minimal_speed , minimal_speed_time);
%             clear pos camaraR posx posy
%             load('laps2.mat','posx','posy');
%             [camaraA,~] = find((camara(:,1)-aversiveTS_run(1)/1000)>0,1,'first'); %TimeStamp of the begining of aversive
%             pos = [camara(camaraA : camaraA+length(posx)-1),posx,posy];
%             %interpolation of dropped frames
%             ejeX = pos(~isnan(pos(:,2)),1); dX = pos(~isnan(pos(:,2)),2); dX_int = interp1(ejeX , dX , pos(:,1));
%             ejeY = pos(~isnan(pos(:,3)),1); dY = pos(~isnan(pos(:,3)),3); dY_int = interp1(ejeY , dY , pos(:,1));
%             pos(:,2) =dX_int; pos(:,3) =dY_int;%saving corrected pos
%             behavior.pos.aversive = [pos];
%             behavior.speed.aversive = LinearVelocity(pos,0);
%             behavior.speed.aversive(:,2) = smoothdata(behavior.speed.aversive(:,2),'gaussian',1,'SamplePoints',behavior.speed.aversive(:,1));
%             behavior.quiet.aversive = QuietPeriods(behavior.speed.aversive , minimal_speed , minimal_speed_time);
%             clear pos camaraR
%         else
%             load('laps2.mat','posx','posy');
%             [camaraR,~] = find((camara(:,1)-rewardTS_run(1)/1000)>0,1,'first'); %TimeStamp of the begining of aversive
%             %         [camaraR2,~] = find((camara(:,1)-rewardTS_run(2)/1000)<0,1,'last'); %TimeStamp of the ending of aversive
%             pos = [camara(camaraR : camaraR+length(posx)-1),posx,posy];
%             %interpolation of dropped frames
%             ejeX = pos(~isnan(pos(:,2)),1); dX = pos(~isnan(pos(:,2)),2); dX_int = interp1(ejeX , dX , pos(:,1));
%             ejeY = pos(~isnan(pos(:,3)),1); dY = pos(~isnan(pos(:,3)),3); dY_int = interp1(ejeY , dY , pos(:,1));
%             pos(:,2) =dX_int; pos(:,3) =dY_int;%saving corrected pos
%             behavior.pos.reward = [pos];
%             behavior.speed.reward = LinearVelocity(pos,0);
%             behavior.speed.reward(:,2) = smoothdata(behavior.speed.reward(:,2),'gaussian',1,'SamplePoints',behavior.speed.reward(:,1));
%             behavior.quiet.reward = QuietPeriods( behavior.speed.reward , minimal_speed , minimal_speed_time);
%             clear pos camaraR posx posy
%             load('laps1.mat','posx','posy');
%             [camaraA,~] = find((camara(:,1)-aversiveTS_run(1)/1000)>0,1,'first'); %TimeStamp of the begining of aversive
%             pos = [camara(camaraA : camaraA+length(posx)-1),posx,posy];
%             %interpolation of dropped frames
%             ejeX = pos(~isnan(pos(:,2)),1); dX = pos(~isnan(pos(:,2)),2); dX_int = interp1(ejeX , dX , pos(:,1));
%             ejeY = pos(~isnan(pos(:,3)),1); dY = pos(~isnan(pos(:,3)),3); dY_int = interp1(ejeY , dY , pos(:,1));
%             pos(:,2) =dX_int; pos(:,3) =dY_int; %saving corrected pos
%             behavior.pos.aversive = [pos];
%             behavior.speed.aversive = LinearVelocity(pos,0);
%             behavior.speed.aversive(:,2) = smoothdata(behavior.speed.aversive(:,2),'gaussian',1,'SamplePoints',behavior.speed.aversive(:,1));
%             behavior.quiet.aversive = QuietPeriods(behavior.speed.aversive , minimal_speed , minimal_speed_time);
%             clear pos camaraA posx posy
%         end
%         
%         % Generation of no-movements periods
%         % Reward
%         start = behavior.speed.reward(1,1);   stop = behavior.speed.reward(end,1);
%         movement.reward = InvertIntervals(behavior.quiet.reward , start , stop); %keep only those higher than criteria
%         %         movement.reward(movement.reward(:,2) - movement.reward(:,1) <1,:)=[]; %eliminate 1sec segments
%         clear tmp start stop
%         % Aversive
%         start = behavior.speed.aversive(1,1);   stop = behavior.speed.aversive(end,1);
%         movement.aversive = InvertIntervals(behavior.quiet.aversive , start , stop);%keep only those higher than criteria
%         %         movement.aversive(movement.aversive(:,2) - movement.aversive(:,1) <1,:)=[];
%         clear tmp start stop
%         
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
        
        %         % to save the clusters I used for further analysis
        %         save([cd,'\SUclusters.mat'],'clusters')
        
        %% Assemblies detection
        if or(numberD > 3 , numberV > 3)
            % --- Aversive ---
            disp('Lets go for the assemblies')
            if isfile('dorsalventral_assemblies_aversive.mat')
                disp('Loading Aversive template')
                load('dorsalventral_assemblies_aversive.mat')
            end
            Thresholded.aversive.all = Th;
            patterns.all.aversive = pat;
            clear cond Th pat
            
            % Detection of members
            if not(isempty(Thresholded.aversive.all))
                if numberD>0
                    cond1 =  sum(Thresholded.aversive.all(1:size(clusters.dHPC,1),:),1)>0; %checking of dHPC SU
                    cond2 =  sum(Thresholded.aversive.all(size(clusters.dHPC,1)+1:end,:),1)>0; %checking of vHPC SU
                    cond.dHPC.aversive = and(cond1 , not(cond2));
                    cond.vHPC.aversive = and(cond2 , not(cond1));
                    cond.both.aversive = and(cond1 , cond2); clear cond1 cond2
                else
                    cond1 =  logical(zeros(1,size(Thresholded.aversive.all,2))); %checking of dHPC SU
                    cond2 =  logical(ones(1,size(Thresholded.aversive.all,2))); %checking of vHPC SU
                    cond.dHPC.aversive = and(cond1 , not(cond2));
                    cond.vHPC.aversive = and(cond2 , not(cond1));
                    cond.both.aversive = and(cond1 , cond2); clear cond1 cond2
                end
            else
                cond1 =  false; %checking of dHPC SU
                cond2 =  logical(0); %checking of vHPC SU
                cond.dHPC.aversive = and(cond1 , not(cond2));
                cond.vHPC.aversive = and(cond2 , not(cond1));
                cond.both.aversive = and(cond1 , cond2); clear cond1 cond2
            end
            
            % --- Reward ---
            disp('Loading Reward template')
            if isfile('dorsalventral_assemblies_reward.mat')
                load('dorsalventral_assemblies_reward.mat')
             end
            Thresholded.reward.all = Th;
            patterns.all.reward = pat;
            clear Th pat
            
            % Detection of members using
            if not(isempty(Thresholded.reward.all))
                if numberD>0
                    cond1 =  sum(Thresholded.reward.all(1:size(clusters.dHPC,1),:),1)>0; %checking of dHPC SU
                    cond2 =  sum(Thresholded.reward.all(size(clusters.dHPC,1)+1:end,:),1)>0; %checking of vHPC SU
                    cond.dHPC.reward = and(cond1 , not(cond2));
                    cond.vHPC.reward = and(cond2 , not(cond1));
                    cond.both.reward = and(cond1 , cond2); clear cond1 cond2
                else
                    cond1 =  logical(zeros(1,size(Thresholded.reward.all,2))); %checking of dHPC SU
                    cond2 =  logical(ones(1,size(Thresholded.reward.all,2))); %checking of vHPC SU
                    cond.dHPC.reward = and(cond1 , not(cond2));
                    cond.vHPC.reward = and(cond2 , not(cond1));
                    cond.both.reward = and(cond1 , cond2); clear cond1 cond2
                end
            else
                cond1 =  logical(0); %checking of dHPC SU
                cond2 =  logical(0); %checking of vHPC SU
                cond.dHPC.reward = and(cond1 , not(cond2));
                cond.vHPC.reward = and(cond2 , not(cond1));
                cond.both.reward = and(cond1 , cond2); clear cond1 cond2
            end
            
            
            %% save members id
            % Aversive joint assemblies
            if sum(cond.both.aversive)>0
                clusters.members.aversive.dHPC = cell(1,sum(cond.both.aversive));
                clusters.members.aversive.vHPC = cell(1,sum(cond.both.aversive));
                
                indD = Thresholded.aversive.all(1:size(clusters.dHPC,1),cond.both.aversive);
                indV = Thresholded.aversive.all(size(clusters.dHPC,1)+1:end,cond.both.aversive);
                
                for i = 1 : sum(cond.both.aversive)
                    clusters.members.aversive.dHPC{i} = clusters.dHPC(indD(:,i));
                    clusters.members.aversive.vHPC{i} = clusters.vHPC(indV(:,i));
                end
                clear indD indV
            end
            
            % Reward joint assemblies
            if sum(cond.both.reward)>0
                clusters.members.reward.dHPC = cell(1,sum(cond.both.reward));
                clusters.members.reward.vHPC = cell(1,sum(cond.both.reward));
                
                indD = Thresholded.reward.all(1:size(clusters.dHPC,1),cond.both.reward);
                indV = Thresholded.reward.all(size(clusters.dHPC,1)+1:end,cond.both.reward);
                
                for i = 1 : sum(cond.both.reward)
                    clusters.members.reward.dHPC{i} = clusters.dHPC(indD(:,i));
                    clusters.members.reward.vHPC{i} = clusters.vHPC(indV(:,i));
                end
                clear indD indV
            end
            
            %% Constructing Spiketrains
            freq = 1/binSize;
            limits = [0 segments.Var1(end)/1000];
            spiketrains_dHPC.aversive = cell(1,sum(cond.both.aversive));
            spiketrains_vHPC.aversive = cell(1,sum(cond.both.aversive));
            for i = 1:sum(cond.both.aversive)
                [Spikes , bins , Clusters] = spike_train_construction(spks_dHPC, clusters.members.aversive.dHPC{i}, cellulartype, binSize, limits, [], false, true);
                spiketrains_dHPC.aversive{i} = Spikes;    
                clear Spikes Clusters
                
                [Spikes , bins , Clusters] = spike_train_construction(spks_vHPC, clusters.members.aversive.vHPC{i}, cellulartype, binSize, limits, [], false, true);
                spiketrains_vHPC.aversive{i} = Spikes;    
                clear Spikes Clusters                
            end
            
            spiketrains_dHPC.reward = cell(1,sum(cond.both.reward));
            spiketrains_vHPC.reward = cell(1,sum(cond.both.reward));
            for i = 1:sum(cond.both.reward)
                [Spikes , bins , Clusters] = spike_train_construction(spks_dHPC, clusters.members.reward.dHPC{i}, cellulartype, binSize, limits, [], false, true);
                spiketrains_dHPC.reward{i} = Spikes;    
                clear Spikes Clusters
                
                [Spikes , bins , Clusters] = spike_train_construction(spks_vHPC, clusters.members.reward.vHPC{i}, cellulartype, binSize, limits, [], false, true);
                spiketrains_vHPC.reward{i} = Spikes;    
                clear Spikes Clusters                
            end   
            
            %% Correlation calculation
            if and(size(spiketrains_vHPC.pyr,2) >= criteria_n(1),size(spiketrains_dHPC.pyr,2) >= criteria_n(2))
                %% Binning of sessions in 100 setps sleep sessions
                % construction of time vector
                s = [0 : win : segments.Var1(end)/1000]; clear dt 
                % Restriction of time vector to conditions and NREM
                segmentation = s(InIntervals(s,sort([NREM.all;aversiveTS_run./1000;rewardTS_run./1000])));
                % creation of iterator
                tmp = [];
                for i = 2 : size(segmentation,2)-1
                    tmp = [tmp , InIntervals(bins,[segmentation(i-1) segmentation(i+1)])];
                end
                segmentation = logical(tmp);
                clear tmp i
                
                % Aversive Joint assemblies
                correlations.aversive = cell(sum(cond.both.aversive),size(segmentation,2));
                matrix.aversive = cell(sum(cond.both.aversive),1);
                for ii = 1 : sum(cond.both.aversive)
                    for i = 1 : size(segmentation,2)
                        x = spiketrains_dHPC.aversive{ii}(segmentation(:,i),:);
                        y = spiketrains_vHPC.aversive{ii}(segmentation(:,i),:);
                        correlations.aversive{ii,i} = corr(x,y);
                        clear x y
                    end
                    clear i
                    
                    for ind1 = 1 : size(segmentation,2)
                        tmpC = [];
                        for ind2 = 1 : size(segmentation,2)
                            x = correlations.aversive{ii,ind1};
                            y = correlations.aversive{ii,ind2};
                            c = corrcoef(x,y,'rows','complete');
                            tmpC = [tmpC c(1,2)];
                            clear x y c
                        end
                        matrix.aversive{ii} = [matrix.aversive{ii} ; tmpC];
                        clear tmpC
                    end
                end
                
    
                
                % Reward Joint assemblies
                correlations.reward = cell(sum(cond.both.reward),size(segmentation,2));
                matrix.reward = cell(sum(cond.both.reward),1);
                for ii = 1 : sum(cond.both.reward)
                    for i = 1 : size(segmentation,2)
                        x = spiketrains_dHPC.reward{ii}(segmentation(:,i),:);
                        y = spiketrains_vHPC.reward{ii}(segmentation(:,i),:);
                        correlations.reward{ii,i} = corr(x,y);
                        clear x y
                    end
                    clear i
                    
                    for ind1 = 1 : size(segmentation,2)
                        tmpC = [];
                        for ind2 = 1 : size(segmentation,2)
                            x = correlations.reward{ii,ind1};
                            y = correlations.reward{ii,ind2};
                            c = corrcoef(x,y,'rows','complete');
                            tmpC = [tmpC c(1,2)];
                            clear x y c
                        end
                        matrix.reward{ii} = [matrix.reward{ii} ; tmpC];
                        clear tmpC
                    end
                end                
                
                
                
                
                
                
                
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
                
                
                
            end
        end
    end
