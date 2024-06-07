clear
clc
close all
%% PARAMETERS - run this section before run any other
path = {'\\Maryjackson\e\Rat127\Ephys\pyr';'\\Maryjackson\e\Rat128\Ephys\in_pyr\ready';'\\Maryjackson\e\Rat103\usable';'\\Maryjackson\e\Rat132\recordings\in_pyr'; '\\Maryjackson\e\Rat165\in_pyr'; ...
    '\\Maryjackson\e\Rat126\Ephys\in_Pyr'};%List of folders from the path
% path = {'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat103\usable';'E:\Rat132\recordings\in_pyr'};%List of folders from the path
% for SU
criteria_fr = 0.01; %criteria to include or not a SU into the analysis
criteria_n = 6; % minimal number of neurons from each structure
criteria_type = 0; %criteria for celltype (0:pyr, 1:int, 2:all)
n_SU_V = [];
n_SU_D = [];
Xedges = 60; %number of bins for RateMap construction - 2.5cm bin (because of the extremes cutting)
sigma = 2;%round(15/(180/Xedges)); %defined for gauss kernel of 15cm
binSize = 0.001; % bin size for replay events detection
bin_size = 1; % to bin pos ans spks in between/within  
% Behavior
minimal_speed = 2.5;% minimal speed to detect quite periods
minimal_speed_time = 2; % minimal time to detect quite periods

%% MAIN LOOP ORIGINAL, to iterate across sessions  - within/between by bins 
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
        
        %Session output
        dHPC = {};% one cell per pc
        vHPC = {}; 
        
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
        % periods of movment during each condition
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
        %Video sampling rate 
        dt = (mean(diff(behavior.pos.aversive(:,1)))); 
        
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
        
        for ii=1:size(group_dHPC,1)

            cluster = group_dHPC(ii,1);
            celltype = logical(cellulartype(cellulartype(:,1) == cluster,2));
            
            if celltype % check if pyr
                spks = spks_dHPC(spks_dHPC(:,1)==cluster,2); % select tspk from cluster

                % --- Aversive ---
                spks_ave = spks; 
                pos_ave = behavior.pos.aversive(:,1:2); 
                in_maze = ToIntervals(pos_ave(:,1),and(pos_ave(:,2)>30 , pos_ave(:,2)<180));% eliminating the extrems of the maze
                    
                %Check lenght of in_maze segments to define real laps
                in_lapA = []; % keeps only full laps (maybe too restrictive)
                    for ix=1:size(in_maze,1)
                        ini = in_maze(ix,1);
                        fin = in_maze(ix,2);
                        ini_pos = pos_ave(pos_ave(:,1)==ini,2); 
                        fin_pos = pos_ave(pos_ave(:,1)==fin,2); 
                        if abs(fin_pos-ini_pos)>100
                             in_lapA = [in_lapA;in_maze(ix,:)]; 
                        end
                    end
                    
                %Restrict spk and position to laps and movement periods
                spks_ave = Restrict(spks_ave,in_lapA);spks_ave = Restrict(spks_ave , movement.aversive); % Restrict spikes to movement periods and laps
                    
                pos_ave = Restrict(pos_ave,in_lapA);pos_ave = Restrict(pos_ave, movement.aversive); % Restrict pos to movement periods and laps
                pos_ave(:,2) = pos_ave(:,2)-min(pos_ave(:,2)); pos_ave(:,2) = pos_ave(:,2)/max(pos_ave(:,2)); %normalization of position
               
                % --- Reward ---
                spks_rew = spks; 
                pos_rew = behavior.pos.reward(:,1:2);
                in_maze = ToIntervals(pos_rew(:,1),and(pos_rew(:,2)>30 , pos_rew(:,2)<180));% eliminating the extrems of the maze
               
                %Check lenght of in_maze segments to define real laps
                in_lapR = []; 
                    for ix=1:size(in_maze,1)
                        ini = in_maze(ix,1);
                        fin = in_maze(ix,2);
                        ini_pos = pos_rew(pos_rew(:,1)==ini,2); 
                        fin_pos = pos_rew(pos_rew(:,1)==fin,2); 
                        if abs(fin_pos-ini_pos)>100
                             in_lapR = [in_lapR;in_maze(ix,:)]; 
                        end
                    end 
                    
                %Restrict spk and position to laps and movement periods
                spks_rew = Restrict(spks_rew,in_lapR); spks_rew = Restrict(spks_rew, movement.reward); % Restrict to movement periods
                pos_rew = Restrict(pos_rew,in_lapR); pos_rew = Restrict(pos_rew , movement.reward);
                pos_rew(:,2) = pos_rew(:,2)-min(pos_rew(:,2)); pos_rew(:,2) = pos_rew(:,2)/max(pos_rew(:,2)); %normalization of position 
                    
                    %%%%%%%%Control plot %%%%%%%%%%
%                     figure(1);clf; hold on;subplot(1,2,1); plot(pos_ave(:,2),pos_ave(:,1),'color', [.3 .3 .3]);hold on; 
%                     Find interpolated position of each spike:
%                     xs = interp1(pos_ave(:,1),pos_ave(:,2),spks_ave);scatter(xs,spks_ave,'*'); 
%                     title(['Aversive - Mean fr:',num2str(fr_ave)]);
%                    
%                     subplot(1,2,2); plot(pos_rew(:,2),pos_rew(:,1),'color', [.3 .3 .3]);hold on; 
%                     Find interpolated position of each spike:
%                     xs = interp1(pos_rew(:,1),pos_rew(:,2),spks_rew); scatter(xs,spks_rew,'*');
%                     title(['Reward- Mean fr:',num2str(fr_rew)]);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
   
                %Checking pc criteria
                % 1. Mean firing rate above tresh; 2. At least 1 pf > 4 bins in one of the cond  
     
                [curveA , statsA] = FiringCurve(pos_ave , spks_ave , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                [curveR , statsR] = FiringCurve(pos_rew , spks_rew , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    
                fr_ave= size(spks_ave,1)/sum(in_lapA(:,2)-in_lapA(:,1)); 
                fr_rew= size(spks_rew,1)/sum(in_lapR(:,2)-in_lapR(:,1)); 
                tresh=0.1; % mean firng rate treshold 
                    
                %if there is no PF, assigne 0 to correct assesment in the next if
                if isempty(statsA.field)
                      statsA.field=0;
                end 
                if isempty(statsR.field)
                      statsR.field=0;
                end 
                    
                    if and(or(sum(statsA.field(:,:,1))>= 4 , sum(statsR.field(:,:,1))>= 4),or(fr_ave>tresh,fr_rew>tresh)) % you are a PC
                        %Within aversive
                        [withinA,withinA_tresh] = Within_pc(pos_ave,spks_ave,1,sigma,Xedges);
                        %Within reward 
                        [withinR,withinR_tresh] = Within_pc(pos_rew,spks_rew,1,sigma,Xedges);
                        %Between
                        [between] = Between_pc(pos_ave,spks_ave,pos_rew,spks_rew,bin_size,sigma,Xedges);
                        %Subsampled
                        [sub_pos_ave,sub_spk_ave,sub_pos_rew,sub_spk_rew] = Subsampling_1d(pos_ave, spks_ave,pos_rew, spks_rew,dt,Xedges,sigma);
                        [between_sub] = Between_pc(sub_pos_ave,sub_spk_ave,sub_pos_rew,sub_spk_rew,bin_size,sigma,Xedges);
                        [withinA_sub,~] = Within_pc(sub_pos_ave,sub_spk_ave,1,sigma,Xedges);
                        [withinR_sub,~] = Within_pc(sub_pos_rew,sub_spk_rew,1,sigma,Xedges);
                        
                    
                        % Save PC varaiables 
                        n.id = cluster; 
                        n.nlap_ave = size(in_lapA,1);
                        n.nlap_rew = size(in_lapR,1);
                        n.frMap_ave = curveA.rate;
                        n.frMap_rew = curveR.rate;
                        n.stats_ave = statsA;
                        n.stats_rew = statsR;
                        n.between = between;
                        n.within_ave = withinA; 
                        n.within_rew = withinR;
                        n.subsampled.between = between_sub;
                        n.subsampled.within_ave = withinA_sub; 
                        n.subsampled.within_rew = withinR_sub; 
                        
                        dHPC{ii}= n;
                   end 
          
            end
         end
           
        disp('vHPC Firing rate map calculation')
     
        for ii=1:size(group_vHPC,1)

            cluster = group_vHPC(ii,1);
            celltype = logical(cellulartype(cellulartype(:,1) == cluster,2));
            
            if celltype % check if pyr
                spks = spks_vHPC(spks_vHPC(:,1)==cluster,2); % select tspk from cluster

                % --- Aversive ---
                spks_ave = spks; 
                pos_ave = behavior.pos.aversive(:,1:2); 
                in_maze = ToIntervals(pos_ave(:,1),and(pos_ave(:,2)>30 , pos_ave(:,2)<180));% eliminating the extrems of the maze
                    
                %Check lenght of in_maze segments to define real laps
                in_lapA = []; % keeps only full laps (maybe too restrictive)
                    for ix=1:size(in_maze,1)
                        ini = in_maze(ix,1);
                        fin = in_maze(ix,2);
                        ini_pos = pos_ave(pos_ave(:,1)==ini,2); 
                        fin_pos = pos_ave(pos_ave(:,1)==fin,2); 
                        if abs(fin_pos-ini_pos)>100
                             in_lapA = [in_lapA;in_maze(ix,:)]; 
                        end
                    end
                    
                %Restrict spk and position to laps and movement periods
                spks_ave = Restrict(spks_ave,in_lapA);spks_ave = Restrict(spks_ave , movement.aversive); % Restrict spikes to movement periods and laps
                    
                pos_ave = Restrict(pos_ave,in_lapA);pos_ave = Restrict(pos_ave, movement.aversive); % Restrict pos to movement periods and laps
                pos_ave(:,2) = pos_ave(:,2)-min(pos_ave(:,2)); pos_ave(:,2) = pos_ave(:,2)/max(pos_ave(:,2)); %normalization of position
               
                % --- Reward ---
                spks_rew = spks; 
                pos_rew = behavior.pos.reward(:,1:2);
                in_maze = ToIntervals(pos_rew(:,1),and(pos_rew(:,2)>30 , pos_rew(:,2)<180));% eliminating the extrems of the maze
               
                %Check lenght of in_maze segments to define real laps
                in_lapR = []; 
                    for ix=1:size(in_maze,1)
                        ini = in_maze(ix,1);
                        fin = in_maze(ix,2);
                        ini_pos = pos_rew(pos_rew(:,1)==ini,2); 
                        fin_pos = pos_rew(pos_rew(:,1)==fin,2); 
                        if abs(fin_pos-ini_pos)>100
                             in_lapR = [in_lapR;in_maze(ix,:)]; 
                        end
                    end 
                    
                %Restrict spk and position to laps and movement periods
                spks_rew = Restrict(spks_rew,in_lapR); spks_rew = Restrict(spks_rew, movement.reward); % Restrict to movement periods
                pos_rew = Restrict(pos_rew,in_lapR); pos_rew = Restrict(pos_rew , movement.reward);
                pos_rew(:,2) = pos_rew(:,2)-min(pos_rew(:,2)); pos_rew(:,2) = pos_rew(:,2)/max(pos_rew(:,2)); %normalization of position 
                    
                    %%%%%%%%Control plot %%%%%%%%%%
%                     figure(1);clf; hold on;subplot(1,2,1); plot(pos_ave(:,2),pos_ave(:,1),'color', [.3 .3 .3]);hold on; 
%                     Find interpolated position of each spike:
%                     xs = interp1(pos_ave(:,1),pos_ave(:,2),spks_ave);scatter(xs,spks_ave,'*'); 
%                     title(['Aversive - Mean fr:',num2str(fr_ave)]);
%                    
%                     subplot(1,2,2); plot(pos_rew(:,2),pos_rew(:,1),'color', [.3 .3 .3]);hold on; 
%                     Find interpolated position of each spike:
%                     xs = interp1(pos_rew(:,1),pos_rew(:,2),spks_rew); scatter(xs,spks_rew,'*');
%                     title(['Reward- Mean fr:',num2str(fr_rew)]);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
   
                %Checking pc criteria
                % 1. Mean firing rate above tresh; 2. At least 1 pf > 4 bins in one of the cond  
     
                [curveA , statsA] = FiringCurve(pos_ave , spks_ave , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                [curveR , statsR] = FiringCurve(pos_rew , spks_rew , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    
                fr_ave= size(spks_ave,1)/sum(in_lapA(:,2)-in_lapA(:,1)); 
                fr_rew= size(spks_rew,1)/sum(in_lapR(:,2)-in_lapR(:,1)); 
                tresh=0.1; % mean firng rate treshold 
                    
                %if there is no PF, assigne 0 to correct assesment in the next if
                if isempty(statsA.field)
                      statsA.field=0;
                end 
                if isempty(statsR.field)
                      statsR.field=0;
                end 
                    
                    if and(or(sum(statsA.field(:,:,1))>= 4 , sum(statsR.field(:,:,1))>= 4),or(fr_ave>tresh,fr_rew>tresh)) % you are a PC
                        %Within aversive
                        [withinA,withinA_tresh] = Within_pc(pos_ave,spks_ave,1,sigma,Xedges);
                        %Within reward 
                        [withinR,withinR_tresh] = Within_pc(pos_rew,spks_rew,1,sigma,Xedges);
                        %Between
                        [between] = Between_pc(pos_ave,spks_ave,pos_rew,spks_rew,bin_size,sigma,Xedges);
                        %Subsampled
                        [sub_pos_ave,sub_spk_ave,sub_pos_rew,sub_spk_rew] = Subsampling_1d(pos_ave, spks_ave,pos_rew, spks_rew,dt,Xedges,sigma);
                        [between_sub] = Between_pc(sub_pos_ave,sub_spk_ave,sub_pos_rew,sub_spk_rew,bin_size,sigma,Xedges);
                        [withinA_sub,~] = Within_pc(sub_pos_ave,sub_spk_ave,1,sigma,Xedges);
                        [withinR_sub,~] = Within_pc(sub_pos_rew,sub_spk_rew,1,sigma,Xedges);
                        
                    
                        % Save PC varaiables 
                        n.id = cluster; 
                        n.nlap_ave = size(in_lapA,1);
                        n.nlap_rew = size(in_lapR,1);
                        n.frMap_ave = curveA.rate;
                        n.frMap_rew = curveR.rate;
                        n.stats_ave = statsA;
                        n.stats_rew = statsR;
                        n.between = between;
                        n.within_ave = withinA; 
                        n.within_rew = withinR;
                        n.subsampled.between = between_sub;
                        n.subsampled.within_ave = withinA_sub; 
                        n.subsampled.within_rew = withinR_sub; 
                        
                        vHPC{ii}= n;
                   end 
          
            end
         end

    %% Saveing PC INFO 
    if ~isempty(dHPC)
        dHPC = dHPC(~cellfun('isempty',dHPC));
        save([cd,'\dHPC_pc.mat'],'dHPC'); 
    end
    if ~isempty(vHPC)
        vHPC = vHPC(~cellfun('isempty',vHPC));
        save([cd,'\vHPC_pc.mat'],'vHPC'); 
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

%% COPY of main loop to re-run only some parameters 
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
        
        %Session output
        dHPC = {};% one cell per pc
        vHPC = {}; 
        
        if exist('Spikesorting/dHPC_pc.mat', 'file')~= 0
            load('Spikesorting/dHPC_pc.mat')
        end 
        
        if exist('Spikesorting/vHPC_pc.mat', 'file')~= 0
            load('Spikesorting/vHPC_pc.mat')
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
        % periods of movment during each condition
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
        %Video sampling rate 
        dt = (mean(diff(behavior.pos.aversive(:,1)))); 
        
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
        count=1;  
        for ii=1:size(group_dHPC,1)

            cluster = group_dHPC(ii,1);
            celltype = logical(cellulartype(cellulartype(:,1) == cluster,2));
            
            if celltype % check if pyr
                
                spks = spks_dHPC(spks_dHPC(:,1)==cluster,2); % select tspk from cluster
              
                % Select neurons with fr greater than tresh in aversive or reward - PC criteria 1   
                tresh = 0.1; % hz
                spks_tmp = Restrict(spks,movement.reward);
                fr_rew= size(spks_tmp,1)/sum(movement.reward(:,2)-movement.reward(:,1));
                
                spks_tmp = Restrict(spks,movement.aversive);
                fr_ave= size(spks_tmp,1)/sum(movement.aversive(:,2)-movement.aversive(:,1));
                
                if fr_rew >= tresh || fr_ave >= tresh
                    m = 1;
                else 
                    m = nan;
                end 
                
                if ~isnan(m)
                    % --- Aversive ---
                    spks_tmp = spks; 
                    pos_tmp = behavior.pos.aversive(:,1:2); 
                    in_maze = ToIntervals(pos_tmp(:,1),and(pos_tmp(:,2)>30 , pos_tmp(:,2)<180));% eliminating the extrems of the maze
                    
                    %Check lenght of in_maze segments to define real laps
                    in_lapA = []; % keeps only full laps (maybe too restrictive)
                    for ix=1:size(in_maze,1)
                        ini = in_maze(ix,1);
                        fin = in_maze(ix,2);
                        ini_pos = pos_tmp(pos_tmp(:,1)==ini,2); 
                        fin_pos = pos_tmp(pos_tmp(:,1)==fin,2); 
                        if abs(fin_pos-ini_pos)>100
                             in_lapA = [in_lapA;in_maze(ix,:)]; 
                        end
                    end
                    
%                     %%%%%%%%Control plot
%                     figure(2);clf; hold on; subplot(1,2,1);plot(pos_tmp(:,2),pos_tmp(:,1));hold on; 
%                     ini_lap = [];
%                     fin_lap = []; 
%                     for ix=1:size(in_maze,1)
%                         ini = in_maze(ix,1);
%                         fin = in_maze(ix,2);
%                         ini_lap = [ini_lap; ini, pos_tmp(pos_tmp(:,1)==ini,2)]; 
%                         fin_lap = [fin_lap; fin, pos_tmp(pos_tmp(:,1)==fin,2)]; 
%                     end 
%                     scatter(ini_lap(:,2),ini_lap(:,1),'r');scatter(fin_lap(:,2),fin_lap(:,1),'b')
%                     title('In Maze');
%                     subplot(1,2,2);hold on;plot(pos_tmp(:,2),pos_tmp(:,1));
%                     ini_lap = [];
%                     fin_lap = []; 
%                     for ix=1:size(in_lapA,1)
%                         ini = in_lapA(ix,1);
%                         fin = in_lapA(ix,2);
%                         ini_lap = [ini_lap; ini, pos_tmp(pos_tmp(:,1)==ini,2)]; 
%                         fin_lap = [fin_lap; fin, pos_tmp(pos_tmp(:,1)==fin,2)]; 
%                     end 
%                     scatter(ini_lap(:,2),ini_lap(:,1),'r');scatter(fin_lap(:,2),fin_lap(:,1),'b')
%                     title('In Lap');
%                     %%%%%%%%
                    
                    %Restrict spk and position to laps and movement periods
                    spks_tmp = Restrict(spks_tmp,in_lapA);spks_tmp = Restrict(spks_tmp , movement.aversive); % Restrict spikes to movement periods
                    pos_tmp = Restrict(pos_tmp,in_lapA);pos_tmp = Restrict(pos_tmp, movement.aversive); % Restrict pos to movement periods

                    pos_tmp(:,2) = pos_tmp(:,2)-min(pos_tmp(:,2)); pos_tmp(:,2) = pos_tmp(:,2)/max(pos_tmp(:,2)); %normalization of position
                    
                    %Store par for between comp
                    spks_ave = spks_tmp;
                    pos_ave = pos_tmp;
                    
                    %%%%%%%%Control plot %%%%%%%%%%
%                     figure(3);clf; hold on; plot(pos_tmp(:,2),pos_tmp(:,1),'color', [.3 .3 .3]);hold on; 
%                     %Find interpolated position of each spike:
%                     xs = interp1(pos_tmp(:,1),pos_tmp(:,2),spks_tmp);
%                     scatter(xs,spks_tmp,'*'); 
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    

                    %Firing curve construction
                    [curveA , statsA] = FiringCurve(pos_tmp , spks_tmp , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    curveA_dhpc = curveA.rate; 
                    
                    %%%%%%%% Control plot %%%%%%%%
%                     figure(1);clf;hold on; 
%                     subplot(2,1,1);imagesc(curveA.rate), colormap 'jet'
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    %%%Within-trial pc parameters or Lap parameters
%                     if contains(comp,'whitin') 
                        [withinA,withinA_tresh] = Within_pc(pos_tmp,spks_tmp,1,sigma,Xedges);
%                     elseif contains(comp,'laps')
                        %Compare laps 
                        [withinA_lap] = Within_lap(pos_tmp,spks_tmp,in_lapA,sigma,Xedges);
                        [withinA_lap_r] = Within_lap_random(pos_tmp,spks_tmp,in_lapA,sigma,Xedges);
%                     else
%                         disp('No comparison method selected')
%                     end 
                    
                    % --- Reward ---
                    spks_tmp = spks; 
                    pos_tmp = behavior.pos.reward(:,1:2);
                    in_maze = ToIntervals(pos_tmp(:,1),and(pos_tmp(:,2)>30 , pos_tmp(:,2)<180));% eliminating the extrems of the maze
                    
                    %Check lenght of in_maze segments to define real laps
                    in_lapR = []; 
                    for ix=1:size(in_maze,1)
                        ini = in_maze(ix,1);
                        fin = in_maze(ix,2);
                        ini_pos = pos_tmp(pos_tmp(:,1)==ini,2); 
                        fin_pos = pos_tmp(pos_tmp(:,1)==fin,2); 
                        if abs(fin_pos-ini_pos)>100
                             in_lapR = [in_lapR;in_maze(ix,:)]; 
                        end
                    end 
                    
%                     %%%%%%%%Control plot%%%%%%%%
%                     figure(3);clf; hold on; subplot(1,2,1);plot(pos_tmp(:,2),pos_tmp(:,1));hold on; 
%                     ini_lap = [];
%                     fin_lap = []; 
%                     for ix=1:size(in_maze,1)
%                         ini = in_maze(ix,1);
%                         fin = in_maze(ix,2);
%                         ini_lap = [ini_lap; ini, pos_tmp(pos_tmp(:,1)==ini,2)]; 
%                         fin_lap = [fin_lap; fin, pos_tmp(pos_tmp(:,1)==fin,2)]; 
%                     end 
%                     scatter(ini_lap(:,2),ini_lap(:,1),'r');scatter(fin_lap(:,2),fin_lap(:,1),'b')
%                     title('In Maze');
%                     subplot(1,2,2);hold on;plot(pos_tmp(:,2),pos_tmp(:,1));
%                     ini_lap = [];
%                     fin_lap = []; 
%                     for ix=1:size(in_lapR,1)
%                         ini = in_lapR(ix,1);
%                         fin = in_lapR(ix,2);
%                         ini_lap = [ini_lap; ini, pos_tmp(pos_tmp(:,1)==ini,2)]; 
%                         fin_lap = [fin_lap; fin, pos_tmp(pos_tmp(:,1)==fin,2)]; 
%                     end 
%                     scatter(ini_lap(:,2),ini_lap(:,1),'r');scatter(fin_lap(:,2),fin_lap(:,1),'b')
%                     title('In Lap');
%                     %%%%%%%%%%%%%%%%%

                    %Restrict spks and pos to valid laps 
                    spks_tmp = Restrict(spks_tmp,in_lapR); spks_tmp = Restrict(spks_tmp, movement.reward); % Restrict to movement periods
                    pos_tmp = Restrict(pos_tmp,in_lapR); pos_tmp = Restrict(pos_tmp , movement.reward);
                    pos_tmp(:,2) = pos_tmp(:,2)-min(pos_tmp(:,2)); pos_tmp(:,2) = pos_tmp(:,2)/max(pos_tmp(:,2)); %normalization of position
                    %Store par for between comp
                    spks_rew = spks_tmp;
                    pos_rew = pos_tmp;
                    
                    %%%%%%%Control plot %%%%%%%%%%
%                     figure(4);clf; hold on; plot(pos_tmp(:,2),pos_tmp(:,1),'color', [.3 .3 .3]);hold on; 
%                     %Find interpolated position of each spike:
%                     xs = interp1(pos_tmp(:,1),pos_tmp(:,2),spks_tmp);
%                     scatter(xs,spks_tmp,'*'); 
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    %Firing curve construction
                    [curveR , statsR] = FiringCurve(pos_tmp , spks_tmp , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4 , 'minPeak' , 0.2);
                    curveR_dhpc = curveR.rate;
                    
%                    subplot(2,1,2);imagesc(curveR.rate), colormap 'jet'
%                    sgtitle('Firing curvs original');
                    
                    %%%Within-trial pc parameters or Lap parameters
%                     if contains(comp,'whitin') 
                        [withinR,withinR_tresh] = Within_pc(pos_tmp,spks_tmp,1,sigma,Xedges);
%                     elseif contains(comp,'laps')
                        %Compare laps 
                        [withinR_lap] = Within_lap(pos_tmp,spks_tmp,in_lapR,sigma,Xedges);
                        [withinR_lap_r] = Within_lap_random(pos_tmp,spks_tmp,in_lapR,sigma,Xedges);
%                     else
%                         disp('No comparison method selected')
%                     end 
                    
                    %Subsampling pos and spk for both conditions
                    [sub_pos_ave,sub_spk_ave,sub_pos_rew,sub_spk_rew] = Subsampling_1d(pos_ave, spks_ave,pos_rew, spks_rew,dt,Xedges,sigma);
                    
                    %%%PC criteria: have a PF bigger or equal to 4 bins 
                    %if there is no PF, assigne 0 to correct assesment in
                    %the next if
                    if isempty(statsA.field)
                        statsA.field=0;
                    end 
                     if isempty(statsR.field)
                        statsR.field=0;
                    end 
                    
                    if or(sum(statsA.field(:,:,1))>= 4 , sum(statsR.field(:,:,1))>= 4) 
              
                        %Subsampling remapping variables   
                        [between_sub] = Between_pc(sub_pos_ave,sub_spk_ave,sub_pos_rew,sub_spk_rew,bin_size,sigma,Xedges);
                        [withinA_sub,~] = Within_pc(sub_pos_ave,sub_spk_ave,1,sigma,Xedges);
                        [withinR_sub,~] = Within_pc(sub_pos_rew,sub_spk_rew,1,sigma,Xedges);

                        %%%Between pc parameters or Between Lap parameters
%                         if contains(comp,'whitin') 
                            [between] = Between_pc(pos_ave,spks_ave,pos_rew,spks_rew,bin_size,sigma,Xedges);
%                         elseif contains(comp,'laps')
                        %Compare laps 
                            [between_lap] = Between_lap(pos_ave,spks_ave,in_lapA,pos_rew,spks_rew,in_lapR,sigma,Xedges);
                            [between_lap_r] = Between_lap_random(pos_ave,spks_ave,in_lapA,pos_rew,spks_rew,in_lapR,sigma,Xedges);
%                         else
%                             disp('No comparison method selected')
%                         end 
                        

                        %Store pc info 
                        n.ave = withinA_lap; 
                        n.rew = withinR_lap ; 
                        n.between = between_lap; 
                        n.ave_r = withinA_lap_r; 
                        n.rew_r = withinR_lap_r ; 
                        n.between_r = between_lap_r;
                        dHPC{count}.lap= n;
                        dHPC{count}.nlap_rew= size(in_lapR,1);
                        dHPC{count}.nlap_ave= size(in_lapA,1);
                        
                        count = count+1; 
                    end
                    clear s q curve OccMap Nspikes spks_tmp pos_tmp m curve stats sp curve1
                    clear curveA curveR qA qR curveSA curveSR statsSA statsSR
                end
            end
            clear celltype tmp b cluster
            
        end
        
        disp('vHPC Firing rate map calculation')
    count = 1; 
    for ii=1:size(group_vHPC,1)
            cluster = group_vHPC(ii,1);
            celltype = logical(cellulartype(cellulartype(:,1) == cluster,2));
            
            if celltype % check if pyr
                
                spks = spks_vHPC(spks_vHPC(:,1)==cluster,2); % select tspk from cluster
              
                % Select neurons with fr greater than tresh in aversive or reward - PC criteria 1   
                tresh = 0.1; % hz
                spks_tmp = Restrict(spks,movement.reward);
                fr_rew= size(spks_tmp,1)/sum(movement.reward(:,2)-movement.reward(:,1));
                
                spks_tmp = Restrict(spks,movement.aversive);
                fr_ave= size(spks_tmp,1)/sum(movement.aversive(:,2)-movement.aversive(:,1));
                
                if fr_rew >= tresh || fr_ave >= tresh
                    m = 1;
                else 
                    m = nan;
                end 
                
                if ~isnan(m)
                    % --- Aversive ---
                    spks_tmp = spks; 
                    pos_tmp = behavior.pos.aversive(:,1:2); 
                    in_maze = ToIntervals(pos_tmp(:,1),and(pos_tmp(:,2)>30 , pos_tmp(:,2)<180));% eliminating the extrems of the maze
                    
                    %Check lenght of in_maze segments to define real laps
                    in_lapA = []; % keeps only full laps (maybe too restrictive)
                    for ix=1:size(in_maze,1)
                        ini = in_maze(ix,1);
                        fin = in_maze(ix,2);
                        ini_pos = pos_tmp(pos_tmp(:,1)==ini,2); 
                        fin_pos = pos_tmp(pos_tmp(:,1)==fin,2); 
                        if abs(fin_pos-ini_pos)>100
                             in_lapA = [in_lapA;in_maze(ix,:)]; 
                        end
                    end
                    
                    %%%%%%%%Control plot
%                     figure(2);clf; hold on; subplot(1,2,1);plot(pos_tmp(:,2),pos_tmp(:,1));hold on; 
%                     ini_lap = [];
%                     fin_lap = []; 
%                     for ix=1:size(in_maze,1)
%                         ini = in_maze(ix,1);
%                         fin = in_maze(ix,2);
%                         ini_lap = [ini_lap; ini, pos_tmp(pos_tmp(:,1)==ini,2)]; 
%                         fin_lap = [fin_lap; fin, pos_tmp(pos_tmp(:,1)==fin,2)]; 
%                     end 
%                     scatter(ini_lap(:,2),ini_lap(:,1),'r');scatter(fin_lap(:,2),fin_lap(:,1),'b')
%                     title('In Maze');
%                     subplot(1,2,2);hold on;plot(pos_tmp(:,2),pos_tmp(:,1));
%                     ini_lap = [];
%                     fin_lap = []; 
%                     for ix=1:size(in_lapA,1)
%                         ini = in_lapA(ix,1);
%                         fin = in_lapA(ix,2);
%                         ini_lap = [ini_lap; ini, pos_tmp(pos_tmp(:,1)==ini,2)]; 
%                         fin_lap = [fin_lap; fin, pos_tmp(pos_tmp(:,1)==fin,2)]; 
%                     end 
%                     scatter(ini_lap(:,2),ini_lap(:,1),'r');scatter(fin_lap(:,2),fin_lap(:,1),'b')
%                     title('In Lap');
                    %%%%%%%%
                    
                    %Restrict spk and position to laps and movement periods
                    spks_tmp = Restrict(spks_tmp,in_lapA);spks_tmp = Restrict(spks_tmp , movement.aversive); % Restrict spikes to movement periods
                    pos_tmp = Restrict(pos_tmp,in_lapA);pos_tmp = Restrict(pos_tmp, movement.aversive); % Restrict pos to movement periods

                    pos_tmp(:,2) = pos_tmp(:,2)-min(pos_tmp(:,2)); pos_tmp(:,2) = pos_tmp(:,2)/max(pos_tmp(:,2)); %normalization of position
                    
                    %Store par for between comp
                    spks_ave = spks_tmp;
                    pos_ave = pos_tmp;
                    
                    %%%%%%%%Control plot %%%%%%%%%%
%                     figure(3);clf; hold on; plot(pos_tmp(:,2),pos_tmp(:,1),'color', [.3 .3 .3]);hold on; 
%                     %Find interpolated position of each spike:
%                     xs = interp1(pos_tmp(:,1),pos_tmp(:,2),spks_tmp);
%                     scatter(xs,spks_tmp,'*'); 
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    

                    %Firing curve construction
                    [curveA , statsA] = FiringCurve(pos_tmp , spks_tmp , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    %%%%%%%% Control plot %%%%%%%%
%                     figure(1);clf;hold on; 
%                     subplot(2,1,1);imagesc(curveA.rate), colormap 'jet'
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    %%%Within-trial pc parameters or Lap parameters
%                     if contains(comp,'whitin') 
%                         [withinA,withinA_tresh] = Within_pc(pos_tmp,spks_tmp,1,sigma,Xedges);
%                     elseif contains(comp,'laps')
                        %Compare laps 
%                         [withinA_lap] = Within_lap(pos_tmp,spks_tmp,in_lapA,sigma,Xedges);
%                         [withinA_lap_r] = Within_lap_random(pos_tmp,spks_tmp,in_lapA,sigma,Xedges);
%                     else
%                         disp('No comparison method selected')
%                     end 
                    
                    
                     % --- Reward ---
                    spks_tmp = spks; 
                    pos_tmp = behavior.pos.reward(:,1:2);
                    in_maze = ToIntervals(pos_tmp(:,1),and(pos_tmp(:,2)>30 , pos_tmp(:,2)<180));% eliminating the extrems of the maze
                    %Check lenght of in_maze segments to define real laps
                    in_lapR = []; 
                    for ix=1:size(in_maze,1)
                        ini = in_maze(ix,1);
                        fin = in_maze(ix,2);
                        ini_pos = pos_tmp(pos_tmp(:,1)==ini,2); 
                        fin_pos = pos_tmp(pos_tmp(:,1)==fin,2); 
                        if abs(fin_pos-ini_pos)>100
                             in_lapR = [in_lapR;in_maze(ix,:)]; 
                        end
                    end 
                    
%                     %%%%%%%%Control plot%%%%%%%%
%                     figure(3);clf; hold on; subplot(1,2,1);plot(pos_tmp(:,2),pos_tmp(:,1));hold on; 
%                     ini_lap = [];
%                     fin_lap = []; 
%                     for ix=1:size(in_maze,1)
%                         ini = in_maze(ix,1);
%                         fin = in_maze(ix,2);
%                         ini_lap = [ini_lap; ini, pos_tmp(pos_tmp(:,1)==ini,2)]; 
%                         fin_lap = [fin_lap; fin, pos_tmp(pos_tmp(:,1)==fin,2)]; 
%                     end 
%                     scatter(ini_lap(:,2),ini_lap(:,1),'r');scatter(fin_lap(:,2),fin_lap(:,1),'b')
%                     title('In Maze');
%                     subplot(1,2,2);hold on;plot(pos_tmp(:,2),pos_tmp(:,1));
%                     ini_lap = [];
%                     fin_lap = []; 
%                     for ix=1:size(in_lapR,1)
%                         ini = in_lapR(ix,1);
%                         fin = in_lapR(ix,2);
%                         ini_lap = [ini_lap; ini, pos_tmp(pos_tmp(:,1)==ini,2)]; 
%                         fin_lap = [fin_lap; fin, pos_tmp(pos_tmp(:,1)==fin,2)]; 
%                     end 
%                     scatter(ini_lap(:,2),ini_lap(:,1),'r');scatter(fin_lap(:,2),fin_lap(:,1),'b')
%                     title('In Lap');
%                     %%%%%%%%%%%%%%%%%

                    %Restrict spks and pos to valid laps 
                    spks_tmp = Restrict(spks_tmp,in_lapR); spks_tmp = Restrict(spks_tmp, movement.reward); % Restrict to movement periods
                    pos_tmp = Restrict(pos_tmp,in_lapR); pos_tmp = Restrict(pos_tmp , movement.reward);
                    pos_tmp(:,2) = pos_tmp(:,2)-min(pos_tmp(:,2)); pos_tmp(:,2) = pos_tmp(:,2)/max(pos_tmp(:,2)); %normalization of position
                    %Store par for between comp
                    spks_rew = spks_tmp;
                    pos_rew = pos_tmp;
                    
                    %%%%%%%Control plot %%%%%%%%%%
%                     figure(4);clf; hold on; plot(pos_tmp(:,2),pos_tmp(:,1),'color', [.3 .3 .3]);hold on; 
%                     %Find interpolated position of each spike:
%                     xs = interp1(pos_tmp(:,1),pos_tmp(:,2),spks_tmp);
%                     scatter(xs,spks_tmp,'*'); 
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    %Firing curve construction
                    [curveR , statsR] = FiringCurve(pos_tmp , spks_tmp , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4 , 'minPeak' , 0.2);

%                    subplot(2,1,2);imagesc(curveR.rate), colormap 'jet'
%                    sgtitle('Firing curvs original');
                    
                    %%%Within-trial pc parameters or Lap parameters
%                     if contains(comp,'whitin') 
%                         [withinR,withinR_tresh] = Within_pc(pos_tmp,spks_tmp,1,sigma,Xedges);
%                     elseif contains(comp,'laps')
                        %Compare laps 
%                         [withinR_lap] = Within_lap(pos_tmp,spks_tmp,in_lapR,sigma,Xedges);
%                         [withinR_lap_r] = Within_lap_random(pos_tmp,spks_tmp,in_lapR,sigma,Xedges);
%                     else
%                         disp('No comparison method selected')
%                     end 
                    
                    %Subsamplin pos and spk 
%                     [sub_pos_ave,sub_spk_ave,sub_pos_rew,sub_spk_rew] = Subsampling_1d(pos_ave, spks_ave,pos_rew, spks_rew,dt,Xedges,sigma); 
                    
                    %PC cirteria: have a PF of at least 4 bins 
                    %if there is no PF, assigne 0 to correct assesment in
                    %the next if
                    if isempty(statsA.field)
                        statsA.field=0;
                    end 
                     if isempty(statsR.field)
                        statsR.field=0;
                    end 
                    
                    if or(sum(statsA.field(:,:,1))>= 4 , sum(statsR.field(:,:,1))>= 4) 
                        
                        %Subsampling
%                         [between_sub] = Between_pc(sub_pos_ave,sub_spk_ave,sub_pos_rew,sub_spk_rew,bin_size,sigma,Xedges);
%                         [withinA_sub,~] = Within_pc(sub_pos_ave,sub_spk_ave,1,sigma,Xedges);
%                         [withinR_sub,~] = Within_pc(sub_pos_rew,sub_spk_rew,1,sigma,Xedges);

                        %%%Between pc parameters or Between Lap parameters
%                         if contains(comp,'whitin') 
%                             [between] = Between_pc(pos_ave,spks_ave,pos_rew,spks_rew,bin_size,sigma,Xedges);
%                         elseif contains(comp,'laps')
                        %Compare laps 
%                             [between_lap] = Between_lap(pos_ave,spks_ave,in_lapA,pos_rew,spks_rew,in_lapR,sigma,Xedges);
%                             [between_lap_r] = Between_lap_random(pos_ave,spks_ave,in_lapA,pos_rew,spks_rew,in_lapR,sigma,Xedges);
%                         else
%                             disp('No comparison method selected')
%                         end 

                        %Save pc info 
%                         n.ave = withinA_lap; 
%                         n.rew = withinR_lap; 
%                         n.between = between_lap; 
%                         n.ave_r = withinA_lap_r; 
%                         n.rew_r = withinR_lap_r ; 
%                         n.between_r = between_lap_r;
%                         vHPC{count}.lap = n; 
                        vHPC{count}.nlap_rew= size(in_lapR,1);
                        vHPC{count}.nlap_ave= size(in_lapA,1);
                        
                        count = count +1;
                        
                        clear fr_ave fr_rew map_ave map_rew
                    end
                    clear s q curve OccMap Nspikes spks_tmp pos_tmp m curve stats sp curve1
                    clear curveA curveR qA qR curveSA curveSR statsSA statsSR
                end
            end
            clear celltype tmp b cluster
            
    end
        

    %% Saveing PC INFO 
    if ~isempty(dHPC)
        dHPC = dHPC(~cellfun('isempty',dHPC));
        save([cd,'\dHPC_pc.mat'],'dHPC'); 
    end
    if ~isempty(vHPC)
        vHPC = vHPC(~cellfun('isempty',vHPC));
        save([cd,'\vHPC_pc.mat'],'vHPC'); 
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

%% MAIN LOOP LAPS RANDOM spliting by random laps subsampling laps  
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
        
        %Session output
        dHPC_lap = {};% one cell per pc
        vHPC_lap = {}; 
        
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
        % periods of movment during each condition
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
        %Video sampling rate 
        dt = (mean(diff(behavior.pos.aversive(:,1)))); 
        
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
        
        for ii=1:size(group_dHPC,1)

            cluster = group_dHPC(ii,1);
            celltype = logical(cellulartype(cellulartype(:,1) == cluster,2));
            
            if celltype % check if pyr
                
                spks = spks_dHPC(spks_dHPC(:,1)==cluster,2); % select tspk from cluster

                % --- Aversive ---
                spks_ave = spks; 
                pos_ave = behavior.pos.aversive(:,1:2); 
                in_maze = ToIntervals(pos_ave(:,1),and(pos_ave(:,2)>30 , pos_ave(:,2)<180));% eliminating the extrems of the maze
                    
                %Check lenght of in_maze segments to define real laps
                in_lapA = []; % keeps only full laps (maybe too restrictive)
                    for ix=1:size(in_maze,1)
                        ini = in_maze(ix,1);
                        fin = in_maze(ix,2);
                        ini_pos = pos_ave(pos_ave(:,1)==ini,2); 
                        fin_pos = pos_ave(pos_ave(:,1)==fin,2); 
                        if abs(fin_pos-ini_pos)>100
                             in_lapA = [in_lapA;in_maze(ix,:)]; 
                        end
                    end
                    
                %Restrict spk and position to laps and movement periods
                spks_ave = Restrict(spks_ave,in_lapA);spks_ave = Restrict(spks_ave , movement.aversive); % Restrict spikes to movement periods and laps
                    
                pos_ave = Restrict(pos_ave,in_lapA);pos_ave = Restrict(pos_ave, movement.aversive); % Restrict pos to movement periods and laps
                pos_ave(:,2) = pos_ave(:,2)-min(pos_ave(:,2)); pos_ave(:,2) = pos_ave(:,2)/max(pos_ave(:,2)); %normalization of position
               
                % --- Reward ---
                spks_rew = spks; 
                pos_rew = behavior.pos.reward(:,1:2);
                in_maze = ToIntervals(pos_rew(:,1),and(pos_rew(:,2)>30 , pos_rew(:,2)<180));% eliminating the extrems of the maze
               
                %Check lenght of in_maze segments to define real laps
                in_lapR = []; 
                    for ix=1:size(in_maze,1)
                        ini = in_maze(ix,1);
                        fin = in_maze(ix,2);
                        ini_pos = pos_rew(pos_rew(:,1)==ini,2); 
                        fin_pos = pos_rew(pos_rew(:,1)==fin,2); 
                        if abs(fin_pos-ini_pos)>100
                             in_lapR = [in_lapR;in_maze(ix,:)]; 
                        end
                    end 
                    
                %Restrict spk and position to laps and movement periods
                spks_rew = Restrict(spks_rew,in_lapR); spks_rew = Restrict(spks_rew, movement.reward); % Restrict to movement periods
                pos_rew = Restrict(pos_rew,in_lapR); pos_rew = Restrict(pos_rew , movement.reward);
                pos_rew(:,2) = pos_rew(:,2)-min(pos_rew(:,2)); pos_rew(:,2) = pos_rew(:,2)/max(pos_rew(:,2)); %normalization of position 
                    
                    %%%%%%%%Control plot %%%%%%%%%%
%                     figure(1);clf; hold on;subplot(1,2,1); plot(pos_ave(:,2),pos_ave(:,1),'color', [.3 .3 .3]);hold on; 
%                     Find interpolated position of each spike:
%                     xs = interp1(pos_ave(:,1),pos_ave(:,2),spks_ave);scatter(xs,spks_ave,'*'); 
%                     title(['Aversive - Mean fr:',num2str(fr_ave)]);
%                    
%                     subplot(1,2,2); plot(pos_rew(:,2),pos_rew(:,1),'color', [.3 .3 .3]);hold on; 
%                     Find interpolated position of each spike:
%                     xs = interp1(pos_rew(:,1),pos_rew(:,2),spks_rew); scatter(xs,spks_rew,'*');
%                     title(['Reward- Mean fr:',num2str(fr_rew)]);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
   
                %Checking pc criteria
                % 1. Mean firing rate above tresh; 2. At least 1 pf > 4 bins in one of the cond  
     
                [curveA , statsA] = FiringCurve(pos_ave , spks_ave , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                [curveR , statsR] = FiringCurve(pos_rew , spks_rew , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    
                fr_ave= size(spks_ave,1)/sum(in_lapA(:,2)-in_lapA(:,1)); 
                fr_rew= size(spks_rew,1)/sum(in_lapR(:,2)-in_lapR(:,1)); 
                tresh=0.1; % mean firng rate treshold 
                    
                %if there is no PF, assigne 0 to correct assesment in the next if
                if isempty(statsA.field)
                      statsA.field=0;
                end 
                if isempty(statsR.field)
                      statsR.field=0;
                end 
                    
                    if and(or(sum(statsA.field(:,:,1))>= 4 , sum(statsR.field(:,:,1))>= 4),or(fr_ave>tresh,fr_rew>tresh)) % you are a PC  
                    
                    %Subsampling laps - keep the same number of laps in each condition. 
                    % Selects the n_lap of the cond with less laps, and randomly select *100 times n_lap from the cond with more laps 
                    % and caluclate within/between comparisons 
                    
                        if size(in_lapR,1)> size(in_lapA,1) % if more laps in rew
                            n_lap = size(in_lapA,1); 
                            [within_mean1,within_mean2,between] = pc_parameters_laps(n_lap,pos_rew,spks_rew,in_lapR,pos_ave,spks_ave,in_lapA,sigma,Xedges);
                            withinA = within_mean2;
                            withinR = within_mean1;                     
                        elseif size(in_lapR,1)< size(in_lapA,1)  %if more laps in ave
                            n_lap = size(in_lapR,1); 
                            [within_mean1,within_mean2,between] = pc_parameters_laps(n_lap,pos_ave,spks_ave,in_lapA,pos_rew,spks_rew,in_lapR,sigma,Xedges);                            
                            withinA = within_mean1;
                            withinR = within_mean2;
                        else % equal # laps in cond 
                            [withinA] = Within_lap_random(pos_ave,spks_ave,in_lapA,sigma,Xedges);
                            [withinR] = Within_lap_random(pos_rew,spks_rew,in_lapR,sigma,Xedges);
                            [between] = Between_lap_random(pos_ave,spks_ave,in_lapA,pos_rew,spks_rew,in_lapR,sigma,Xedges);
                        end 
                    
                    % Save PC varaiables 
                       
                    n.id = cluster; 
                    n.n_lap =n_lap; 
                    n.fr_mean.ave= fr_ave;
                    n.fr_mean.rew= fr_rew;
                    n.frMap_ave = curveA.rate;
                    n.frMap_rew = curveR.rate;
                    n.stats_ave = statsA;
                    n.stats_rew = statsR;
                    n.between = between;
                    n.within.ave = withinA; 
                    n.within.rew = withinR;
       
                    dHPC_lap{ii}= n;
                    end 
            end
        end

        disp('vHPC Firing rate map calculation')
        for ii=1:size(group_vHPC,1)

            cluster = group_vHPC(ii,1);
            celltype = logical(cellulartype(cellulartype(:,1) == cluster,2));
            
            if celltype % check if pyr
                
                spks = spks_vHPC(spks_vHPC(:,1)==cluster,2); % select tspk from cluster

                % --- Aversive ---
                spks_ave = spks; 
                pos_ave = behavior.pos.aversive(:,1:2); 
                in_maze = ToIntervals(pos_ave(:,1),and(pos_ave(:,2)>30 , pos_ave(:,2)<180));% eliminating the extrems of the maze
                    
                %Check lenght of in_maze segments to define real laps
                in_lapA = []; % keeps only full laps (maybe too restrictive)
                    for ix=1:size(in_maze,1)
                        ini = in_maze(ix,1);
                        fin = in_maze(ix,2);
                        ini_pos = pos_ave(pos_ave(:,1)==ini,2); 
                        fin_pos = pos_ave(pos_ave(:,1)==fin,2); 
                        if abs(fin_pos-ini_pos)>100
                             in_lapA = [in_lapA;in_maze(ix,:)]; 
                        end
                    end
                    
                %Restrict spk and position to laps and movement periods
                spks_ave = Restrict(spks_ave,in_lapA);spks_ave = Restrict(spks_ave , movement.aversive); % Restrict spikes to movement periods and laps
                    
                pos_ave = Restrict(pos_ave,in_lapA);pos_ave = Restrict(pos_ave, movement.aversive); % Restrict pos to movement periods and laps
                pos_ave(:,2) = pos_ave(:,2)-min(pos_ave(:,2)); pos_ave(:,2) = pos_ave(:,2)/max(pos_ave(:,2)); %normalization of position
               
                % --- Reward ---
                spks_rew = spks; 
                pos_rew = behavior.pos.reward(:,1:2);
                in_maze = ToIntervals(pos_rew(:,1),and(pos_rew(:,2)>30 , pos_rew(:,2)<180));% eliminating the extrems of the maze
               
                %Check lenght of in_maze segments to define real laps
                in_lapR = []; 
                    for ix=1:size(in_maze,1)
                        ini = in_maze(ix,1);
                        fin = in_maze(ix,2);
                        ini_pos = pos_rew(pos_rew(:,1)==ini,2); 
                        fin_pos = pos_rew(pos_rew(:,1)==fin,2); 
                        if abs(fin_pos-ini_pos)>100
                             in_lapR = [in_lapR;in_maze(ix,:)]; 
                        end
                    end 
                    
                %Restrict spk and position to laps and movement periods
                spks_rew = Restrict(spks_rew,in_lapR); spks_rew = Restrict(spks_rew, movement.reward); % Restrict to movement periods
                pos_rew = Restrict(pos_rew,in_lapR); pos_rew = Restrict(pos_rew , movement.reward);
                pos_rew(:,2) = pos_rew(:,2)-min(pos_rew(:,2)); pos_rew(:,2) = pos_rew(:,2)/max(pos_rew(:,2)); %normalization of position 
                    
                    %%%%%%%%Control plot %%%%%%%%%%
%                     figure(1);clf; hold on;subplot(1,2,1); plot(pos_ave(:,2),pos_ave(:,1),'color', [.3 .3 .3]);hold on; 
%                     %Find interpolated position of each spike:
%                     xs = interp1(pos_ave(:,1),pos_ave(:,2),spks_ave);scatter(xs,spks_ave,'*'); 
%                     title('Aversive');
%                    
%                     subplot(1,2,2); plot(pos_rew(:,2),pos_rew(:,1),'color', [.3 .3 .3]);hold on; 
%                     %Find interpolated position of each spike:
%                     xs = interp1(pos_rew(:,1),pos_rew(:,2),spks_rew); scatter(xs,spks_rew,'*');
%                     title('Reward');
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
   
                %Checking pc criteria
                % 1. Mean firing rate above tresh; 2. At least 1 pf > 4 bins in one of the cond  
     
                [curveA , statsA] = FiringCurve(pos_ave , spks_ave , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                [curveR , statsR] = FiringCurve(pos_rew , spks_rew , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    
                fr_ave= size(spks_ave,1)/sum(in_lapA(:,2)-in_lapA(:,1)); 
                fr_rew= size(spks_rew,1)/sum(in_lapR(:,2)-in_lapR(:,1)); 
                tresh=0.1; % mean firng rate treshold 
                    
                %if there is no PF, assigne 0 to correct assesment in the next if
                if isempty(statsA.field)
                      statsA.field=0;
                end 
                if isempty(statsR.field)
                      statsR.field=0;
                end 
                    
                    if and(or(sum(statsA.field(:,:,1))>= 4 , sum(statsR.field(:,:,1))>= 4),or(fr_ave>tresh,fr_rew>tresh)) % you are a PC  
                    
                    %Subsampling laps - keep the same number of laps in each condition. 
                    % Selects the n_lap of the cond with less laps, and randomly select *100 times n_lap from the cond with more laps 
                    % and caluclate within/between comparisons 
                    
                        if size(in_lapR,1)> size(in_lapA,1) % if more laps in rew
                            n_lap = size(in_lapA,1); 
                            [within_mean1,within_mean2,between] = pc_parameters_laps(n_lap,pos_rew,spks_rew,in_lapR,pos_ave,spks_ave,in_lapA,sigma,Xedges);
                            withinA = within_mean2;
                            withinR = within_mean1;                     
                        elseif size(in_lapR,1)< size(in_lapA,1)  %if more laps in ave
                            n_lap = size(in_lapR,1); 
                            [within_mean1,within_mean2,between] = pc_parameters_laps(n_lap,pos_ave,spks_ave,in_lapA,pos_rew,spks_rew,in_lapR,sigma,Xedges);                            
                            withinA = within_mean1;
                            withinR = within_mean2;
                        else % equal # laps in cond 
                            [withinA] = Within_lap_random(pos_ave,spks_ave,in_lapA,sigma,Xedges);
                            [withinR] = Within_lap_random(pos_rew,spks_rew,in_lapR,sigma,Xedges);
                            [between] = Between_lap_random(pos_ave,spks_ave,in_lapA,pos_rew,spks_rew,in_lapR,sigma,Xedges);
                        end 
                    
                    % Save PC varaiables 
                       
                    n.id = cluster; 
                    n.n_lap =n_lap; 
                    n.frMap_ave = curveA.rate;
                    n.frMap_rew = curveR.rate;
                    n.stats_ave = statsA;
                    n.stats_rew = statsR;
                    n.between = between;
                    n.within.ave = withinA; 
                    n.within.rew = withinR;
       
                    vHPC_lap{ii}= n;
                    end 
            end
        end     

    %% Saveing PC INFO 
    if ~isempty(dHPC_lap)
        dHPC_lap = dHPC_lap(~cellfun('isempty',dHPC_lap));
        save([cd,'\dHPC_pc_lap.mat'],'dHPC_lap'); 
    end
    if ~isempty(vHPC_lap)
        vHPC_lap = vHPC_lap(~cellfun('isempty',vHPC_lap));
        save([cd,'\vHPC_pc_lap.mat'],'vHPC_lap'); 
    end
  

        
        clear aversiveTS aversiveTS_run rewardTS rewardTS_run baselineTS
        clear behavior camara Cell_type_classification cellulartype cluster
        clear group_dHPC group_vHPC K Kinfo leftvalve rightvalve movement R
        clear Rewards_filt Shocks_filt spks_dHPC spks_vHPC r rZ
        clear cluster cluster_dHPC cluster_vHPC coordinated coordinatedV coordinatedV_refined
        clear camaraA count dX dX_int dY dY_int gc grps i ii MUA segments tmp WAKE REM NREM
        disp(['-- Finished folder #' , num2str(t) , ' from rat #' , num2str(tt) , ' ---'])
        disp('  ')
    end
    disp(['-------- Finished rat#' , num2str(tt) , ' --------'])
    disp(' ')
end

%% REMAPPING PARAMETERS PLOTS and STATS -  LAPS
% Load remapping parameters and creates matrix to plot and stats.
% Prev to run this section, you should have run the main loop. So, this
% section assumes that inside each session folder, you have a dHPC_pc and
% vHPC_pc file.

%c1 = spatial c2= fr_change c3=  overlap c4: pf shift c5:  1 = between 2 = within aversive 3= within reward 
clear
clc
close all
% Parameters
path = {'\\Maryjackson\e\Rat127\Ephys\pyr';'\\Maryjackson\e\Rat128\Ephys\in_pyr\ready';'\\Maryjackson\e\Rat103\usable';'\\Maryjackson\e\Rat132\recordings\in_pyr'; '\\Maryjackson\e\Rat165\in_pyr'; ...
    '\\Maryjackson\e\Rat126\Ephys\in_Pyr'};%List of folders from the path

%Output matrix 
dhpc_sub = [];
vhpc_sub = [];

for tt = 1:length(path)
    %List of folders from the path
    files = dir(path{tt});
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
    clear files dirFlags
    dhpc_temp = [];
    vhpc_temp = [];
    for t = 1 : length(subFolders)-2
        disp(['-- Initiating analysis of folder #' , num2str(t) , ' from rat #',num2str(tt) , ' --'])
        session = [subFolders(t+2).folder,'\',subFolders(t+2).name];
        cd([session,'\Spikesorting'])
        %Loading pc matrix of the sessions
        disp('Uploading session pc matrix');
        try
            load('dHPC_pc.mat');
             %Save session parameters 
        for d=1:size(dHPC,2)
            pc = dHPC{d}; 
            bet = [pc.subsampled.between,1]; 
            wa = [pc.subsampled.within_ave,2];
            wr = [pc.subsampled.within_rew,3];
            dhpc_temp = [dhpc_temp;bet;wa;wr]; 
            
        end
        catch 
            dHPC = []; disp(['No dHPC_pc_lap.mat file in ',session]);
        end 
        
        try
            load('vHPC_pc.mat');
            for d=1:size(vHPC,2)
                pc = vHPC{d}; 
                bet = [pc.subsampled.between,1]; 
                wa = [pc.subsampled.within_ave,2];
                wr = [pc.subsampled.within_rew,3];
                vhpc_temp = [vhpc_temp;bet;wa;wr]; 
            end 
     
        catch
            vHPC= []; disp(['No vHPC_pc_lap.mat file in ',session]);
        end
        
    end 
    %Save sub parameters 
    dhpc_sub = [dhpc_sub;dhpc_temp];
    vhpc_sub = [vhpc_sub;vhpc_temp];
     
end 

% Plots + stats 
%dHP
ylabels ={'Spatial correlation', 'Fr change', 'Overlap', 'Pf shift', 'Mean fr ratio'}; 
xlabels = {'Between', 'Within Ave', 'Within Rew'}; 
figure(1);clf;hold on, 
sgtitle('dHPC')
for p=1:size(ylabels,2)
    subplot(1,5,p); hold on; 
    x = dhpc_sub(:,6);
    y= dhpc_sub(:,p);%Change the column number (1-4) to choose which variable to plot 
    c = [.3,.3,.3];
    scatter(x,y,[],c,"filled",'jitter','on', 'jitterAmount',0.1); xlim([0 4]);
    ylabel(ylabels{p});
    xticks([1 2 3])
    xticklabels({'Bet', 'WA', 'WR'});
    hold on
    s1=scatter(1,nanmedian(dhpc_sub(dhpc_sub(:,6)==1,p)), "filled");s1.MarkerFaceColor = [0 0 0];
    s2=scatter(2,nanmedian(dhpc_sub(dhpc_sub(:,6)==2,p)), "filled");s2.MarkerFaceColor = [1 0.1 0.2];
    s3=scatter(3,nanmedian(dhpc_sub(dhpc_sub(:,6)==3,p)), "filled");s3.MarkerFaceColor = [0.1 0.3 1];
end 


%Only 3 plots
ylabels ={'Spatial correlation', 'Overlap', 'Pf shift'};
idx = [1,3,4]; 
xlabels = {'Between', 'Within Ave', 'Within Rew'}; 
figure(1);clf;hold on, 
sgtitle('dHPC')
for a=1:size(ylabels,2)
    p = idx(a);
    subplot(1,3,a); hold on; 
    x = dhpc_sub(:,6);
    y= dhpc_sub(:,p);%Change the column number (1-4) to choose which variable to plot 
    c = [.3,.3,.3];
    scatter(x,y,30,c,"filled",'jitter','on', 'jitterAmount',0.3); xlim([0 4]);
    ylabel(ylabels{a});
    xticks([1 2 3])
    xticklabels({'Bet', 'WA', 'WR'});
    hold on
    s1=scatter(1,nanmedian(dhpc_sub(dhpc_sub(:,6)==1,p)), "filled");s1.MarkerFaceColor = [0 0 0];
    s2=scatter(2,nanmedian(dhpc_sub(dhpc_sub(:,6)==2,p)), "filled");s2.MarkerFaceColor = [1 0.1 0.2];
    s3=scatter(3,nanmedian(dhpc_sub(dhpc_sub(:,6)==3,p)), "filled");s3.MarkerFaceColor = [0.1 0.3 1];
end

%Stats
[P,ANOVATAB,STATS] = kruskalwallis(dhpc_sub(:,5),dhpc_sub(:,6));
c = multcompare(STATS)

tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])

% Paired non-parametric anova - Friedman test
%Prepare data to test
p = 4; % 1=sp corr 2=fr change 3=overlap 4=pf shift
data=[dhpc_sub(dhpc_sub(:,6)==1,p), dhpc_sub(dhpc_sub(:,6)==2,p), ...
    dhpc_sub(dhpc_sub(:,6)==3,p), dhpc_sub(dhpc_sub(:,6)==4,p),dhpc_sub(dhpc_sub(:,6)==5,p)];
data(any(isnan(data), 2), :) = [];

[p,~,stats] =  friedman(data,1); 
c = multcompare(stats);
tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])


%vHP
ylabels ={'Spatial correlation', 'Fr change', 'Overlap', 'Pf shift', 'Mean fr ratio'}; 
xlabels = {'Between', 'Within Ave', 'Within Rew'}; 
figure(2);clf;hold on, 
sgtitle('vHPC')
for p=1:size(ylabels,2)
    subplot(1,5,p); hold on; 
    x = vhpc_sub(:,6);
    y=  vhpc_sub(:,p);%Change the column number (1-4) to choose which variable to plot 
    c = [.3,.3,.3];
    scatter(x,y,[],c,"filled",'jitter','on', 'jitterAmount',0.1); xlim([0 4]);
    ylabel(ylabels{p});
    xticks([1 2 3])
    xticklabels({'Bet', 'WA', 'WR'});
    hold on
    s1=scatter(1,nanmedian(vhpc_sub(vhpc_sub(:,6)==1,p)), "filled");s1.MarkerFaceColor = [0 0 0];
    s2=scatter(2,nanmedian(vhpc_sub(vhpc_sub(:,6)==2,p)), "filled");s2.MarkerFaceColor = [1 0.1 0.2];
    s3=scatter(3,nanmedian(vhpc_sub(vhpc_sub(:,6)==3,p)), "filled");s3.MarkerFaceColor = [0.1 0.3 1];
end 


%Only 3 plots
ylabels ={'Spatial correlation', 'Overlap', 'Pf shift'};
idx = [1,3,4]; 
xlabels = {'Between', 'Within Ave', 'Within Rew'}; 
figure(1);clf;hold on, 
sgtitle('vHPC')
for a=1:size(ylabels,2)
    p = idx(a);
    subplot(1,3,a); hold on; 
    x = vhpc_sub(:,6);
    y= vhpc_sub(:,p);%Change the column number (1-4) to choose which variable to plot 
    c = [.3,.3,.3];
    scatter(x,y,30,c,"filled",'jitter','on', 'jitterAmount',0.3); xlim([0 4]);
    ylabel(ylabels{a});
    xticks([1 2 3])
    xticklabels({'Bet', 'WA', 'WR'});
    hold on
    s1=scatter(1,nanmedian(vhpc_sub(vhpc_sub(:,6)==1,p)), "filled");s1.MarkerFaceColor = [0 0 0];
    s2=scatter(2,nanmedian(vhpc_sub(vhpc_sub(:,6)==2,p)), "filled");s2.MarkerFaceColor = [1 0.1 0.2];
    s3=scatter(3,nanmedian(vhpc_sub(vhpc_sub(:,6)==3,p)), "filled");s3.MarkerFaceColor = [0.1 0.3 1];
end



%Stats
[P,ANOVATAB,STATS] = kruskalwallis(vhpc_sub(:,1),vhpc_sub(:,6));
c = multcompare(STATS)

tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])

% Paired non-parametric anova - Friedman test
%Prepare data to test
p = 4; % 1=sp corr 2=fr change 3=overlap 4=pf shift
data=[vhpc_sub(vhpc_sub(:,6)==1,p), vhpc_sub(vhpc_sub(:,6)==2,p), ...
    vhpc_sub(vhpc_sub(:,6)==3,p), vhpc_sub(vhpc_sub(:,6)==4,p),vhpc_sub(vhpc_sub(:,6)==5,p)];
data(any(isnan(data), 2), :) = [];

[p,~,stats] =  friedman(data,1); 
c = multcompare(stats)
tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])


%Stats
[P,ANOVATAB,STATS] = kruskalwallis(vhpc_sub(:,1),vhpc_sub(:,6));
c = multcompare(STATS)

tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])

%% REMAPPING PARAMETERS PLOTS AND STATS
% Load remapping parameters and creates matrix to plot and stats.
% Prev to run this section, you should have run the main loop. So, this
% section assumes that inside each session folder, you have a dHPC_pc and
% vHPC_pc file.

%c1 = spatial c2= fr_change c3=  overlap c4: pf shift c5:  1 = between 2 = within aversive 3= within reward 
clear
clc
close all
% Parameters
path = {'\\Maryjackson\e\Rat127\Ephys\pyr';'\\Maryjackson\e\Rat128\Ephys\in_pyr\ready';'\\Maryjackson\e\Rat103\usable';'\\Maryjackson\e\Rat132\recordings\in_pyr'; '\\Maryjackson\e\Rat165\in_pyr'; ...
    '\\Maryjackson\e\Rat126\Ephys\in_Pyr'};%List of folders from the path

%Output matrix 
dhpc_sub = [];
vhpc_sub = [];

for tt = 1:length(path)
    %List of folders from the path
    files = dir(path{tt});
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
    clear files dirFlags
    dhpc_temp = [];
    vhpc_temp = [];
    for t = 1 : length(subFolders)-2
        disp(['-- Initiating analysis of folder #' , num2str(t) , ' from rat #',num2str(tt) , ' --'])
        session = [subFolders(t+2).folder,'\',subFolders(t+2).name];
        cd([session,'\Spikesorting'])
        %Loading pc matrix of the sessions
        disp('Uploading session pc matrix');
        try
            load('dHPC_pc.mat');
             %Save session parameters 
            for d=1:size(dHPC,2)
                pc = dHPC{d}; 
                bet = [pc.between,1]; 
                wa = [pc.within_ave,2];
                wr = [pc.within_rew,3];
                dhpc_temp = [dhpc_temp;bet;wa;wr]; 
            
            end
        catch 
            dHPC = []; disp(['No dHPC_pc.mat file in ',session]);
        end 
        
        try
            load('vHPC_pc.mat');
            for d=1:size(vHPC,2)
                pc = vHPC{d}; 
                bet = [pc.between,1]; 
                wa = [pc.within_ave,2];
                wr = [pc.within_rew,3];
                vhpc_temp = [vhpc_temp;bet;wa;wr]; 
            end 
     
        catch
            vHPC= []; disp(['No vHPC_pc.mat file in ',session]);
        end
        
    end 
    %Save sub parameters 
    dhpc_sub = [dhpc_sub;dhpc_temp];
    vhpc_sub = [vhpc_sub;vhpc_temp];
     
end 

% Plots + stats 
%dHP
ylabels ={'Spatial correlation', 'Fr change', 'Overlap', 'Pf shift', 'Mean fr ratio'}; 
xlabels = {'Between', 'Within Ave', 'Within Rew'}; 
figure(1);clf;hold on, 
sgtitle('dHPC')
for p=1:size(ylabels,2)
    subplot(1,5,p); hold on; 
    x = dhpc_sub(:,6);
    y= dhpc_sub(:,p);%Change the column number (1-4) to choose which variable to plot 
    c = [.3,.3,.3];
    scatter(x,y,[],c,"filled",'jitter','on', 'jitterAmount',0.1); xlim([0 4]);
    ylabel(ylabels{p});
    xticks([1 2 3])
    xticklabels({'Bet', 'WA', 'WR'});
    hold on
    s1=scatter(1,nanmedian(dhpc_sub(dhpc_sub(:,6)==1,p)), "filled");s1.MarkerFaceColor = [0 0 0];
    s2=scatter(2,nanmedian(dhpc_sub(dhpc_sub(:,6)==2,p)), "filled");s2.MarkerFaceColor = [1 0.1 0.2];
    s3=scatter(3,nanmedian(dhpc_sub(dhpc_sub(:,6)==3,p)), "filled");s3.MarkerFaceColor = [0.1 0.3 1];
end 

%Stats
[P,ANOVATAB,STATS] = kruskalwallis(dhpc_sub(:,5),dhpc_sub(:,6));
c = multcompare(STATS)

tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])

% Paired non-parametric anova - Friedman test
%Prepare data to test
p = 1; % 1=sp corr 2=fr change 3=overlap 4=pf shift
data=[dhpc_sub(dhpc_sub(:,6)==1,p), dhpc_sub(dhpc_sub(:,6)==2,p), ...
    dhpc_sub(dhpc_sub(:,6)==3,p), dhpc_sub(dhpc_sub(:,6)==4,p),dhpc_sub(dhpc_sub(:,6)==5,p)];
data(any(isnan(data), 2), :) = [];

[p,~,stats] =  friedman(data,1); 
c = multcompare(stats);
tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])


%vHP
ylabels ={'Spatial correlation', 'Fr change', 'Overlap', 'Pf shift', 'Mean fr ratio'}; 
xlabels = {'Between', 'Within Ave', 'Within Rew'}; 
figure(2);clf;hold on, 
sgtitle('vHPC')
for p=1:size(ylabels,2)
    subplot(1,5,p); hold on; 
    x = vhpc_sub(:,6);
    y=  vhpc_sub(:,p);%Change the column number (1-4) to choose which variable to plot 
    c = [.3,.3,.3];
    scatter(x,y,[],c,"filled",'jitter','on', 'jitterAmount',0.1); xlim([0 4]);
    ylabel(ylabels{p});
    xticks([1 2 3])
    xticklabels({'Bet', 'WA', 'WR'});
    hold on
    s1=scatter(1,nanmedian(vhpc_sub(vhpc_sub(:,6)==1,p)), "filled");s1.MarkerFaceColor = [0 0 0];
    s2=scatter(2,nanmedian(vhpc_sub(vhpc_sub(:,6)==2,p)), "filled");s2.MarkerFaceColor = [1 0.1 0.2];
    s3=scatter(3,nanmedian(vhpc_sub(vhpc_sub(:,6)==3,p)), "filled");s3.MarkerFaceColor = [0.1 0.3 1];
end 


%Stats
[P,ANOVATAB,STATS] = kruskalwallis(vhpc_sub(:,5),vhpc_sub(:,6));
c = multcompare(STATS)

tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])

% Paired non-parametric anova - Friedman test
%Prepare data to test
p = 1; % 1=sp corr 2=fr change 3=overlap 4=pf shift
data=[vhpc_sub(vhpc_sub(:,6)==1,p), vhpc_sub(vhpc_sub(:,6)==2,p), ...
    vhpc_sub(vhpc_sub(:,6)==3,p), vhpc_sub(vhpc_sub(:,6)==4,p),vhpc_sub(vhpc_sub(:,6)==5,p)];
data(any(isnan(data), 2), :) = [];

[p,~,stats] =  friedman(data,1); 
c = multcompare(stats)
tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])


%Stats
[P,ANOVATAB,STATS] = kruskalwallis(vhpc_sub(:,1),vhpc_sub(:,6));
c = multcompare(STATS)

tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])
%% VELOCITY AND OCCUPANCY MAPS
%Main loop
occ_ave=[]; 
occ_rew=[]; 

v_ave=[]; 
v_rew=[]; 

for tt = 1:length(path)
    %List of folders from the path
    files = dir(path{tt});
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
    clear files dirFlags
    
    % Output figures
    shock_fig=figure(3);clf;hold on % shock plot
    
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
        % periods of movment during each condition
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
        %Video sampling rate 
        dt = (mean(diff(behavior.pos.aversive(:,1)))); 
        
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
        
        %% Shock plot
        %Find shock position
        pos_ave = behavior.pos.aversive(:,1:2);
        x_shock = interp1(pos_ave(:,1),pos_ave(:,2),Shocks_filt);
        y=ones(size(x_shock,1),1);
        scatter(x_shock,y,'*');
        figure(shock_fig)
       
        %% Velocity and occupancy plot 
        pos_ave(:,2) =pos_ave(:,2)-min(pos_ave(:,2)); pos_ave(:,2) = pos_ave(:,2)/max(pos_ave(:,2)); %normalization of position
        
        pos_rew = behavior.pos.reward(:,1:2);
        pos_rew(:,2) =pos_rew(:,2)-min(pos_rew(:,2)); pos_rew(:,2) = pos_rew(:,2)/max(pos_rew(:,2));
        
        %Occupancy 
        [Na, edges_a,bin_ave] = histcounts(pos_ave(:,2),Xedges); 
        occ_ave=[occ_ave;Na*dt]; 
        
        [Nr, edges_r,bin_rew] = histcounts(pos_rew(:,2),Xedges); 
        occ_rew=[occ_rew;Nr*dt]; 
        
        clear Na Nr
        %Velocity 
        vel_ave = [behavior.speed.aversive(:,2),bin_ave];
        Na=zeros(1,60); 
        for v=1:60
          Na(v)=nanmean(vel_ave(vel_ave(:,2)==v,1),1);
        end
        v_ave=[v_ave;Na]; 
        
        vel_rew = [behavior.speed.reward(:,2),bin_rew];
        Nr=zeros(1,60); 
        for v=1:60
          Nr(v)=nanmean(vel_rew(vel_rew(:,2)==v,1),1);
        end
        v_rew=[v_rew;Nr]; 


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

%Plot
fig_ave=figure(1);clf;hold on

mocc = nanmean(occ_ave,1);
semocc = nansem(occ_ave);
subplot(2,1,1);plot(1:60,Smooth(mocc,1), 'LineWidth',2,'Color',[1 0.2 0.4]);hold on;
ciplot(Smooth(mocc-semocc,1),Smooth(mocc+semocc,1),1:60,'r'); alpha 0.1; hold on;
title('Occupancy');  
mv = nanmean(v_ave,1);
semv = nansem(v_ave);
subplot(2,1,2);plot(1:60,Smooth(mv,1), 'LineWidth',2,'Color',[1 0.4 0.5]);hold on;
ciplot(Smooth(mv-semv,1),Smooth(mv+semv,1),1:60,'r'); alpha 0.1; hold on;
title('Velocity');  
sgtitle('Aversive')

fig_rew=figure(2);clf;hold on 

mocc = nanmean(occ_rew,1);
semocc = nansem(occ_rew);
subplot(2,1,1);plot(1:60,Smooth(mocc,1), 'LineWidth',2,'Color',[0.2 0.2 1]);hold on;
ciplot(Smooth(mocc-semocc,1),Smooth(mocc+semocc,1),1:60,[0.2 0.2 1]); alpha 0.1; hold on;
title('Occupancy');  
mv = nanmean(v_rew,1);
semv = nansem(v_rew);
subplot(2,1,2);plot(1:60,Smooth(mv,1), 'LineWidth',2,'Color',[0.2 0.2 1]);hold on;
ciplot(Smooth(mv-semv,1),Smooth(mv+semv,1),1:60,[0.2 0.2 1]); alpha 0.1; hold on;
title('Velocity');  
sgtitle('Reward')
%% TOTAL NUMBER OF RECORDED NEURONS - next time put this inside the main loop
clear
clc
close all
% Parameters
path = {'\\Maryjackson\e\Rat127\Ephys\pyr';'\\Maryjackson\e\Rat128\Ephys\in_pyr\ready';'\\Maryjackson\e\Rat103\usable';'\\Maryjackson\e\Rat132\recordings\in_pyr'; '\\Maryjackson\e\Rat165\in_pyr'; ...
    '\\Maryjackson\e\Rat126\Ephys\in_Pyr'};%List of folders from the path

% c1: rat c2: session c3: #dHPC neurons c4: dHPC #pyr c5: #vHPCneurons  c6:VHPC#pyr 
n_total = nan(2000,6); 

ind_global = 1; 

for tt = 1:length(path)
    %List of folders from the path
    files = dir(path{tt});
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
    
    for t = 1 : length(subFolders)-2
        session = [subFolders(t+2).folder,'\',subFolders(t+2).name];
        cd(session)
        rat = str2num(files(3).name(4:6)); 
        %Output
        n_total(ind_global,1) = rat; % rat
        n_total(ind_global,2) = str2num(session(end-7:end)); %session
        
        %% Spikes

        cd 'Spikesorting'
        spks = double([readNPY('spike_clusters.npy') readNPY('spike_times.npy')]);
        K = tdfread('cluster_group.tsv'); % import clusters ID and groups (MUA,Noise,Good)
        Kinfo = tdfread('cluster_info.tsv'); % import information of clusters
        K = [K.cluster_id(K.group(:,1) == 'g') , Kinfo.ch(K.group(:,1) == 'g'),]; % defining good clusters

        % Load neuronal classification
        load('Cell_type_classification')
        K = [K , Cell_type_classification(:,7:8)];
        group_dHPC = K(K(:,2) > 63,:);
        group_vHPC = K(K(:,2) <= 63,:);
      
        %Output
        n_total(ind_global,3) = nansum([n_total(ind_global,3),size(group_dHPC,1)]); %# neurons per session dhpc
        n_total(ind_global,4) = nansum([n_total(ind_global,4),sum(group_dHPC(:,3))]); %# pyr per session dhpc
        
        n_total(ind_global,5) = nansum([n_total(ind_global,5),size(group_vHPC,1)]); %# neurons per session vhpc
        n_total(ind_global,6) = nansum([n_total(ind_global,6),sum(group_vHPC(:,3))]); %# pyr per session vhpc 
        
       ind_global = ind_global+1;
          
    end
    disp(['-------- Finished rat#' , num2str(tt) , ' --------'])
    disp(' ')
end

n_total(any(isnan(n_total), 2), :) = [];

% Count - change manually
size(unique(n_total(n_total(:,1)==165,2)))

sum(n_total(n_total(:,1)==132,6))

%% PC TOTAL NUMBER
clear
clc
close all
% Parameters
path = {'\\Maryjackson\e\Rat127\Ephys\pyr';'\\Maryjackson\e\Rat128\Ephys\in_pyr\ready';'\\Maryjackson\e\Rat103\usable';'\\Maryjackson\e\Rat132\recordings\in_pyr'; '\\Maryjackson\e\Rat165\in_pyr'; ...
    '\\Maryjackson\e\Rat126\Ephys\in_Pyr'};%List of folders from the path

%Output matrix 
dhpc_pc = [];
vhpc_pc = [];

for tt = 1:length(path)
    %List of folders from the path
    files = dir(path{tt});
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
    clear files dirFlags
    dhpc_temp = [];
    vhpc_temp = [];
    for t = 1 : length(subFolders)-2
        disp(['-- Initiating analysis of folder #' , num2str(t) , ' from rat #',num2str(tt) , ' --'])
        session = [subFolders(t+2).folder,'\',subFolders(t+2).name];
        cd([session,'\Spikesorting'])
        %Loading pc matrix of the sessions
        disp('Uploading session pc matrix');
        try
            load('dHPC_pc_lap.mat');
        catch 
            dHPC = []; disp(['No dHPC_pc_lap.mat file in ',session]);
        end 
        try
            load('vHPC_pc_la.mat');
        catch
            vHPC= []; disp(['No vHPC_pc.mat file in ',session]);
        end
        
        %Save session sub parameters 
        for d=1:size(dHPC,2)
            pc = dHPC{d}; 
            bet = [pc.lap.between_r,1]; 
            wa = [pc.lap.ave_r,2];
            wr = [pc.lap.rew_r,3];
            dhpc_temp = [dhpc_temp;bet;wa;wr]; 
            
        end

        for d=1:size(vHPC,2)
            pc = vHPC{d}; 
            bet = [pc.lap.between_r,1]; 
            wa = [pc.lap.ave_r,2];
            wr = [pc.lap.rew_r,3];
            vhpc_temp = [vhpc_temp;bet;wa;wr]; 
            
        end
        
    end
    
    %Save sub parameters 
    dhpc_lap = [dhpc_lap;dhpc_temp];
    vhpc_lap = [vhpc_lap;vhpc_temp];
     
end 


%% PC STABILITY - WITHIN 1st AND 2nd 
% 1.Compares remapping parameters between the 1st and the 2nd half of the
% session 
% Prev to run this section, you should have run the main loop. So, this
% section assumes that inside each session folder, you have a dHPC_pc and
% vHPC_pc
% file. 
% 
%Temporal output
all_stability_ave_dHPC = []; 
all_stability_rew_dHPC = []; 

all_stability_ave_vHPC = []; 
all_stability_rew_vHPC = []; 


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
        
        %Session output
        dHPC = {};% one cell per pc
        vHPC = {}; 
        
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
        %Video sampling rate 
        dt = (mean(diff(behavior.pos.aversive(:,1)))); 
        
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

%         
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
            cellulartype = [K(:,1) , ones(length(K),1)];
        end
        % Definition of the emotional transition
        if aversiveTS_run(1)<rewardTS_run(1)
            cond = 1; % 1 if the order was Aversive -> Reward
        else
            cond = 2;% 1 if the order was Reward -> Aversive
        end
        
        %% Firing maps calculation
        disp('dHPC Stability calculation')

        load('dHPC_pc.mat'); 
        
        temp_stability_ave_dHPC = []; 
        temp_stability_rew_dHPC = []; 
        
        for ii=1:size(dHPC,2)
            
            cluster = dHPC{ii}.id;
  
            spks = spks_dHPC(spks_dHPC(:,1)==cluster,2); % select tspk from cluster
              
            % --- Aversive ---
             spks_tmp = Restrict(spks , movement.aversive); % Restrict spikes to movement periods
             pos_tmp = Restrict(behavior.pos.aversive(:,1:2) , movement.aversive); % Restrict pos to movement periods
             in_maze = ToIntervals(pos_tmp(:,1),and(pos_tmp(:,2)>30 , pos_tmp(:,2)<170));% eliminating the extrems of the maze
             spks_tmp = Restrict(spks_tmp,in_maze);
             pos_tmp = Restrict(pos_tmp,in_maze); 
             pos_tmp(:,2) = pos_tmp(:,2)-min(pos_tmp(:,2)); pos_tmp(:,2) = pos_tmp(:,2)/max(pos_tmp(:,2)); %normalization of position
                    
             % 1st half vs 2nd half 
             half_tresh = round((size(pos_tmp,1))/2);
             
             % Dividing spks and pos by half_tresh
             pos_1=pos_tmp(1:half_tresh,:); 
             pos_2= pos_tmp(half_tresh+1:end,:); 
             spks_1 = spks_tmp(spks_tmp<pos_tmp(half_tresh,1));
             spks_2 = spks_tmp(spks_tmp>pos_tmp(half_tresh,1));
             
             %Dividing by time -> not valid if unequal exploration
%              total_time = max(pos_tmp(:,1))-min(pos_tmp(:,1)); 
%              half_time = total_time/2;
             %Divide spikes time in two by time 
%              spks_1 = spks_tmp(spks_tmp<(min(pos_tmp(:,1)+half_time)));
%              pos_1 = pos_tmp(pos_tmp(:,1)<(min(pos_tmp(:,1)+half_time)),:); 
%              spks_2 = spks_tmp(spks_tmp>=(min(pos_tmp(:,1)+half_time)));
%              pos_2 =  pos_tmp(pos_tmp(:,1)>(min(pos_tmp(:,1)+half_time)),:); 
    
            %Calculate remapping parameters 
            [curve1,stats1] = FiringCurve(pos_1, spks_1 , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 2 , 'minPeak' , 0.1);
            [curve2,stats2] = FiringCurve(pos_2, spks_2, 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 2 , 'minPeak' , 0.1);
                    
            fr_1= nanmean(curve1.rate);
            fr_2= nanmean(curve2.rate);
                    
            %Fr change
            fr_change = abs((fr_1 - fr_2)/(fr_1 + fr_2));

            %Rate overlap
            if fr_1<=fr_2 
                overlap = fr_1/fr_2;
            else 
                overlap = fr_2/fr_1;
            end
    
            %Spatial  corr
            s = corrcoef(curve1.rate, curve2.rate);
            spatial = s(1,2);
    
            %Peak shift 
            shift = abs(stats1.x(1) - stats2.x(1));
            
            within1vs2 = [spatial, fr_change, overlap, shift]; 
            
            dHPC{ii}.within1vs2.ave =  within1vs2; 
            temp_stability_ave_dHPC = [temp_stability_ave_dHPC;within1vs2]; 
        
            
            % --- Reward ---
            spks_tmp = Restrict(spks , movement.reward); % Restrict to movement periods
            pos_tmp = Restrict(behavior.pos.reward(:,1:2), movement.reward);
            in_maze = ToIntervals(pos_tmp(:,1),and(pos_tmp(:,2)>30 , pos_tmp(:,2)<170));% eliminating the extrems of the maze
            spks_tmp = Restrict(spks_tmp,in_maze);
            pos_tmp = Restrict(pos_tmp,in_maze); 
            pos_tmp(:,2) = pos_tmp(:,2)-min(pos_tmp(:,2)); pos_tmp(:,2) = pos_tmp(:,2)/max(pos_tmp(:,2)); %normalization of position
                    
%             % 1st half vs 2nd half 
%             total_time = max(pos_tmp(:,1))-min(pos_tmp(:,1)); 
%             half_time = total_time/2;
%             %Divide spikes time in two
%             spks_1 = spks_tmp(spks_tmp<(min(pos_tmp(:,1)+half_time)));
%             pos_1 = pos_tmp(pos_tmp(:,1)<(min(pos_tmp(:,1)+half_time)),:); 
%             spks_2 = spks_tmp(spks_tmp>=(min(pos_tmp(:,1)+half_time)));
%             pos_2 =  pos_tmp(pos_tmp(:,1)>(min(pos_tmp(:,1)+half_time)),:); 

            % 1st half vs 2nd half 
             half_tresh = round((size(pos_tmp,1))/2);
             
             % Dividing spks and pos by half_tresh
             pos_1=pos_tmp(1:half_tresh,:); 
             pos_2= pos_tmp(half_tresh+1:end,:); 
             spks_1 = spks_tmp(spks_tmp<pos_tmp(half_tresh,1));
             spks_2 = spks_tmp(spks_tmp>pos_tmp(half_tresh,1));
                  
            %Calculate remapping parameters 
            [curve1,stats1] = FiringCurve(pos_1, spks_1 , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 2 , 'minPeak' , 0.1);
            [curve2,stats2] = FiringCurve(pos_2, spks_2, 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 2 , 'minPeak' , 0.1);
                    
            fr_1= nanmean(curve1.rate);
            fr_2= nanmean(curve2.rate);
                    
            %Fr change
            fr_change = abs((fr_1 - fr_2)/(fr_1 + fr_2));

            %Rate overlap
            if fr_1<=fr_2 
                overlap = fr_1/fr_2;
            else 
                overlap = fr_2/fr_1;
            end
    
            %Spatial  corr
            s = corrcoef(curve1.rate, curve2.rate);
            spatial = s(1,2);
    
            %Peak shift 
            shift = abs(stats1.x(1) - stats2.x(1));
            
            within1vs2 = [spatial, fr_change, overlap, shift]; 
            dHPC{ii}.within1vs2.rew =  within1vs2; 
            
            temp_stability_rew_dHPC = [temp_stability_rew_dHPC;within1vs2]; 
          
             
        end
        
        %Store 
        all_stability_ave_dHPC = [all_stability_ave_dHPC; temp_stability_ave_dHPC];  
        all_stability_rew_dHPC = [all_stability_rew_dHPC; temp_stability_rew_dHPC]; 
        
        disp('vHPC Stability calculation')
        
        load('vHPC_pc.mat'); 
       
        temp_stability_ave_vHPC = []; 
        temp_stability_rew_vHPC = []; 
        
        for ii=1:size(vHPC,2)
            
            cluster = vHPC{ii}.id;
  
            spks = spks_vHPC(spks_vHPC(:,1)==cluster,2); % select tspk from cluster
              
            % --- Aversive ---
             spks_tmp = Restrict(spks , movement.aversive); % Restrict spikes to movement periods
             pos_tmp = Restrict(behavior.pos.aversive(:,1:2) , movement.aversive); % Restrict pos to movement periods
             in_maze = ToIntervals(pos_tmp(:,1),and(pos_tmp(:,2)>30 , pos_tmp(:,2)<170));% eliminating the extrems of the maze
             spks_tmp = Restrict(spks_tmp,in_maze);
             pos_tmp = Restrict(pos_tmp,in_maze); 
             pos_tmp(:,2) = pos_tmp(:,2)-min(pos_tmp(:,2)); pos_tmp(:,2) = pos_tmp(:,2)/max(pos_tmp(:,2)); %normalization of position
              
             % 1st half vs 2nd half 
             half_tresh = round((size(pos_tmp,1))/2);
             
             % Dividing spks and pos by half_tresh
             pos_1=pos_tmp(1:half_tresh,:); 
             pos_2= pos_tmp(half_tresh+1:end,:); 
             spks_1 = spks_tmp(spks_tmp<pos_tmp(half_tresh,1));
             spks_2 = spks_tmp(spks_tmp>pos_tmp(half_tresh,1));
                   
             
%              % 1st half vs 2nd half 
%              total_time = max(pos_tmp(:,1))-min(pos_tmp(:,1)); 
%              half_time = total_time/2;
%              
%              %Divide spikes time in two
%              spks_1 = spks_tmp(spks_tmp<(min(pos_tmp(:,1)+half_time)));
%              pos_1 = pos_tmp(pos_tmp(:,1)<(min(pos_tmp(:,1)+half_time)),:); 
%              
%              spks_2 = spks_tmp(spks_tmp>=(min(pos_tmp(:,1)+half_time)));
%              pos_2 =  pos_tmp(pos_tmp(:,1)>(min(pos_tmp(:,1)+half_time)),:); 
    
            %Calculate remapping parameters 
            [curve1,stats1] = FiringCurve(pos_1, spks_1 , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 2 , 'minPeak' , 0.1);
            [curve2,stats2] = FiringCurve(pos_2, spks_2, 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 2 , 'minPeak' , 0.1);
                    
            fr_1= nanmean(curve1.rate);
            fr_2= nanmean(curve2.rate);
                    
            %Fr change
            fr_change = abs((fr_1 - fr_2)/(fr_1 + fr_2));

            %Rate overlap
            if fr_1<=fr_2 
                overlap = fr_1/fr_2;
            else 
                overlap = fr_2/fr_1;
            end
    
            %Spatial  corr
            s = corrcoef(curve1.rate, curve2.rate);
            spatial = s(1,2);
    
            %Peak shift 
            shift = abs(stats1.x(1) - stats2.x(1));
            
            within1vs2 = [spatial, fr_change, overlap, shift]; 
            
            vHPC{ii}.within1vs2.ave =  within1vs2; 
            temp_stability_ave_vHPC = [temp_stability_ave_vHPC;within1vs2];
            % --- Reward ---
            spks_tmp = Restrict(spks , movement.reward); % Restrict to movement periods
            pos_tmp = Restrict(behavior.pos.reward(:,1:2), movement.reward);
            in_maze = ToIntervals(pos_tmp(:,1),and(pos_tmp(:,2)>30 , pos_tmp(:,2)<170));% eliminating the extrems of the maze
            spks_tmp = Restrict(spks_tmp,in_maze);
            pos_tmp = Restrict(pos_tmp,in_maze); 
            pos_tmp(:,2) = pos_tmp(:,2)-min(pos_tmp(:,2)); pos_tmp(:,2) = pos_tmp(:,2)/max(pos_tmp(:,2)); %normalization of position
             
             % 1st half vs 2nd half 
             half_tresh = round((size(pos_tmp,1))/2);
             
             % Dividing spks and pos by half_tresh
             pos_1=pos_tmp(1:half_tresh,:); 
             pos_2= pos_tmp(half_tresh+1:end,:); 
             spks_1 = spks_tmp(spks_tmp<pos_tmp(half_tresh,1));
             spks_2 = spks_tmp(spks_tmp>pos_tmp(half_tresh,1));
             
%             % 1st half vs 2nd half 
%             total_time = max(pos_tmp(:,1))-min(pos_tmp(:,1)); 
%             half_time = total_time/2;
%              
%              %Divide spikes time in two
%              spks_1 = spks_tmp(spks_tmp<(min(pos_tmp(:,1)+half_time)));
%              pos_1 = pos_tmp(pos_tmp(:,1)<(min(pos_tmp(:,1)+half_time)),:); 
%              
%              spks_2 = spks_tmp(spks_tmp>=(min(pos_tmp(:,1)+half_time)));
%              pos_2 =  pos_tmp(pos_tmp(:,1)>(min(pos_tmp(:,1)+half_time)),:); 
    
            %Calculate remapping parameters 
            [curve1,stats1] = FiringCurve(pos_1, spks_1 , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 2 , 'minPeak' , 0.1);
            [curve2,stats2] = FiringCurve(pos_2, spks_2, 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 2 , 'minPeak' , 0.1);
                    
            fr_1= nanmean(curve1.rate);
            fr_2= nanmean(curve2.rate);
                    
            %Fr change
            fr_change = abs((fr_1 - fr_2)/(fr_1 + fr_2));

            %Rate overlap
            if fr_1<=fr_2 
                overlap = fr_1/fr_2;
            else 
                overlap = fr_2/fr_1;
            end
    
            %Spatial  corr
            s = corrcoef(curve1.rate, curve2.rate);
            spatial = s(1,2);
    
            %Peak shift 
            shift = abs(stats1.x(1) - stats2.x(1));
            
            within1vs2 = [spatial, fr_change, overlap, shift]; 
            vHPC{ii}.within1vs2.rew =  within1vs2; 
            temp_stability_rew_vHPC = [temp_stability_rew_vHPC;within1vs2];
             
        end 
        
         
        %Store 
        all_stability_ave_vHPC = [all_stability_ave_vHPC; temp_stability_ave_vHPC];  
        all_stability_rew_vHPC = [all_stability_rew_vHPC; temp_stability_rew_vHPC]; 
        
    end 
    
 end      

save('W:\Remapping-analysis-Facu\pc_all_stability_within.mat', 'all_stability_ave_dHPC', 'all_stability_rew_dHPC', 'all_stability_ave_vHPC', 'all_stability_rew_vHPC');

%Plot + stats stability 1st vs 2nd 
%c1 = spatial c2= fr_change c3=  overlap c4: pf shift 

data = [all_stability_ave_dHPC,ones(size(all_stability_ave_dHPC,1),1);...
    all_stability_rew_dHPC,ones(size(all_stability_rew_dHPC,1),1)*2];

data = [all_stability_ave_vHPC,ones(size(all_stability_ave_vHPC,1),1);...
    all_stability_rew_vHPC,ones(size(all_stability_rew_vHPC,1),1)*2];

figure, 
subplot(1,2,2); hold on; 
x = data(:,5);
y= data(:,4);%Change the column number (1-4) to choose which variable to plot 
scatter(x,y,"filled",'jitter','on', 'jitterAmount',0.1) , xlim([0 2]), % ylim([-0.8 1])
title('Stability 1st vs 2nd - vHPC')
ylabel('PF shift')
hold on
x = [1 2];
y = [nanmedian(data(data(:,5)==1,4)) , nanmedian(data(data(:,5)==2,4))]; % select the same c than y 
scatter(x,y, "filled") , xlim([0 3]), hold on

[h p] = ranksum(data(data(:,5)==1,4) , data(data(:,5)==2,4))



%% FIRING MAP PLOTS 

clear
clc
close all
% Parameters
path = {'\\Maryjackson\e\Rat127\Ephys\pyr';'\\Maryjackson\e\Rat128\Ephys\in_pyr\ready';'\\Maryjackson\e\Rat103\usable';'\\Maryjackson\e\Rat132\recordings\in_pyr'; '\\Maryjackson\e\Rat165\in_pyr'; ...
    '\\Maryjackson\e\Rat126\Ephys\in_Pyr'};%List of folders from the path

%Output matrix 
dhpc_ave = [];dhpc_rew = [];
vhpc_ave = [];vhpc_rew = [];

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
        cd([session,'\Spikesorting'])
        %Loading pc matrix of the sessions
        disp('Uploading session pc matrix');
        try
            load('dHPC_pc.mat');
            %Save session sub parameters 
        for d=1:size(dHPC,2)
            pc = dHPC{d}; 
            A = pc.frMap_ave - min(pc.frMap_ave);A = A ./ max(A);
            dhpc_ave=[dhpc_ave; A];
            clear A
          
            A = pc.frMap_rew - min(pc.frMap_rew);A = A ./ max(A);
            dhpc_rew=[dhpc_rew; A];
            clear A
        end
        
        catch 
            dHPC = []; disp(['No dHPC file in ',session]);
            
        end 
        
        
        try
            load('vHPC_pc.mat');
            for d=1:size(vHPC,2)
                pc = vHPC{d}; 
                A = pc.frMap_ave - min(pc.frMap_ave);A = A ./ max(A);
                vhpc_ave=[vhpc_ave; A];
                clear A
            
                A = pc.frMap_rew - min(pc.frMap_rew);A = A ./ max(A);
                vhpc_rew=[vhpc_rew; A];
                clear A
            end
        catch
            vHPC= []; disp(['No vHPCfile in ',session]);
        end
        
    end 
     
end 

%%%%% dHPC Plots %%%%%

[h idx] = max (dhpc_ave, [],2);
[m mm] = sort(idx); 

figure(1);clf;hold on;
fr = dhpc_ave(mm,:);  fr(~any(~isnan(fr), 2),:)=[];
subplot(1,2,1);imagesc([0:50:150], [1:1:size(dhpc_ave,1)],fr), colormap 'gray'; title('Aversive');
fr = dhpc_rew(mm,:);  fr(~any(~isnan(fr), 2),:)=[];
subplot(1,2,2);imagesc([0:50:150], [1:1:size(dhpc_rew,1)],fr), colormap 'gray'; title('Reward');
sgtitle('dHPC firing maps');

%%%%% vHPC Plots %%%%%

[h idx] = max (vhpc_ave, [],2);
[m mm] = sort(idx); 


figure(2);clf;
fr = vhpc_ave(mm,:);  fr(~any(~isnan(fr), 2),:)=[];
subplot(1,2,1); imagesc([0:50:150], [1:1:size(vhpc_ave,1)],fr), caxis([0 1]),colormap 'gray'; axis tight; title('Aversive');
fr = vhpc_rew(mm,:);  fr(~any(~isnan(fr), 2),:)=[];
subplot(1,2,2); imagesc([0:50:150], [1:1:size(vhpc_rew,1)],fr), caxis([0 1]), colormap 'gray'; title('Reward');axis tight; 


%% PC PARAMETERS:  pf size - in progress
% Firts run Parameters section
%c1 = spatial c2= fr_change c3=  overlap c4: pf shift c5:  1 = between 2 = within aversive 3= within reward 
clear
clc
close all
% Parameters
path = {'\\Maryjackson\e\Rat127\Ephys\pyr';'\\Maryjackson\e\Rat128\Ephys\in_pyr\ready';'\\Maryjackson\e\Rat103\usable';'\\Maryjackson\e\Rat132\recordings\in_pyr'; '\\Maryjackson\e\Rat165\in_pyr'; ...
    '\\Maryjackson\e\Rat126\Ephys\in_Pyr'};%List of folders from the path

%Output matrix 
dhpc_size = []; dhpc_skaggs = []; dhpc_n=0; 
vhpc_size = []; vhpc_skaggs = []; vhpc_n=0;

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
        cd([session,'\Spikesorting'])
        %Loading pc matrix of the sessions
        disp('Uploading session pc matrix');
        try
            load('dHPC_pc.mat');
             %Save session parameters 
        for d=1:size(dHPC,2)
            pc = dHPC{d}; 
            dhpc_size = [dhpc_size;pc.stats_ave.size(1),1;pc.stats_rew.size(1),2]; 
            dhpc_skaggs = [dhpc_skaggs;pc.stats_ave.specificity,1;pc.stats_rew.specificity,2];
            dhpc_n = dhpc_n+1;
        end
        catch 
            dHPC = []; disp(['No dHPC_pc_lap.mat file in ',session]);
        end 
        
        try
            load('vHPC_shock_remap.mat');
            for d=1:size(vHPC,2)
                pc = vHPC{d}; 
                vhpc_size = [vhpc_size;pc.stats_ave.size(1),1;pc.stats_rew.size(1),2]; 
                vhpc_skaggs = [vhpc_skaggs;pc.stats_ave.specificity,1;pc.stats_rew.specificity,2];
                vhpc_n = vhpc_n+1;
            end 
     
        catch
            vHPC= []; disp(['No vHPC_pc_lap.mat file in ',session]);
        end
        
    end 

     
end 


% Plots pf size
ylabels ={'Pf size(bins)'}; 
xlabels = {'Ave', 'Rew'}; 
figure(1);clf;hold on, 
sgtitle('Place field size')

    subplot(1,2,1); hold on; 
    x = dhpc_size(:,2);
    y= dhpc_size(:,1);
    c = [.3,.3,.3];
    scatter(x,y,[],c,"filled",'jitter','on', 'jitterAmount',0.1); xlim([0 3]);
    ylabel(ylabels);
    xticks([1 2])
    xticklabels({'Ave','Rew'});
    hold on
    s1=scatter(1,nanmedian(dhpc_size(dhpc_size(:,2)==1,1)), "filled");s1.MarkerFaceColor = [1 0.1 0.2];
    s3=scatter(2,nanmedian(dhpc_size(dhpc_size(:,2)==2,1)), "filled");s3.MarkerFaceColor = [0.1 0.3 1];
    title('dhpc'); 
    
    subplot(1,2,2); hold on; 
    x = vhpc_size(:,2);
    y= vhpc_size(:,1);
    c = [.3,.3,.3];
    scatter(x,y,[],c,"filled",'jitter','on', 'jitterAmount',0.1); xlim([0 3]);
    ylabel(ylabels);
    xticks([1 2])
    xticklabels({'Ave','Rew'});
    hold on
    s1=scatter(1,nanmedian(vhpc_size(vhpc_size(:,2)==1,1)), "filled");s1.MarkerFaceColor = [1 0.1 0.2];
    s3=scatter(2,nanmedian(vhpc_size(vhpc_size(:,2)==2,1)), "filled");s3.MarkerFaceColor = [0.1 0.3 1];
    title('vhpc'); 
 
% Plots skaggs
 
ylabels ={'Skaggs'}; 
xlabels = {'Ave', 'Rew'}; 
figure(1);clf;hold on, 
sgtitle('Skaggs')

    subplot(1,2,1); hold on; 
    x = dhpc_skaggs(:,2);
    y= dhpc_skaggs(:,1);
    c = [.3,.3,.3];
    scatter(x,y,[],c,"filled",'jitter','on', 'jitterAmount',0.1); xlim([0 3]);
    ylabel(ylabels);
    xticks([1 2])
    xticklabels({'Ave','Rew'});
    hold on
    s1=scatter(1,nanmedian(dhpc_skaggs(dhpc_skaggs(:,2)==1,1)), "filled");s1.MarkerFaceColor = [1 0.1 0.2];
    s3=scatter(2,nanmedian(dhpc_skaggs(dhpc_skaggs(:,2)==2,1)), "filled");s3.MarkerFaceColor = [0.1 0.3 1];
    title('dhpc'); 
    
    subplot(1,2,2); hold on; 
    x = vhpc_skaggs(:,2);
    y= vhpc_skaggs(:,1);
    c = [.3,.3,.3];
    scatter(x,y,[],c,"filled",'jitter','on', 'jitterAmount',0.1); xlim([0 3]);
    ylabel(ylabels);
    xticks([1 2])
    xticklabels({'Ave','Rew'});
    hold on
    s1=scatter(1,nanmedian(vhpc_skaggs(vhpc_skaggs(:,2)==1,1)), "filled");s1.MarkerFaceColor = [1 0.1 0.2];
    s3=scatter(2,nanmedian(vhpc_skaggs(vhpc_skaggs(:,2)==2,1)), "filled");s3.MarkerFaceColor = [0.1 0.3 1];
    title('vhpc'); 



%% SUBSAMPLED PLOTS
% Load subsampled parameters and creates matrix to plot and stats.
% Prev to run this section, you should have run the main loop. So, this
% section assumes that inside each session folder, you have a dHPC_pc and
% vHPC_pc file.
%c1 = spatial c2= fr_change c3=  overlap c4: pf shift c5:  1 = between 2 = within aversive 3= within reward 
clear
clc
close all
% Parameters
path = {'\\Maryjackson\e\Rat127\Ephys\pyr';'\\Maryjackson\e\Rat128\Ephys\in_pyr\ready';'\\Maryjackson\e\Rat103\usable';'\\Maryjackson\e\Rat132\recordings\in_pyr'; '\\Maryjackson\e\Rat165\in_pyr'; ...
    '\\Maryjackson\e\Rat126\Ephys\in_Pyr'};%List of folders from the path

%Output matrix 
dhpc_sub = [];
vhpc_sub = [];

for tt = 1:length(path)
    %List of folders from the path
    files = dir(path{tt});
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
    clear files dirFlags
    dhpc_temp = [];
    vhpc_temp = [];
    for t = 1 : length(subFolders)-2
        disp(['-- Initiating analysis of folder #' , num2str(t) , ' from rat #',num2str(tt) , ' --'])
        session = [subFolders(t+2).folder,'\',subFolders(t+2).name];
        cd([session,'\Spikesorting'])
        %Loading pc matrix of the sessions
        disp('Uploading session pc matrix');
        try
            load('dHPC_pc.mat');
        catch 
            dHPC = []; disp(['No dHPC_pc.mat file in ',session]);
        end 
        try
            load('vHPC_pc.mat');
        catch
            vHPC= []; disp(['No vHPC_pc.mat file in ',session]);
        end
        
        %Save session sub parameters 
        for d=1:size(dHPC,2)
            pc = dHPC{d}; 
            bet = [pc.between_sub,1]; 
            wa = [pc.within.ave_sub,2];
            wr = [pc.within.rew_sub,3];
            dhpc_temp = [dhpc_temp;bet;wa;wr]; 
            
        end

        for d=1:size(vHPC,2)
            pc = vHPC{d}; 
            bet = [pc.between_sub,1]; 
            wa = [pc.within.ave_sub,2];
            wr = [pc.within.rew_sub,3];
            vhpc_temp = [vhpc_temp;bet;wa;wr]; 
            
        end
        
    end
    
    %Save sub parameters 
    dhpc_sub = [dhpc_sub;dhpc_temp];
    vhpc_sub = [vhpc_sub;vhpc_temp];
     
end 

% Plots + stats 
%dHP
ylabels ={'Spatial correlation', 'Fr change', 'Overlap', 'Pf shift'}; 
xlabels = {'Between', 'Within Ave', 'Within Rew'}; 
figure(1);clf;hold on, 
sgtitle('dHPC')
for p=1:size(ylabels,2)
    subplot(1,4,p); hold on; 
    x = dhpc_sub(:,5);
    y= dhpc_sub(:,p);%Change the column number (1-4) to choose which variable to plot 
    c = [.3,.3,.3];
    scatter(x,y,[],c,"filled",'jitter','on', 'jitterAmount',0.1); xlim([0 4]);
    ylabel(ylabels{p});
    xticks([1 2 3])
    xticklabels({'Between', 'Within Ave', 'Within Rew'});
    hold on
    x = [1 2 3];
    y = [nanmedian(dhpc_sub(dhpc_sub(:,5)==1,p)), nanmedian(dhpc_sub(dhpc_sub(:,5)==2,p)),nanmedian(dhpc_sub(dhpc_sub(:,5)==3,p))]; 
    scatter(x,y, "filled");
end 

%Stats
[P,ANOVATAB,STATS] = kruskalwallis(dhpc_sub(:,1),dhpc_sub(:,5));
c = multcompare(STATS)

tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])

% Paired non-parametric anova - Friedman test
%Prepare data to test
p = 4; % 1=sp corr 2=fr change 3=overlap 4=pf shift
data=[dhpc_sub(dhpc_sub(:,5)==1,p), dhpc_sub(dhpc_sub(:,5)==2,p), ...
    dhpc_sub(dhpc_sub(:,5)==3,p)];
data(any(isnan(data), 2), :) = [];

[p,~,stats] =  friedman(data,1); 
c = multcompare(stats)
tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])


%vHP
ylabels ={'Spatial correlation', 'Fr change', 'Overlap', 'Pf shift'}; 
xlabels = {'Between', 'Within Ave', 'Within Rew'}; 
figure(2);clf;hold on, 
sgtitle('vHPC')
for p=1:size(ylabels,2)
    subplot(1,4,p); hold on; 
    x = vhpc_sub(:,5);
    y= vhpc_sub(:,p);%Change the column number (1-4) to choose which variable to plot 
    c = [.3,.3,.3];
    scatter(x,y,[],c,"filled",'jitter','on', 'jitterAmount',0.1); xlim([0 4]);
    ylabel(ylabels{p});
    xticks([1 2 3])
    xticklabels({'Between', 'Within Ave', 'Within Rew'});
    hold on
    x = [1 2 3];
    y = [nanmedian(vhpc_sub(vhpc_sub(:,5)==1,p)), nanmedian(vhpc_sub(vhpc_sub(:,5)==2,p)),nanmedian(vhpc_sub(vhpc_sub(:,5)==3,p))]; 
    scatter(x,y, "filled");
end 

%Stats
% [P,ANOVATAB,STATS] = kruskalwallis(vhpc_sub(:,1),vhpc_sub(:,5));
% c = multcompare(STATS)
% 
% tbl = array2table(c,"VariableNames", ...
%     ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])

% Paired non-parametric anova - Friedman test
%Prepare data to test
p = 4; % 1=sp corr 2=fr change 3=overlap 4=pf shift
data=[vhpc_sub(vhpc_sub(:,5)==1,p), vhpc_sub(vhpc_sub(:,5)==2,p), ...
    vhpc_sub(vhpc_sub(:,5)==3,p)];
data(any(isnan(data), 2), :) = [];

[p,~,stats] =  friedman(data,1); 
c = multcompare(stats)
tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])

%% LAP PLOTS
% Load lap parameters and creates matrix to plot and stats.
% Prev to run this section, you should have run the main loop. So, this
% section assumes that inside each session folder, you have a dHPC_pc and
% vHPC_pc file.
%c1 = spatial c2= fr_change c3=  overlap c4: pf shift c5:  1 = between 2 = within aversive 3= within reward 
clear
clc
close all
% Parameters
path = {'\\Maryjackson\e\Rat127\Ephys\pyr';'\\Maryjackson\e\Rat128\Ephys\in_pyr\ready';'\\Maryjackson\e\Rat103\usable';'\\Maryjackson\e\Rat132\recordings\in_pyr'; '\\Maryjackson\e\Rat165\in_pyr'; ...
    '\\Maryjackson\e\Rat126\Ephys\in_Pyr'};%List of folders from the path

%Output matrix 
dhpc_lap = [];
vhpc_lap = [];

for tt = 1:length(path)
    %List of folders from the path
    files = dir(path{tt});
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
    clear files dirFlags
    dhpc_temp = [];
    vhpc_temp = [];
    for t = 1 : length(subFolders)-2
        disp(['-- Initiating analysis of folder #' , num2str(t) , ' from rat #',num2str(tt) , ' --'])
        session = [subFolders(t+2).folder,'\',subFolders(t+2).name];
        cd([session,'\Spikesorting'])
        %Loading pc matrix of the sessions
        disp('Uploading session pc matrix');
        try
            load('dHPC_pc.mat');
        catch 
            dHPC = []; disp(['No dHPC_pc.mat file in ',session]);
        end 
        try
            load('vHPC_pc.mat');
        catch
            vHPC= []; disp(['No vHPC_pc.mat file in ',session]);
        end
        
        %Save session sub parameters 
        for d=1:size(dHPC,2)
            pc = dHPC{d}; 
            bet = [pc.lap.between_r,1]; 
            wa = [pc.lap.ave_r,2];
            wr = [pc.lap.rew_r,3];
            dhpc_temp = [dhpc_temp;bet;wa;wr]; 
            
        end

        for d=1:size(vHPC,2)
            pc = vHPC{d}; 
            bet = [pc.lap.between_r,1]; 
            wa = [pc.lap.ave_r,2];
            wr = [pc.lap.rew_r,3];
            vhpc_temp = [vhpc_temp;bet;wa;wr]; 
            
        end
        
    end
    
    %Save sub parameters 
    dhpc_lap = [dhpc_lap;dhpc_temp];
    vhpc_lap = [vhpc_lap;vhpc_temp];
     
end 

% Plots + stats 
%dHP
ylabels ={'Spatial correlation', 'Fr change', 'Overlap', 'Pf shift'}; 
xlabels = {'Between', 'Within Ave', 'Within Rew'}; 
figure(1);clf;hold on, 
sgtitle('dHPC random laps')
for p=1:size(ylabels,2)
    subplot(1,4,p); hold on; 
    x = dhpc_lap(:,5);
    y= dhpc_lap(:,p);%Change the column number (1-4) to choose which variable to plot 
    c = [.3,.3,.3];
    scatter(x,y,[],c,"filled",'jitter','on', 'jitterAmount',0.1); xlim([0 4]);
    ylabel(ylabels{p});
    xticks([1 2 3])
    xticklabels({'Between', 'Within Ave', 'Within Rew'})
    hold on
    x = [1 2 3];
    y = [nanmedian(dhpc_lap(dhpc_lap(:,5)==1,p)), nanmedian(dhpc_lap(dhpc_lap(:,5)==2,p)),nanmedian(dhpc_lap(dhpc_lap(:,5)==3,p))]; 
    scatter(x,y, "filled");
end 


% Paired non-parametric anova - Friedman test
%Prepare data to test
p = 1; % 1=sp corr 2=fr change 3=overlap 4=pf shift
data=[dhpc_lap(dhpc_lap(:,5)==1,p), dhpc_lap(dhpc_lap(:,5)==2,p), ...
    dhpc_lap(dhpc_lap(:,5)==3,p)];
data(any(isnan(data), 2), :) = [];

[p,~,stats] =  friedman(data,1); 
c = multcompare(stats)
tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])


%vHP
ylabels ={'Spatial correlation', 'Fr change', 'Overlap', 'Pf shift'}; 
xlabels = {'Between', 'Within Ave', 'Within Rew'}; 
figure(2);clf;hold on, 
sgtitle('vHPC random laps')
for p=1:size(ylabels,2)
    subplot(1,4,p); hold on; 
    x = vhpc_lap(:,5);
    y= vhpc_lap(:,p);%Change the column number (1-4) to choose which variable to plot 
    c = [.3,.3,.3];
    scatter(x,y,[],c,"filled",'jitter','on', 'jitterAmount',0.1); xlim([0 4]);
    ylabel(ylabels{p});
    xticks([1 2 3])
    xticklabels({'Between', 'Within Ave', 'Within Rew'});
    hold on
    x = [1 2 3];
    y = [nanmedian(vhpc_lap(vhpc_lap(:,5)==1,p)), nanmedian(vhpc_lap(vhpc_lap(:,5)==2,p)),nanmedian(vhpc_lap(vhpc_lap(:,5)==3,p))]; 
    scatter(x,y, "filled");
end 

%Stats
% [P,ANOVATAB,STATS] = kruskalwallis(vhpc_sub(:,1),vhpc_sub(:,5));
% c = multcompare(STATS)
% 
% tbl = array2table(c,"VariableNames", ...
%     ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])

% Paired non-parametric anova - Friedman test
%Prepare data to test
p = 1; % 1=sp corr 2=fr change 3=overlap 4=pf shift
data=[vhpc_lap(vhpc_lap(:,5)==1,p), vhpc_lap(vhpc_lap(:,5)==2,p), ...
    vhpc_lap(vhpc_lap(:,5)==3,p)];
data(any(isnan(data), 2), :) = [];

[p,~,stats] =  friedman(data,1); 
c = multcompare(stats)
tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])

%% N° TYPE LAP
clear
clc
close all
% Parameters
path = {'\\Maryjackson\e\Rat127\Ephys\pyr';'\\Maryjackson\e\Rat128\Ephys\in_pyr\ready';'\\Maryjackson\e\Rat103\usable';'\\Maryjackson\e\Rat132\recordings\in_pyr'; '\\Maryjackson\e\Rat165\in_pyr'; ...
    '\\Maryjackson\e\Rat126\Ephys\in_Pyr'};%List of folders from the path

%Output matrix 
dhpc_nlap = [];
vhpc_nlap = [];

for tt = 1:length(path)
    %List of folders from the path
    files = dir(path{tt});
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
    clear files dirFlags
    dhpc_temp = [];
    vhpc_temp = [];
    for t = 1 : length(subFolders)-2
        disp(['-- Initiating analysis of folder #' , num2str(t) , ' from rat #',num2str(tt) , ' --'])
        session = [subFolders(t+2).folder,'\',subFolders(t+2).name];
        cd([session,'\Spikesorting'])
        %Loading pc matrix of the sessions
        disp('Uploading session pc matrix');
        try
            load('dHPC_pc.mat');
            dhpc_temp = [dhpc_temp;dHPC{1}.nlap_ave,dHPC{1}.nlap_rew]; 
        catch 
            dHPC = []; disp(['No dHPC_pc.mat file in ',session]);
        end 
        try
            load('vHPC_pc.mat');
            vhpc_temp = [vhpc_temp;vHPC{1}.nlap_ave,vHPC{1}.nlap_rew]; 
        catch
            vHPC= []; disp(['No vHPC_pc.mat file in ',session]);
        end

    end
    
    %Save sub parameters 
    dhpc_nlap = [dhpc_nlap;dhpc_temp];
    vhpc_nlap = [vhpc_nlap;vhpc_temp];
     
end 


figure(3);clf;hold on;
sgtitle('dHPC laps')
histogram(dhpc_nlap(:,1),'FaceColor','r','EdgeColor','none','BinWidth',5);
histogram(dhpc_nlap(:,2),'FaceColor','b','EdgeColor','none','BinWidth',5);
ylim([0 12]);

figure(4);clf;hold on;
sgtitle('vHPC laps')
histogram(vhpc_nlap(:,1),'FaceColor','r','EdgeColor','none','BinWidth',5);
histogram(vhpc_nlap(:,2),'FaceColor','b','EdgeColor','none','BinWidth',5);


%% PLOT PC EXAMPLES

clear
clc
close all
path = {'\\Maryjackson\e\Rat127\Ephys\pyr';'\\Maryjackson\e\Rat128\Ephys\in_pyr\ready';'\\Maryjackson\e\Rat103\usable';'\\Maryjackson\e\Rat132\recordings\in_pyr'; '\\Maryjackson\e\Rat165\in_pyr'; ...
    '\\Maryjackson\e\Rat126\Ephys\in_Pyr'};%List of folders from the path
n_SU_V = [];
n_SU_D = [];
Xedges = 60; %number of bins for RateMap construction - 3cm bin
sigma = 2;%round(15/(180/Xedges)); %defined for gauss kernel of 15cm
binSize = 0.001; % bin size for replay events detection
bin_size = 1; % to bin pos ans spks in between/within  
% Behavior
minimal_speed = 2.5;% minimal speed to detect quite periods
minimal_speed_time = 2; % minimal time to detect quite periods


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
        % periods of movment during each condition
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
        %Video sampling rate 
        dt = (mean(diff(behavior.pos.aversive(:,1)))); 
        
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
            cellulartype = [K(:,1) , ones(length(K),1)];
        end
       
        
        %% Firing maps calculation
        disp('Plot dHPC Firing rate map')

        % Select pc from already run
        try
            load('dHPC_pc.mat');
        catch 
            dHPC = []; disp(['No dHPC_pc.mat file in ',session]);
        end 
        try
            load('vHPC_pc.mat');
        catch
            vHPC= []; disp(['No vHPC_pc.mat file in ',session]);
        end
        dhpc_pc = [];
        for d=1:size(dHPC,2)
            pc = dHPC{d}; 
            dhpc_pc = [dhpc_pc;pc.id]; 
            
        end

        vhpc_pc = [];
        for d=1:size(vHPC,2)
            pc = vHPC{d}; 
            vhpc_pc = [vhpc_pc;pc.id]; 
            
        end
        
        group_dHPC = group_dHPC(ismember(group_dHPC(:,1),dhpc_pc),:);
        group_vHPC = group_vHPC(ismember(group_vHPC(:,1),vhpc_pc),:);
        
        % Plotting: 
        for ii=1:size(group_dHPC,1)
            cluster = group_dHPC(ii,1);

            spks = spks_dHPC(spks_dHPC(:,1)==cluster,2); % select tspk from cluster
              
            
            % --- Aversive ---
            spks_tmp = spks; 
            pos_tmp = behavior.pos.aversive(:,1:2); 
            in_maze = ToIntervals(pos_tmp(:,1),and(pos_tmp(:,2)>30 , pos_tmp(:,2)<180));% eliminating the extrems of the maze
                    
            %Check lenght of in_maze segments to define real laps
             in_lapA = []; % keeps only full laps (maybe too restrictive)
             for ix=1:size(in_maze,1)
                        ini = in_maze(ix,1);
                        fin = in_maze(ix,2);
                        ini_pos = pos_tmp(pos_tmp(:,1)==ini,2); 
                        fin_pos = pos_tmp(pos_tmp(:,1)==fin,2); 
                        if abs(fin_pos-ini_pos)>100
                             in_lapA = [in_lapA;in_maze(ix,:)]; 
                        end
              end
 
            %Restrict spk and position to laps and movement periods
            spks_tmp = Restrict(spks_tmp,in_lapA);spks_tmp = Restrict(spks_tmp , movement.aversive); % Restrict spikes to movement periods
            pos_tmp = Restrict(pos_tmp,in_lapA);
            pos_tmp = Restrict(pos_tmp, movement.aversive); % Restrict pos to movement periods

            pos_tmp(:,2) = pos_tmp(:,2)-min(pos_tmp(:,2)); pos_tmp(:,2) = pos_tmp(:,2)/max(pos_tmp(:,2)); %normalization of position
  
            %Firing curve construction
             [curveA , statsA] = FiringCurve(pos_tmp , spks_tmp , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
             curveA_dhpc = curveA.rate; 

                    
          %%%%%%% PLOT
          figure; hold on; 
          sgtitle(['dHPC ',session(end-14:end),' id ',num2str(cluster),' iteration ',num2str(ii)])
          subplot(4,2,[1 5]); hold on; plot(pos_tmp(:,2),pos_tmp(:,1),'color','red');hold on; axis tight;
          title('Aversive');
          %Find interpolated position of each spike:
          xs = interp1(pos_tmp(:,1),pos_tmp(:,2),spks_tmp); scatter(xs,spks_tmp,'filled', 'MarkerFaceColor',[.3 .3 .3]); 
          %Fr map 
          subplot(4,2,7); hold on;imagesc(curveA.rate), colormap 'jet';
  
          % --- Reward ---
          spks_tmp = spks; 
          pos_tmp = behavior.pos.reward(:,1:2);
          in_maze = ToIntervals(pos_tmp(:,1),and(pos_tmp(:,2)>30 , pos_tmp(:,2)<180));% eliminating the extrems of the maze
                    
          %Check lenght of in_maze segments to define real laps
          in_lapR = []; 
          for ix=1:size(in_maze,1)
                        ini = in_maze(ix,1);
                        fin = in_maze(ix,2);
                        ini_pos = pos_tmp(pos_tmp(:,1)==ini,2); 
                        fin_pos = pos_tmp(pos_tmp(:,1)==fin,2); 
                        if abs(fin_pos-ini_pos)>100
                             in_lapR = [in_lapR;in_maze(ix,:)]; 
                        end
          end 
                    

          %Restrict spks and pos to valid laps 
          spks_tmp = Restrict(spks_tmp,in_lapR); spks_tmp = Restrict(spks_tmp, movement.reward); % Restrict to movement periods
          pos_tmp = Restrict(pos_tmp,in_lapR); pos_tmp = Restrict(pos_tmp , movement.reward);
          pos_tmp(:,2) = pos_tmp(:,2)-min(pos_tmp(:,2)); pos_tmp(:,2) = pos_tmp(:,2)/max(pos_tmp(:,2)); %normalization of position

          %Firing curve construction
          [curveR , statsR] = FiringCurve(pos_tmp , spks_tmp , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4 , 'minPeak' , 0.2);
          curveR_dhpc = curveR.rate;
                    
          subplot(4,2,[2 6]); hold on; plot(pos_tmp(:,2),pos_tmp(:,1),'color',[0.175 0.54 0.60] );hold on; axis tight;
          title('Reward');
          %Find interpolated position of each spike:
          xs = interp1(pos_tmp(:,1),pos_tmp(:,2),spks_tmp); scatter(xs,spks_tmp,'filled', 'MarkerFaceColor',[.3 .3 .3]); 
          %Fr map 
          subplot(4,2,8); hold on;imagesc(curveR.rate), colormap 'jet';

        end
        
        disp('Plot vHPC Firing rate map')
          
        for ii=1:size(group_vHPC,1)
            cluster = group_vHPC(ii,1);

            spks = spks_vHPC(spks_vHPC(:,1)==cluster,2); % select tspk from cluster
              
            
            % --- Aversive ---
            spks_tmp = spks; 
            pos_tmp = behavior.pos.aversive(:,1:2); 
            in_maze = ToIntervals(pos_tmp(:,1),and(pos_tmp(:,2)>30 , pos_tmp(:,2)<180));% eliminating the extrems of the maze
                    
            %Check lenght of in_maze segments to define real laps
             in_lapA = []; % keeps only full laps (maybe too restrictive)
             for ix=1:size(in_maze,1)
                        ini = in_maze(ix,1);
                        fin = in_maze(ix,2);
                        ini_pos = pos_tmp(pos_tmp(:,1)==ini,2); 
                        fin_pos = pos_tmp(pos_tmp(:,1)==fin,2); 
                        if abs(fin_pos-ini_pos)>100
                             in_lapA = [in_lapA;in_maze(ix,:)]; 
                        end
              end
 
            %Restrict spk and position to laps and movement periods
            spks_tmp = Restrict(spks_tmp,in_lapA);spks_tmp = Restrict(spks_tmp , movement.aversive); % Restrict spikes to movement periods
            pos_tmp = Restrict(pos_tmp,in_lapA);
            pos_tmp = Restrict(pos_tmp, movement.aversive); % Restrict pos to movement periods

            pos_tmp(:,2) = pos_tmp(:,2)-min(pos_tmp(:,2)); pos_tmp(:,2) = pos_tmp(:,2)/max(pos_tmp(:,2)); %normalization of position
  
            %Firing curve construction
             [curveA , statsA] = FiringCurve(pos_tmp , spks_tmp , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
             curveA_dhpc = curveA.rate; 

                    
          %%%%%%% PLOT
          figure; hold on; 
          sgtitle(['vHPC ',session(end-14:end),' id ',num2str(cluster),' iteration ',num2str(ii)])
          subplot(4,2,[1 5]); hold on; plot(pos_tmp(:,2),pos_tmp(:,1),'color', 'red');hold on; axis tight;
          title('Aversive');
          %Find interpolated position of each spike:
          xs = interp1(pos_tmp(:,1),pos_tmp(:,2),spks_tmp); scatter(xs,spks_tmp,'filled', 'MarkerFaceColor',[.3 .3 .3] ); 
          %Fr map 
          subplot(4,2,7); hold on;imagesc(curveA.rate), colormap 'jet';
  
          % --- Reward ---
          spks_tmp = spks; 
          pos_tmp = behavior.pos.reward(:,1:2);
          in_maze = ToIntervals(pos_tmp(:,1),and(pos_tmp(:,2)>30 , pos_tmp(:,2)<180));% eliminating the extrems of the maze
                    
          %Check lenght of in_maze segments to define real laps
          in_lapR = []; 
          for ix=1:size(in_maze,1)
                        ini = in_maze(ix,1);
                        fin = in_maze(ix,2);
                        ini_pos = pos_tmp(pos_tmp(:,1)==ini,2); 
                        fin_pos = pos_tmp(pos_tmp(:,1)==fin,2); 
                        if abs(fin_pos-ini_pos)>100
                             in_lapR = [in_lapR;in_maze(ix,:)]; 
                        end
          end 
                    

          %Restrict spks and pos to valid laps 
          spks_tmp = Restrict(spks_tmp,in_lapR); spks_tmp = Restrict(spks_tmp, movement.reward); % Restrict to movement periods
          pos_tmp = Restrict(pos_tmp,in_lapR); pos_tmp = Restrict(pos_tmp , movement.reward);
          pos_tmp(:,2) = pos_tmp(:,2)-min(pos_tmp(:,2)); pos_tmp(:,2) = pos_tmp(:,2)/max(pos_tmp(:,2)); %normalization of position

          %Firing curve construction
          [curveR , statsR] = FiringCurve(pos_tmp , spks_tmp , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4 , 'minPeak' , 0.2);
          curveR_dhpc = curveR.rate;
                    
          subplot(4,2,[2 6]); hold on; plot(pos_tmp(:,2),pos_tmp(:,1),'color',[0.175 0.54 0.60]);hold on; axis tight;
          title('Reward');
          %Find interpolated position of each spike:
          xs = interp1(pos_tmp(:,1),pos_tmp(:,2),spks_tmp); scatter(xs,spks_tmp,'filled', 'MarkerFaceColor',[.3 .3 .3]); 
          %Fr map 
          subplot(4,2,8); hold on;imagesc(curveR.rate), colormap 'jet';

       end
        

        disp(['-- Finished folder #' , num2str(t) , ' from rat #' , num2str(tt) , ' ---'])
        disp('  ')
    end
    disp(['-------- Finished rat#' , num2str(tt) , ' --------'])
    disp(' ')
end


%% FIRING CURVES OF ALL NEURONS 
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
        
        %Session output
        dHPC = {};% one cell per pc
        vHPC = {}; 
        
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
        % periods of movment during each condition
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
        %Video sampling rate 
        dt = (mean(diff(behavior.pos.aversive(:,1)))); 
        
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
        
        for ii=1:size(group_dHPC,1)

            cluster = group_dHPC(ii,1);
            celltype = logical(cellulartype(cellulartype(:,1) == cluster,2));
            
            if celltype % check if pyr
            spks = spks_dHPC(spks_dHPC(:,1)==cluster,2); % select tspk from cluster

            % --- Aversive ---
            spks_ave = spks; 
            pos_ave = behavior.pos.aversive(:,1:2); 
            in_maze = ToIntervals(pos_ave(:,1),and(pos_ave(:,2)>30 , pos_ave(:,2)<180));% eliminating the extrems of the maze
                    
            %Check lenght of in_maze segments to define real laps
             in_lapA = []; % keeps only full laps (maybe too restrictive)
                    for ix=1:size(in_maze,1)
                        ini = in_maze(ix,1);
                        fin = in_maze(ix,2);
                        ini_pos = pos_ave(pos_ave(:,1)==ini,2); 
                        fin_pos = pos_ave(pos_ave(:,1)==fin,2); 
                        if abs(fin_pos-ini_pos)>100
                             in_lapA = [in_lapA;in_maze(ix,:)]; 
                        end
                    end
                    
             %Restrict spk and position to laps and movement periods
             spks_ave = Restrict(spks_ave,in_lapA);spks_ave = Restrict(spks_ave , movement.aversive); % Restrict spikes to movement periods and laps
                    
             pos_ave = Restrict(pos_ave,in_lapA);pos_ave = Restrict(pos_ave, movement.aversive); % Restrict pos to movement periods and laps
             pos_ave(:,2) = pos_ave(:,2)-min(pos_ave(:,2)); pos_ave(:,2) = pos_ave(:,2)/max(pos_ave(:,2)); %normalization of position
               
             % --- Reward ---
             spks_rew = spks; 
             pos_rew = behavior.pos.reward(:,1:2);
             in_maze = ToIntervals(pos_rew(:,1),and(pos_rew(:,2)>30 , pos_rew(:,2)<180));% eliminating the extrems of the maze
               
             %Check lenght of in_maze segments to define real laps
             in_lapR = []; 
                    for ix=1:size(in_maze,1)
                        ini = in_maze(ix,1);
                        fin = in_maze(ix,2);
                        ini_pos = pos_rew(pos_rew(:,1)==ini,2); 
                        fin_pos = pos_rew(pos_rew(:,1)==fin,2); 
                        if abs(fin_pos-ini_pos)>100
                             in_lapR = [in_lapR;in_maze(ix,:)]; 
                        end
                    end 
                    
             %Restrict spk and position to laps and movement periods
             spks_rew = Restrict(spks_rew,in_lapR); spks_rew = Restrict(spks_rew, movement.reward); % Restrict to movement periods
             pos_rew = Restrict(pos_rew,in_lapR); pos_rew = Restrict(pos_rew , movement.reward);
             pos_rew(:,2) = pos_rew(:,2)-min(pos_rew(:,2)); pos_rew(:,2) = pos_rew(:,2)/max(pos_rew(:,2)); %normalization of position 
 
             [curveA , statsA] = FiringCurve(pos_ave , spks_ave , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
             [curveR , statsR] = FiringCurve(pos_rew , spks_rew , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
  
             % Save PC varaiables 
                        n.id = cluster;  
                        n.nlap_ave = size(in_lapA,1);
                        n.nlap_rew = size(in_lapR,1);
                        n.frMap_ave = curveA.rate;
                        n.frMap_rew = curveR.rate;
                        n.stats_ave = statsA;
                        n.stats_rew = statsR;
                        dHPC{ii}= n;
                  
          
            end 
         end
           
        disp('vHPC Firing rate map calculation')
     
        for ii=1:size(group_vHPC,1)

            cluster = group_vHPC(ii,1);
            celltype = logical(cellulartype(cellulartype(:,1) == cluster,2));
            
             if celltype % check if pyr
                spks = spks_vHPC(spks_vHPC(:,1)==cluster,2); % select tspk from cluster

                % --- Aversive ---
                spks_ave = spks; 
                pos_ave = behavior.pos.aversive(:,1:2); 
                in_maze = ToIntervals(pos_ave(:,1),and(pos_ave(:,2)>30 , pos_ave(:,2)<180));% eliminating the extrems of the maze
                    
                %Check lenght of in_maze segments to define real laps
                in_lapA = []; % keeps only full laps (maybe too restrictive)
                    for ix=1:size(in_maze,1)
                        ini = in_maze(ix,1);
                        fin = in_maze(ix,2);
                        ini_pos = pos_ave(pos_ave(:,1)==ini,2); 
                        fin_pos = pos_ave(pos_ave(:,1)==fin,2); 
                        if abs(fin_pos-ini_pos)>100
                             in_lapA = [in_lapA;in_maze(ix,:)]; 
                        end
                    end
                    
                %Restrict spk and position to laps and movement periods
                spks_ave = Restrict(spks_ave,in_lapA);spks_ave = Restrict(spks_ave , movement.aversive); % Restrict spikes to movement periods and laps
                    
                pos_ave = Restrict(pos_ave,in_lapA);pos_ave = Restrict(pos_ave, movement.aversive); % Restrict pos to movement periods and laps
                pos_ave(:,2) = pos_ave(:,2)-min(pos_ave(:,2)); pos_ave(:,2) = pos_ave(:,2)/max(pos_ave(:,2)); %normalization of position
               
                % --- Reward ---
                spks_rew = spks; 
                pos_rew = behavior.pos.reward(:,1:2);
                in_maze = ToIntervals(pos_rew(:,1),and(pos_rew(:,2)>30 , pos_rew(:,2)<180));% eliminating the extrems of the maze
               
                %Check lenght of in_maze segments to define real laps
                in_lapR = []; 
                    for ix=1:size(in_maze,1)
                        ini = in_maze(ix,1);
                        fin = in_maze(ix,2);
                        ini_pos = pos_rew(pos_rew(:,1)==ini,2); 
                        fin_pos = pos_rew(pos_rew(:,1)==fin,2); 
                        if abs(fin_pos-ini_pos)>100
                             in_lapR = [in_lapR;in_maze(ix,:)]; 
                        end
                    end 
                    
                %Restrict spk and position to laps and movement periods
                spks_rew = Restrict(spks_rew,in_lapR); spks_rew = Restrict(spks_rew, movement.reward); % Restrict to movement periods
                pos_rew = Restrict(pos_rew,in_lapR); pos_rew = Restrict(pos_rew , movement.reward);
                pos_rew(:,2) = pos_rew(:,2)-min(pos_rew(:,2)); pos_rew(:,2) = pos_rew(:,2)/max(pos_rew(:,2)); %normalization of position 
                    

     
                [curveA , statsA] = FiringCurve(pos_ave , spks_ave , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                [curveR , statsR] = FiringCurve(pos_rew , spks_rew , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    
                
                        
                    
                        % Save PC varaiables 
                        n.id = cluster; 
                        n.celltype = celltype; 
                        n.nlap_ave = size(in_lapA,1);
                        n.nlap_rew = size(in_lapR,1);
                        n.frMap_ave = curveA.rate;
                        n.frMap_rew = curveR.rate;
                        n.stats_ave = statsA;
                        n.stats_rew = statsR;
                       
                        vHPC{ii}= n;

             end 
         end

    %% Saveing PC INFO 
    if ~isempty(dHPC)
        dHPC = dHPC(~cellfun('isempty',dHPC));
        save([cd,'\dHPC_all.mat'],'dHPC'); 
    end
    if ~isempty(vHPC)
        vHPC = vHPC(~cellfun('isempty',vHPC));
        save([cd,'\vHPC_all.mat'],'vHPC'); 
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


