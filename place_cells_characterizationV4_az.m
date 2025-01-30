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

binSize = 0.001; % bin size for replay events detection
bin_size = 1; % to bin pos ans spks in between/within  
% Behavior
minimal_speed = 2.5;% minimal speed to detect quite periods
minimal_speed_time = 2; % minimal time to detect quite periods
%Place cells map
%Parameters for FiringCurve.m
min_time = 0.15; %min time in a given bin to be consider
min_peak = 0.2; %
min_size = 4;% min number of bins for the place fields
Xedges = 60; %number of bins for RateMap construction - 2.5cm bin (because of the extremes cutting)
sigma = 2;%round(15/(180/Xedges)); %defined for gauss kernel of 15cm
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
        
        %Restrict pos to inside maze and only movement periods
        %---Aversive----
        pos_ave = behavior.pos.aversive(:,1:2); 
        pos_ave(:,2) = pos_ave(:,2)-min(pos_ave(:,2)); pos_ave(:,2) = pos_ave(:,2)/max(pos_ave(:,2)); %normalization of position 
        in_maze = ToIntervals(pos_ave(:,1),and(pos_ave(:,2)>0.09 , pos_ave(:,2)<1-0.09));% eliminating the extrems of the maze
                    
       %Check lenght of in_maze segments to define real laps
        in_lapA = []; % keeps only full laps (maybe too restrictive)
        for ix=1:size(in_maze,1)
                        ini = in_maze(ix,1);
                        fin = in_maze(ix,2);
                        ini_pos = pos_ave(pos_ave(:,1)==ini,2); 
                        fin_pos = pos_ave(pos_ave(:,1)==fin,2); 
                        if abs(fin_pos-ini_pos)>0.6
                             in_lapA = [in_lapA;in_maze(ix,:)]; 
                        end
        end
        
        pos_ave = Restrict(pos_ave,in_lapA);pos_ave = Restrict(pos_ave, movement.aversive); % Restrict pos to movement periods and laps
       
         %---Rewarded----
        pos_rew = behavior.pos.reward(:,1:2);
        pos_rew(:,2) = pos_rew(:,2)-min(pos_rew(:,2)); pos_rew(:,2) = pos_rew(:,2)/max(pos_rew(:,2)); %normalization of position 
        in_maze = ToIntervals(pos_rew(:,1),and(pos_rew(:,2)>0.09 , pos_rew(:,2)<1-0.09));% eliminating the extrems of the maze (10cm-reward box)
               
        %Check lenght of in_maze segments to define real laps
        in_lapR = []; 
        for ix=1:size(in_maze,1)
                        ini = in_maze(ix,1);
                        fin = in_maze(ix,2);
                        ini_pos = pos_rew(pos_rew(:,1)==ini,2); 
                        fin_pos = pos_rew(pos_rew(:,1)==fin,2); 
                        if abs(fin_pos-ini_pos)>0.6
                             in_lapR = [in_lapR;in_maze(ix,:)]; 
                        end
        end 
        pos_rew = Restrict(pos_rew,in_lapR); pos_rew = Restrict(pos_rew , movement.reward);
        
        
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
                %Restrict spk and position to laps and movement periods
                spks_ave = Restrict(spks_ave,in_lapA);spks_ave = Restrict(spks_ave , movement.aversive); % Restrict spikes to movement periods and laps
               
                % --- Reward ---
                spks_rew = spks; 
                %Restrict spk  to laps and movement periods
                spks_rew = Restrict(spks_rew,in_lapR); spks_rew = Restrict(spks_rew, movement.reward); % Restrict to movement periods
               
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
                % 1. Mean firing rate above tresh.
                % 2.Skaggs > skaggs random. 
                % 3. At least 1 pf > 4 bins in one of the cond.  
     
                [curveA , statsA] = FiringCurve(pos_ave , spks_ave , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 6, 'minPeak' , 0.2);
                [curveR , statsR] = FiringCurve(pos_rew , spks_rew , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 6, 'minPeak' , 0.2);
                    
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
                
                
                %Compute Q90 skaggs
                qa = SkaggsRandomFMT_1D(spks_ave, pos_ave,sigma, Xedges, 0.90);
                qr = SkaggsRandomFMT_1D(spks_rew, pos_rew,sigma, Xedges, 0.90);
%                 disp(['qa= ', num2str(qa), 'skaggs=' num2str(statsA.specificity)]);
%                 disp(['qr= ', num2str(qr), 'skaggs=' num2str(statsR.specificity)]);

                if or(statsA.specificity>qa, statsR.specificity>qr) % 2. criteria
                    if and(or(sum(statsA.field(:,:,1))>= 4 , sum(statsR.field(:,:,1))>= 4),or(fr_ave>tresh,fr_rew>tresh)) % 1. and 3. criteria
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
                        
                        %Place field centero of mass COM
                        % Aversive
                        comA=nan; 
                        if sum(sum(~isnan(statsA.fieldX)))>0
                            pf_lim = statsA.fieldX; 
                            pc_frmap = curveA.rate; 
                            field = pc_frmap(pf_lim(1,1):pf_lim(1,2));
                
                            %Center of mass of the field 
                            c = 1:size(field,2);                             
                            com = sum(c .* field) / sum(field);
                            comA = com + pf_lim(1,1); % back in general scale bins 
                            clear pf_lim pc_frmap field c com
                        end 
                        
                        % Reward
                        comR=[];
                        if sum(sum(~isnan(statsR.fieldX)))>0
                            pf_lim = statsR.fieldX; 
                            pc_frmap = curveR.rate; 
                            field = pc_frmap(pf_lim(1,1):pf_lim(1,2));
                
                            %Center of mass of the field 
                            c = 1:size(field,2);                             
                            com = sum(c .* field) / sum(field);
                            comR = com + pf_lim(1,1); % back in general scale bins 
                            clear pf_lim pc_frmap field c com
                        end 
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
                        n.remap_tresh.ave=withinA_tresh;
                        n.remap_tresh.rew=withinR_tresh;
                        n.subsampled.between = between_sub;
                        n.subsampled.within_ave = withinA_sub; 
                        n.subsampled.within_rew = withinR_sub; 
                        n.com.ave = comA;
                        n.com.rew = comR;
                        
                        dHPC{ii}= n;
                    end 
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
                %Restrict spk and position to laps and movement periods
                spks_ave = Restrict(spks_ave,in_lapA);spks_ave = Restrict(spks_ave , movement.aversive); % Restrict spikes to movement periods and laps

                % --- Reward ---
                spks_rew = spks; 
                %Restrict spk and position to laps and movement periods
                spks_rew = Restrict(spks_rew,in_lapR); spks_rew = Restrict(spks_rew, movement.reward); % Restrict to movement periods
                
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
                % 1. Mean firing rate above tresh.
                % 2.Skaggs > skaggs random. 
                % 3. At least 1 pf > 4 bins in one of the cond.  
     
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
                 
                %Compute Q90 skaggs
                qa = SkaggsRandomFMT_1D(spks_ave, pos_ave,sigma, Xedges, 0.90);
                qr = SkaggsRandomFMT_1D(spks_rew, pos_rew,sigma, Xedges, 0.90);


                if or(statsA.specificity>qa, statsR.specificity>qr) % 2. criteria
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
                        
                        %Place field centero of mass COM
                        % Aversive
                        comA=nan; 
                        if sum(sum(~isnan(statsA.fieldX)))>0
                            pf_lim = statsA.fieldX; 
                            pc_frmap = curveA.rate; 
                            field = pc_frmap(pf_lim(1,1):pf_lim(1,2));
                
                            %Center of mass of the field 
                            c = 1:size(field,2);                             
                            com = sum(c .* field) / sum(field);
                            comA = com + pf_lim(1,1); % back in general scale bins 
                            clear pf_lim pc_frmap field c com
                        end
                        
                        %Reward
                        comR=nan; 
                        if sum(sum(~isnan(statsR.fieldX)))>0
                            pf_lim = statsR.fieldX; 
                            pc_frmap = curveR.rate; 
                            field = pc_frmap(pf_lim(1,1):pf_lim(1,2));
                
                            %Center of mass of the field 
                            c = 1:size(field,2);                             
                            com = sum(c .* field) / sum(field);
                            comR = com + pf_lim(1,1); % back in general scale bins 
                            clear pf_lim pc_frmap field c com
                        end 
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
                        n.remap_tresh.ave=withinA_tresh;
                        n.remap_tresh.rew=withinR_tresh;
                        n.subsampled.between = between_sub;
                        n.subsampled.within_ave = withinA_sub; 
                        n.subsampled.within_rew = withinR_sub; 
                        n.com.ave = comA;
                        n.com.rew = comR;
                        
                        vHPC{ii}= n;
                    end 
                end
            end
         end

    %% Saveing PC INFO 
    if ~isempty(dHPC)
        dHPC = dHPC(~cellfun('isempty',dHPC));
        save([cd,'\dHPC_pc_skaggs_circular.mat'],'dHPC'); 
    end
    if ~isempty(vHPC)
        vHPC = vHPC(~cellfun('isempty',vHPC));
        save([cd,'\vHPC_pc_skaggs_circular.mat'],'vHPC'); 
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

%% # PC PER CONDITION
count_dhpc = 0; count_pyr_dhpc=[]; 
pc_ave_dhpc = 0;pc_rew_dhpc = 0;
count_vhpc = 0; count_pyr_vhpc=[]; 
pc_ave_vhpc = 0;pc_rew_vhpc = 0;

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
        
        %Restrict pos to inside maze and only movement periods
        %---Aversive----
        pos_ave = behavior.pos.aversive(:,1:2); 
        pos_ave(:,2) = pos_ave(:,2)-min(pos_ave(:,2)); pos_ave(:,2) = pos_ave(:,2)/max(pos_ave(:,2)); %normalization of position 
        in_maze = ToIntervals(pos_ave(:,1),and(pos_ave(:,2)>0.09 , pos_ave(:,2)<1-0.09));% eliminating the extrems of the maze
                    
       %Check lenght of in_maze segments to define real laps
        in_lapA = []; % keeps only full laps (maybe too restrictive)
        for ix=1:size(in_maze,1)
                        ini = in_maze(ix,1);
                        fin = in_maze(ix,2);
                        ini_pos = pos_ave(pos_ave(:,1)==ini,2); 
                        fin_pos = pos_ave(pos_ave(:,1)==fin,2); 
                        if abs(fin_pos-ini_pos)>0.6
                             in_lapA = [in_lapA;in_maze(ix,:)]; 
                        end
        end
        
        pos_ave = Restrict(pos_ave,in_lapA);pos_ave = Restrict(pos_ave, movement.aversive); % Restrict pos to movement periods and laps
       
         %---Rewarded----
        pos_rew = behavior.pos.reward(:,1:2);
        pos_rew(:,2) = pos_rew(:,2)-min(pos_rew(:,2)); pos_rew(:,2) = pos_rew(:,2)/max(pos_rew(:,2)); %normalization of position 
        in_maze = ToIntervals(pos_rew(:,1),and(pos_rew(:,2)>0.09 , pos_rew(:,2)<1-0.09));% eliminating the extrems of the maze (10cm-reward box)
               
        %Check lenght of in_maze segments to define real laps
        in_lapR = []; 
        for ix=1:size(in_maze,1)
                        ini = in_maze(ix,1);
                        fin = in_maze(ix,2);
                        ini_pos = pos_rew(pos_rew(:,1)==ini,2); 
                        fin_pos = pos_rew(pos_rew(:,1)==fin,2); 
                        if abs(fin_pos-ini_pos)>0.6
                             in_lapR = [in_lapR;in_maze(ix,:)]; 
                        end
        end 
        pos_rew = Restrict(pos_rew,in_lapR); pos_rew = Restrict(pos_rew , movement.reward);
        
        
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
                
                
                count_pyr_dhpc=count_pyr_dhpc+1; 
                spks = spks_dHPC(spks_dHPC(:,1)==cluster,2); % select tspk from cluster

                % --- Aversive ---
                spks_ave = spks; 
                %Restrict spk and position to laps and movement periods
                spks_ave = Restrict(spks_ave,in_lapA);spks_ave = Restrict(spks_ave , movement.aversive); % Restrict spikes to movement periods and laps
               
                % --- Reward ---
                spks_rew = spks; 
                %Restrict spk  to laps and movement periods
                spks_rew = Restrict(spks_rew,in_lapR); spks_rew = Restrict(spks_rew, movement.reward); % Restrict to movement periods
               
                                 
   
                %Checking pc criteria
                % 1. Mean firing rate above tresh.
                % 2.Skaggs > skaggs random. 
                % 3. At least 1 pf > 4 bins in one of the cond.  
     
                [curveA , statsA] = FiringCurve(pos_ave , spks_ave , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 6, 'minPeak' , 0.2);
                [curveR , statsR] = FiringCurve(pos_rew , spks_rew , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 6, 'minPeak' , 0.2);
                    
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
                
                
                %Compute Q90 skaggs
                qa = SkaggsRandomFMT_1D(spks_ave, pos_ave,movement.aversive(1,1),movement.aversive(end,1),sigma, Xedges, 0.90);
                qr = SkaggsRandomFMT_1D(spks_rew, pos_rew,movement.reward(1,1),movement.reward(end,1), sigma, Xedges, 0.90);
%                 disp(['qa= ', num2str(qa), 'skaggs=' num2str(statsA.specificity)]);
%                 disp(['qr= ', num2str(qr), 'skaggs=' num2str(statsR.specificity)]);

                % PC AVE 
                if statsA.specificity>qa % 2. criteria
                    if and(sum(statsA.field(:,:,1))>= 4 ,fr_ave>tresh) % 1. and 3. criteria
                        pc_ave_dhpc = pc_ave_dhpc+1;
                    end 
                end

                %PC REWARD 
                if statsR.specificity>qr % 2. criteria
                    if and(sum(statsR.field(:,:,1))>= 4 ,fr_rew>tresh) % 1. and 3. criteria
                        pc_rew_dhpc = pc_rew_dhpc+1;
                    end 
                end
                
            end
        end
         
        count_dhpc= count_dhpc+size(group_dHPC,1);
        
        disp('vHPC Firing rate map calculation')
     
        for ii=1:size(group_vHPC,1)

            cluster = group_vHPC(ii,1);
            celltype = logical(cellulartype(cellulartype(:,1) == cluster,2));
            
            if celltype % check if pyr
                
                
                count_pyr_vhpc=count_pyr_vhpc+1; 
                spks = spks_vHPC(spks_vHPC(:,1)==cluster,2); % select tspk from cluster

                % --- Aversive ---
                spks_ave = spks; 
                %Restrict spk and position to laps and movement periods
                spks_ave = Restrict(spks_ave,in_lapA);spks_ave = Restrict(spks_ave , movement.aversive); % Restrict spikes to movement periods and laps
               
                % --- Reward ---
                spks_rew = spks; 
                %Restrict spk  to laps and movement periods
                spks_rew = Restrict(spks_rew,in_lapR); spks_rew = Restrict(spks_rew, movement.reward); % Restrict to movement periods
               
                                 
   
                %Checking pc criteria
                % 1. Mean firing rate above tresh.
                % 2.Skaggs > skaggs random. 
                % 3. At least 1 pf > 4 bins in one of the cond.  
     
                [curveA , statsA] = FiringCurve(pos_ave , spks_ave , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 6, 'minPeak' , 0.2);
                [curveR , statsR] = FiringCurve(pos_rew , spks_rew , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 6, 'minPeak' , 0.2);
                    
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
                
                
                %Compute Q90 skaggs
                qa = SkaggsRandomFMT_1D(spks_ave, pos_ave,sigma, Xedges, 0.90);
                qr = SkaggsRandomFMT_1D(spks_rew, pos_rew, sigma, Xedges, 0.90);
%                 disp(['qa= ', num2str(qa), 'skaggs=' num2str(statsA.specificity)]);
%                 disp(['qr= ', num2str(qr), 'skaggs=' num2str(statsR.specificity)]);

                % PC AVE 
                if statsA.specificity>qa % 2. criteria
                    if and(sum(statsA.field(:,:,1))>= 4 ,fr_ave>tresh) % 1. and 3. criteria
                        pc_ave_vhpc = pc_ave_vhpc+1;
                    end 
                end

                %PC REWARD 
                if statsR.specificity>qr % 2. criteria
                    if and(sum(statsR.field(:,:,1))>= 4 ,fr_rew>tresh) % 1. and 3. criteria
                        pc_rew_vhpc = pc_rew_vhpc+1;
                    end 
                end
                
            end
        end
         
        count_vhpc= count_vhpc+size(group_vHPC,1); 

    %% Saveing PC INFO 
    if ~isempty(dHPC)
        dHPC = dHPC(~cellfun('isempty',dHPC));
        save([cd,'\dHPC_pc_skaggs_pos.mat'],'dHPC'); 
    end
    if ~isempty(vHPC)
        vHPC = vHPC(~cellfun('isempty',vHPC));
        save([cd,'\vHPC_pc_skaggs_pos.mat'],'vHPC'); 
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
        
        %Restrict pos to inside maze and only movement periods
        %---Aversive----
        pos_ave = behavior.pos.aversive(:,1:2); 
        in_maze = ToIntervals(pos_ave(:,1),and(pos_ave(:,2)>35 , pos_ave(:,2)<175));% eliminating the extrems of the maze
                    
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
        
        pos_ave = Restrict(pos_ave,in_lapA);pos_ave = Restrict(pos_ave, movement.aversive); % Restrict pos to movement periods and laps
        pos_ave(:,2) = pos_ave(:,2)-min(pos_ave(:,2)); pos_ave(:,2) = pos_ave(:,2)/max(pos_ave(:,2)); %normalization of position           
       
         %---Rewarded----
        pos_rew = behavior.pos.reward(:,1:2);
        in_maze = ToIntervals(pos_rew(:,1),and(pos_rew(:,2)>35 , pos_rew(:,2)<175));% eliminating the extrems of the maze (15cm-reward box)
               
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
        pos_rew = Restrict(pos_rew,in_lapR); pos_rew = Restrict(pos_rew , movement.reward);
        pos_rew(:,2) = pos_rew(:,2)-min(pos_rew(:,2)); pos_rew(:,2) = pos_rew(:,2)/max(pos_rew(:,2)); %normalization of position 
        
        
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
        %Keep only PC
        try
            load('dHPC_pc_skaggs_circular.mat');
            dhpc_sub =[];
            for d=1:size(dHPC,2)
                dhpc_sub = [dhpc_sub;dHPC{d}.id];
            end 
            temp_id= ismember(group_dHPC(:,1),dhpc_sub); group_dHPC=group_dHPC(temp_id,:);
            
             for ii=1:size(group_dHPC,1)

            cluster = group_dHPC(ii,1);

                spks = spks_dHPC(spks_dHPC(:,1)==cluster,2); % select tspk from cluster

                % --- Aversive ---
                spks_ave = spks; 
                %Restrict spk and position to laps and movement periods
                spks_ave = Restrict(spks_ave,in_lapA);spks_ave = Restrict(spks_ave , movement.aversive); % Restrict spikes to movement periods and laps
                    
               
                % --- Reward ---
                spks_rew = spks; 
                %Restrict spk and position to laps and movement periods
                spks_rew = Restrict(spks_rew,in_lapR); spks_rew = Restrict(spks_rew, movement.reward); % Restrict to movement periods
               
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
                % Firirng curve whole session 
                [curveA , statsA] = FiringCurve(pos_ave , spks_ave , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                [curveR , statsR] = FiringCurve(pos_rew , spks_rew , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
 
                %Firing curve lap by lap 
                frMap_ave_lap = [];
                for ix=1:size(in_lapA,1)
                    spks_lap = Restrict(spks_ave,in_lapA(ix,:));
                    pos_lap = Restrict(pos_ave,in_lapA(ix,:));
                    [curveA , statsA] = FiringCurve(pos_lap, spks_lap , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    frMap_ave_lap = [frMap_ave_lap;curveA.rate];
                end
                
                frMap_rew_lap = [];
                for ix=1:size(in_lapR,1)
                    spks_lap = Restrict(spks_rew,in_lapR(ix,:));
                    pos_lap = Restrict(pos_rew,in_lapR(ix,:));
                    [curveR , statsR] = FiringCurve(pos_lap, spks_lap , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    frMap_rew_lap = [frMap_rew_lap;curveR.rate];
                end
                
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
                    n.frMap_lap.ave = frMap_ave_lap;
                    n.frMap_lap.rew = frMap_ave_rew;
                    n.stats_ave = statsA;
                    n.stats_rew = statsR;
                    n.between = between;
                    n.within.ave = withinA; 
                    n.within.rew = withinR;
       
                    dHPC_lap{ii}= n;
                   
           
            end
            
        catch 
            dHPC = []; disp(['No PC file in ',session]);
        end

        disp('vHPC Firing rate map calculation')
        
        %Keep only PC
        try
            load('vHPC_pc_skaggs_circular.mat');
            vhpc_sub =[];
            for d=1:size(vHPC,2)
                vhpc_sub = [vhpc_sub;vHPC{d}.id];
            end 
            temp_id= ismember(group_vHPC(:,1),vhpc_sub); group_vHPC=group_vHPC(temp_id,:);
            
             for ii=1:size(group_vHPC,1)

            cluster = group_vHPC(ii,1);
            
                
                spks = spks_vHPC(spks_vHPC(:,1)==cluster,2); % select tspk from cluster

                % --- Aversive ---
                spks_ave = spks; 
                pos_ave = behavior.pos.aversive(:,1:2); 
                in_maze = ToIntervals(pos_ave(:,1),and(pos_ave(:,2)>35 , pos_ave(:,2)<175));% eliminating the extrems of the maze
                    
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

     
                [curveA , statsA] = FiringCurve(pos_ave , spks_ave , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                [curveR , statsR] = FiringCurve(pos_rew , spks_rew , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    
               %Firing curve lap by lap 
                frMap_ave_lap = [];
                for ix=1:size(in_lapA,1)
                    spks_lap = Restrict(spks_ave,in_lapA(ix,:));
                    pos_lap = Restrict(pos_ave,in_lapA(ix,:));
                    [curveA , statsA] = FiringCurve(pos_lap, spks_lap , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    frMap_ave_lap = [frMap_ave_lap;curveA.rate];
                end
                
                frMap_rew_lap = [];
                for ix=1:size(in_lapR,1)
                    spks_lap = Restrict(spks_rew,in_lapR(ix,:));
                    pos_lap = Restrict(pos_rew,in_lapR(ix,:));
                    [curveR , statsR] = FiringCurve(pos_lap, spks_lap , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    frMap_rew_lap = [frMap_rew_lap;curveR.rate];
                end
                    
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
                    n.frMap_lap.ave = frMap_ave_lap;
                    n.frMap_lap.rew = frMap_ave_rew;
                    n.stats_ave = statsA;
                    n.stats_rew = statsR;
                    n.between = between;
                    n.within.ave = withinA; 
                    n.within.rew = withinR;
       
                    vHPC_lap{ii}= n;
                
         
            end     

        catch 
            vHPC = []; disp(['No PC file in ',session]);
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

 %% MAIN LOOP ODD/EVEN LAPS  

dhpc = []; 
vhpc = [];

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
        
        %Restrict pos to inside maze and only movement periods
        %---Aversive----
        pos_ave = behavior.pos.aversive(:,1:2); 
        in_maze = ToIntervals(pos_ave(:,1),and(pos_ave(:,2)>35 , pos_ave(:,2)<175));% eliminating the extrems of the maze
                    
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
        
        pos_ave = Restrict(pos_ave,in_lapA);pos_ave = Restrict(pos_ave, movement.aversive); % Restrict pos to movement periods and laps
        pos_ave(:,2) = pos_ave(:,2)-min(pos_ave(:,2)); pos_ave(:,2) = pos_ave(:,2)/max(pos_ave(:,2)); %normalization of position           
       
        laps_odd_ave =in_lapA(1:2:end,:);
        laps_even_ave =in_lapA(2:2:end,:);
        
        %---Rewarded----
        pos_rew = behavior.pos.reward(:,1:2);
        in_maze = ToIntervals(pos_rew(:,1),and(pos_rew(:,2)>35 , pos_rew(:,2)<175));% eliminating the extrems of the maze (15cm-reward box)
               
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
        pos_rew = Restrict(pos_rew,in_lapR); pos_rew = Restrict(pos_rew , movement.reward);
        pos_rew(:,2) = pos_rew(:,2)-min(pos_rew(:,2)); pos_rew(:,2) = pos_rew(:,2)/max(pos_rew(:,2)); %normalization of position 
        
        laps_odd_rew =in_lapR(1:2:end,:);
        laps_even_rew =in_lapR(2:2:end,:);
        
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
        %Keep only PC
        try
            load('dHPC_pc_skaggs_circular.mat');
            dhpc_sub =[];
            for d=1:size(dHPC,2)
                dhpc_sub = [dhpc_sub;dHPC{d}.id];
            end 
            temp_id= ismember(group_dHPC(:,1),dhpc_sub); group_dHPC=group_dHPC(temp_id,:);
            
             for ii=1:size(group_dHPC,1)

                cluster = group_dHPC(ii,1);
                spks = spks_dHPC(spks_dHPC(:,1)==cluster,2); % select tspk from cluster
                % --- Aversive ---
                spks_ave = spks; 
                %Restrict spk and position to laps and movement periods
                spks_ave = Restrict(spks_ave,in_lapA);spks_ave = Restrict(spks_ave , movement.aversive); % Restrict spikes to movement periods and laps
                % --- Reward ---
                spks_rew = spks; 
                %Restrict spk and position to laps and movement periods
                spks_rew = Restrict(spks_rew,in_lapR); spks_rew = Restrict(spks_rew, movement.reward); % Restrict to movement periods


                %%%%%%%%% Between %%%%%%%%%%%% 
                [curve1,stats1] = FiringCurve(pos_ave, spks_ave , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4 , 'minPeak' , 0.2,'minTime',0.15);
                [curve2,stats2] = FiringCurve(pos_rew, spks_rew , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4 , 'minPeak' , 0.2,'minTime',0.15);
                %%%
%                 figure(1);clf;hold on; subplot(1,2,1);imagesc(curve1.rate);subplot(1,2,2);imagesc(curve2.rate);
                %%%
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
                
                dhpc = [dhpc;spatial,fr_change,overlap, shift,1];
                
                %%%%%%%%%%% WITHIN Aversive %%%%%%%%%%%%%%
                spks_odd = Restrict(spks_ave,laps_odd_ave);
                pos_odd = Restrict(pos_ave,laps_odd_ave); 
                
                spks_even = Restrict(spks_ave,laps_even_ave);
                pos_even = Restrict(pos_ave,laps_even_ave); 
                
                [curve1,stats1] = FiringCurve(pos_odd, spks_odd , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4 , 'minPeak' , 0.2,'minTime',0.15);
                [curve2,stats2] = FiringCurve(pos_even, spks_even , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4 , 'minPeak' , 0.2,'minTime',0.15);
                %%%
%                 figure(1);clf;hold on; subplot(1,2,1);imagesc(curve1.rate);subplot(1,2,2);imagesc(curve2.rate);
                %%%
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
                
                dhpc = [dhpc;spatial,fr_change,overlap, shift,2];
                
                %%%%%%%%% WITHIN Reward %%%%%%%%%%
                spks_odd = Restrict(spks_rew,laps_odd_rew);
                pos_odd = Restrict(pos_rew,laps_odd_rew); 
                
                spks_even = Restrict(spks_rew,laps_even_rew);
                pos_even = Restrict(pos_rew,laps_even_rew); 
                
                [curve1,stats1] = FiringCurve(pos_odd, spks_odd , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4 , 'minPeak' , 0.2,'minTime',0.15);
                [curve2,stats2] = FiringCurve(pos_even, spks_even , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4 , 'minPeak' , 0.2,'minTime',0.15);
                %%%
%                 figure(1);clf;hold on; subplot(1,2,1);imagesc(curve1.rate);subplot(1,2,2);imagesc(curve2.rate);
                %%%
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
                
                dhpc = [dhpc;spatial,fr_change,overlap, shift,3];
               
            end
            
        catch 
            dHPC = []; disp(['No PC file in ',session]);
        end

        
        disp('vHPC Firing rate map calculation')
        %Keep only PC
        try
            load('vHPC_pc_skaggs_circular.mat');
            vhpc_sub =[];
            for d=1:size(vHPC,2)
                vhpc_sub = [vhpc_sub;vHPC{d}.id];
            end 
            temp_id= ismember(group_vHPC(:,1),vhpc_sub); group_vHPC=group_vHPC(temp_id,:);
            
             for ii=1:size(group_vHPC,1)

                cluster = group_vHPC(ii,1);
                spks = spks_vHPC(spks_vHPC(:,1)==cluster,2); % select tspk from cluster
                % --- Aversive ---
                spks_ave = spks; 
                %Restrict spk and position to laps and movement periods
                spks_ave = Restrict(spks_ave,in_lapA);spks_ave = Restrict(spks_ave , movement.aversive); % Restrict spikes to movement periods and laps
                % --- Reward ---
                spks_rew = spks; 
                %Restrict spk and position to laps and movement periods
                spks_rew = Restrict(spks_rew,in_lapR); spks_rew = Restrict(spks_rew, movement.reward); % Restrict to movement periods
  
                 % Between 
                [curve1,stats1] = FiringCurve(pos_ave, spks_ave , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4 , 'minPeak' , 0.2,'minTime',0.15);
                [curve2,stats2] = FiringCurve(pos_rew, spks_rew , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4 , 'minPeak' , 0.2,'minTime',0.15);
                %%%
%                 figure(1);clf;hold on; subplot(1,2,1);imagesc(curve1.rate);subplot(1,2,2);imagesc(curve2.rate);
                %%%
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
                
                vhpc = [vhpc;spatial,fr_change,overlap, shift,1];
                
                %Aversive
                spks_odd = Restrict(spks_ave,laps_odd_ave);
                pos_odd = Restrict(pos_ave,laps_odd_ave); 
                
                spks_even = Restrict(spks_ave,laps_even_ave);
                pos_even = Restrict(pos_ave,laps_even_ave); 
                
                [curve1,stats1] = FiringCurve(pos_odd, spks_odd , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4 , 'minPeak' , 0.2,'minTime',0.15);
                [curve2,stats2] = FiringCurve(pos_even, spks_even , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4 , 'minPeak' , 0.2,'minTime',0.15);
                %%%
%                 figure(1);clf;hold on; subplot(1,2,1);imagesc(curve1.rate);subplot(1,2,2);imagesc(curve2.rate);
                %%%
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
                
                vhpc = [vhpc;spatial,fr_change,overlap, shift,2];
                
                %Reward
                spks_odd = Restrict(spks_rew,laps_odd_rew);
                pos_odd = Restrict(pos_rew,laps_odd_rew); 
                
                spks_even = Restrict(spks_rew,laps_even_rew);
                pos_even = Restrict(pos_rew,laps_even_rew); 
                
                [curve1,stats1] = FiringCurve(pos_odd, spks_odd , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4 , 'minPeak' , 0.2,'minTime',0.15);
                [curve2,stats2] = FiringCurve(pos_even, spks_even , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4 , 'minPeak' , 0.2,'minTime',0.15);
                %%%
%                 figure(1);clf;hold on; subplot(1,2,1);imagesc(curve1.rate);subplot(1,2,2);imagesc(curve2.rate);
                %%%
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
                
                vhpc = [vhpc;spatial,fr_change,overlap, shift,3];
               
         
            end     

        catch 
            vHPC = []; disp(['No PC file in ',session]);
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


%Only 3 plots
ylabels ={'Spatial correlation', 'Overlap', 'Pf shift'};
idx = [1,3,4]; 
xlabels = {'Between', 'Within Ave', 'Within Rew'}; 
figure(1);clf;hold on, 
sgtitle('dHPC')
for a=1:size(ylabels,2)
    p = idx(a);
    subplot(1,3,a); hold on; 
    x = dhpc(:,5);
    y= dhpc(:,p);%Change the column number (1-4) to choose which variable to plot 
    c = [.3,.3,.3];
    scatter(x,y,30,c,"filled",'jitter','on', 'jitterAmount',0.3); xlim([0 4]);
    ylabel(ylabels{a});
    xticks([1 2 3])
    xticklabels({'Bet', 'WA', 'WR'});
    hold on
    s1=scatter(1,nanmedian(dhpc(dhpc(:,5)==1,p)), "filled");s1.MarkerFaceColor = [0 0 0];
    s2=scatter(2,nanmedian(dhpc(dhpc(:,5)==2,p)), "filled");s2.MarkerFaceColor = [1 0.1 0.2];
    s3=scatter(3,nanmedian(dhpc(dhpc(:,5)==3,p)), "filled");s3.MarkerFaceColor = [0.1 0.3 1];
end

%Stats
[P,ANOVATAB,STATS] = kruskalwallis(dhpc_sub(:,5),dhpc_sub(:,6));
c = multcompare(STATS)

tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])

% Paired non-parametric anova - Friedman test
%Prepare data to test
p = 3; % 1=sp corr 2=fr change 3=overlap 4=pf shift
data=[dhpc(dhpc(:,5)==1,p), dhpc(dhpc(:,5)==2,p), ...
    dhpc(dhpc(:,5)==3,p), dhpc(dhpc(:,5)==4,p)];
data(any(isnan(data), 2), :) = [];

[p,~,stats] =  friedman(data,1); 
c = multcompare(stats);
tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])

%Only 3 plots
ylabels ={'Spatial correlation', 'Overlap', 'Pf shift'};
idx = [1,3,4]; 
xlabels = {'Between', 'Within Ave', 'Within Rew'}; 
figure(2);clf;hold on, 
sgtitle('vHPC')
for a=1:size(ylabels,2)
    p = idx(a);
    subplot(1,3,a); hold on; 
    x = vhpc(:,5);
    y= vhpc(:,p);%Change the column number (1-4) to choose which variable to plot 
    c = [.3,.3,.3];
    scatter(x,y,30,c,"filled",'jitter','on', 'jitterAmount',0.3); xlim([0 4]);
    ylabel(ylabels{a});
    xticks([1 2 3])
    xticklabels({'Bet', 'WA', 'WR'});
    hold on
    s1=scatter(1,nanmedian(vhpc(vhpc(:,5)==1,p)), "filled");s1.MarkerFaceColor = [0 0 0];
    s2=scatter(2,nanmedian(vhpc(vhpc(:,5)==2,p)), "filled");s2.MarkerFaceColor = [1 0.1 0.2];
    s3=scatter(3,nanmedian(vhpc(vhpc(:,5)==3,p)), "filled");s3.MarkerFaceColor = [0.1 0.3 1];
end

%Stats
[P,ANOVATAB,STATS] = kruskalwallis(vhpc(:,1),vhpc(:,5));
c = multcompare(STATS);tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])

% Paired non-parametric anova - Friedman test
%Prepare data to test
p = 3; % 1=sp corr 2=fr change 3=overlap 4=pf shift
data=[vhpc(vhpc(:,5)==1,p), vhpc(vhpc(:,5)==2,p), ...
    vhpc(vhpc(:,5)==3,p), vhpc(vhpc(:,5)==4,p)];
data(any(isnan(data), 2), :) = [];

[p,~,stats] =  friedman(data,1); 
c = multcompare(stats);tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])
%% REMAPPING SPLITING ODD/EVEN
%Outputs
dhpc_odd = []; dhpc_even = []; 
vhpc_odd = []; vhpc_even= [];

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
        
        %Restrict pos to inside maze and only movement periods
        %---Aversive----
        pos_ave = behavior.pos.aversive(:,1:2); 
        in_maze = ToIntervals(pos_ave(:,1),and(pos_ave(:,2)>35 , pos_ave(:,2)<175));% eliminating the extrems of the maze
                    
       %Check lenght of in_maze segments to define real laps
        in_lapA = []; % keeps only full laps (maybe too restrictive)
        for ix=1:size(in_maze,1)
                        ini = in_maze(ix,1);
                        fin = in_maze(ix,2);
                        ini_pos = pos_ave(pos_ave(:,1)==ini,2); 
                        fin_pos = pos_ave(pos_ave(:,1)==fin,2); 
                        if abs(fin_pos-ini_pos)>100
                             if fin_pos-ini_pos>0
                                 in_lapA = [in_lapA;in_maze(ix,:),1]; 
                             else
                                 in_lapA = [in_lapA;in_maze(ix,:),2]; 
                             end 
                        end
        end
        
        pos_ave = Restrict(pos_ave,in_lapA(:,1:2));pos_ave = Restrict(pos_ave, movement.aversive); % Restrict pos to movement periods and laps
        pos_ave(:,2) = pos_ave(:,2)-min(pos_ave(:,2)); pos_ave(:,2) = pos_ave(:,2)/max(pos_ave(:,2)); %normalization of position           
       
        
        laps_odd_ave =in_lapA(in_lapA(:,3)==1,1:2);
        laps_even_ave =in_lapA(in_lapA(:,3)==2,1:2);
        
        %---Rewarded----
        pos_rew = behavior.pos.reward(:,1:2);
        in_maze = ToIntervals(pos_rew(:,1),and(pos_rew(:,2)>35 , pos_rew(:,2)<175));% eliminating the extrems of the maze (15cm-reward box)
               
        %Check lenght of in_maze segments to define real laps
        in_lapR = []; 
        for ix=1:size(in_maze,1)
                        ini = in_maze(ix,1);
                        fin = in_maze(ix,2);
                        ini_pos = pos_rew(pos_rew(:,1)==ini,2); 
                        fin_pos = pos_rew(pos_rew(:,1)==fin,2); 
                        if abs(fin_pos-ini_pos)>100
                            if fin_pos-ini_pos>0
                                 in_lapR = [in_lapR;in_maze(ix,:),1]; 
                             else
                                 in_lapR = [in_lapR;in_maze(ix,:),2]; 
                             end  
                        end
        end 
        pos_rew = Restrict(pos_rew,in_lapR(:,1:2)); pos_rew = Restrict(pos_rew , movement.reward);
        pos_rew(:,2) = pos_rew(:,2)-min(pos_rew(:,2)); pos_rew(:,2) = pos_rew(:,2)/max(pos_rew(:,2)); %normalization of position 
        
        laps_odd_rew =in_lapR(in_lapR(:,3)==1,1:2);
        laps_even_rew =in_lapR(in_lapR(:,3)==2,1:2);
        
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
        if  and(size(in_lapR,1)>=10,size(in_lapA,1)>=10)
        %Keep only PC
        try
            load('dHPC_pc_skaggs_circular.mat');
            
            dhpc_sub =[];
            for d=1:size(dHPC,2)
                dhpc_sub = [dhpc_sub;dHPC{d}.id];
            end 
            temp_id= ismember(group_dHPC(:,1),dhpc_sub); group_dHPC=group_dHPC(temp_id,:);
          
             for ii=1:size(group_dHPC,1)

                cluster = group_dHPC(ii,1);
                spks = spks_dHPC(spks_dHPC(:,1)==cluster,2); % select tspk from cluster
                % --- Aversive ---
                spks_ave = spks; 
                %Restrict spk and position to laps and movement periods
                spks_ave = Restrict(spks_ave,in_lapA(:,1:2));spks_ave = Restrict(spks_ave , movement.aversive); % Restrict spikes to movement periods and laps
                % --- Reward ---
                spks_rew = spks; 
                %Restrict spk and position to laps and movement periods
                spks_rew = Restrict(spks_rew,in_lapR(:,1:2)); spks_rew = Restrict(spks_rew, movement.reward); % Restrict to movement periods
                
                [between_corr] = remapping_laps_between(laps_odd_ave,pos_ave,spks_ave,laps_odd_rew, pos_rew,spks_rew, sigma, Xedges,min_size, min_peak, min_time);
                bet = nanmean(between_corr,"all");
                dhpc_odd = [dhpc_odd;bet,1];
                dhpc_even = [dhpc_even;bet,1];
                
                [within_odd_ave] = remapping_laps_within(laps_odd_ave ,pos_ave,spks_ave, sigma, Xedges,min_size, min_peak, min_time);
                 within = tril(within_odd_ave,-1);within(within==0) = NaN;within=nanmean(within,'all');
                 dhpc_odd = [dhpc_odd; within,2];
                 
                [within_even_ave] = remapping_laps_within(laps_even_ave ,pos_ave,spks_ave, sigma, Xedges,min_size, min_peak, min_time);
                within = tril(within_even_ave,-1);within(within==0) = NaN;within=nanmean(within,'all');
                dhpc_even = [dhpc_even;within,2];
                
                [within_odd_rew] = remapping_laps_within(laps_odd_rew ,pos_rew,spks_rew,sigma,Xedges,min_size,min_peak,min_time);
                within = tril(within_odd_rew,-1);within(within==0) = NaN;within=nanmean(within,'all');
                 dhpc_odd = [dhpc_odd; within,3];
                 
                [within_even_rew] = remapping_laps_within(laps_even_rew ,pos_rew,spks_rew,sigma,Xedges,min_size,min_peak,min_time);
                within = tril(within_even_rew,-1);within(within==0) = NaN;within=nanmean(within,'all');
                dhpc_even = [dhpc_even;within,3];
           
            end
            
        catch 
             disp(['No PC file in ',session]);
        end

        disp('vHPC Firing rate map calculation')
        
        %Keep only PC
        try
            load('vHPC_pc_skaggs_circular.mat');
            vhpc_sub =[];
            for d=1:size(vHPC,2)
                vhpc_sub = [vhpc_sub;vHPC{d}.id];
            end 
            temp_id= ismember(group_vHPC(:,1),vhpc_sub); group_vHPC=group_vHPC(temp_id,:);
            
             for ii=1:size(group_vHPC,1)

                cluster = group_vHPC(ii,1);
                spks = spks_vHPC(spks_vHPC(:,1)==cluster,2); % select tspk from cluster
                % --- Aversive ---
                spks_ave = spks; 
                %Restrict spk and position to laps and movement periods
                spks_ave = Restrict(spks_ave,in_lapA(:,1:2));spks_ave = Restrict(spks_ave , movement.aversive); % Restrict spikes to movement periods and laps
                % --- Reward ---
                spks_rew = spks; 
                %Restrict spk and position to laps and movement periods
                spks_rew = Restrict(spks_rew,in_lapR(:,1:2)); spks_rew = Restrict(spks_rew, movement.reward); % Restrict to movement periods
  
                [between_corr] = remapping_laps_between(laps_odd_ave,pos_ave,spks_ave,laps_odd_rew, pos_rew,spks_rew, sigma, Xedges,min_size, min_peak, min_time);
                bet = nanmean(between_corr,"all");
                vhpc_odd = [vhpc_odd;bet,1];
                vhpc_even = [vhpc_even;bet,1];
                
                [within_odd_ave] = remapping_laps_within(laps_odd_ave ,pos_ave,spks_ave, sigma, Xedges,min_size, min_peak, min_time);
                 within = tril(within_odd_ave,-1);within(within==0) = NaN;within=nanmean(within,'all');
                 vhpc_odd = [vhpc_odd; within,2];
                 
                [within_even_ave] = remapping_laps_within(laps_even_ave ,pos_ave,spks_ave, sigma, Xedges,min_size, min_peak, min_time);
                within = tril(within_even_ave,-1);within(within==0) = NaN;within=nanmean(within,'all');
                vhpc_even = [vhpc_even;within,2];
                
                [within_odd_rew] = remapping_laps_within(laps_odd_rew ,pos_rew,spks_rew,sigma,Xedges,min_size,min_peak,min_time);
                within = tril(within_odd_rew,-1);within(within==0) = NaN;within=nanmean(within,'all');
                vhpc_odd = [vhpc_odd; within,3];
                 
                [within_even_rew] = remapping_laps_within(laps_even_rew ,pos_rew,spks_rew,sigma,Xedges,min_size,min_peak,min_time);
                within = tril(within_even_rew,-1);within(within==0) = NaN;within=nanmean(within,'all');
                vhpc_even = [vhpc_even;within,3];

            end     

        catch 
            vHPC = []; disp(['No PC file in ',session]);
        end
        
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

%Plots
odd = vhpc_odd;
even = vhpc_even;
c = [.3,.3,.3];
figure(1);clf;hold on;sgtitle('dHPC')
subplot(1,2,1); hold on;title('Odd'); 
scatter(odd(:,2),odd(:,1),30,c,"filled",'jitter','on', 'jitterAmount',0.3,'MarkerFaceAlpha',.5); xlim([0 4]); ylim([-0.2 1]);
ylabel('Spatial correlation');xticks([1 2 3]);xticklabels({'Bet', 'WA', 'WR'});hold on;
s1=scatter(1,nanmedian(odd(odd(:,2)==1,1)), "filled");s1.MarkerFaceColor = [0 0 0];
s2=scatter(2,nanmedian(odd(odd(:,2)==2,1)), "filled");s2.MarkerFaceColor = [1 0.1 0.2];
s3=scatter(3,nanmedian(odd(odd(:,2)==3,1)), "filled");s3.MarkerFaceColor = [0.1 0.3 1];

subplot(1,2,2); hold on;title('Even'); 
scatter(even(:,2),even(:,1),30,c,"filled",'jitter','on', 'jitterAmount',0.3,'MarkerFaceAlpha',.5); xlim([0 4]); ylim([-0.2 1]);
ylabel('Spatial correlation');xticks([1 2 3]);xticklabels({'Bet', 'WA', 'WR'});hold on;
s1=scatter(1,nanmedian(even(even(:,2)==1,1)), "filled");s1.MarkerFaceColor = [0 0 0];
s2=scatter(2,nanmedian(even(even(:,2)==2,1)), "filled");s2.MarkerFaceColor = [1 0.1 0.2];
s3=scatter(3,nanmedian(even(even(:,2)==3,1)), "filled");s3.MarkerFaceColor = [0.1 0.3 1];

%Stats
[P,ANOVATAB,STATS] = kruskalwallis(dhpc_sub(:,5),dhpc_sub(:,6));
c = multcompare(STATS);tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])
% Paired non-parametric anova - Friedman test
dhpc= vhpc_even; %change for the output you want 
data=[dhpc(dhpc(:,2)==1,1), dhpc(dhpc(:,2)==2,1), ...
    dhpc(dhpc(:,2)==3,1)];data(any(isnan(data), 2), :) = [];
[p,~,stats] =  friedman(data,1);c = multcompare(stats);tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])

%% RATEMAP DIRECTIONALITY PLOTS
dhpc_odd_ave = []; dhpc_odd_rew = [];dhpc_even_ave = [];dhpc_even_rew = [];  
vhpc_odd_ave = []; vhpc_odd_rew = [];vhpc_even_ave = [];vhpc_even_rew = [];  

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
        
        %Restrict pos to inside maze and only movement periods
        %---Aversive----
        pos_ave = behavior.pos.aversive(:,1:2); 
        in_maze = ToIntervals(pos_ave(:,1),and(pos_ave(:,2)>35 , pos_ave(:,2)<175));% eliminating the extrems of the maze
                    
       %Check lenght of in_maze segments to define real laps
        in_lapA = []; % keeps only full laps (maybe too restrictive)
        for ix=1:size(in_maze,1)
                        ini = in_maze(ix,1);
                        fin = in_maze(ix,2);
                        ini_pos = pos_ave(pos_ave(:,1)==ini,2); 
                        fin_pos = pos_ave(pos_ave(:,1)==fin,2); 
                        if abs(fin_pos-ini_pos)>100
                             if fin_pos-ini_pos>0
                                 in_lapA = [in_lapA;in_maze(ix,:),1]; 
                             else
                                 in_lapA = [in_lapA;in_maze(ix,:),2]; 
                             end 
                        end
        end
        
        pos_ave = Restrict(pos_ave,in_lapA(:,1:2));pos_ave = Restrict(pos_ave, movement.aversive); % Restrict pos to movement periods and laps
        pos_ave(:,2) = pos_ave(:,2)-min(pos_ave(:,2)); pos_ave(:,2) = pos_ave(:,2)/max(pos_ave(:,2)); %normalization of position           
       
        
        laps_odd_ave =in_lapA(in_lapA(:,3)==1,1:2);
        laps_even_ave =in_lapA(in_lapA(:,3)==2,1:2);
        
        %---Rewarded----
        pos_rew = behavior.pos.reward(:,1:2);
        in_maze = ToIntervals(pos_rew(:,1),and(pos_rew(:,2)>35 , pos_rew(:,2)<175));% eliminating the extrems of the maze (15cm-reward box)
               
        %Check lenght of in_maze segments to define real laps
        in_lapR = []; 
        for ix=1:size(in_maze,1)
                        ini = in_maze(ix,1);
                        fin = in_maze(ix,2);
                        ini_pos = pos_rew(pos_rew(:,1)==ini,2); 
                        fin_pos = pos_rew(pos_rew(:,1)==fin,2); 
                        if abs(fin_pos-ini_pos)>100
                            if fin_pos-ini_pos>0
                                 in_lapR = [in_lapR;in_maze(ix,:),1]; 
                             else
                                 in_lapR = [in_lapR;in_maze(ix,:),2]; 
                             end  
                        end
        end 
        pos_rew = Restrict(pos_rew,in_lapR(:,1:2)); pos_rew = Restrict(pos_rew , movement.reward);
        pos_rew(:,2) = pos_rew(:,2)-min(pos_rew(:,2)); pos_rew(:,2) = pos_rew(:,2)/max(pos_rew(:,2)); %normalization of position 
        
        laps_odd_rew =in_lapR(in_lapR(:,3)==1,1:2);
        laps_even_rew =in_lapR(in_lapR(:,3)==2,1:2);
        
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
       
        %Keep only PC
        try
            load('dHPC_pc_skaggs_circular.mat');
            
            dhpc_sub =[];
            for d=1:size(dHPC,2)
                dhpc_sub = [dhpc_sub;dHPC{d}.id];
            end 
            temp_id= ismember(group_dHPC(:,1),dhpc_sub); group_dHPC=group_dHPC(temp_id,:);
          
             for ii=1:size(group_dHPC,1)

                cluster = group_dHPC(ii,1);
                spks = spks_dHPC(spks_dHPC(:,1)==cluster,2); % select tspk from cluster
                % --- Aversive ---
                spks_ave = spks; 
                %Restrict spk and position to laps and movement periods
                spks_ave = Restrict(spks_ave,in_lapA(:,1:2));spks_ave = Restrict(spks_ave , movement.aversive); % Restrict spikes to movement periods and laps
                % --- Reward ---
                spks_rew = spks; 
                %Restrict spk and position to laps and movement periods
                spks_rew = Restrict(spks_rew,in_lapR(:,1:2)); spks_rew = Restrict(spks_rew, movement.reward); % Restrict to movement periods
                
                %ODD
                spks = Restrict(spks_ave, laps_odd_ave);
                pos = Restrict(pos_ave, laps_odd_ave); 
                [curve1] = FiringCurve(pos, spks , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , min_size , 'minPeak' , min_peak,'minTime',min_time);
                A = curve1.rate - min(curve1.rate);A = A ./ max(A);
                dhpc_odd_ave = [dhpc_odd_ave; A];% saving stack of lap rate maps 
                
                spks = Restrict(spks_rew, laps_odd_rew);
                pos = Restrict(pos_rew, laps_odd_rew); 
                [curve1] = FiringCurve(pos, spks , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , min_size , 'minPeak' , min_peak,'minTime',min_time);
                A = curve1.rate - min(curve1.rate);A = A ./ max(A);
                dhpc_odd_rew = [dhpc_odd_rew; A];% saving stack of lap rate maps 
                
                %EVEN
                spks = Restrict(spks_ave, laps_even_ave);
                pos = Restrict(pos_ave, laps_even_ave); 
                [curve1] = FiringCurve(pos, spks , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , min_size , 'minPeak' , min_peak,'minTime',min_time);
                A = curve1.rate - min(curve1.rate);A = A ./ max(A);
                dhpc_even_ave = [dhpc_even_ave; A];% saving stack of lap rate maps 
                
                spks = Restrict(spks_rew, laps_even_rew);
                pos = Restrict(pos_rew, laps_even_rew); 
                [curve1] = FiringCurve(pos, spks , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , min_size , 'minPeak' , min_peak,'minTime',min_time);
                A = curve1.rate - min(curve1.rate);A = A ./ max(A);
                dhpc_even_rew = [dhpc_even_rew; A];% saving stack of lap rate maps 
                
            end
            
        catch 
             disp(['No PC file in ',session]);
        end

        disp('vHPC Firing rate map calculation')
        
        %Keep only PC
        try
            load('vHPC_pc_skaggs_circular.mat');
            vhpc_sub =[];
            for d=1:size(vHPC,2)
                vhpc_sub = [vhpc_sub;vHPC{d}.id];
            end 
            temp_id= ismember(group_vHPC(:,1),vhpc_sub); group_vHPC=group_vHPC(temp_id,:);
            
             for ii=1:size(group_vHPC,1)

                cluster = group_vHPC(ii,1);
                spks = spks_vHPC(spks_vHPC(:,1)==cluster,2); % select tspk from cluster
                % --- Aversive ---
                spks_ave = spks; 
                %Restrict spk and position to laps and movement periods
                spks_ave = Restrict(spks_ave,in_lapA(:,1:2));spks_ave = Restrict(spks_ave , movement.aversive); % Restrict spikes to movement periods and laps
                % --- Reward ---
                spks_rew = spks; 
                %Restrict spk and position to laps and movement periods
                spks_rew = Restrict(spks_rew,in_lapR(:,1:2)); spks_rew = Restrict(spks_rew, movement.reward); % Restrict to movement periods
  
                %ODD
                spks = Restrict(spks_ave, laps_odd_ave);
                pos = Restrict(pos_ave, laps_odd_ave); 
                [curve1] = FiringCurve(pos, spks , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , min_size , 'minPeak' , min_peak,'minTime',min_time);
                A = curve1.rate - min(curve1.rate);A = A ./ max(A);
                vhpc_odd_ave = [vhpc_odd_ave; A];% saving stack of lap rate maps 
                
                spks = Restrict(spks_rew, laps_odd_rew);
                pos = Restrict(pos_rew, laps_odd_rew); 
                [curve1] = FiringCurve(pos, spks , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , min_size , 'minPeak' , min_peak,'minTime',min_time);
                A = curve1.rate - min(curve1.rate);A = A ./ max(A);
                vhpc_odd_rew = [vhpc_odd_rew; A];% saving stack of lap rate maps 
                
                %EVEN
                spks = Restrict(spks_ave, laps_even_ave);
                pos = Restrict(pos_ave, laps_even_ave); 
                [curve1] = FiringCurve(pos, spks , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , min_size , 'minPeak' , min_peak,'minTime',min_time);
                A = curve1.rate - min(curve1.rate);A = A ./ max(A);
                vhpc_even_ave = [vhpc_even_ave; A];% saving stack of lap rate maps 
                
                spks = Restrict(spks_rew, laps_even_rew);
                pos = Restrict(pos_rew, laps_even_rew); 
                [curve1] = FiringCurve(pos, spks , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , min_size , 'minPeak' , min_peak,'minTime',min_time);
                A = curve1.rate - min(curve1.rate);A = A ./ max(A);
                vhpc_even_rew = [vhpc_even_rew; A];% saving stack of lap rate maps 

            end     

        catch 
            vHPC = []; disp(['No PC file in ',session]);
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


%%%%% Ratemap Plots %%%%%
% chaange vhpc or dhpc 
[h idx] = max (dhpc_odd_ave, [],2);
[m mm] = sort(idx); 

figure(2);clf;hold on;
fr = dhpc_odd_ave(mm,:);  fr(~any(~isnan(fr), 2),:)=[];
subplot(1,4,1);imagesc([0:20:140], [1:1:size(dhpc_odd_ave,1)],fr), colormap 'gray'; title('Aversive');
fr = dhpc_odd_rew(mm,:);  fr(~any(~isnan(fr), 2),:)=[];
subplot(1,4,2);imagesc([0:20:140], [1:1:size(dhpc_odd_rew,1)],fr), colormap 'gray'; title('Reward');
sgtitle('dHPC firing maps ODD - EVEN sorted by ODD AVE');

[h idx] = max (vhpc_even_ave, [],2);
[m mm] = sort(idx); 

fr = dhpc_even_ave(mm,:);  fr(~any(~isnan(fr), 2),:)=[];
subplot(1,4,3);imagesc([0:20:140], [1:1:size(dhpc_even_ave,1)],fr), colormap 'gray'; title('Aversive');
fr = dhpc_even_rew(mm,:);  fr(~any(~isnan(fr), 2),:)=[];
subplot(1,4,4);imagesc([0:20:140], [1:1:size(dhpc_even_rew,1)],fr), colormap 'gray'; title('Reward');


%%%%% Spatial correlation plot %%%%%

corr_ave_dhpc = []; corr_rew_dhpc = []; 
corr_ave_vhpc = []; corr_rew_vhpc = []; 

for i=1:size(dhpc_odd_ave,1)
    sp = corrcoef(dhpc_odd_ave(i,:),dhpc_even_ave(i,:)); 
    corr_ave_dhpc= [ corr_ave_dhpc;sp(1,2),1]; 
    
    sp = corrcoef(dhpc_odd_rew(i,:),dhpc_even_rew(i,:)); 
    corr_rew_dhpc= [ corr_rew_dhpc;sp(1,2),2]; 
end


for i=1:size(vhpc_odd_ave,1)
    sp = corrcoef(vhpc_odd_ave(i,:),vhpc_even_ave(i,:)); 
    corr_ave_vhpc= [ corr_ave_vhpc;sp(1,2),1]; 
    
    sp = corrcoef(vhpc_odd_rew(i,:),vhpc_even_rew(i,:)); 
    corr_rew_vhpc= [ corr_rew_vhpc;sp(1,2),2]; 
end

%%% dhpc %%
figure(3);clf;hold on
subplot(1,2,1); hold on; 
data = [corr_ave_dhpc;corr_rew_dhpc];
x = data(:,2);
y = data(:,1); 
c = [.3,.3,.3];
scatter(x,y,[],c,"filled",'jitter','on', 'jitterAmount',0.1); xlim([0.5 2.5]);
ylabel('Spatial corr even-odd');
xticks([1 2])
xticklabels({'Ave', 'Rew'});
hold on
s1=scatter(1,nanmedian(corr_ave_dhpc(:,1)), "filled");s1.MarkerFaceColor = [1 0.1 0.2];
s2=scatter(2,nanmedian(corr_rew_dhpc(:,1)), "filled");s2.MarkerFaceColor = [0 0.3 1];

%Stats
[h,p] = ttest(corr_ave_dhpc(:,1),corr_rew_dhpc(:,1)); 

subplot(1,2,2); hold on; 
data = [corr_ave_vhpc;corr_rew_vhpc];
x = data(:,2);
y = data(:,1); 
c = [.3,.3,.3];
scatter(x,y,[],c,"filled",'jitter','on', 'jitterAmount',0.1); xlim([0.5 2.5]);
ylabel('Spatial corr even-odd');
xticks([1 2])
xticklabels({'Ave', 'Rew'});
hold on
s1=scatter(1,nanmedian(corr_ave_vhpc(:,1)), "filled");s1.MarkerFaceColor = [1 0.1 0.2];
s2=scatter(2,nanmedian(corr_rew_vhpc(:,1)), "filled");s2.MarkerFaceColor = [0 0.3 1];

[h,p] = ttest(corr_ave_vhpc(:,1),corr_rew_vhpc(:,1))
   
%% REMAPPING PARAMETERS PLOTS and STATS 
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
            load('dHPC_pc_skaggs_circular.mat');%'dHPC_pc_skaggs_pos.mat'
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
            load('vHPC_pc_skaggs_circular.mat');%'vHPC_pc_skaggs_pos.mat'
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




%%%%% vhpc vs dhpc comparison 

%Only 3 plots
ylabels ={'Spatial correlation', 'Overlap', 'Pf shift'};
figure(1);clf;hold on, 
sgtitle('dHPC vs vHPC')
subplot(1,3,1); hold on; 
x =[ones(size(dhpc_sub(dhpc_sub(:,6)==1,1),1),1);2*ones(size(vhpc_sub(vhpc_sub(:,6)==1,1),1),1)]; 
y= [dhpc_sub(dhpc_sub(:,6)==1,1);vhpc_sub(vhpc_sub(:,6)==1,1)]; 
c = [.3,.3,.3];
scatter(x,y,30,c,"filled",'jitter','on', 'jitterAmount',0.3); xlim([0 3]);
ylabel('Spatial correlation');
xticks([1 2]);
xticklabels({'dHPC', 'vHPC'});
hold on
s1=scatter(1,nanmedian(dhpc_sub(dhpc_sub(:,6)==1,1)), "filled");s1.MarkerFaceColor = [0 1 0.5];
s2=scatter(2,nanmedian(vhpc_sub(vhpc_sub(:,6)==1,1)), "filled");s2.MarkerFaceColor = [0 0.5 1];
  

[p,tbl,stats] = kruskalwallis(y,x);
c = multcompare(stats);
tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])


subplot(1,3,2); hold on; 
x =[ones(size(dhpc_sub(dhpc_sub(:,6)==1,1),1),1);2*ones(size(vhpc_sub(vhpc_sub(:,6)==1,1),1),1)]; 
y= [dhpc_sub(dhpc_sub(:,6)==1,3);vhpc_sub(vhpc_sub(:,6)==1,3)]; 
c = [.3,.3,.3];
scatter(x,y,30,c,"filled",'jitter','on', 'jitterAmount',0.3); xlim([0 3]);
ylabel('Overlap');
xticks([1 2]);
xticklabels({'dHPC', 'vHPC'});
hold on
s1=scatter(1,nanmedian(dhpc_sub(dhpc_sub(:,6)==1,3)), "filled");s1.MarkerFaceColor = [0 1 0.5];
s2=scatter(2,nanmedian(vhpc_sub(vhpc_sub(:,6)==1,3)), "filled");s2.MarkerFaceColor = [0 0.5 1];

[p,tbl,stats] = kruskalwallis(y,x);
c = multcompare(stats);
tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])

subplot(1,3,3); hold on; 
x =[ones(size(dhpc_sub(dhpc_sub(:,6)==1,1),1),1);2*ones(size(vhpc_sub(vhpc_sub(:,6)==1,1),1),1)]; 
y= [dhpc_sub(dhpc_sub(:,6)==1,4);vhpc_sub(vhpc_sub(:,6)==1,4)]; 
c = [.3,.3,.3];
scatter(x,y,30,c,"filled",'jitter','on', 'jitterAmount',0.3); xlim([0 3]);
ylabel('Pf shift');
xticks([1 2]);
xticklabels({'dHPC', 'vHPC'});
hold on
s1=scatter(1,nanmedian(dhpc_sub(dhpc_sub(:,6)==1,4)), "filled");s1.MarkerFaceColor = [0 1 0.5];
s2=scatter(2,nanmedian(vhpc_sub(vhpc_sub(:,6)==1,4)), "filled");s2.MarkerFaceColor = [0 0.5 1];

[p,tbl,stats] = kruskalwallis(y,x);
c = multcompare(stats);
tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])


%% REMAPPING PARAMETERS PLOTS AND STATS - LAPS
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
            load('dHPC_pc_lap.mat');%'dHPC_pc_skaggs_pos.mat'
             %Save session parameters 
        for d=1:size(dHPC_lap,2)
            pc = dHPC_lap{d}; 
            bet = [pc.between,1]; 
            wa = [pc.within.ave,2];
            wr = [pc.within.rew,3];
            dhpc_temp = [dhpc_temp;bet;wa;wr]; 
            
        end
        catch 
            disp(['No dHPC_pc_lap.mat file in ',session]);
        end 
        
        try
            load('vHPC_pc_lap.mat');%'vHPC_pc_skaggs_pos.mat'
            for d=1:size(vHPC_lap,2)
                pc = vHPC_lap{d}; 
                bet = [pc.between,1]; 
                wa = [pc.within.ave,2];
                wr = [pc.within.rew,3];
                vhpc_temp = [vhpc_temp;bet;wa;wr]; 
            end 
     
        catch
            disp(['No vHPC_pc_lap.mat file in ',session]);
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
figure(2);clf;hold on, 
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

%%%%%%%% STATS VHPC
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
[P,ANOVATAB,STATS] = kruskalwallis(vhpc_sub(:,4),vhpc_sub(:,6));
c = multcompare(STATS)

tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])




%%%%% vhpc vs dhpc comparison 

%Only 3 plots
ylabels ={'Spatial correlation', 'Overlap', 'Pf shift'};
figure(1);clf;hold on, 
sgtitle('dHPC vs vHPC')
subplot(1,3,1); hold on; 
x =[ones(size(dhpc_sub(dhpc_sub(:,6)==1,1),1),1);2*ones(size(vhpc_sub(vhpc_sub(:,6)==1,1),1),1)]; 
y= [dhpc_sub(dhpc_sub(:,6)==1,1);vhpc_sub(vhpc_sub(:,6)==1,1)]; 
c = [.3,.3,.3];
scatter(x,y,30,c,"filled",'jitter','on', 'jitterAmount',0.3); xlim([0 3]);
ylabel('Spatial correlation');
xticks([1 2]);
xticklabels({'dHPC', 'vHPC'});
hold on
s1=scatter(1,nanmedian(dhpc_sub(dhpc_sub(:,6)==1,1)), "filled");s1.MarkerFaceColor = [0 1 0.5];
s2=scatter(2,nanmedian(vhpc_sub(vhpc_sub(:,6)==1,1)), "filled");s2.MarkerFaceColor = [0 0.5 1];
  

[p,tbl,stats] = kruskalwallis(y,x);
c = multcompare(stats);
tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])


subplot(1,3,2); hold on; 
x =[ones(size(dhpc_sub(dhpc_sub(:,6)==1,1),1),1);2*ones(size(vhpc_sub(vhpc_sub(:,6)==1,1),1),1)]; 
y= [dhpc_sub(dhpc_sub(:,6)==1,3);vhpc_sub(vhpc_sub(:,6)==1,3)]; 
c = [.3,.3,.3];
scatter(x,y,30,c,"filled",'jitter','on', 'jitterAmount',0.3); xlim([0 3]);
ylabel('Overlap');
xticks([1 2]);
xticklabels({'dHPC', 'vHPC'});
hold on
s1=scatter(1,nanmedian(dhpc_sub(dhpc_sub(:,6)==1,3)), "filled");s1.MarkerFaceColor = [0 1 0.5];
s2=scatter(2,nanmedian(vhpc_sub(vhpc_sub(:,6)==1,3)), "filled");s2.MarkerFaceColor = [0 0.5 1];

[p,tbl,stats] = kruskalwallis(y,x);
c = multcompare(stats);
tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])

subplot(1,3,3); hold on; 
x =[ones(size(dhpc_sub(dhpc_sub(:,6)==1,1),1),1);2*ones(size(vhpc_sub(vhpc_sub(:,6)==1,1),1),1)]; 
y= [dhpc_sub(dhpc_sub(:,6)==1,4);vhpc_sub(vhpc_sub(:,6)==1,4)]; 
c = [.3,.3,.3];
scatter(x,y,30,c,"filled",'jitter','on', 'jitterAmount',0.3); xlim([0 3]);
ylabel('Pf shift');
xticks([1 2]);
xticklabels({'dHPC', 'vHPC'});
hold on
s1=scatter(1,nanmedian(dhpc_sub(dhpc_sub(:,6)==1,4)), "filled");s1.MarkerFaceColor = [0 1 0.5];
s2=scatter(2,nanmedian(vhpc_sub(vhpc_sub(:,6)==1,4)), "filled");s2.MarkerFaceColor = [0 0.5 1];

[p,tbl,stats] = kruskalwallis(y,x);
c = multcompare(stats);
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
ylim([0 100])
title('Occupancy');  

mv = nanmean(v_ave,1);
semv = nansem(v_ave);
subplot(2,1,2);plot(1:60,Smooth(mv,1), 'LineWidth',2,'Color',[1 0.4 0.5]);hold on;
ciplot(Smooth(mv-semv,1),Smooth(mv+semv,1),1:60,'r'); alpha 0.1; hold on;
ylim([0 60])
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
            load('dHPC_pc_skaggs_circular.mat');
        catch 
            dHPC = []; disp(['No dHPC_pc_lap.mat file in ',session]);
        end 
        try
            load('vHPC_pc_skaggs_circular.mat');
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
        
        %%%%%%% Restrict pos to inside maze and only movement periods %%%%%
        %---Aversive----
        pos_ave = behavior.pos.aversive(:,1:2); 
        pos_ave(:,2) = pos_ave(:,2)-min(pos_ave(:,2)); pos_ave(:,2) = pos_ave(:,2)/max(pos_ave(:,2)); %normalization of position 
        in_maze = ToIntervals(pos_ave(:,1),and(pos_ave(:,2)>0.09 , pos_ave(:,2)<1-0.09));% eliminating the extrems of the maze
                    
       %Check lenght of in_maze segments to define real laps
        in_lapA = []; % keeps only full laps (maybe too restrictive)
        for ix=1:size(in_maze,1)
                        ini = in_maze(ix,1);
                        fin = in_maze(ix,2);
                        ini_pos = pos_ave(pos_ave(:,1)==ini,2); 
                        fin_pos = pos_ave(pos_ave(:,1)==fin,2); 
                        if abs(fin_pos-ini_pos)>0.6
                             in_lapA = [in_lapA;in_maze(ix,:)]; 
                        end
        end
        
        pos_ave = Restrict(pos_ave,in_lapA);pos_ave = Restrict(pos_ave, movement.aversive); % Restrict pos to movement periods and laps
       
         %---Rewarded----
        pos_rew = behavior.pos.reward(:,1:2);
        pos_rew(:,2) = pos_rew(:,2)-min(pos_rew(:,2)); pos_rew(:,2) = pos_rew(:,2)/max(pos_rew(:,2)); %normalization of position 
        in_maze = ToIntervals(pos_rew(:,1),and(pos_rew(:,2)>0.09 , pos_rew(:,2)<1-0.09));% eliminating the extrems of the maze (10cm-reward box)
               
        %Check lenght of in_maze segments to define real laps
        in_lapR = []; 
        for ix=1:size(in_maze,1)
                        ini = in_maze(ix,1);
                        fin = in_maze(ix,2);
                        ini_pos = pos_rew(pos_rew(:,1)==ini,2); 
                        fin_pos = pos_rew(pos_rew(:,1)==fin,2); 
                        if abs(fin_pos-ini_pos)>0.6
                             in_lapR = [in_lapR;in_maze(ix,:)]; 
                        end
        end 
        pos_rew = Restrict(pos_rew,in_lapR); pos_rew = Restrict(pos_rew , movement.reward);
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

        try
            load('dHPC_pc_skaggs_pos.mat'); 
        catch 
            dHPC = []; disp(['No dHPC file in ',session]);
        end 

    
        for ii=1:size(dHPC,2)
            
            cluster = dHPC{ii}.id;
  
            spks = spks_dHPC(spks_dHPC(:,1)==cluster,2); % select tspk from cluster
              
            % --- Aversive ---
            spks_ave = spks; 
            %Restrict spk and position to laps and movement periods
            spks_ave = Restrict(spks_ave,in_lapA);spks_ave = Restrict(spks_ave , movement.aversive); % Restrict spikes to movement periods and laps
        
            % --- Reward ---
             spks_rew = spks; 
             %Restrict spk  to laps and movement periods
             spks_rew = Restrict(spks_rew,in_lapR); spks_rew = Restrict(spks_rew, movement.reward); % Restrict to movement periods
               
             
                    
             [sub_stab_ave, sub_stab_rew] = Subsamplin_stability_1d(pos_ave,pos_rew,spks_ave,spks_rew,Xedges,sigma);
              
             % Save PC varaiables 
              n.id = cluster;
              n.withinStab_ave = sub_stab_ave;
              n.withinStab_rew = sub_stab_rew;
              dHPC{ii}= n;
              
             
        end 
        
        disp('vHPC Stability calculation')
        
        try
            load('vHPC_pc_skaggs_pos.mat'); 
        catch 
            vHPC = []; disp(['No vHPC file in ',session]);
        end 

    
        for ii=1:size(vHPC,2)
            
            cluster = vHPC{ii}.id;
  
            spks = spks_vHPC(spks_vHPC(:,1)==cluster,2); % select tspk from cluster
              
            % --- Aversive ---
            spks_ave = spks; 
            %Restrict spk and position to laps and movement periods
            spks_ave = Restrict(spks_ave,in_lapA);spks_ave = Restrict(spks_ave , movement.aversive); % Restrict spikes to movement periods and laps
        
            % --- Reward ---
             spks_rew = spks; 
             %Restrict spk  to laps and movement periods
             spks_rew = Restrict(spks_rew,in_lapR); spks_rew = Restrict(spks_rew, movement.reward); % Restrict to movement periods
               
             
                    
             [sub_stab_ave, sub_stab_rew] = Subsamplin_stability_1d(pos_ave,pos_rew,spks_ave,spks_rew,Xedges,sigma);
              
             % Save PC varaiables 
              n.id = cluster;
              n.withinStab_ave = sub_stab_ave;
              n.withinStab_rew = sub_stab_rew;
              vHPC{ii}= n;
              
             
        end 
        
        
        %% Saveing PC INFO 
    if ~isempty(dHPC)
        dHPC = dHPC(~cellfun('isempty',dHPC));
        save([cd,'\dHPC_pc_within_stab.mat'],'dHPC'); 
    end
    if ~isempty(vHPC)
        vHPC = vHPC(~cellfun('isempty',vHPC));
        save([cd,'\vHPC_pc_within_stab.mat'],'vHPC'); 
    end
  

    end 
    
end


%%%%%%%%%% Plots and stats %%%%%%%%%%

%c1 = spatial c2= fr_change c3=  overlap c4: pf shift c5:  1 = between 2 = within aversive 3= within reward 
clear
clc
close all
% Parameters
path = {'\\Maryjackson\e\Rat127\Ephys\pyr';'\\Maryjackson\e\Rat128\Ephys\in_pyr\ready';'\\Maryjackson\e\Rat103\usable';'\\Maryjackson\e\Rat132\recordings\in_pyr'; '\\Maryjackson\e\Rat165\in_pyr'; ...
    '\\Maryjackson\e\Rat126\Ephys\in_Pyr'};%List of folders from the path

%Output matrix 
dhpc_sub = []; dave = []; drew=[];
vhpc_sub = []; vave = []; vrew=[];

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
            load('dHPC_pc_within_stab.mat');
             %Save session parameters 
        for d=1:size(dHPC,2)
            pc = dHPC{d}; 
            dhpc_sub = [dhpc_sub;pc.withinStab_ave,1; pc.withinStab_rew,2];
            dave = [dave;pc.withinStab_ave];
            drew = [drew; pc.withinStab_rew];
            
        end
        catch 
            dHPC = []; disp(['No dHPC file in ',session]);
        end 
        
        try
           load('vHPC_pc_within_stab.mat');
            for d=1:size(vHPC,2)
                pc = vHPC{d}; 
                vhpc_sub = [vhpc_sub;pc.withinStab_ave,1; pc.withinStab_rew,2];
                vave = [vave;pc.withinStab_ave];
                vrew = [vrew; pc.withinStab_rew];
            end 
     
        catch
            vHPC= []; disp(['No vHPC file in ',session]);
        end
        
    end 
  
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
xlabels = {'Aversive', 'Reward'}; 
figure(1);clf;hold on, 
sgtitle('dHPC stability 1st vs 2nd')
for a=1:size(ylabels,2)
    p = idx(a);
    subplot(1,3,a); hold on; 
    x = dhpc_sub(:,5);
    y= dhpc_sub(:,p);%Change the column number (1-4) to choose which variable to plot 
    c = [.3,.3,.3];
    scatter(x,y,30,c,"filled",'jitter','on', 'jitterAmount',0.3); xlim([0 3]);
    ylabel(ylabels{a});
    xticks([1 2])
    xticklabels({'Ave', 'Rew'});
    hold on
    s2=scatter(1,nanmedian(dhpc_sub(dhpc_sub(:,5)==1,p)), "filled");s2.MarkerFaceColor = [1 0.1 0.2];
    s3=scatter(2,nanmedian(dhpc_sub(dhpc_sub(:,5)==2,p)), "filled");s3.MarkerFaceColor = [0.1 0.3 1];
end

x=3;
[p,h] = ranksum(dave(:,x),drew(:,x)) 

%vhpc
%Only 3 plots
ylabels ={'Spatial correlation', 'Overlap', 'Pf shift'};
idx = [1,3,4]; 
xlabels = {'Aversive', 'Reward'}; 
figure(2);clf;hold on, 
sgtitle('vHPC stability 1st vs 2nd')
for a=1:size(ylabels,2)
    p = idx(a);
    subplot(1,3,a); hold on; 
    x = vhpc_sub(:,5);
    y= vhpc_sub(:,p);%Change the column number (1-4) to choose which variable to plot 
    c = [.3,.3,.3];
    scatter(x,y,30,c,"filled",'jitter','on', 'jitterAmount',0.3); xlim([0 3]);
    ylabel(ylabels{a});
    xticks([1 2])
    xticklabels({'Ave', 'Rew'});
    hold on
    s2=scatter(1,nanmedian(vhpc_sub(vhpc_sub(:,5)==1,p)), "filled");s2.MarkerFaceColor = [1 0.1 0.2];
    s3=scatter(2,nanmedian(vhpc_sub(vhpc_sub(:,5)==2,p)), "filled");s3.MarkerFaceColor = [0.1 0.3 1];
end

x=3;
[p,h] = ranksum(vave(:,x),vrew(:,x)) 


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
n_dhpc = 0; n_vhpc=0; 
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
            load('dHPC_pc_skaggs_circular.mat');
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
            n_dhpc=n_dhpc+size(dHPC,2);
        catch 
            dHPC = []; disp(['No dHPC file in ',session]);
            
        end 
        
        
        try
            load('vHPC_pc_skaggs_circular.mat');
            for d=1:size(vHPC,2)
                pc = vHPC{d}; 
                A = pc.frMap_ave - min(pc.frMap_ave);A = A ./ max(A);
                vhpc_ave=[vhpc_ave; A];
                clear A
            
                A = pc.frMap_rew - min(pc.frMap_rew);A = A ./ max(A);
                vhpc_rew=[vhpc_rew; A];
                clear A
            end
            n_vhpc=n_vhpc+size(vHPC,2);
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
subplot(1,2,1);imagesc([0:20:140], [1:1:size(dhpc_ave,1)],fr), colormap 'gray'; title('Aversive');
fr = dhpc_rew(mm,:);  fr(~any(~isnan(fr), 2),:)=[];
subplot(1,2,2);imagesc([0:20:140], [1:1:size(dhpc_rew,1)],fr), colormap 'gray'; title('Reward');
sgtitle('dHPC firing maps');


% sorted by reward
[h idx] = max (dhpc_rew, [],2);
[m mm] = sort(idx); 

figure(1);clf;hold on;
fr = dhpc_rew(mm,:);  fr(~any(~isnan(fr), 2),:)=[];
subplot(1,2,1);imagesc([0:20:140], [1:1:size(dhpc_rew,1)],fr), colormap 'gray'; title('Reward');
fr = dhpc_ave(mm,:);  fr(~any(~isnan(fr), 2),:)=[];
subplot(1,2,2);imagesc([0:20:140], [1:1:size(dhpc_ave,1)],fr), colormap 'gray'; title('Aversive');
sgtitle('dHPC firing maps');


%%%%% vHPC Plots %%%%%

[h idx] = max (vhpc_ave, [],2);
[m mm] = sort(idx); 


figure(2);clf;
fr = vhpc_ave(mm,:);  fr(~any(~isnan(fr), 2),:)=[];
subplot(1,2,1); imagesc([0:20:140], [1:1:size(vhpc_ave,1)],fr), caxis([0 1]),colormap 'gray'; axis tight; title('Aversive');
fr = vhpc_rew(mm,:);  fr(~any(~isnan(fr), 2),:)=[];
subplot(1,2,2); imagesc([0:20:140], [1:1:size(vhpc_rew,1)],fr), caxis([0 1]), colormap 'gray'; title('Reward');axis tight; 

%reward sorted
[h idx] = max (vhpc_rew, [],2);
[m mm] = sort(idx); 


figure(2);clf;
fr = vhpc_rew(mm,:);  fr(~any(~isnan(fr), 2),:)=[];
subplot(1,2,1); imagesc([0:20:140], [1:1:size(vhpc_rew,1)],fr), caxis([0 1]),colormap 'gray'; axis tight; title('Rewarded');
fr = vhpc_rew(mm,:);  fr(~any(~isnan(fr), 2),:)=[];
subplot(1,2,2); imagesc([0:20:140], [1:1:size(vhpc_rew,1)],fr), caxis([0 1]), colormap 'gray'; title('Reward');axis tight; 



%% PC PARAMETERS: skaggs and pf size 
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
            load('dHPC_pc_skaggs_circular.mat');
             %Save session parameters 
        for d=1:size(dHPC,2)
            pc = dHPC{d}; 
            dhpc_size = [dhpc_size;(pc.stats_ave.size(1)/60)*100,1;(pc.stats_rew.size(1)/60)*100,2]; 
            dhpc_skaggs = [dhpc_skaggs;pc.stats_ave.specificity,1;pc.stats_rew.specificity,2];
            dhpc_n = dhpc_n+1;
        end
        catch 
            dHPC = []; disp(['No dHPC_pc_skaggs_circular.mat file in',session]);
        end 
        
        try
            load('vHPC_pc_skaggs_circular.mat');
            for d=1:size(vHPC,2)
                pc = vHPC{d}; 
                vhpc_size = [vhpc_size;(pc.stats_ave.size(1)/60)*100,1;(pc.stats_rew.size(1)/60)*100,2]; 
                vhpc_skaggs = [vhpc_skaggs;pc.stats_ave.specificity,1;pc.stats_rew.specificity,2];
                vhpc_n = vhpc_n+1;
            end 
     
        catch
            vHPC= []; disp(['No vHPC_pc_skaggs_circular.matfile in ',session]);
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
    s1=scatter(1,nanmean(dhpc_size(dhpc_size(:,2)==1,1)), "filled");s1.MarkerFaceColor = [1 0.1 0.2];
    s3=scatter(2,nanmean(dhpc_size(dhpc_size(:,2)==2,1)), "filled");s3.MarkerFaceColor = [0.1 0.3 1];
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
    s1=scatter(1,nanmean(vhpc_size(vhpc_size(:,2)==1,1)), "filled");s1.MarkerFaceColor = [1 0.1 0.2];
    s3=scatter(2,nanmean(vhpc_size(vhpc_size(:,2)==2,1)), "filled");s3.MarkerFaceColor = [0.1 0.3 1];
    title('vhpc'); 
    
    
    
    figure(1);
    boxplot(dhpc_size(:,1), dhpc_size(:,2))
    xticks([1 2])
    xticklabels({'Ave','Rew'});
    ylabel('Place field size (% area)');
     title('dhpc');
     
     
     figure(2);
    boxplot(vhpc_size(:,1), vhpc_size(:,2))
    xticks([1 2])
    xticklabels({'Ave','Rew'});
    ylabel('Place field size (% area)');
     title('vhpc');
 

data=[vhpc_size(vhpc_size(:,2)==1,1), vhpc_size(vhpc_size(:,2)==2,1)];
data(any(isnan(data), 2), :) = [];

[p,~,stats] =  friedman(data,1); 
c = multcompare(stats);
tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])    
     
     
     
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
    s1=scatter(1,nanmean(dhpc_skaggs(dhpc_skaggs(:,2)==1,1)), "filled");s1.MarkerFaceColor = [1 0.1 0.2];
    s3=scatter(2,nanmean(dhpc_skaggs(dhpc_skaggs(:,2)==2,1)), "filled");s3.MarkerFaceColor = [0.1 0.3 1];
    title('dhpc'); 
    boxplot(dhpc_skaggs(:,1), dhpc_skaggs(:,2))
    xticks([1 2])
    xticklabels({'Ave','Rew'});
    
    
    
    subplot(1,2,2); hold on; 
    x = vhpc_skaggs(:,2);
    y= vhpc_skaggs(:,1);
    c = [.3,.3,.3];
    scatter(x,y,[],c,"filled",'jitter','on', 'jitterAmount',0.1); xlim([0 3]);
    ylabel(ylabels);
    xticks([1 2])
    xticklabels({'Ave','Rew'});
    hold on
    s1=scatter(1,nanmean(vhpc_skaggs(vhpc_skaggs(:,2)==1,1)), "filled");s1.MarkerFaceColor = [1 0.1 0.2];
    s3=scatter(2,nanmean(vhpc_skaggs(vhpc_skaggs(:,2)==2,1)), "filled");s3.MarkerFaceColor = [0.1 0.3 1];
    title('vhpc'); 
    boxplot(vhpc_skaggs(:,1), vhpc_skaggs(:,2))
    xticklabels({'Ave','Rew'});
    
    
    
    data=[dhpc_skaggs(dhpc_skaggs(:,2)==1,1), dhpc_skaggs(dhpc_skaggs(:,2)==2,1)];
    data(any(isnan(data), 2), :) = [];

[p,~,stats] =  friedman(data,1); 
c = multcompare(stats);
tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])    

% PLOT skaggs general

% Plots skaggs
 
ylabels ={'Skaggs'}; 
xlabels = {'dHPC', 'vHPC'}; 
figure(3);clf;hold on, 
title('Skaggs')

   
y = [dhpc_skaggs(:,1);vhpc_skaggs(:,1)];
x=[ones(size(dhpc_skaggs,1),1);2*ones(size(vhpc_skaggs,1),1)];
c = [.3,.3,.3];
scatter(x,y,[],c,"filled",'jitter','on', 'jitterAmount',0.1); xlim([0.5 2.5]);
ylabel(ylabels);
xticks([1 2])
xticklabels({'dHPC', 'vHPC'});
hold on
s1=scatter(1,nanmean(dhpc_skaggs(:,1)), "filled");s1.MarkerFaceColor = [1 0.1 0.2];
s3=scatter(2,nanmean(vhpc_skaggs(:,1)), "filled");s3.MarkerFaceColor = [0.1 0.3 1];
   

[p,h] = ranksum(dhpc_skaggs(:,1),vhpc_skaggs(:,1))   

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

figure(3);clf;hold on;
sgtitle('vHPC laps')
histogram(vhpc_nlap(:,1),'FaceColor','r','EdgeColor','none','BinWidth',5);
histogram(vhpc_nlap(:,2),'FaceColor','b','EdgeColor','none','BinWidth',5);


%% PLOT PC EXAMPLES

clear
clc
% close all
path = {'\\Maryjackson\e\Rat127\Ephys\pyr';'\\Maryjackson\e\Rat128\Ephys\in_pyr\ready';'\\Maryjackson\e\Rat103\usable';'\\Maryjackson\e\Rat132\recordings\in_pyr'; '\\Maryjackson\e\Rat165\in_pyr'; ...
    '\\Maryjackson\e\Rat126\Ephys\in_Pyr'};%List of folders from the path
n_SU_V = [];
n_SU_D = [];
Xedges = 60; %number of bins for RateMap construction - 3cm bin
sigma = 2;%round(15/(180/Xedges)); %defined for gauss kernel of 15cm
binSize = 0.001; % bin size for replay events detection
bin_size = 1; % to bin pos ans spks in between/within  
criteria_type = 0; %criteria for celltype (0:pyr, 1:int, 2:all)
% Behavior
minimal_speed = 2.5;% minimal speed to detect quite periods
minimal_speed_time = 2; % minimal time to detect quite periods


pc_dhpc_all=[];
pc_vhpc_all=[];

for tt = 4:length(path)
    %List of folders from the path
    files = dir(path{tt});
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
   
    clear files dirFlags
    
    for t =7 : length(subFolders)-2
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
        
        %Restrict pos to inside maze and only movement periods
        %---Aversive----
        pos_ave = behavior.pos.aversive(:,1:2); 
        pos_ave(:,2) = pos_ave(:,2)-min(pos_ave(:,2)); pos_ave(:,2) = pos_ave(:,2)/max(pos_ave(:,2)); %normalization of position 
        in_maze = ToIntervals(pos_ave(:,1),and(pos_ave(:,2)>0.09 , pos_ave(:,2)<1-0.09));% eliminating the extrems of the maze
                    
       %Check lenght of in_maze segments to define real laps
        in_lapA = []; % keeps only full laps (maybe too restrictive)
        for ix=1:size(in_maze,1)
                        ini = in_maze(ix,1);
                        fin = in_maze(ix,2);
                        ini_pos = pos_ave(pos_ave(:,1)==ini,2); 
                        fin_pos = pos_ave(pos_ave(:,1)==fin,2); 
                        if abs(fin_pos-ini_pos)>0.6
                             in_lapA = [in_lapA;in_maze(ix,:)]; 
                        end
        end
        
        pos_ave = Restrict(pos_ave,in_lapA);pos_ave = Restrict(pos_ave, movement.aversive); % Restrict pos to movement periods and laps
       
         %---Rewarded----
        pos_rew = behavior.pos.reward(:,1:2);
        pos_rew(:,2) = pos_rew(:,2)-min(pos_rew(:,2)); pos_rew(:,2) = pos_rew(:,2)/max(pos_rew(:,2)); %normalization of position 
        in_maze = ToIntervals(pos_rew(:,1),and(pos_rew(:,2)>0.09 , pos_rew(:,2)<1-0.09));% eliminating the extrems of the maze (10cm-reward box)
               
        %Check lenght of in_maze segments to define real laps
        in_lapR = []; 
        for ix=1:size(in_maze,1)
                        ini = in_maze(ix,1);
                        fin = in_maze(ix,2);
                        ini_pos = pos_rew(pos_rew(:,1)==ini,2); 
                        fin_pos = pos_rew(pos_rew(:,1)==fin,2); 
                        if abs(fin_pos-ini_pos)>0.6
                             in_lapR = [in_lapR;in_maze(ix,:)]; 
                        end
        end 
        pos_rew = Restrict(pos_rew,in_lapR); pos_rew = Restrict(pos_rew , movement.reward);

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
            load('dHPC_pc_skaggs_circular.mat');
            dhpc_pc = [];
            for d=1:size(dHPC,2)
                pc = dHPC{d}; 
                dhpc_pc = [dhpc_pc;pc.id]; 
            
            end
            group_dHPC = group_dHPC(ismember(group_dHPC(:,1),dhpc_pc),:);
            
            pc_dhpc_session = [];
            
            % Plotting: 
            for ii=12:size(group_dHPC,1)
                cluster = group_dHPC(ii,1);

                spks = spks_dHPC(spks_dHPC(:,1)==cluster,2); % select tspk from cluster
                
                %Save pc info 
                pc_dhpc_session(ii,1)=cluster; 
                
                pc_dhpc_session(ii,2)=dHPC{1,ii}.stats_ave.specificity;% skaggs ave
                pc_dhpc_session(ii,3)=dHPC{1,ii}.stats_rew.specificity;% skaggs rew
                
                pc_dhpc_session(ii,4)=dHPC{1,ii}.between(1);% between spatial
                pc_dhpc_session(ii,5)=dHPC{1,ii}.between(3);% between overlap
                
                pc_dhpc_session(ii,6)=dHPC{1,ii}.within_ave(1);%  spatial
                pc_dhpc_session(ii,7)=dHPC{1,ii}.within_ave(3);% n overlap
                
                pc_dhpc_session(ii,8)=dHPC{1,ii}.within_rew(1);%  spatial
                pc_dhpc_session(ii,9)=dHPC{1,ii}.within_rew(3);% n overlap
                
                 
               % --- Aversive ---
                spks_ave = spks; 
                %Restrict spk and position to laps and movement periods
                spks_ave = Restrict(spks_ave,in_lapA);spks_ave = Restrict(spks_ave , movement.aversive); % Restrict spikes to movement periods and laps
               
                % --- Reward ---
                spks_rew = spks; 
                %Restrict spk  to laps and movement periods
                spks_rew = Restrict(spks_rew,in_lapR); spks_rew = Restrict(spks_rew, movement.reward); % Restrict to movement periods
                     
   
                %Checking pc criteria
                % 1. Mean firing rate above tresh.
                % 2.Skaggs > skaggs random. 
                % 3. At least 1 pf > 4 bins in one of the cond.  
     
                [curveA , statsA] = FiringCurve(pos_ave , spks_ave , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 6, 'minPeak' , 0.2);
                [curveR , statsR] = FiringCurve(pos_rew , spks_rew , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 6, 'minPeak' , 0.2);

            %%%%%%% PLOT 
            figure(2);clf; hold on; 
            set(gcf,'position',[100,100,800,400])
            sgtitle(['dHPC ',session(end-14:end),' id ',num2str(cluster),' iteration ',num2str(ii)])
            subplot(5,6,[1 15]); hold on; plot(pos_ave(:,2),pos_ave(:,1),'color','red','LineWidth',1);hold on; axis tight; xline(0.25);set(gca,'yticklabel',[],'TickLength',[0 0],'YColor','none');
            title('Aversive'); ylabel('Laps');
            %Find interpolated position of each spike:
            xs = interp1(pos_ave(:,1),pos_ave(:,2),spks_ave); s= scatter(xs,spks_ave,100,'filled','Marker','.', 'MarkerFaceColor',[.3 .3 .3],'MarkerEdgeColor',[.3 .3 .3]); 
            %Fr map 
            subplot(5,6,[19 21]); hold on;plot(curveA.rate,'Color',[.3 .3 .3],'LineWidth',1.5);ylabel('Fr(Hz)');xlabel('Spatial Bins');
            subplot(5,6,[25 27]);imagesc(curveA.rate),colormap 'jet';
           
            subplot(5,6,[4 18]); hold on; plot(pos_rew(:,2),pos_rew(:,1),'color',[0.175 0.54 0.60],'LineWidth',1);hold on; axis tight;
            xline(0.25);set(gca,'yticklabel',[],'TickLength',[0 0],'YColor','none');
            title('Reward'); ylabel('Laps');
            %Find interpolated position of each spike:
            xs = interp1(pos_rew(:,1),pos_rew(:,2),spks_rew); scatter(xs,spks_rew,100,'filled','Marker','.', 'MarkerFaceColor',[.3 .3 .3],'MarkerEdgeColor',[.3 .3 .3]); 
            %Fr map 
            subplot(5,6,[22 24]); hold on;plot(curveR.rate,'Color',[.3 .3 .3],'LineWidth',1.5);ylabel('Fr(Hz)');xlabel('Spatial Bins');
            subplot(5,6,[28 30]);imagesc(curveR.rate),colormap 'jet';
            
            % Save figure:
            saveas(figure(2),['W:\Remapping-analysis-Facu\Place cells examples\all\',session(end-14:end),'_dhpc_',num2str(cluster)], 'jpg');
               

            end
       
            
        catch 
            dHPC = []; disp(['No dHPC_pc.mat file in ',session]);

        end
        
        
        try
            load('vHPC_pc_skaggs_circular.mat');
            vhpc_pc = [];
            for d=1:size(vHPC,2)
                pc = vHPC{d}; 
                vhpc_pc = [vhpc_pc;pc.id]; 
            
            end
            group_vHPC = group_vHPC(ismember(group_vHPC(:,1),vhpc_pc),:);
            
             pc_vhpc_session=[]; 
                
            disp('Plot vHPC Firing rate map')
          
            for ii=13:size(group_vHPC,1)
                
                cluster = group_vHPC(ii,1);
                spks = spks_vHPC(spks_vHPC(:,1)==cluster,2); % select tspk from cluster
              
                pc_vhpc_session(ii,1)=cluster; 
                
                pc_vhpc_session(ii,2)=vHPC{1,ii}.stats_ave.specificity;% skaggs ave
                pc_vhpc_session(ii,3)=vHPC{1,ii}.stats_rew.specificity;% skaggs rew
                
                pc_vhpc_session(ii,4)=vHPC{1,ii}.between(1);% between spatial
                pc_vhpc_session(ii,5)=vHPC{1,ii}.between(3);% between overlap
                
                pc_vhpc_session(ii,6)=vHPC{1,ii}.within_ave(1);%  spatial
                pc_vhpc_session(ii,7)=vHPC{1,ii}.within_ave(3);% n overlap
                
                pc_vhpc_session(ii,8)=vHPC{1,ii}.within_rew(1);%  spatial
                pc_vhpc_session(ii,9)=vHPC{1,ii}.within_rew(3);% n overlap
                
               % --- Aversive ---
                spks_ave = spks; 
                %Restrict spk and position to laps and movement periods
                spks_ave = Restrict(spks_ave,in_lapA);spks_ave = Restrict(spks_ave , movement.aversive); % Restrict spikes to movement periods and laps
               
                % --- Reward ---
                spks_rew = spks; 
                %Restrict spk  to laps and movement periods
                spks_rew = Restrict(spks_rew,in_lapR); spks_rew = Restrict(spks_rew, movement.reward); % Restrict to movement periods
                     
   
                %Checking pc criteria
                % 1. Mean firing rate above tresh.
                % 2.Skaggs > skaggs random. 
                % 3. At least 1 pf > 4 bins in one of the cond.  
     
                [curveA , statsA] = FiringCurve(pos_ave , spks_ave , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 6, 'minPeak' , 0.2);
                [curveR , statsR] = FiringCurve(pos_rew , spks_rew , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 6, 'minPeak' , 0.2);

                    
           %%%%%%% PLOT 
            figure(3);clf; hold on; 
            set(gcf,'position',[100,100,800,400])
            sgtitle(['vHPC ',session(end-14:end),' id ',num2str(cluster),' iteration ',num2str(ii)])
            subplot(5,6,[1 15]); hold on; plot(pos_ave(:,2),pos_ave(:,1),'color','red','LineWidth',1);hold on; axis tight; xline(0.25);set(gca,'yticklabel',[],'TickLength',[0 0],'YColor','none');
            title('Aversive'); ylabel('Laps');
            %Find interpolated position of each spike:
            xs = interp1(pos_ave(:,1),pos_ave(:,2),spks_ave); s= scatter(xs,spks_ave,100,'filled','Marker','.', 'MarkerFaceColor',[.3 .3 .3],'MarkerEdgeColor',[.3 .3 .3]); 
            %Fr map 
            subplot(5,6,[19 21]); hold on;plot(curveA.rate,'Color',[.3 .3 .3],'LineWidth',1.5);ylabel('Fr(Hz)');xlabel('Spatial Bins');
            subplot(5,6,[25 27]);imagesc(curveA.rate),colormap 'jet';
           
            subplot(5,6,[4 18]); hold on; plot(pos_rew(:,2),pos_rew(:,1),'color',[0.175 0.54 0.60],'LineWidth',1);hold on; axis tight;
            xline(0.25);set(gca,'yticklabel',[],'TickLength',[0 0],'YColor','none');
            title('Reward'); ylabel('Laps');
            %Find interpolated position of each spike:
            xs = interp1(pos_rew(:,1),pos_rew(:,2),spks_rew); scatter(xs,spks_rew,100,'filled','Marker','.', 'MarkerFaceColor',[.3 .3 .3],'MarkerEdgeColor',[.3 .3 .3]); 
            %Fr map 
            subplot(5,6,[22 24]); hold on;plot(curveR.rate,'Color',[.3 .3 .3],'LineWidth',1.5);ylabel('Fr(Hz)');xlabel('Spatial Bins');
            subplot(5,6,[28 30]);imagesc(curveR.rate),colormap 'jet';
            
            % Save figure:
            saveas(figure(3),['W:\Remapping-analysis-Facu\Place cells examples\all\',session(end-14:end),'_vhpc_',num2str(cluster)], 'pdf');
            end

        catch
            vHPC= []; disp(['No vHPC_pc.mat file in ',session]);
        end
        
        pc_dhpc_all = [pc_dhpc_all;pc_dhpc_session]; 
        pc_vhpc_all = [pc_vhpc_all;pc_vhpc_session ]; 
        
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

%% PC per session 
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

    for t = 1 : length(subFolders)-2
        disp(['-- Initiating analysis of folder #' , num2str(t) , ' from rat #',num2str(tt) , ' --'])
        session = [subFolders(t+2).folder,'\',subFolders(t+2).name];
        cd([session,'\Spikesorting'])
        
        dhpc_sub =[]; 
        
        
        %Loading pc matrix of the sessions
        disp('Uploading session pc matrix');
        try
            load('dHPC_pc_skaggs_pos.mat');
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
            load('vHPC_pc_skaggs_pos.mat');
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

%% # PC per session 
clear
clc
close all
% Parameters
path = {'\\Maryjackson\e\Rat127\Ephys\pyr';'\\Maryjackson\e\Rat128\Ephys\in_pyr\ready';'\\Maryjackson\e\Rat103\usable';'\\Maryjackson\e\Rat132\recordings\in_pyr'; '\\Maryjackson\e\Rat165\in_pyr'; ...
    '\\Maryjackson\e\Rat126\Ephys\in_Pyr'};%List of folders from the path

%Output matrix 
pc_session= [];
ind_global=1;

for tt = 1:length(path)
    %List of folders from the path
    files = dir(path{tt});
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
    clear files dirFlags
    rat = str2num(subFolders(3).name(4:6));
    
    for t = 1 : length(subFolders)-2
        disp(['-- Initiating analysis of folder #' , num2str(t) , ' from rat #',num2str(tt) , ' --'])
        session = [subFolders(t+2).folder,'\',subFolders(t+2).name];
        cd([session,'\Spikesorting'])
        %Loading pc matrix of the session
        s= str2num(session(end-7:end));
        
        disp('Uploading session pc matrix');
        pc_session(ind_global,1) =rat; 
        pc_session(ind_global,2) =s; 
        try
            load('dHPC_pc_skaggs_circular.mat');%'dHPC_pc_skaggs_pos.mat'
             %Save session parameters 
           
            pc_session(ind_global,3) =size(dHPC,2); 
            
        catch 
            dHPC = []; disp(['No dHPC_pc_lap.mat file in ',session]);
        end 
        
        try
            load('vHPC_pc_skaggs_circular.mat');%'vHPC_pc_skaggs_pos.mat'
            pc_session(ind_global,4) =size(vHPC,2); 
        catch
            vHPC= []; disp(['No vHPC_pc_lap.mat file in ',session]);
        end
        ind_global = ind_global+1; 
    end 
    
     
end 


T = array2table(pc_session);
% Default heading for the columns will be A1, A2 and so on. 
% You can assign the specific headings to your table in the following manner
T.Properties.VariableNames(1:4) = {'Rat','Session','dHPC pc','vHPC pc'}

%% TRAJECTORY AND VELOCITY PLOTS
%First run initial segment 
% tt chose rat
% t chose session

for tt = 1:length(path)
    %List of folders from the path
    files = dir(path{tt});
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
   
    clear files dirFlags
    %t=9
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
            clear pos camaraR posx posy22
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
        
        %Restrict pos to inside maze and only movement periods
        %---Aversive----
        pos_ave = behavior.pos.aversive(:,1:2); 
        pos_ave(:,2) = pos_ave(:,2)-min(pos_ave(:,2)); pos_ave(:,2) = pos_ave(:,2)/max(pos_ave(:,2)); %normalization of position 
        in_maze = ToIntervals(pos_ave(:,1),and(pos_ave(:,2)>0.09 , pos_ave(:,2)<1-0.09));% eliminating the extrems of the maze
                    
       %Check lenght of in_maze segments to define real laps
        in_lapA = []; % keeps only full laps (maybe too restrictive)
        for ix=1:size(in_maze,1)
                        ini = in_maze(ix,1);
                        fin = in_maze(ix,2);
                        ini_pos = pos_ave(pos_ave(:,1)==ini,2); 
                        fin_pos = pos_ave(pos_ave(:,1)==fin,2); 
                        if abs(fin_pos-ini_pos)>0.6
                             in_lapA = [in_lapA;in_maze(ix,:)]; 
                        end
        end
        
        pos_ave = Restrict(pos_ave,in_lapA);pos_ave = Restrict(pos_ave, movement.aversive); % Restrict pos to movement periods and laps
       
         %---Rewarded----
        pos_rew = behavior.pos.reward(:,1:2);
        pos_rew(:,2) = pos_rew(:,2)-min(pos_rew(:,2)); pos_rew(:,2) = pos_rew(:,2)/max(pos_rew(:,2)); %normalization of position 
        in_maze = ToIntervals(pos_rew(:,1),and(pos_rew(:,2)>0.09 , pos_rew(:,2)<1-0.09));% eliminating the extrems of the maze (10cm-reward box)
               
        %Check lenght of in_maze segments to define real laps
        in_lapR = []; 
        for ix=1:size(in_maze,1)
                        ini = in_maze(ix,1);
                        fin = in_maze(ix,2);
                        ini_pos = pos_rew(pos_rew(:,1)==ini,2); 
                        fin_pos = pos_rew(pos_rew(:,1)==fin,2); 
                        if abs(fin_pos-ini_pos)>0.6
                             in_lapR = [in_lapR;in_maze(ix,:)]; 
                        end
        end 
        pos_rew = Restrict(pos_rew,in_lapR); pos_rew = Restrict(pos_rew , movement.reward);

       
        
        %% Plot
       
          
%             %%%%%%% PLOT op 1
%             figure(4);clf; hold on; 
%             set(gcf,'position',[100,100,400,800])
%             sgtitle(['Behavior ',session(end-14:end)])
%             %%%Ave 
%             %Position
%             subplot(2,1,1); hold on; plot(behavior.pos.aversive(:,2),behavior.pos.aversive(:,1)/.60000,'color','red','LineWidth',1);hold on; axis tight;%set(gca,'yticklabel',[],'TickLength',[0 0],'YColor','none');
% %             line([50 50], [1.865452772500000e+04 (1.865452772500000e+04 +60000)]);
%             title('Aversive run'); 
%             %Find interpolated position of each shock:
%             xs = interp1(behavior.pos.aversive(:,1),behavior.pos.aversive(:,2),Shocks_filt); 
%             s= scatter(xs,Shocks_filt,300,'filled','Marker','.', 'MarkerFaceColor','yellow','MarkerEdgeColor','yellow'); 
%             %Velocity 
%              plot(behavior.speed.aversive(:,2),behavior.speed.aversive(:,1),'color',[.3 .3 .3],'LineWidth',1)
%             
%              
%             %%%Rew
%             %Position
%             subplot(2,1,2); hold on; plot(behavior.pos.reward(:,2),behavior.pos.reward(:,1),'color',[0.175 0.54 0.60],'LineWidth',1);hold on; axis tight;set(gca,'yticklabel',[],'TickLength',[0 0],'YColor','none');
%             title('Reward run'); ylabel('Laps');
%             %Water
%             xs = interp1(behavior.pos.reward(:,1),behavior.pos.reward(:,2), Rewards_filt(:,1)); 
%             s= scatter(xs, Rewards_filt(:,1),300,'filled','Marker','.', 'MarkerFaceColor','blue','MarkerEdgeColor','blue'); 
%             %Velocity 
%             plot(behavior.speed.reward(:,2),behavior.speed.reward(:,1),'color',[.3 .3 .3],'LineWidth',1)
    
            
            %%%%%%% PLOT OP 2
            
            figure();clf; hold on; 
            set(gcf,'position',[100,100,400,800])
            sgtitle(['Behavior ',session(end-14:end)])
            %%%Ave 
    
            subplot(2,1,1); hold on; 
%             scatter(behavior.pos.aversive(:,2),behavior.pos.aversive(:,1)/.60000, 30, behavior.speed.aversive(:,2), 'filled');% velocity
            plot(behavior.pos.aversive(:,2),behavior.pos.aversive(:,1)/.60000);%line
            hold on; axis tight;set(gca,'yticklabel',[],'TickLength',[0 0],'YColor','none');
            
            title('Aversive run'); 
            %Find interpolated position of each shock:
%             xs = interp1(behavior.pos.aversive(:,1),behavior.pos.aversive(:,2),Shocks_filt); 
%             s= scatter(xs,Shocks_filt,300,'filled','Marker','.', 'MarkerFaceColor','yellow','MarkerEdgeColor','yellow'); 
%            
             
            %%%Rew
           
            %Position
            subplot(2,1,2); hold on; %scatter(behavior.pos.reward(:,2),behavior.pos.reward(:,1)/.60000, 30, behavior.speed.reward(:,2), 'filled');
            plot(behavior.pos.reward(:,2),behavior.pos.reward(:,1)/.60000);
            hold on; axis tight;set(gca,'yticklabel',[],'TickLength',[0 0],'YColor','none');
            title('Reward run'); ylabel('Laps');
            %Water
%             xs = interp1(behavior.pos.reward(:,1),behavior.pos.reward(:,2), Rewards_filt(:,1)); 
%             s= scatter(xs, Rewards_filt(:,1),300,'filled','Marker','.', 'MarkerFaceColor','blue','MarkerEdgeColor','blue'); 
%          
     
            %%%%%%%%%%%
            
            % Save figure:
%             saveas(figure(2),['W:\Remapping-analysis-Facu\Behavior plots\','trajectory_velocity',session(end-14:end)], 'pdf');
               

        disp(['-- Finished folder #' , num2str(t) , ' from rat #' , num2str(tt) , ' ---'])
        disp('  ')
    end
    disp(['-------- Finished rat#' , num2str(tt) , ' --------'])
    disp(' ')
end


%% BEHAVIOR PLOTS: elpased time and velocity 

 mov_rew_time = [];mov_ave_time = []; no_mov_rew_time= []; no_mov_ave_time= [];
 
 mov_rew_speed = [];mov_ave_speed = [];  mov_rew_speed_res = [];mov_ave_speed_res = [];


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
        time_rew = stop-start; 
        movement.reward = InvertIntervals(behavior.quiet.reward , start , stop); %keep only those higher than criteria
        movement.reward(movement.reward(:,2) - movement.reward(:,1) <1,:)=[];
        clear tmp start stop
        % Aversive
        start = behavior.speed.aversive(1,1);   stop = behavior.speed.aversive(end,1);
        time_ave = stop-start; 
        movement.aversive = InvertIntervals(behavior.quiet.aversive , start , stop);%keep only those higher than criteria
        movement.aversive(movement.aversive(:,2) - movement.aversive(:,1) <1,:)=[];
        clear tmp start stop
        
        %Restrict pos to inside maze and only movement periods
        %---Aversive----
        pos_ave = behavior.pos.aversive(:,1:2); 
        pos_ave(:,2) = pos_ave(:,2)-min(pos_ave(:,2)); pos_ave(:,2) = pos_ave(:,2)/max(pos_ave(:,2)); %normalization of position 
        in_maze = ToIntervals(pos_ave(:,1),and(pos_ave(:,2)>0.09 , pos_ave(:,2)<1-0.09));% eliminating the extrems of the maze
                    
       %Check lenght of in_maze segments to define real laps
        in_lapA = []; % keeps only full laps (maybe too restrictive)
        for ix=1:size(in_maze,1)
                        ini = in_maze(ix,1);
                        fin = in_maze(ix,2);
                        ini_pos = pos_ave(pos_ave(:,1)==ini,2); 
                        fin_pos = pos_ave(pos_ave(:,1)==fin,2); 
                        if abs(fin_pos-ini_pos)>0.6
                             in_lapA = [in_lapA;in_maze(ix,:)]; 
                        end
        end
        
        pos_ave = Restrict(pos_ave,in_lapA);pos_ave = Restrict(pos_ave, movement.aversive); % Restrict pos to movement periods and laps
       
         %---Rewarded----
        pos_rew = behavior.pos.reward(:,1:2);
        pos_rew(:,2) = pos_rew(:,2)-min(pos_rew(:,2)); pos_rew(:,2) = pos_rew(:,2)/max(pos_rew(:,2)); %normalization of position 
        in_maze = ToIntervals(pos_rew(:,1),and(pos_rew(:,2)>0.09 , pos_rew(:,2)<1-0.09));% eliminating the extrems of the maze (10cm-reward box)
               
        %Check lenght of in_maze segments to define real laps
        in_lapR = []; 
        for ix=1:size(in_maze,1)
                        ini = in_maze(ix,1);
                        fin = in_maze(ix,2);
                        ini_pos = pos_rew(pos_rew(:,1)==ini,2); 
                        fin_pos = pos_rew(pos_rew(:,1)==fin,2); 
                        if abs(fin_pos-ini_pos)>0.6
                             in_lapR = [in_lapR;in_maze(ix,:)]; 
                        end
        end 
        pos_rew = Restrict(pos_rew,in_lapR); pos_rew = Restrict(pos_rew , movement.reward);
        
       %% Elapsed time
       
       mov_rew_time = [mov_rew_time;(sum(movement.reward(:,2)-movement.reward(:,1))/time_rew)*100];
       
       mov_ave_time = [mov_ave_time;(sum(movement.aversive(:,2)-movement.aversive(:,1))/time_ave)*100];
       
       no_mov_rew_time= [no_mov_rew_time; (sum(behavior.quiet.reward(:,2)-behavior.quiet.reward(:,1))/time_rew)*100]; 
       
       no_mov_ave_time= [no_mov_ave_time; (sum(behavior.quiet.aversive(:,2)-behavior.quiet.aversive(:,1))/time_ave)*100]; 
       
       
       
       mov_rew_speed = [mov_rew_speed;nanmean(behavior.speed.reward(:,2))];
       
       mov_ave_speed = [mov_ave_speed;nanmean(behavior.speed.aversive(:,2))];
       
       temp=Restrict(behavior.speed.reward,movement.reward);
       mov_rew_speed_res = [mov_rew_speed_res; nanmean(temp(:,2))];
       
       temp=Restrict(behavior.speed.aversive,movement.aversive);
       mov_ave_speed_res = [mov_ave_speed_res; nanmean(temp(:,2))];
       


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


%Elapsed time 
figure(4);clf;hold on, title('Elapsed time')
x = [ones(47,1);ones(47,1)*2;ones(47,1)*3;ones(47,1)*4];
y= [no_mov_ave_time;no_mov_rew_time;mov_ave_time;mov_rew_time];
c = [.3,.3,.3];
scatter(x,y,20,c,"filled",'jitter','on', 'jitterAmount',0.3); xlim([0.5 4.5]);
ylabel({'Elapsed time(%)'});
xticks([1 2 3 4])
xticklabels({'Ave', 'Rew','Ave', 'Rew'});
hold on

s1=scatter(1,nanmedian(no_mov_ave_time), "filled");s1.MarkerFaceColor = [1 0 0];
s2=scatter(2,nanmedian(no_mov_rew_time), "filled");s2.MarkerFaceColor = [0 0.2 1];

s3=scatter(3,nanmedian(mov_ave_time), "filled");s3.MarkerFaceColor = [1 0 0];
s4=scatter(4,nanmedian(mov_rew_time), "filled");s4.MarkerFaceColor = [0 0.2 1];

p = signrank(mov_ave_time,mov_rew_time)


%Elapsed time plots index 

index_mov = mov_rew_time./mov_ave_time; index_mov =[index_mov, ones(size(index_mov,1),1)];
index_no_mov = no_mov_rew_time./no_mov_ave_time;index_no_mov =[index_no_mov, 2*ones(size(index_no_mov,1),1)];
index= [index_mov;index_no_mov];

figure(1);clf;hold on, 
subplot(1,3,1);hold on; title('Elapsed time')
x = index(:,2);
y= index(:,1);

ylabel({'Elapsed time index rew/ave'});
yline(1);
xticks([1 2])
xticklabels({'Acrive', 'Quiet'});
hold on
s1=scatter(1,nanmedian(index_mov(:,1)), "filled");s1.MarkerFaceColor = [0 0 1];
s2=scatter(2,nanmedian(index_no_mov(:,1)), "filled");s2.MarkerFaceColor = [0 0.2 1];

[h,p] = kstest(index_mov(:,1))

p = signrank(index_mov(:,1),index_no_mov(:,1))


h = ttest(index_mov(:,1),1,'Alpha',0.05)

%Speed plots

index_speed =  mov_rew_speed./ mov_ave_speed; index_speed =[index_speed, ones(size(index_speed,1),1)];
index_speed_res =  mov_rew_speed_res./ mov_ave_speed_res; 
index_speed_res =[index_speed_res,ones(size(index_speed_res,1),1)*2];
index= [index_speed_res;index_speed];


% figure(3);clf;hold on, 
% sgtitle('Speed')
% 
subplot(1,3,2);hold on; title('All')
y = [mov_ave_speed;mov_rew_speed];
x= [ones(size(mov_ave_speed,1),1);2*ones(size(mov_rew_speed,1),1)];
c = [.3,.3,.3];
scatter(x,y,20,c,"filled",'jitter','on', 'jitterAmount',0.3); xlim([0.5 2.5]);ylim([2 25]);
ylabel({'Speed cm/sec'});
xticks([1 2])
xticklabels({'Ave', 'Rew'});
hold on
s1=scatter(1,nanmedian(mov_ave_speed), "filled");s1.MarkerFaceColor = [1 0 0];
s2=scatter(2,nanmedian(mov_rew_speed), "filled");s2.MarkerFaceColor = [0 0 1];


subplot(1,3,3);hold on; title('Active')
y = [mov_ave_speed_res;mov_rew_speed_res];
x= [ones(size(mov_ave_speed_res,1),1);2*ones(size(mov_rew_speed_res,1),1)];
c = [.3,.3,.3];
scatter(x,y,20,c,"filled",'jitter','on', 'jitterAmount',0.3); xlim([0.5 2.5]); ylim([2 25]);
ylabel({'Speed cm/sec'});
xticks([1 2])
xticklabels({'Ave', 'Rew'});
hold on
s1=scatter(1,nanmedian(mov_ave_speed_res), "filled");s1.MarkerFaceColor = [1 0 0];
s2=scatter(2,nanmedian(mov_rew_speed_res), "filled");s2.MarkerFaceColor = [0 0 1];


[h,p] = kstest(mov_ave_speed(:,1))
[h,p] = kstest(mov_ave_speed_res(:,1))
[h,p] = ttest(mov_rew_speed_res(:,1),mov_ave_speed_res(:,1))
[h,p] = ttest(mov_rew_speed(:,1),mov_ave_speed(:,1))
% plots 

figure(2);clf;hold on, 
sgtitle('Speed index')
x = index(:,2);
y= index(:,1);
c = [.3,.3,.3];
scatter(x,y,30,c,"filled",'jitter','on', 'jitterAmount',0.3); xlim([0.5 2.5]);
ylabel({'Speed index rew/ave'});
xticks([1 2])
xticklabels({'Active', 'All'});
hold on
s1=scatter(1,nanmedian(index_speed_res(:,1)), "filled");s1.MarkerFaceColor = [0 0 1];
s2=scatter(2,nanmedian(index_speed(:,1)), "filled");s2.MarkerFaceColor = [0 0.2 1];
yline(1);


p = ranksum(index_speed_res(:,1),index_speed(:,1))

p = ranksum(index_speed(:,1))

%% VELOCITY / FR CORRELATION 
%Spesific parameters
 bin_size = 1; %sec

%Output

dhpc_r_all= []; 
vhpc_r_all= []; 


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
        %Rewards selection -c1:open valve c2: close valve
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
        
        %%%%%%%%% Binned velocity 
        
        %Aversive
        start = behavior.speed.aversive(1,1);   stop = behavior.speed.aversive(end,1);
        bin_limits = [start:bin_size:stop]; 
        [mov_ave] = InIntervals(bin_limits,behavior.quiet.aversive); 
        mov_ave=~mov_ave;
        bin_velocity_ave= nan(size(bin_limits,2),1); 
        for v = 1:size(bin_velocity_ave,1)
            
          if v == size(bin_velocity_ave,1) % for the last bin
              ini = bin_limits(v); fin = stop; 
              [status,~,~] = InIntervals(behavior.speed.aversive(:,1),[ini,fin]);
              bin_velocity_ave(v)=nanmean(behavior.speed.aversive(status,2));
          else 
              
          ini = bin_limits(v); fin = bin_limits(v+1); 
          [status,~,~] = InIntervals(behavior.speed.aversive(:,1),[ini,fin]);
          bin_velocity_ave(v)=nanmean(behavior.speed.aversive(status,2));
          end
         end
        
        bin_velocity_ave=[bin_limits', bin_velocity_ave];
        mean_velocity_ave = nanmean(behavior.speed.aversive(:,2));
        
        %Reward
        start = behavior.speed.reward(1,1);   stop = behavior.speed.reward(end,1);
        bin_limits = [start:bin_size:stop]; 
        [mov_rew] = InIntervals(bin_limits,behavior.quiet.reward); 
        mov_rew=~mov_rew;
        bin_velocity_rew= nan(size(bin_limits,2),1); 
        for v = 1:size(bin_velocity_rew,1)
            
          if v == size(bin_velocity_rew,1)
              ini = bin_limits(v); fin = stop; 
              [status,~,~] = InIntervals(behavior.speed.reward(:,1),[ini,fin]);
              bin_velocity_rew(v)=nanmean(behavior.speed.reward(status,2));
          else 
          ini = bin_limits(v); fin = bin_limits(v+1); 
          [status,~,~] = InIntervals(behavior.speed.reward(:,1),[ini,fin]);
          bin_velocity_rew(v)=nanmean(behavior.speed.reward(status,2));
          end
        end
        
        bin_velocity_rew=[bin_limits', bin_velocity_rew];
        mean_velocity_rew = nanmean(behavior.speed.reward(:,2));
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

        %% Fr/velocity correlation 
        
        %%%%%%% PC
        try
            load('dHPC_pc_skaggs_circular.mat');
            dhpc_pc = [];
            for d=1:size(dHPC,2)
                pc = dHPC{d}; 
                dhpc_pc = [dhpc_pc;pc.id]; 
            end
            
            group_dHPC = group_dHPC(ismember(group_dHPC(:,1),dhpc_pc),:); % keep only pc 
             
            for d=1:size(group_dHPC,1)
                
                cluster = group_dHPC(d,1);
                spks = spks_dHPC(spks_dHPC(:,1)==cluster,2); % select tspk from cluster
                
                %Binned spikes
                [dN,~]=binspikes(spks,1/bin_size,[behavior.speed.aversive(1,1) behavior.speed.aversive(end,1)]); 
                binned_fr = dN./bin_size; binned_fr =binned_fr(1:end-1);
                
%                 [R_ave,~] = corrcoef(bin_velocity_ave(mov_ave,2),binned_fr(mov_ave)); R_ave=R_ave(1,2);
%                 
%                 figure(2);clf; hold on; scatter(bin_velocity_ave(mov_ave,2),binned_fr(mov_ave), 'filled'); xlabel('Velocity'); ylabel('Fr')
                dhpc_table = array2table([bin_velocity_ave(mov_ave,2),binned_fr(mov_ave)],'VariableNames',{'bin velocity','bin fr'}); 
                mdl = fitlm(dhpc_table);R_ave= mdl.Rsquared.Ordinary; 
%                 figure(1);plot(mdl)

                %Binned spikes
                [dN,~]=binspikes(spks,1/bin_size,[behavior.speed.reward(1,1) behavior.speed.reward(end,1)]); 
                binned_fr = dN./bin_size; binned_fr =binned_fr(1:end-1);
                
%                 [R_rew,~] = corrcoef(bin_velocity_rew(mov_rew,2),binned_fr(mov_rew));R_rew= R_rew(1,2);
%               figure(); hold on; scatter(bin_velocity_rew(mov_rew,2),binned_fr(mov_rew), 'filled'); xlabel('Velocity'); ylabel('Fr');
                

                dhpc_table = array2table([bin_velocity_rew(mov_rew,2),binned_fr(mov_rew)],'VariableNames',{'bin velocity','bin fr'}); 
                mdl = fitlm(dhpc_table);R_rew= mdl.Rsquared.Ordinary;
%                 figure; plot(mdl)
                
                dhpc_r_all= [dhpc_r_all;R_ave, R_rew]; 
                
            end

        catch 
            dHPC = []; disp(['No dHPC_pc.mat file in ',session]);
        end 
        
        try
            load('vHPC_pc_skaggs_circular.mat');
            vhpc_pc = [];
            for d=1:size(vHPC,2)
               pc = vHPC{d}; 
               vhpc_pc = [vhpc_pc;pc.id];
            end
            group_vHPC = group_vHPC(ismember(group_vHPC(:,1),vhpc_pc),:); % keep only pc 
            for d=1:size(group_vHPC,1)
                
                cluster = group_vHPC(d,1);
                spks = spks_vHPC(spks_vHPC(:,1)==cluster,2); % select tspk from cluster
                
                %Binned spikes
                [dN,~]=binspikes(spks,1/bin_size,[behavior.speed.aversive(1,1) behavior.speed.aversive(end,1)]); 
                binned_fr = dN./bin_size; binned_fr =binned_fr(1:end-1);
                
%                 [R_ave,~] = corrcoef(bin_velocity_ave(mov_ave,2),binned_fr(mov_ave)); R_ave=R_ave(1,2);
      
%                 figure; hold on; scatter(bin_velocity',binned_fr(1:end-1) , 'filled'); xlabel('Velocity'); ylabel('Fr')
               
                table = array2table([bin_velocity_ave(mov_ave,2),binned_fr(mov_ave)],'VariableNames',{'bin velocity','bin fr'}); 
                mdl = fitlm(table);R_ave= mdl.Rsquared.Ordinary; 
                
                %Binned spikes
                [dN,~]=binspikes(spks,1/bin_size,[behavior.speed.reward(1,1) behavior.speed.reward(end,1)]); 
                binned_fr = dN./bin_size;binned_fr =binned_fr(1:end-1);
                
%                 [R_rew,~] = corrcoef(bin_velocity_rew(mov_rew,2),binned_fr(mov_rew)); R_rew=R_rew(1,2);
%                 figure; hold on; scatter(bin_velocity_rew',binned_fr(1:end-1) , 'filled'); xlabel('Velocity'); ylabel('Fr')

                table = array2table([bin_velocity_rew(mov_rew,2),binned_fr(mov_rew)],'VariableNames',{'bin velocity','bin fr'}); 
                mdl = fitlm(table);R_rew= mdl.Rsquared.Ordinary;
                
                vhpc_r_all= [vhpc_r_all;R_ave,R_rew]; 
                
            end
            
        catch
            vHPC= []; disp(['No vHPC_pc.mat file in ',session]);
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


%%%%%Plot

%Cumilative 
figure(6);clf;hold on; cdfplot([dhpc_r_all(:,1);dhpc_r_all(:,2)]); 
cdfplot([vhpc_r_all(:,1);vhpc_r_all(:,2)]); 

%Correlational plot 
figure(1);clf; hold on;title('dhpc-green vhpc-blue');
scatter(dhpc_r_all(:,1),dhpc_r_all(:,2),[],[0 1 0],'filled','g'); xlabel('Rave'); ylabel('Rrew'); ylim([-0.2 0.8]);xlim([-0.2 0.6]);
scatter(vhpc_r_all(:,1),vhpc_r_all(:,2),[],[0 0 1],'filled','b'); xlabel('Rave'); ylabel('Rrew'); ylim([-0.2 0.8]);xlim([-0.2 0.6]);

% plot([-0.2 0.6],[-0.2 0.8]);

figure(2);clf; hold on; 
subplot(1,2,1);mdl_dhpc = fitlm(array2table([dhpc_r_all(:,1),dhpc_r_all(:,2)])); plot(mdl_dhpc); title('dHPC');xlabel('Rave'); ylabel('Rrew');
subplot(1,2,2);mdl_vhpc = fitlm(array2table([vhpc_r_all(:,1),vhpc_r_all(:,2)]));plot(mdl_vhpc); title('vHPC');xlabel('Rave'); ylabel('Rrew');

%%%%% Stats
[h,p] = kstest([dhpc_r_all(:,1);dhpc_r_all(:,2)]) %  no normal data 
[h,p] = kstest([vhpc_r_all(:,1);vhpc_r_all(:,2)]) %  no normal data 

[p,h] = signrank(dhpc_r_all(:,1),dhpc_r_all(:,2))
[p, h] = signrank(dhpc_r_all(:,1)); % h=1, p sig, reject null hypothesis -> data differ from 0
[p, h] = signrank(dhpc_r_all(:,2)); % h=1, p sig, reject null hypothesis -> data differ from 0

[p,h] = signrank(vhpc_r_all(:,1),vhpc_r_all(:,2))
[p, h] = signrank(vhpc_r_all(:,1)); % h=0,p no sig, data do not differ from 0
[p, h] = signrank(vhpc_r_all(:,2)); % h=1, p  sig, data  differ from 0

% Scatter plots 
figure(2);clf;hold on
sgtitle('Firing rate - velocity')
subplot(1,2,1); title('dHPC');
x = [ones(size(dhpc_r_all,1),1);ones(size(dhpc_r_all,1),1)*2];
y= [dhpc_r_all(:,1);dhpc_r_all(:,2)];

scatter(x,y,30,"filled",'jitter','on', 'jitterAmount',0.3); xlim([0.5 2.5]); ylim([-0.2 0.8]);
ylabel('r2');
xticks([1 2]); xticklabels({'Ave', 'Rew'});
hold on
s1=scatter(1,nanmedian(dhpc_r_all(:,1)), "filled");s1.MarkerFaceColor = [1 0 0];
s2=scatter(2,nanmedian(dhpc_r_all(:,2)), "filled");s2.MarkerFaceColor = [0 0 1];
 
subplot(1,2,2); title('vHPC');
x = [ones(size(vhpc_r_all,1),1);ones(size(vhpc_r_all,1),1)*2];
y= [vhpc_r_all(:,1);vhpc_r_all(:,2)];

scatter(x,y,30,"filled",'jitter','on', 'jitterAmount',0.3); xlim([0.5 2.5]); ylim([-0.2 0.8]);
ylabel('r2');
xticks([1 2]); xticklabels({'Ave', 'Rew'});
hold on
s1=scatter(1,nanmedian(vhpc_r_all(:,1)), "filled");s1.MarkerFaceColor = [1 0 0];
s2=scatter(2,nanmedian(vhpc_r_all(:,2)), "filled");s2.MarkerFaceColor = [0 0 1];

%%%%%Ratio r2 ave/ r2 rew
figure(3);clf;hold on
sgtitle('R2 ratio')
x = [ones(size(dhpc_r_all,1),1);ones(size(vhpc_r_all,1),1)*2];
y= [dhpc_r_all(:,1)./dhpc_r_all(:,2);vhpc_r_all(:,1)./vhpc_r_all(:,2)];

scatter(x,y,30,"filled",'jitter','on', 'jitterAmount',0.3); xlim([0.5 2.5]); ylim([-100 100]);
ylabel('r2 ratio');
xticks([1 2]); xticklabels({'dHPC', 'vHPC'});
hold on
s1=scatter(1,nanmedian(dhpc_r_all(:,1)./dhpc_r_all(:,2)), "filled");s1.MarkerFaceColor = [1 0 0];
s2=scatter(2,nanmedian(vhpc_r_all(:,1)./vhpc_r_all(:,2)), "filled");s2.MarkerFaceColor = [0 0 1];
yline(0); 

%stats 
[p,h] = ranksum(dhpc_r_all(:,1)./dhpc_r_all(:,2),vhpc_r_all(:,1)./vhpc_r_all(:,2))
[p, h] = signrank(dhpc_r_all(:,1)./dhpc_r_all(:,2)); % h=1, p sig, reject null hypothesis -> data differ from 0
[p, h] = signrank(vhpc_r_all(:,1)./vhpc_r_all(:,2)); % h=1, p sig, reject null hypothesis -> data differ from 0

%% SUBSTRACT NEURONS WITH HIGH VELOCITY/FR CORRELATION 
% using output of previous section
 figure(1);clf; 
 sgtitle('dHPC')
 Q = quantile(dhpc_r_all(:,1),[0.025 0.25 0.5 0.75 0.975]); 
 subplot(1,2,1);h=histogram(dhpc_r_all(:,1));h.EdgeColor = 'none';h.FaceColor=[.3,.3,.3];xlim([0 0.4]);ylim([0 450]);
 xl = xline(Q(4),'r');
 xlabel('r2 aversive'); ylabel('Counts');legend(xl,'Q 0.75')
 subplot(1,2,2);h=histogram(dhpc_r_all(:,2));xlim([0 0.4]);ylim([0 450]);h.EdgeColor = 'none';h.FaceColor=[.3,.3,.3]; 
 xline(Q(4),'r'); xlabel('r2 reward'); ylabel('Counts');
 tresh_dhpc = Q(4); 
 
 figure(1);clf; 
 sgtitle('vHPC')
 subplot(1,2,1);h= histogram(vhpc_r_all(:,1));h.EdgeColor = 'none';h.FaceColor=[.3,.3,.3];xlim([0 0.2]);ylim([0 130])
 Q = quantile(vhpc_r_all(:,1),[0.025 0.25 0.5 0.75 0.975]); 
 xl = xline(Q(4),'r'); xlabel('r2 aversive'); ylabel('Counts');legend(xl,'Q 0.75')
 subplot(1,2,2);h= histogram(vhpc_r_all(:,2));h.EdgeColor = 'none';h.FaceColor=[.3,.3,.3];xlim([0 0.2]);ylim([0 130]);
 xline(Q(4),'r'); xlabel('r2 reward'); ylabel('Counts');
 tresh_vhpc = Q(4); 
 

%Spesific parameters
 bin_size = 1; %sec
%Output
dhpc_params= []; 
vhpc_params= []; 

% not the fastest way of doing this 
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
        %Rewards selection -c1:open valve c2: close valve
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
        
        %%%%%%%%% Binned velocity 
        
        %Aversive
        start = behavior.speed.aversive(1,1);   stop = behavior.speed.aversive(end,1);
        bin_limits = [start:bin_size:stop]; 
        [mov_ave] = InIntervals(bin_limits,behavior.quiet.aversive); 
        
        bin_velocity_ave= nan(size(bin_limits,2),1); 
        for v = 1:size(bin_velocity_ave,1)
            
          if v == size(bin_velocity_ave,1)
              ini = bin_limits(v); fin = stop; 
              [status,~,~] = InIntervals(behavior.speed.aversive(:,1),[ini,fin]);
              bin_velocity_ave(v)=nanmean(behavior.speed.aversive(status,2));
          else 
          ini = bin_limits(v); fin = bin_limits(v+1); 
          [status,~,~] = InIntervals(behavior.speed.aversive(:,1),[ini,fin]);
          bin_velocity_ave(v)=nanmean(behavior.speed.aversive(status,2));
          end
         end
        
        bin_velocity_ave=[bin_limits', bin_velocity_ave];
        mean_velocity_ave = nanmean(behavior.speed.aversive(:,2));
        
        %Reward
        start = behavior.speed.reward(1,1);   stop = behavior.speed.reward(end,1);
        bin_limits = [start:bin_size:stop]; 
        [mov_rew] = InIntervals(bin_limits,behavior.quiet.reward); 
        bin_velocity_rew= nan(size(bin_limits,2),1); 
        for v = 1:size(bin_velocity_rew,1)
            
          if v == size(bin_velocity_rew,1)
              ini = bin_limits(v); fin = stop; 
              [status,~,~] = InIntervals(behavior.speed.reward(:,1),[ini,fin]);
              bin_velocity_rew(v)=nanmean(behavior.speed.reward(status,2));
          else 
          ini = bin_limits(v); fin = bin_limits(v+1); 
          [status,~,~] = InIntervals(behavior.speed.reward(:,1),[ini,fin]);
          bin_velocity_rew(v)=nanmean(behavior.speed.reward(status,2));
          end
        end
        
        bin_velocity_rew=[bin_limits', bin_velocity_rew];
        mean_velocity_rew = nanmean(behavior.speed.reward(:,2));
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

        %% Fr/velocity correlation 
        
        %%%%%%% PC
        try
            load('dHPC_pc_skaggs_circular.mat');
            dhpc_pc = [];
            for d=1:size(dHPC,2)
                pc = dHPC{d}; 
                dhpc_pc = [dhpc_pc;pc.id]; 
            end
            
            group_dHPC = group_dHPC(ismember(group_dHPC(:,1),dhpc_pc),:); % keep only pc 
             
            for d=1:size(group_dHPC,1)
                
                cluster = group_dHPC(d,1);
                spks = spks_dHPC(spks_dHPC(:,1)==cluster,2); % select tspk from cluster
                
                %Binned spikes
                [dN,~]=binspikes(spks,1/bin_size,[behavior.speed.aversive(1,1) behavior.speed.aversive(end,1)]); 
                binned_fr = dN./bin_size; binned_fr =binned_fr(1:end-1);
                
                [R_ave,~] = corrcoef(bin_velocity_ave(mov_ave,2),binned_fr(mov_ave)); 
      
%                 figure(4); hold on; scatter(bin_velocity_ave(mov_ave,2),binned_fr(mov_ave), 'filled'); xlabel('Velocity'); ylabel('Fr')
               
                %Binned spikes
                [dN,~]=binspikes(spks,1/bin_size,[behavior.speed.reward(1,1) behavior.speed.reward(end,1)]); 
                binned_fr = dN./bin_size; binned_fr =binned_fr(1:end-1);
                
                [R_rew,~] = corrcoef(bin_velocity_rew(mov_rew,2),binned_fr(mov_rew)); 
%                 figure; hold on; scatter(bin_velocity_rew',binned_fr(1:end-1) , 'filled'); xlabel('Velocity'); ylabel('Fr')
                
                if R_ave(1,2) && R_rew(1,2) < tresh_dhpc
                    dhpc_params= [dhpc_params;[dHPC{d}.between,1;dHPC{d}.within_ave,2;dHPC{d}.within_rew,3]]; 
                end 
            end

        catch 
            dHPC = []; disp(['No dHPC_pc.mat file in ',session]);
        end 
        
        try
            load('vHPC_pc_skaggs_circular.mat');
            vhpc_pc = [];
            for d=1:size(vHPC,2)
               pc = vHPC{d}; 
               vhpc_pc = [vhpc_pc;pc.id];
            end
            group_vHPC = group_vHPC(ismember(group_vHPC(:,1),vhpc_pc),:); % keep only pc 
            for d=1:size(group_vHPC,1)
                
                cluster = group_vHPC(d,1);
                spks = spks_vHPC(spks_vHPC(:,1)==cluster,2); % select tspk from cluster
                
                %Binned spikes
                [dN,~]=binspikes(spks,1/bin_size,[behavior.speed.aversive(1,1) behavior.speed.aversive(end,1)]); 
                binned_fr = dN./bin_size; binned_fr =binned_fr(1:end-1);
                
                [R_ave,~] = corrcoef(bin_velocity_ave(mov_ave,2),binned_fr(mov_ave));  
      
%                 figure; hold on; scatter(bin_velocity',binned_fr(1:end-1) , 'filled'); xlabel('Velocity'); ylabel('Fr')
               
                %Binned spikes
                [dN,~]=binspikes(spks,1/bin_size,[behavior.speed.reward(1,1) behavior.speed.reward(end,1)]); 
                binned_fr = dN./bin_size;binned_fr =binned_fr(1:end-1);
                
                [R_rew,~] = corrcoef(bin_velocity_rew(mov_rew,2),binned_fr(mov_rew)); 
%                 figure; hold on; scatter(bin_velocity_rew',binned_fr(1:end-1) , 'filled'); xlabel('Velocity'); ylabel('Fr')
                
               if R_ave(1,2) && R_rew(1,2) < tresh_vhpc
                    vhpc_params= [vhpc_params;[vHPC{d}.between,1;vHPC{d}.within_ave,2;vHPC{d}.within_rew,3]]; 
                end  
                
            end
            
        catch
            vHPC= []; disp(['No vHPC_pc.mat file in ',session]);
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
 

%Only 3 plots
ylabels ={'Spatial correlation', 'Overlap', 'Pf shift'};
idx = [1,3,4]; 
xlabels = {'Between', 'Within Ave', 'Within Rew'}; 
figure(1);clf;hold on, 
sgtitle('dHPC')
for a=1:size(ylabels,2)
    p = idx(a);
    subplot(1,3,a); hold on; 
    x = dhpc_params(:,6);
    y= dhpc_params(:,p);%Change the column number (1-4) to choose which variable to plot 
    c = [.3,.3,.3];
    scatter(x,y,30,c,"filled",'jitter','on', 'jitterAmount',0.3); xlim([0 4]);
    ylabel(ylabels{a});
    xticks([1 2 3])
    xticklabels({'Bet', 'WA', 'WR'});
    hold on
    s1=scatter(1,nanmedian(dhpc_params(dhpc_params(:,6)==1,p)), "filled");s1.MarkerFaceColor = [0 0 0];
    s2=scatter(2,nanmedian(dhpc_params(dhpc_params(:,6)==2,p)), "filled");s2.MarkerFaceColor = [1 0.1 0.2];
    s3=scatter(3,nanmedian(dhpc_params(dhpc_params(:,6)==3,p)), "filled");s3.MarkerFaceColor = [0.1 0.3 1];
end

%Stats
[P,ANOVATAB,STATS] = kruskalwallis(dhpc_sub(:,5),dhpc_sub(:,6));
c = multcompare(STATS)

tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])

% Paired non-parametric anova - Friedman test
%Prepare data to test
p = 1; % 1=sp corr 2=fr change 3=overlap 4=pf shift
data=[dhpc_params(dhpc_params(:,6)==1,p), dhpc_params(dhpc_params(:,6)==2,p), ...
    dhpc_params(dhpc_params(:,6)==3,p), dhpc_params(dhpc_params(:,6)==4,p),dhpc_params(dhpc_params(:,6)==5,p)];
data(any(isnan(data), 2), :) = [];

[p,~,stats] =  friedman(data,1); 
c = multcompare(stats);
tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])

%Only 3 plots
ylabels ={'Spatial correlation', 'Overlap', 'Pf shift'};
idx = [1,3,4]; 
xlabels = {'Between', 'Within Ave', 'Within Rew'}; 
figure(1);clf;hold on, 
sgtitle('vHPC')
for a=1:size(ylabels,2)
    p = idx(a);
    subplot(1,3,a); hold on; 
    x = vhpc_params(:,6);
    y= vhpc_params(:,p);%Change the column number (1-4) to choose which variable to plot 
    c = [.3,.3,.3];
    scatter(x,y,30,c,"filled",'jitter','on', 'jitterAmount',0.3); xlim([0 4]);
    ylabel(ylabels{a});
    xticks([1 2 3])
    xticklabels({'Bet', 'WA', 'WR'});
    hold on
    s1=scatter(1,nanmedian(vhpc_params(vhpc_params(:,6)==1,p)), "filled");s1.MarkerFaceColor = [0 0 0];
    s2=scatter(2,nanmedian(vhpc_params(vhpc_params(:,6)==2,p)), "filled");s2.MarkerFaceColor = [1 0.1 0.2];
    s3=scatter(3,nanmedian(vhpc_params(vhpc_params(:,6)==3,p)), "filled");s3.MarkerFaceColor = [0.1 0.3 1];
end

%Stats
[P,ANOVATAB,STATS] = kruskalwallis(dhpc_sub(:,5),dhpc_sub(:,6));
c = multcompare(STATS)

tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])

% Paired non-parametric anova - Friedman test
%Prepare data to test
p = 1; % 1=sp corr 2=fr change 3=overlap 4=pf shift
data=[vhpc_params(vhpc_params(:,6)==1,p), vhpc_params(vhpc_params(:,6)==2,p), ...
    vhpc_params(vhpc_params(:,6)==3,p), vhpc_params(vhpc_params(:,6)==4,p),vhpc_params(vhpc_params(:,6)==5,p)];
data(any(isnan(data), 2), :) = [];

[p,~,stats] =  friedman(data,1); 
c = multcompare(stats);
tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])

%% VELOCITY - REMAPPING CORRELATION
%Spesific parameters
 bin_size = 1; %sec

%Output

dhpc_params= []; 
vhpc_params= []; 


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
        %Rewards selection -c1:open valve c2: close valve
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
        
        %%%%%%%%% Binned velocity 
        
        %Aversive
        start = behavior.speed.aversive(1,1);   stop = behavior.speed.aversive(end,1);
        bin_limits = [start:bin_size:stop]; 
        [mov_ave] = InIntervals(bin_limits,behavior.quiet.aversive); 
        
        bin_velocity_ave= nan(size(bin_limits,2),1); 
        for v = 1:size(bin_velocity_ave,1)
            
          if v == size(bin_velocity_ave,1)
              ini = bin_limits(v); fin = stop; 
              [status,~,~] = InIntervals(behavior.speed.aversive(:,1),[ini,fin]);
              bin_velocity_ave(v)=nanmean(behavior.speed.aversive(status,2));
          else 
          ini = bin_limits(v); fin = bin_limits(v+1); 
          [status,~,~] = InIntervals(behavior.speed.aversive(:,1),[ini,fin]);
          bin_velocity_ave(v)=nanmean(behavior.speed.aversive(status,2));
          end
         end
        
        bin_velocity_ave=[bin_limits', bin_velocity_ave];
        mean_velocity_ave = nanmean(behavior.speed.aversive(:,2));
        
        %Reward
        start = behavior.speed.reward(1,1);   stop = behavior.speed.reward(end,1);
        bin_limits = [start:bin_size:stop]; 
        [mov_rew] = InIntervals(bin_limits,behavior.quiet.reward); 
        bin_velocity_rew= nan(size(bin_limits,2),1); 
        for v = 1:size(bin_velocity_rew,1)
            
          if v == size(bin_velocity_rew,1)
              ini = bin_limits(v); fin = stop; 
              [status,~,~] = InIntervals(behavior.speed.reward(:,1),[ini,fin]);
              bin_velocity_rew(v)=nanmean(behavior.speed.reward(status,2));
          else 
          ini = bin_limits(v); fin = bin_limits(v+1); 
          [status,~,~] = InIntervals(behavior.speed.reward(:,1),[ini,fin]);
          bin_velocity_rew(v)=nanmean(behavior.speed.reward(status,2));
          end
        end
        
        bin_velocity_rew=[bin_limits', bin_velocity_rew];
        mean_velocity_rew = nanmean(behavior.speed.reward(:,2));
        %% Saving remapping parameters 
        cd 'Spikesorting'
      
        try
            load('dHPC_pc_skaggs_circular.mat');

             
            for d=1:size(dHPC,2)
                dhpc_params= [dhpc_params;[dHPC{d}.between,mean_velocity_ave, mean_velocity_rew]]; 
            end

        catch 
            dHPC = []; disp(['No dHPC_pc.mat file in ',session]);
        end 
        
        try
            load('vHPC_pc_skaggs_circular.mat');
            
            for d=1:size(vHPC,2)
                vhpc_params= [vhpc_params;[vHPC{d}.between,mean_velocity_ave, mean_velocity_rew]];
  
            end
            
        catch
            vHPC= []; disp(['No vHPC_pc.mat file in ',session]);
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

dhpc_table = array2table([(dhpc_params(:,7)-dhpc_params(:,6))./(dhpc_params(:,7)+dhpc_params(:,6)),dhpc_params(:,3)],'VariableNames',{'Velocity change','overlap'}); 
mdl = fitlm(dhpc_table)

vhpc_table = array2table([(vhpc_params(:,7)-vhpc_params(:,6))./(vhpc_params(:,7)+vhpc_params(:,6)),vhpc_params(:,1)],'VariableNames',{'Velocity change','overlap'}); 
mdl = fitlm(vhpc_table)

plot(mdl)
 %plot  
%% REMAPPING DORSAL vs VENTRAL - MEAN WITHIN 
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
            load('dHPC_pc_skaggs_circular.mat');%'dHPC_pc_skaggs_pos.mat'
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
            load('vHPC_pc_skaggs_circular.mat');%'vHPC_pc_skaggs_pos.mat'
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

%%%%%%%%% Dividing by genereal mean 
%DHPC 
%Compute general mean within 
within_all = dhpc_sub(or(dhpc_sub(:,6)== 2,dhpc_sub(:,6)== 3),1:5); 
within_mean = nanmean(within_all);
%Divide by the mean 
between_all = dhpc_sub(dhpc_sub(:,6)== 1,1:5); 
between_mean_dhpc = between_all ./ within_mean;

%vHPC 
%Compute general mean within 
within_all = vhpc_sub(or(vhpc_sub(:,6)== 2,vhpc_sub(:,6)== 3),1:5); 
within_mean = nanmean(within_all);
%Divide by the mean 
between_all = vhpc_sub(vhpc_sub(:,6)== 1,1:5); 
between_mean_vhpc = between_all ./ within_mean; 

%%%%%% Plot dhpc vs vhpc
%Only 3 plots
ylabels ={'Spatial correlation', 'Overlap', 'Pf shift'};
idx = [1,3,4]; 
xlabels = {'dHPC', 'vHPC'}; 
figure(1);clf;hold on, 
sgtitle('dHPC vs vHPC')
for a=1:size(ylabels,2)
    p = idx(a);
    subplot(1,3,a); hold on; 
    x = [ones(size(between_mean_dhpc(:,p),1),1);2*ones(size(between_mean_vhpc(:,p),1),1)];
    y= [between_mean_dhpc(:,p);between_mean_vhpc(:,p)];%Change the column number (1-4) to choose which variable to plot 
    c = [.3,.3,.3];
    scatter(x,y,30,c,"filled",'jitter','on', 'jitterAmount',0.3); xlim([0.5 2.5]);
    ylabel(ylabels{a});
    xticks([1 2])
    xticklabels(xlabels);
    hold on
    s1=scatter(1,nanmedian(between_mean_dhpc(:,p)), "filled");s1.MarkerFaceColor = [0.3 1 0.2];
    s2=scatter(2,nanmedian(between_mean_vhpc(:,p)), "filled");s2.MarkerFaceColor = [0 0.3 1];
    
end

%Stats
p=4; % 1=sp corr 2=fr change 3=overlap 4=pf shift

x = [ones(size(between_mean_dhpc(:,p),1),1);2*ones(size(between_mean_vhpc(:,p),1),1)];
y= [between_mean_dhpc(:,p);between_mean_vhpc(:,p)];%Change the column number (1-4) to choose which variable to plot 

[P,ANOVATAB,STATS] = kruskalwallis(y,x);
c = multcompare(STATS)

tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])


%%%%%%%%% Dividing each neuron by its mean 

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
            load('dHPC_pc_skaggs_circular.mat');%'dHPC_pc_skaggs_pos.mat'
             %Save session parameters 
        for d=1:size(dHPC,2)
            pc = dHPC{d}; 
            within = nanmean([pc.subsampled.within_ave;pc.subsampled.within_rew]);
            bet = pc.subsampled.between./ within; 
            dhpc_temp = [dhpc_temp;bet]; 
            
        end
        catch 
            dHPC = []; disp(['No dHPC_pc_lap.mat file in ',session]);
        end 
        
        try
            load('vHPC_pc_skaggs_circular.mat');%'vHPC_pc_skaggs_pos.mat'
            for d=1:size(vHPC,2)
                pc = vHPC{d}; 
                within = nanmean([pc.subsampled.within_ave;pc.subsampled.within_rew]);
                bet = pc.subsampled.between./ within; 
                vhpc_temp = [vhpc_temp;bet]; 
            end 
     
        catch
            vHPC= []; disp(['No vHPC_pc_lap.mat file in ',session]);
        end
        
    end 
    %Save sub parameters 
    dhpc_sub = [dhpc_sub;dhpc_temp];
    vhpc_sub = [vhpc_sub;vhpc_temp];
     
end 



%%%%%% Plot dhpc vs vhpc
%Only 3 plots
ylabels ={'Spatial correlation', 'Overlap', 'Pf shift'};
idx = [1,3,4]; 
xlabels = {'dHPC', 'vHPC'}; 
figure(1);clf;hold on, 
sgtitle('dHPC vs vHPC')
for a=1:size(ylabels,2)
    p = idx(a);
    subplot(1,3,a); hold on; 
    x = [ones(size(dhpc_sub(:,p),1),1);2*ones(size(vhpc_sub(:,p),1),1)];
    y= [dhpc_sub(:,p);vhpc_sub(:,p)];%Change the column number (1-4) to choose which variable to plot 
    c = [.3,.3,.3];
    scatter(x,y,30,c,"filled",'jitter','on', 'jitterAmount',0.3); xlim([0.5 2.5]);
    ylabel(ylabels{a});
    xticks([1 2])
    xticklabels(xlabels);
    hold on
    s1=scatter(1,nanmedian(dhpc_sub(:,p)), "filled");s1.MarkerFaceColor = [0.3 1 0.2];
    s2=scatter(2,nanmedian(vhpc_sub(:,p)), "filled");s2.MarkerFaceColor = [0 0.3 1];
    
end

%Stats
p=4; % 1=sp corr 2=fr change 3=overlap 4=pf shift

x = [ones(size(dhpc_sub(:,p),1),1);2*ones(size(vhpc_sub(:,p),1),1)];
y= [dhpc_sub(:,p);vhpc_sub(:,p)];%Change the column number (1-4) to choose which variable to plot 

[P,ANOVATAB,STATS] = kruskalwallis(y,x);
c = multcompare(STATS)

tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])


%% SPATIAL SUBSAMPLING FR 

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
        
        %Restrict pos to inside maze and only movement periods
        %---Aversive----
        pos_ave = behavior.pos.aversive(:,1:2); 
        pos_ave(:,2) = pos_ave(:,2)-min(pos_ave(:,2)); pos_ave(:,2) = pos_ave(:,2)/max(pos_ave(:,2)); %normalization of position 
        in_maze = ToIntervals(pos_ave(:,1),and(pos_ave(:,2)>0.09 , pos_ave(:,2)<1-0.09));% eliminating the extrems of the maze
                    
       %Check lenght of in_maze segments to define real laps
        in_lapA = []; % keeps only full laps (maybe too restrictive)
        for ix=1:size(in_maze,1)
                        ini = in_maze(ix,1);
                        fin = in_maze(ix,2);
                        ini_pos = pos_ave(pos_ave(:,1)==ini,2); 
                        fin_pos = pos_ave(pos_ave(:,1)==fin,2); 
                        if abs(fin_pos-ini_pos)>0.6
                             in_lapA = [in_lapA;in_maze(ix,:)]; 
                        end
        end
        
        pos_ave = Restrict(pos_ave,in_lapA);pos_ave = Restrict(pos_ave, movement.aversive); % Restrict pos to movement periods and laps
       
         %---Rewarded----
        pos_rew = behavior.pos.reward(:,1:2);
        pos_rew(:,2) = pos_rew(:,2)-min(pos_rew(:,2)); pos_rew(:,2) = pos_rew(:,2)/max(pos_rew(:,2)); %normalization of position 
        in_maze = ToIntervals(pos_rew(:,1),and(pos_rew(:,2)>0.09 , pos_rew(:,2)<1-0.09));% eliminating the extrems of the maze (10cm-reward box)
               
        %Check lenght of in_maze segments to define real laps
        in_lapR = []; 
        for ix=1:size(in_maze,1)
                        ini = in_maze(ix,1);
                        fin = in_maze(ix,2);
                        ini_pos = pos_rew(pos_rew(:,1)==ini,2); 
                        fin_pos = pos_rew(pos_rew(:,1)==fin,2); 
                        if abs(fin_pos-ini_pos)>0.6
                             in_lapR = [in_lapR;in_maze(ix,:)]; 
                        end
        end 
        pos_rew = Restrict(pos_rew,in_lapR); pos_rew = Restrict(pos_rew , movement.reward);
        
        
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
        disp('dHPC PC parameters')
        
         %%%%%%%% DHPC %%%%%%%%%%%
        try
            load('dHPC_pc_skaggs_circular.mat');
            dhpc_pc = [];
            for d=1:size(dHPC,2)
                pc = dHPC{d}; 
                dhpc_pc = [dhpc_pc;pc.id]; 
            end
            group_dHPC = group_dHPC(ismember(group_dHPC(:,1),dhpc_pc),:); % keep only pc 
               
            for n=1:size(group_dHPC,1)
                cluster = group_dHPC(n,1);
                spks = spks_dHPC(spks_dHPC(:,1)==cluster,2); % select tspk from cluster

                % --- Aversive ---
                spks_ave = spks; 
                %Restrict spk and position to laps and movement periods
                spks_ave = Restrict(spks_ave,in_lapA);spks_ave = Restrict(spks_ave , movement.aversive); % Restrict spikes to movement periods and laps
               
                % --- Reward ---
                spks_rew = spks; 
                %Restrict spk  to laps and movement periods
                spks_rew = Restrict(spks_rew,in_lapR); spks_rew = Restrict(spks_rew, movement.reward); % Restrict to movement periods
               
                [sub_ratemap_ave, sub_ratemap_rew] = Subsampling_fr_1d(pos_ave, spks_ave,pos_rew, spks_rew,dt,Xedges,sigma);
                
                
                
            end
                
               
        catch 
           disp(['No dHPC_pc.mat file in ',session]);
        end 
        pv_corr.dhpc_bet = [pv_corr.dhpc_bet;temp_dhpc_bet];
        pv_corr.dhpc_within_ave = [pv_corr.dhpc_within_ave;temp_dhpc_within_ave];
        pv_corr.dhpc_within_rew = [pv_corr.dhpc_within_rew;temp_dhpc_within_rew];
           
        disp('vHPC PC parameters')
     
      
    %% Saveing PC INFO 
    if ~isempty(dHPC)
        dHPC = dHPC(~cellfun('isempty',dHPC));
        save([cd,'\dHPC_pc_skaggs_circular.mat'],'dHPC'); 
    end
    if ~isempty(vHPC)
        vHPC = vHPC(~cellfun('isempty',vHPC));
        save([cd,'\vHPC_pc_skaggs_circular.mat'],'vHPC'); 
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





