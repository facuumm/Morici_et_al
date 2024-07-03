%% Stability pipeline -compare lap by lap stability
clear
clc
%% Parameters
path = {'\\Maryjackson\e\Rat127\Ephys\pyr';'\\Maryjackson\e\Rat128\Ephys\in_pyr\ready';'\\Maryjackson\e\Rat103\usable';'\\Maryjackson\e\Rat132\recordings\in_pyr'; '\\Maryjackson\e\Rat165\in_pyr'; ...
    '\\Maryjackson\e\Rat126\Ephys\in_Pyr'};%List of folders from the path
% path = {'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat103\usable';'E:\Rat132\recordings\in_pyr'};%List of folders from the path
% for SU
criteria_fr = 0.01; %criteria to include or not a SU into the analysis
criteria_n = 6; % minimal number of neurons from each structure
criteria_type = 0; %criteria for celltype (0:pyr, 1:int, 2:all)
n_SU_V = [];
n_SU_D = [];
Xedges = 60; %number of bins for RateMap construction - 3cm bin
sigma = 2;%round(15/(180/Xedges)); %defined for gauss kernel of 15cm
binSize = 0.001; % bin size for replay events detection
bin_size = 1; % to bin pos ans spks in between/within  

% Behavior
minimal_speed = 2.5;% minimal speed to detect quite periods
minimal_speed_time = 2; % minimal time to detect quite periods

 
%% Choose the number of laps to compare 
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
histogram(dhpc_nlap(:,1),'FaceColor','r','EdgeColor','none','BinWidth',1);
histogram(dhpc_nlap(:,2),'FaceColor','b','EdgeColor','none','BinWidth',1);
ylim([0 7]);ylabel('Counts'); xlabel('# laps'); 

figure(4);clf;hold on;
sgtitle('vHPC laps')
histogram(vhpc_nlap(:,1),'FaceColor','r','EdgeColor','none','BinWidth',1);
histogram(vhpc_nlap(:,2),'FaceColor','b','EdgeColor','none','BinWidth',1);
ylim([0 7]);ylabel('Counts'); xlabel('# laps'); 


%I will take sessions with 15 laps or more for AVERSIVE and sessions with
%30 laps or more for REWARD
 

%% Main Loop: stability across time
% Stability parameters 
tresh_lap_ave = 15;
tresh_lap_rew = 30; 
template_laps = 5; %# of laps to contruct the tempalte and compare stability. 

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
        
        %Check if the session has the sufficient # laps 
        if isfile('Spikesorting\dHPC_pc.mat')
            load('Spikesorting\dHPC_pc.mat'); 
            nlap_ave= dHPC{1}.nlap_ave; nlap_rew =dHPC{1}.nlap_rew;
        elseif isfile('Spikesorting\vHPC_pc.mat')
            load('Spikesorting\vHPC_pc.mat');
            nlap_ave= vHPC{1}.nlap_ave; nlap_rew =vHPC{1}.nlap_rew;
        else
            nlap_ave= 0; nlap_rew =0;
        end 
        
        if or(nlap_ave >= tresh_lap_ave,nlap_rew >= tresh_lap_rew)
         
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
            dHPC_stab = {};% one cell per pc
            vHPC_stab = {}; 
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
        %% Firing curve 
        disp('dHPC Stability calculation')
        try
           load('dHPC_pc.mat')
           %Store pc identity
            pc_dhpc = [];
            for pc=1:size(dHPC,2)
                pc_dhpc = [pc_dhpc;dHPC{pc}.id];
            end
        catch
           dHPC = []; disp('No dhpc pc')
        end 
        
        
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
                
               % AVERSIVE 
               [curveA ] = FiringCurve(pos_ave , spks_ave , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);

                % condicional si nlap_ave supera el tresh de laps ave
                if nlap_ave>= tresh_lap_ave
                    in_lapA = in_lapA(1:tresh_lap_ave,:); % keep the first tresh laps 
                    %Construct a template with the last # template_laps laps 
                    template = in_lapA(end-template_laps-1:end,:); spks_t = Restrict(spks_ave,template);pos_t = Restrict(pos_ave,template);  
                    [curveT , statsT] = FiringCurve(pos_t , spks_t , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    %%%%%%%%%%%%%%%%%%%%%control plot%%%%%%%%%%%%%%%%%%%%%
%                     if mod(ii,2)% plot only odds 
%                         figure;hold on 
%                         nc = 2+tresh_lap_ave-template_laps;
%                         subplot(nc,1,1);imagesc(curveA.rate), colormap 'jet'
%                         title('Firing curvs original AVERSIVE');
%                         subplot(nc,1,2);imagesc(curveT.rate), colormap 'jet'
%                         title('Firing curvs template');
%                     end 
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    stab_ave=[];
                    for c=1:size(in_lapA,1)-template_laps
                        lap = in_lapA(c,:); 
                        spks_lap = Restrict(spks_ave,lap);pos_lap = Restrict(pos_ave,lap);
                        [curveL , statsL] = FiringCurve(pos_lap , spks_lap , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                        %Spatial correlation
                        s = corrcoef(curveT.rate, curveL.rate);s = s(1,2);
                        stab_ave=[stab_ave,s];
%                         if mod(ii,2)% plot only odds
%                             subplot(nc,1,c+2);imagesc(curveL.rate), colormap 'jet'
%                         end    
                    end 
                end 
                
                % REWARD
                [curveR] = FiringCurve(pos_rew , spks_rew , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);

                % condicional si nlap_ave supera el tresh de laps ave
                if nlap_rew>= tresh_lap_rew
                    in_lapR = in_lapR(1:tresh_lap_rew,:); % keep the first tresh laps 
                    %Construct a template with the last # template_laps laps 
                    template = in_lapR(end-template_laps-1:end,:); spks_t = Restrict(spks_rew,template);pos_t = Restrict(pos_rew,template);  
                    [curveT , ] = FiringCurve(pos_t , spks_t , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    %%%%%%%%%%%%%%%%%%%%%control plot%%%%%%%%%%%%%%%%%%%%%
%                     if mod(ii,2)% plot only odds 
%                         figure;hold on 
%                         nc = 2+tresh_lap_rew-template_laps;
%                         subplot(nc,1,1);imagesc(curveR.rate), colormap 'jet'
%                         title('Firing curvs original REWARD');
%                         subplot(nc,1,2);imagesc(curveT.rate), colormap 'jet'
%                         title('Firing curvs template');
%                     end 
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    stab_rew=[];
                    for c=1:size(in_lapR,1)-template_laps
                        lap = in_lapR(c,:); 
                        spks_lap = Restrict(spks_rew,lap);pos_lap = Restrict(pos_rew,lap);
                        [curveL , statsL] = FiringCurve(pos_lap , spks_lap , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                        %Spatial correlation
                        s = corrcoef(curveT.rate, curveL.rate);s = s(1,2);
                        stab_rew=[stab_rew,s];
%                         if mod(ii,2)% plot only odds
%                             subplot(nc,1,c+2);imagesc(curveL.rate), colormap 'jet'
%                         end    
                    end 
                end
                
                if ~exist('stab_ave', 'var')
                    stab_ave = []; 
                end   
                if ~exist('stab_rew', 'var')
                    stab_rew = []; 
                end    
                n.id = cluster; 
                n.pc =ismember(cluster, pc_dhpc); 
                n.stability.ave= stab_ave;
                n.stability.rew= stab_rew;
                dHPC_stab{ii}= n;
                    
            end
       end
        
        disp('VHPC Stability calculation')
         try
           load('vHPC_pc.mat') 
           pc_vhpc = [];
            for pc=1:size(vHPC,2)
                pc_vhpc = [pc_vhpc;vHPC{pc}.id];
            end
        
        catch
           dHPC = []; disp('No vhpc pc')
        end 
        %Store pc identity
        
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
                
               % AVERSIVE 
               [curveA ] = FiringCurve(pos_ave , spks_ave , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);

                % condicional si nlap_ave supera el tresh de laps ave
                if nlap_ave>= tresh_lap_ave
                    in_lapA = in_lapA(1:tresh_lap_ave,:); % keep the first tresh laps 
                    %Construct a template with the last # template_laps laps 
                    template = in_lapA(end-template_laps-1:end,:); spks_t = Restrict(spks_ave,template);pos_t = Restrict(pos_ave,template);  
                    [curveT , statsT] = FiringCurve(pos_t , spks_t , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    %%%%%%%%%%%%%%%%%%%%%control plot%%%%%%%%%%%%%%%%%%%%%
%                     if mod(ii,2)% plot only odds 
%                         figure;hold on 
%                         nc = 2+tresh_lap_ave-template_laps;
%                         subplot(nc,1,1);imagesc(curveA.rate), colormap 'jet'
%                         title('Firing curvs original AVERSIVE');
%                         subplot(nc,1,2);imagesc(curveT.rate), colormap 'jet'
%                         title('Firing curvs template');
%                     end 
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    stab_ave=[];
                    for c=1:size(in_lapA,1)-template_laps
                        lap = in_lapA(c,:); 
                        spks_lap = Restrict(spks_ave,lap);pos_lap = Restrict(pos_ave,lap);
                        [curveL , statsL] = FiringCurve(pos_lap , spks_lap , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                        %Spatial correlation
                        s = corrcoef(curveT.rate, curveL.rate);s = s(1,2);
                        stab_ave=[stab_ave,s];
%                         if mod(ii,2)% plot only odds
%                             subplot(nc,1,c+2);imagesc(curveL.rate), colormap 'jet'
%                         end    
                    end 
                end 
                
                % REWARD
                [curveR] = FiringCurve(pos_rew , spks_rew , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);

                % condicional si nlap_ave supera el tresh de laps ave
                if nlap_rew>= tresh_lap_rew
                    in_lapR = in_lapR(1:tresh_lap_rew,:); % keep the first tresh laps 
                    %Construct a template with the last # template_laps laps 
                    template = in_lapR(end-template_laps-1:end,:); spks_t = Restrict(spks_rew,template);pos_t = Restrict(pos_rew,template);  
                    [curveT , ] = FiringCurve(pos_t , spks_t , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    %%%%%%%%%%%%%%%%%%%%%control plot%%%%%%%%%%%%%%%%%%%%%
%                     if mod(ii,2)% plot only odds 
%                         figure;hold on 
%                         nc = 2+tresh_lap_rew-template_laps;
%                         subplot(nc,1,1);imagesc(curveR.rate), colormap 'jet'
%                         title('Firing curvs original REWARD');
%                         subplot(nc,1,2);imagesc(curveT.rate), colormap 'jet'
%                         title('Firing curvs template');
%                     end 
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    stab_rew=[];
                    for c=1:size(in_lapR,1)-template_laps
                        lap = in_lapR(c,:); 
                        spks_lap = Restrict(spks_rew,lap);pos_lap = Restrict(pos_rew,lap);
                        [curveL , statsL] = FiringCurve(pos_lap , spks_lap , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                        %Spatial correlation
                        s = corrcoef(curveT.rate, curveL.rate);s = s(1,2);
                        stab_rew=[stab_rew,s];
%                         if mod(ii,2)% plot only odds
%                             subplot(nc,1,c+2);imagesc(curveL.rate), colormap 'jet'
%                         end    
                    end 
                end
                
                if ~exist('stab_ave', 'var')
                    stab_ave = []; 
                end   
                if ~exist('stab_rew', 'var')
                    stab_rew = []; 
                end    
                n.id = cluster; 
                n.pc =ismember(cluster, pc_vhpc); 
                n.stability.ave= stab_ave;
                n.stability.rew= stab_rew;
                vHPC_stab{ii}= n;
                    
            end
       end
        
         %% Saveing PC INFO 
        if ~isempty(dHPC_stab)
            dHPC_stab = dHPC_stab(~cellfun('isempty',dHPC_stab));
            save([cd,'\dHPC_stab.mat'],'dHPC_stab'); 
        end
        if ~isempty(vHPC_stab)
            vHPC_stab = vHPC_stab(~cellfun('isempty',vHPC_stab));
            save([cd,'\vHPC_stab.mat'],'vHPC_stab'); 
        end
        
        end 
       
    end
    disp(['-------- Finished rat#' , num2str(tt) , ' --------'])
    disp(' ')
end

%% Bar plots
clear
clc
close all
% Parameters
path = {'\\Maryjackson\e\Rat127\Ephys\pyr';'\\Maryjackson\e\Rat128\Ephys\in_pyr\ready';'\\Maryjackson\e\Rat103\usable';'\\Maryjackson\e\Rat132\recordings\in_pyr'; '\\Maryjackson\e\Rat165\in_pyr'; ...
    '\\Maryjackson\e\Rat126\Ephys\in_Pyr'};%List of folders from the path

%Preparing the data 
%Output matrix 
dstab_ave = []; 
dstab_rew = [];

vstab_ave = []; 
vstab_rew = [];

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
            load('dHPC_stab.mat');
            %Save session parameters 
            for d=1:size(dHPC_stab,2)
                clu = dHPC_stab{d};
                if clu.pc ==1
                    dstab_ave = [dstab_ave;clu.stability.ave]; 
                    dstab_rew = [dstab_rew;clu.stability.rew];
                end
            end
        catch 
            dHPC = []; disp(['No dHPC_stab.mat file in ',session]);
        end 
        
        try
            load('vHPC_stab.mat');
            for d=1:size(vHPC_stab,2)
                clu = vHPC_stab{d}; 
                if clu.pc==1
                    vstab_ave = [vstab_ave;clu.stability.ave]; 
                    vstab_rew = [vstab_rew;clu.stability.rew];
                end
            end
     
        catch
            vHPC= []; disp(['No vHPC_stab.mat file in ',session]);
        end
        
    end 

end

%Ploting 
%%%%%%%%%%%Dorsal 
figure(1);clf; hold on;
subplot(1,2,1);hold on;
title('Aversive stability dHPC')
for i=1:size(dstab_ave,2)
     c = [.3,.3,.3];
    scatter(ones(size(dstab_ave,1),1)*i,dstab_ave(:,i),[],c,"filled",'jitter','on', 'jitterAmount',0.1,'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3);
end
bar((1:size(dstab_ave,2)),nanmean(dstab_ave),'r');
ylabel('Spatial correlation'); 
xlabel('Laps'); 

subplot(1,2,2);hold on;
title('Reward stability dHPC')
for i=1:size(dstab_rew,2)
     c = [.3,.3,.3];
    scatter(ones(size(dstab_rew,1),1)*i,dstab_rew(:,i),[],c,"filled",'jitter','on', 'jitterAmount',0.1,'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3);
end
bar((1:size(dstab_rew,2)),nanmean(dstab_rew),'b');
ylabel('Spatial correlation'); 
xlabel('Laps'); 
sgtitle('Only PC dHPC');

%%%%%%%%%% Ventral 
figure(2);clf; hold on;
subplot(1,2,1);hold on;
title('Aversive stability vHPC')
for i=1:size(dstab_ave,2)
     c = [.3,.3,.3];
    scatter(ones(size(vstab_ave,1),1)*i,vstab_ave(:,i),[],c,"filled",'jitter','on', 'jitterAmount',0.1,'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3);
end
bar((1:size(vstab_ave,2)),nanmean(vstab_ave),'r');
ylabel('Spatial correlation'); 
xlabel('Laps'); 
subplot(1,2,2);hold on;
title('Reward stability vHPC')
for i=1:size(vstab_rew,2)
     c = [.3,.3,.3];
    scatter(ones(size(vstab_rew,1),1)*i,vstab_rew(:,i),[],c,"filled",'jitter','on', 'jitterAmount',0.1,'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3);
end
bar((1:size(vstab_rew,2)),nanmean(vstab_rew),'b');
ylabel('Spatial correlation'); 
xlabel('Laps'); 
sgtitle('Only PC vHPC');

%% Main Loop: stability across time controling for nlaps in ave/rew -SLOPE
% Stability parameters 
tresh_lap_ave = 10;
tresh_lap_rew = 10; 
template_laps = 5; %# of laps to contruct the tempalte and compare stability. 

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
        
        %Check if the session has the sufficient # laps 
        if isfile('Spikesorting\dHPC_pc_skaggs_pos.mat')
            load('Spikesorting\dHPC_pc_skaggs_pos.mat'); 
            nlap_ave= dHPC{1}.nlap_ave; nlap_rew =dHPC{1}.nlap_rew;
        elseif isfile('Spikesorting\vHPC_pc_skaggs_pos.mat')
            load('Spikesorting\vHPC_pc_skaggs_pos.mat');
            nlap_ave= vHPC{1}.nlap_ave; nlap_rew =vHPC{1}.nlap_rew;
        else
            nlap_ave= 0; nlap_rew =0;
        end 
        
        if and(nlap_ave >= tresh_lap_ave,nlap_rew >= tresh_lap_rew)
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
            dHPC_stab_subsampled = {};% one cell per pc
            vHPC_stab_subsampled = {}; 
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
        %% Firing curve 
        disp('dHPC Stability calculation')
        try
           load('dHPC_pc_skaggs_pos.mat')
            pc_dhpc = [];
        for pc=1:size(dHPC,2)
            pc_dhpc = [pc_dhpc;dHPC{pc}.id];
        end
        catch
           dHPC = []; disp('No dhpc pc')
        end 

        for ii=1:size(pc_dhpc,1)
           
            cluster = pc_dhpc(ii);

                spks = spks_dHPC(spks_dHPC(:,1)==cluster,2); % select tspk from cluster

                % --- Aversive ---
                spks_ave = spks;
                %Restrict spk and position to laps and movement periods
                spks_ave = Restrict(spks_ave,in_lapA);spks_ave = Restrict(spks_ave , movement.aversive); % Restrict spikes to movement periods and laps

                % --- Reward ---
                spks_rew = spks; 
                %Restrict spk and position to laps and movement periods
                spks_rew = Restrict(spks_rew,in_lapR); spks_rew = Restrict(spks_rew, movement.reward); % Restrict to movement periods
               
%                [curveA ] = FiringCurve(pos_ave , spks_ave , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                
               %%%%%%%%Control plot %%%%%%%%%%
%                figure(1);clf; hold on;subplot(1,2,1); plot(pos_ave(:,2),pos_ave(:,1),'color', [.3 .3 .3]);hold on; 
% % %                Find interpolated position of each spike:
%                xs = interp1(pos_ave(:,1),pos_ave(:,2),spks_ave);scatter(xs,spks_ave,'*'); 
%                title('Aversive');
%                    
%                subplot(1,2,2); plot(pos_rew(:,2),pos_rew(:,1),'color', [.3 .3 .3]);hold on; 
% % %            Find interpolated position of each spike:
%                xs = interp1(pos_rew(:,1),pos_rew(:,2),spks_rew); scatter(xs,spks_rew,'*');
%                title('Reward');
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
        
               %STABILITY 
                if nlap_ave < nlap_rew % subsample rew laps 
                    disp('1')
                    %Aversive
                    %Construct a template with the last # template_laps laps 
                    template = in_lapA(end-template_laps+1:end,:); spks_t = Restrict(spks_ave,template);pos_t = Restrict(pos_ave,template);  
                    [curveT] = FiringCurve(pos_t , spks_t , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    %Compute the correlation between each lap and the template 
                    stab_ave=[];
                    for c=1:size(in_lapA,1)-template_laps
                        lap = in_lapA(c,:); 
                        spks_lap = Restrict(spks_ave,lap);pos_lap = Restrict(pos_ave,lap);
                        [curveL , statsL] = FiringCurve(pos_lap , spks_lap , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                        %Spatial correlation
                        s = corrcoef(curveT.rate, curveL.rate);s = s(1,2);
                        stab_ave=[stab_ave,s];
%                         if mod(ii,2)% plot only odds
%                             subplot(nc,1,c+2);imagesc(curveL.rate), colormap 'jet'
%                         end    
                    end
                    
                   
                    slope_ave = fitlm((1:size(stab_ave,2))',stab_ave');
%                     figure, plot(slope_ave)

                    %Reward
                    %Subsampling 
                    f = round(nlap_rew/nlap_ave);
                    in_lapRs = in_lapR(1:f:end,:);
                    %Checking and fixing subsampling  
                    if size(in_lapRs,1)> nlap_ave
                        %Discard the last ones
                        in_lapRs = in_lapRs(1:nlap_ave,:);
                    elseif size(in_lapRs,1)< nlap_ave
                        la = nlap_ave - size(in_lapRs,1); 
                        %Randomly select la laps from the discarded ones 
                        idx =1:size(in_lapR,1); selected =idx(1:f:end); 
                        non_selected = idx(~ismember(idx,selected));
                        mix_idx=randperm(size(non_selected,2)); mix_idx=sort(mix_idx(1:la));% mixed not selected idx
                        extra_laps = in_lapR(non_selected(mix_idx),:);% picking randomly l laps from the non selected
                        in_lapRs = sort([in_lapRs;extra_laps]); 
                    end 
                    %Construct a template with the last # template_laps laps 
                    template = in_lapRs(end-template_laps+1:end,:); spks_t = Restrict(spks_rew,template);pos_t = Restrict(pos_rew,template);  
                    [curveT] = FiringCurve(pos_t , spks_t , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    %Compute the correlation between each lap and the template 
                    stab_rew=[];
                    for c=1:size(in_lapRs,1)-template_laps
                        lap = in_lapRs(c,:); 
                        spks_lap = Restrict(spks_rew,lap);pos_lap = Restrict(pos_rew,lap);
                        [curveL] = FiringCurve(pos_lap , spks_lap , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                        %Spatial correlation
                        s = corrcoef(curveT.rate, curveL.rate);s = s(1,2);
                        stab_rew=[stab_rew,s];
%                        
                    end
                    
                   
                    slope_rew = fitlm((1:size(stab_rew,2))',stab_rew');
%                     figure, plot(slope_rew)
                    
                    
                elseif nlap_ave > nlap_rew % subsample ave laps
                    disp('2');
                    %Reward
                    %Construct a template with the last # template_laps laps 
                    template = in_lapR(end-template_laps+1:end,:); spks_t = Restrict(spks_rew,template);pos_t = Restrict(pos_rew,template);  
                    [curveT] = FiringCurve(pos_t , spks_t , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    %Compute the correlation between each lap and the template 
                    stab_rew=[];
                    for c=1:size(in_lapR,1)-template_laps
                        lap = in_lapR(c,:); 
                        spks_lap = Restrict(spks_rew,lap);pos_lap = Restrict(pos_rew,lap);
                        [curveL] = FiringCurve(pos_lap , spks_lap , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                        %Spatial correlation
                        s = corrcoef(curveT.rate, curveL.rate);s = s(1,2);
                        stab_rew=[stab_rew,s];
%                         if mod(ii,2)% plot only odds
%                             subplot(nc,1,c+2);imagesc(curveL.rate), colormap 'jet'
%                         end    
                    end
                   
                    slope_rew = fitlm((1:size(stab_rew,2))',stab_rew');
%                     figure, plot(slope_rew)
                    
                    %Aversive
                    %Subsampling 
                    f = round(nlap_ave/nlap_rew);
                    in_lapAs = in_lapA(1:f:end,:);
                    %Checking and fixing subsampling  
                    if size(in_lapAs,1)> nlap_rew
                        %Discard the last ones
                        in_lapA = in_lapAs(1:nlap_rew,:);
                    elseif size(in_lapAs,1)< nlap_rew
                        la = nlap_rew - size(in_lapAs,1); 
                        %Randomly select la laps from the discarded ones 
                        idx =1:size(in_lapA,1); selected =idx(1:f:end); 
                        non_selected = idx(~ismember(idx,selected));
                        mix_idx=randperm(size(non_selected,2)); mix_idx=sort(mix_idx(1:la));% mixed not selected idx
                        extra_laps = in_lapA(non_selected(mix_idx),:);% picking randomly l laps from the non selected
                        in_lapAs = sort([in_lapAs;extra_laps]); 
                    end 
                    %Construct a template with the last # template_laps laps 
                    template = in_lapAs(end-template_laps+1:end,:); spks_t = Restrict(spks_ave,template);pos_t = Restrict(pos_ave,template);  
                    [curveT] = FiringCurve(pos_t , spks_t , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    %Compute the correlation between each lap and the template 
                    
                    stab_ave=[];
                    for c=1:size(in_lapAs,1)-template_laps
                        lap = in_lapAs(c,:); 
                        spks_lap = Restrict(spks_ave,lap);pos_lap = Restrict(pos_ave,lap);
                        [curveL] = FiringCurve(pos_lap , spks_lap , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                        %Spatial correlation
                        s = corrcoef(curveT.rate, curveL.rate);s = s(1,2);
                        stab_ave=[stab_ave,s];
                    end
                    
                   
                    slope_ave = fitlm((1:size(stab_ave,2))',stab_ave');
%                     figure(), plot(slope_ave);

                elseif nlap_ave == nlap_rew % do not subsample
                    disp('3');
                    %Aversive
                    %Construct a template with the last # template_laps laps 
                    template = in_lapA(end-template_laps+1:end,:); spks_t = Restrict(spks_ave,template);pos_t = Restrict(pos_ave,template);  
                    [curveT] = FiringCurve(pos_t , spks_t , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    %Compute the correlation between each lap and the template 
                    stab_ave=[];
                    for c=1:size(in_lapA,1)-template_laps
                        lap = in_lapA(c,:); 
                        spks_lap = Restrict(spks_ave,lap);pos_lap = Restrict(pos_ave,lap);
                        [curveL , statsL] = FiringCurve(pos_lap , spks_lap , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                        %Spatial correlation
                        s = corrcoef(curveT.rate, curveL.rate);s = s(1,2);
                        stab_ave=[stab_ave,s];
%                         if mod(ii,2)% plot only odds
%                             subplot(nc,1,c+2);imagesc(curveL.rate), colormap 'jet'
%                         end    
                    end
                   
                    slope_ave = fitlm((1:size(stab_ave,2))',stab_ave');
                    
                    %Reward
                    %Construct a template with the last # template_laps laps 
                    template = in_lapR(end-template_laps+1:end,:); spks_t = Restrict(spks_rew,template);pos_t = Restrict(pos_rew,template);  
                    [curveT] = FiringCurve(pos_t , spks_t , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    %Compute the correlation between each lap and the template 
                    stab_rew=[];
                    for c=1:size(in_lapR,1)-template_laps
                        lap = in_lapR(c,:); 
                        spks_lap = Restrict(spks_rew,lap);pos_lap = Restrict(pos_rew,lap);
                        [curveL] = FiringCurve(pos_lap , spks_lap , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                        %Spatial correlation
                        s = corrcoef(curveT.rate, curveL.rate);s = s(1,2);
                        stab_rew=[stab_rew,s];
%                         if mod(ii,2)% plot only odds
%                             subplot(nc,1,c+2);imagesc(curveL.rate), colormap 'jet'
%                         end    
                    end
                   
                    slope_rew = fitlm((1:size(stab_rew,2))',stab_rew');
%                     figure, plot(slope)
                end
                
                if ~exist('stab_ave', 'var')
                    stab_ave = []; 
                end   
                if ~exist('stab_rew', 'var')
                    stab_rew = []; 
                end    
                
                %Save neuron info
                n.id = cluster; 
                n.nlap_ave=nlap_ave;
                n.nlap_rew=nlap_rew;
                n.stability.ave= stab_ave;
                n.stability.rew= stab_rew;
                n.slope.ave = slope_ave;
                n.slope.rew =slope_rew;
                dHPC_stab_subsampled{ii}= n;
        end
       
        %%%%%%%%%%%%%%%%%% VHPC%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('VHPC Stability calculation')
        try
           load('vHPC_pc_skaggs_pos.mat') 
           pc_vhpc = [];
            for pc=1:size(vHPC,2)
                pc_vhpc = [pc_vhpc;vHPC{pc}.id];
            end
        
        catch
           vHPC = []; disp('No vhpc pc')
        end 
        %Store pc identity
       
        for ii=1:size(pc_vhpc,1)

            cluster = pc_vhpc(ii);

                spks = spks_vHPC(spks_vHPC(:,1)==cluster,2); % select tspk from cluster

                % --- Aversive ---
                spks_ave = spks; 
 
                %Restrict spk and position to laps and movement periods
                spks_ave = Restrict(spks_ave,in_lapA);spks_ave = Restrict(spks_ave , movement.aversive); % Restrict spikes to movement periods and laps

                % --- Reward ---
                spks_rew = spks; 
               
                %Restrict spk and position to laps and movement periods
                spks_rew = Restrict(spks_rew,in_lapR); spks_rew = Restrict(spks_rew, movement.reward); % Restrict to movement periods
            
%                [curveA ] = FiringCurve(pos_ave , spks_ave , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                
               %%%%%%%%Control plot %%%%%%%%%%
%                figure(1);clf; hold on;subplot(1,2,1); plot(pos_ave(:,2),pos_ave(:,1),'color', [.3 .3 .3]);hold on; 
% % %                Find interpolated position of each spike:
%                xs = interp1(pos_ave(:,1),pos_ave(:,2),spks_ave);scatter(xs,spks_ave,'*'); 
%                title('Aversive');
%                    
%                subplot(1,2,2); plot(pos_rew(:,2),pos_rew(:,1),'color', [.3 .3 .3]);hold on; 
% % %            Find interpolated position of each spike:
%                xs = interp1(pos_rew(:,1),pos_rew(:,2),spks_rew); scatter(xs,spks_rew,'*');
%                title('Reward');
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            

               %STABILITY 
                if nlap_ave < nlap_rew % subsample rew laps 
                    %Aversive
                    %Construct a template with the last # template_laps laps 
                    template = in_lapA(end-template_laps+1:end,:); spks_t = Restrict(spks_ave,template);pos_t = Restrict(pos_ave,template);  
                    [curveT] = FiringCurve(pos_t , spks_t , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    %Compute the correlation between each lap and the template 
                    stab_ave=[];
                    for c=1:size(in_lapA,1)-template_laps
                    lap = in_lapA(c,:); 
                    spks_lap = Restrict(spks_ave,lap);pos_lap = Restrict(pos_ave,lap);
                    [curveL , statsL] = FiringCurve(pos_lap , spks_lap , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    %Spatial correlation
                    s = corrcoef(curveT.rate, curveL.rate);s = s(1,2);
                    stab_ave=[stab_ave,s];
%                         if mod(ii,2)% plot only odds
%                             subplot(nc,1,c+2);imagesc(curveL.rate), colormap 'jet'
%                         end    
                    end
                    slope_ave = fitlm((1:size(stab_ave,2))',stab_ave');
%                     figure, plot(slope_ave)

                    %Reward
                    %Subsampling 
                    f = round(nlap_rew/nlap_ave);
                    in_lapRs = in_lapR(1:f:end,:);
                    %Checking and fixing subsampling  
                    if size(in_lapRs,1)> nlap_ave
                        %Discard the last ones
                        in_lapRs = in_lapRs(1:nlap_ave,:);
                    elseif size(in_lapRs,1)< nlap_ave
                        la = nlap_ave - size(in_lapRs,1); 
                        %Randomly select la laps from the discarded ones 
                        idx =1:size(in_lapR,1); selected =idx(1:f:end); 
                        non_selected = idx(~ismember(idx,selected));
                        mix_idx=randperm(size(non_selected,2)); mix_idx=sort(mix_idx(1:la));% mixed not selected idx
                        extra_laps = in_lapR(non_selected(mix_idx),:);% picking randomly l laps from the non selected
                        in_lapRs = sort([in_lapRs;extra_laps]); 
                    end 
                    %Construct a template with the last # template_laps laps 
                    template = in_lapRs(end-template_laps+1:end,:); spks_t = Restrict(spks_rew,template);pos_t = Restrict(pos_rew,template);  
                    [curveT] = FiringCurve(pos_t , spks_t , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    %Compute the correlation between each lap and the template 
                    stab_rew=[];
                    for c=1:size(in_lapRs,1)-template_laps
                        lap = in_lapRs(c,:); 
                        spks_lap = Restrict(spks_rew,lap);pos_lap = Restrict(pos_rew,lap);
                        [curveL] = FiringCurve(pos_lap , spks_lap , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                        %Spatial correlation
                        s = corrcoef(curveT.rate, curveL.rate);s = s(1,2);
                        stab_rew=[stab_rew,s];
%                        
                    end
                    slope_rew = fitlm((1:size(stab_rew,2))',stab_rew');
%                     figure, plot(slope_rew)
                    
                    
                elseif nlap_ave > nlap_rew % subsample ave laps
                    %Reward
                    %Construct a template with the last # template_laps laps 
                    template = in_lapR(end-template_laps+1:end,:); spks_t = Restrict(spks_rew,template);pos_t = Restrict(pos_rew,template);  
                    [curveT] = FiringCurve(pos_t , spks_t , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    %Compute the correlation between each lap and the template 
                    stab_rew=[];
                    for c=1:size(in_lapR,1)-template_laps
                    lap = in_lapR(c,:); 
                    spks_lap = Restrict(spks_rew,lap);pos_lap = Restrict(pos_rew,lap);
                    [curveL] = FiringCurve(pos_lap , spks_lap , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    %Spatial correlation
                    s = corrcoef(curveT.rate, curveL.rate);s = s(1,2);
                    stab_rew=[stab_rew,s];
%                         if mod(ii,2)% plot only odds
%                             subplot(nc,1,c+2);imagesc(curveL.rate), colormap 'jet'
%                         end    
                    end
                    slope_rew = fitlm((1:size(stab_rew,2))',stab_rew');
%                     figure, plot(slope_rew)
                    
                    %Aversive
                    %Subsampling 
                    f = round(nlap_ave/nlap_rew);
                    in_lapAs = in_lapA(1:f:end,:);
                    %Checking and fixing subsampling  
                    if size(in_lapAs,1)> nlap_rew
                        %Discard the last ones
                        in_lapA = in_lapAs(1:nlap_rew,:);
                    elseif size(in_lapAs,1)< nlap_rew
                        la = nlap_rew - size(in_lapAs,1); 
                        %Randomly select la laps from the discarded ones 
                        idx =1:size(in_lapA,1); selected =idx(1:f:end); 
                        non_selected = idx(~ismember(idx,selected));
                        mix_idx=randperm(size(non_selected,2)); mix_idx=sort(mix_idx(1:la));% mixed not selected idx
                        extra_laps = in_lapA(non_selected(mix_idx),:);% picking randomly l laps from the non selected
                        in_lapAs = sort([in_lapAs;extra_laps]); 
                    end 
                    %Construct a template with the last # template_laps laps 
                    template = in_lapAs(end-template_laps+1:end,:); spks_t = Restrict(spks_ave,template);pos_t = Restrict(pos_ave,template);  
                    [curveT] = FiringCurve(pos_t , spks_t , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    %Compute the correlation between each lap and the template 
                    stab_ave=[];
                    for c=1:size(in_lapAs,1)-template_laps
                        lap = in_lapAs(c,:); 
                        spks_lap = Restrict(spks_ave,lap);pos_lap = Restrict(pos_ave,lap);
                        [curveL] = FiringCurve(pos_lap , spks_lap , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                        %Spatial correlation
                        s = corrcoef(curveT.rate, curveL.rate);s = s(1,2);
                        stab_ave=[stab_ave,s];
                    end
                    slope_ave = fitlm((1:size(stab_ave,2))',stab_ave');
%                     figure(), plot(slope_ave);

                elseif nlap_ave == nlap_rew % do not subsample
                    
                    %Aversive
                    %Construct a template with the last # template_laps laps 
                    template = in_lapA(end-template_laps+1:end,:); spks_t = Restrict(spks_ave,template);pos_t = Restrict(pos_ave,template);  
                    [curveT] = FiringCurve(pos_t , spks_t , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    %Compute the correlation between each lap and the template 
                    stab_ave=[];
                    for c=1:size(in_lapA,1)-template_laps
                    lap = in_lapA(c,:); 
                    spks_lap = Restrict(spks_ave,lap);pos_lap = Restrict(pos_ave,lap);
                    [curveL , statsL] = FiringCurve(pos_lap , spks_lap , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    %Spatial correlation
                    s = corrcoef(curveT.rate, curveL.rate);s = s(1,2);
                    stab_ave=[stab_ave,s];
%                         if mod(ii,2)% plot only odds
%                             subplot(nc,1,c+2);imagesc(curveL.rate), colormap 'jet'
%                         end    
                    end
                    slope_ave = fitlm((1:size(stab_ave,2))',stab_ave');
                    
                    %Reward
                    %Construct a template with the last # template_laps laps 
                    template = in_lapR(end-template_laps+1:end,:); spks_t = Restrict(spks_rew,template);pos_t = Restrict(pos_rew,template);  
                    [curveT] = FiringCurve(pos_t , spks_t , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    %Compute the correlation between each lap and the template 
                    stab_rew=[];
                    for c=1:size(in_lapR,1)-template_laps
                    lap = in_lapR(c,:); 
                    spks_lap = Restrict(spks_rew,lap);pos_lap = Restrict(pos_rew,lap);
                    [curveL] = FiringCurve(pos_lap , spks_lap , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    %Spatial correlation
                    s = corrcoef(curveT.rate, curveL.rate);s = s(1,2);
                    stab_rew=[stab_rew,s];
%                         if mod(ii,2)% plot only odds
%                             subplot(nc,1,c+2);imagesc(curveL.rate), colormap 'jet'
%                         end    
                    end
                    slope_rew = fitlm((1:size(stab_rew,2))',stab_rew');
%                     figure, plot(slope)
                end
                 
                if ~exist('stab_ave', 'var')
                    stab_ave = []; 
                end   
                if ~exist('stab_rew', 'var')
                    stab_rew = []; 
                end    
                
                %Save neuron info
                n.id = cluster; 
                n.pc =ismember(cluster, pc_vhpc);
                n.nlap_ave=nlap_ave;
                n.nlap_rew=nlap_rew;
                n.stability.ave= stab_ave;
                n.stability.rew= stab_rew;
                n.slope.ave = slope_ave;
                n.slope.rew =slope_rew;
                vHPC_stab_subsampled{ii}= n;
       
        end
     
        
         %% Saveing PC INFO 
        if ~isempty(dHPC_stab_subsampled)
            dHPC_stab_subsampled = dHPC_stab_subsampled(~cellfun('isempty',dHPC_stab_subsampled));
            save([cd,'\dHPC_stab_subsampled.mat'],'dHPC_stab_subsampled'); 
        end
        if ~isempty(vHPC_stab_subsampled)
            vHPC_stab_subsampled = vHPC_stab_subsampled(~cellfun('isempty',vHPC_stab_subsampled));
            save([cd,'\vHPC_stab_subsampled.mat'],'vHPC_stab_subsampled'); 
        end
        
        
        end 
    end
    disp(['-------- Finished rat#' , num2str(tt) , ' --------'])
    disp(' ')
end


%% Main Loop: stability across time controling for nlaps in ave/rew -COM DISTRIBUTION
% Stability parameters 
tresh_lap_ave = 15;
tresh_lap_rew = 15; 
template_laps = 5; %# of laps to contruct the tempalte and compare stability. 

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
        
        %Check if the session has the sufficient # laps 
        if isfile('Spikesorting\dHPC_pc_skaggs_pos.mat')
            load('Spikesorting\dHPC_pc_skaggs_pos.mat'); 
            nlap_ave= dHPC{1}.nlap_ave; nlap_rew =dHPC{1}.nlap_rew;
        elseif isfile('Spikesorting\vHPC_pc_skaggs_pos.mat')
            load('Spikesorting\vHPC_pc_skaggs_pos.mat');
            nlap_ave= vHPC{1}.nlap_ave; nlap_rew =vHPC{1}.nlap_rew;
        else
            nlap_ave= 0; nlap_rew =0;
        end 
        
        if and(nlap_ave >= tresh_lap_ave,nlap_rew >= tresh_lap_rew)
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
            dHPC_stab_subsampled = {};% one cell per pc
            vHPC_stab_subsampled = {}; 
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
        %% Firing curve 
        disp('dHPC Stability calculation')
        try
           load('dHPC_pc_skaggs_pos.mat')
            pc_dhpc = [];
        for pc=1:size(dHPC,2)
            pc_dhpc = [pc_dhpc;dHPC{pc}.id];
        end
        catch
           dHPC = []; disp('No dhpc pc')
        end 

        for ii=1:size(pc_dhpc,1)
            
            cluster = pc_dhpc(ii);

                spks = spks_dHPC(spks_dHPC(:,1)==cluster,2); % select tspk from cluster

                % --- Aversive ---
                spks_ave = spks;
                %Restrict spk and position to laps and movement periods
                spks_ave = Restrict(spks_ave,in_lapA);spks_ave = Restrict(spks_ave , movement.aversive); % Restrict spikes to movement periods and laps

                % --- Reward ---
                spks_rew = spks; 
                %Restrict spk and position to laps and movement periods
                spks_rew = Restrict(spks_rew,in_lapR); spks_rew = Restrict(spks_rew, movement.reward); % Restrict to movement periods
               
%                [curveA ] = FiringCurve(pos_ave , spks_ave , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                
               %%%%%%%%Control plot %%%%%%%%%%
%                figure(1);clf; hold on;subplot(1,2,1); plot(pos_ave(:,2),pos_ave(:,1),'color', [.3 .3 .3]);hold on; 
% % %                Find interpolated position of each spike:
%                xs = interp1(pos_ave(:,1),pos_ave(:,2),spks_ave);scatter(xs,spks_ave,'*'); 
%                title('Aversive');
%                    
%                subplot(1,2,2); plot(pos_rew(:,2),pos_rew(:,1),'color', [.3 .3 .3]);hold on; 
% % %            Find interpolated position of each spike:
%                xs = interp1(pos_rew(:,1),pos_rew(:,2),spks_rew); scatter(xs,spks_rew,'*');
%                title('Reward');
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
              
               %STABILITY 
                if nlap_ave < nlap_rew % subsample rew laps 
                    %Aversive
                    %Construct a template with the last # template_laps laps 
                    template = in_lapA(end-template_laps+1:end,:); spks_t = Restrict(spks_ave,template);pos_t = Restrict(pos_ave,template);  
                    [curveT] = FiringCurve(pos_t , spks_t , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    %Compute the correlation between each lap and the template 
                    stab_ave=[];
                    for c=1:size(in_lapA,1)-template_laps
                        lap = in_lapA(c,:); 
                        spks_lap = Restrict(spks_ave,lap);pos_lap = Restrict(pos_ave,lap);
                        [curveL , statsL] = FiringCurve(pos_lap , spks_lap , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                        %Spatial correlation
                        s = corrcoef(curveT.rate, curveL.rate);s = s(1,2);
                        stab_ave=[stab_ave,s];
%                         if mod(ii,2)% plot only odds
%                             subplot(nc,1,c+2);imagesc(curveL.rate), colormap 'jet'
%                         end    
                    end
                    
                    stab_ave=stab_ave+ rescale; % rescaling to avoid 'rank deficient to within machine precision'
                    slope_ave = fitlm((1:size(stab_ave,2))',stab_ave'+100);
%                     figure, plot(slope_ave)

                    %Reward
                    %Subsampling 
                    f = round(nlap_rew/nlap_ave);
                    in_lapRs = in_lapR(1:f:end,:);
                    %Checking and fixing subsampling  
                    if size(in_lapRs,1)> nlap_ave
                        %Discard the last ones
                        in_lapRs = in_lapRs(1:nlap_ave,:);
                    elseif size(in_lapRs,1)< nlap_ave
                        la = nlap_ave - size(in_lapRs,1); 
                        %Randomly select la laps from the discarded ones 
                        idx =1:size(in_lapR,1); selected =idx(1:f:end); 
                        non_selected = idx(~ismember(idx,selected));
                        mix_idx=randperm(size(non_selected,2)); mix_idx=sort(mix_idx(1:la));% mixed not selected idx
                        extra_laps = in_lapR(non_selected(mix_idx),:);% picking randomly l laps from the non selected
                        in_lapRs = sort([in_lapRs;extra_laps]); 
                    end 
                    %Construct a template with the last # template_laps laps 
                    template = in_lapRs(end-template_laps+1:end,:); spks_t = Restrict(spks_rew,template);pos_t = Restrict(pos_rew,template);  
                    [curveT] = FiringCurve(pos_t , spks_t , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    %Compute the correlation between each lap and the template 
                    stab_rew=[];
                    for c=1:size(in_lapRs,1)-template_laps
                        lap = in_lapRs(c,:); 
                        spks_lap = Restrict(spks_rew,lap);pos_lap = Restrict(pos_rew,lap);
                        [curveL] = FiringCurve(pos_lap , spks_lap , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                        %Spatial correlation
                        s = corrcoef(curveT.rate, curveL.rate);s = s(1,2);
                        stab_rew=[stab_rew,s];
%                        
                    end
                    
                    stab_rew=stab_rew+ rescale;% rescaling to avoid 'rank deficient to within machine precision'
                    slope_rew = fitlm((1:size(stab_rew,2))',stab_rew');
%                     figure, plot(slope_rew)
                    
                    
                elseif nlap_ave > nlap_rew % subsample ave laps
                    %Reward
                    %Construct a template with the last # template_laps laps 
                    template = in_lapR(end-template_laps+1:end,:); spks_t = Restrict(spks_rew,template);pos_t = Restrict(pos_rew,template);  
                    [curveT] = FiringCurve(pos_t , spks_t , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    %Compute the correlation between each lap and the template 
                    stab_rew=[];
                    for c=1:size(in_lapR,1)-template_laps
                        lap = in_lapR(c,:); 
                        spks_lap = Restrict(spks_rew,lap);pos_lap = Restrict(pos_rew,lap);
                        [curveL] = FiringCurve(pos_lap , spks_lap , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                        %Spatial correlation
                        s = corrcoef(curveT.rate, curveL.rate);s = s(1,2);
                        stab_rew=[stab_rew,s];
%                         if mod(ii,2)% plot only odds
%                             subplot(nc,1,c+2);imagesc(curveL.rate), colormap 'jet'
%                         end    
                    end
                    
                    stab_rew=stab_rew+ rescale;% rescaling to avoid 'rank deficient to within machine precision'
                    slope_rew = fitlm((1:size(stab_rew,2))',stab_rew');
%                     figure, plot(slope_rew)
                    
                    %Aversive
                    %Subsampling 
                    f = round(nlap_ave/nlap_rew);
                    in_lapAs = in_lapA(1:f:end,:);
                    %Checking and fixing subsampling  
                    if size(in_lapAs,1)> nlap_rew
                        %Discard the last ones
                        in_lapA = in_lapAs(1:nlap_rew,:);
                    elseif size(in_lapAs,1)< nlap_rew
                        la = nlap_rew - size(in_lapAs,1); 
                        %Randomly select la laps from the discarded ones 
                        idx =1:size(in_lapA,1); selected =idx(1:f:end); 
                        non_selected = idx(~ismember(idx,selected));
                        mix_idx=randperm(size(non_selected,2)); mix_idx=sort(mix_idx(1:la));% mixed not selected idx
                        extra_laps = in_lapA(non_selected(mix_idx),:);% picking randomly l laps from the non selected
                        in_lapAs = sort([in_lapAs;extra_laps]); 
                    end 
                    %Construct a template with the last # template_laps laps 
                    template = in_lapAs(end-template_laps+1:end,:); spks_t = Restrict(spks_ave,template);pos_t = Restrict(pos_ave,template);  
                    [curveT] = FiringCurve(pos_t , spks_t , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    %Compute the correlation between each lap and the template 
                    
                    stab_ave=[];
                    for c=1:size(in_lapAs,1)-template_laps
                        lap = in_lapAs(c,:); 
                        spks_lap = Restrict(spks_ave,lap);pos_lap = Restrict(pos_ave,lap);
                        [curveL] = FiringCurve(pos_lap , spks_lap , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                        %Spatial correlation
                        s = corrcoef(curveT.rate, curveL.rate);s = s(1,2);
                        stab_ave=[stab_ave,s];
                    end
                    
                    stab_ave=stab_ave+ rescale;% rescaling to avoid 'rank deficient to within machine precision'
                    slope_ave = fitlm((1:size(stab_ave,2))',stab_ave');
%                     figure(), plot(slope_ave);

                elseif nlap_ave == nlap_rew % do not subsample
                    
                    %Aversive
                    %Construct a template with the last # template_laps laps 
                    template = in_lapA(end-template_laps+1:end,:); spks_t = Restrict(spks_ave,template);pos_t = Restrict(pos_ave,template);  
                    [curveT] = FiringCurve(pos_t , spks_t , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    %Compute the correlation between each lap and the template 
                    stab_ave=[];
                    for c=1:size(in_lapA,1)-template_laps
                        lap = in_lapA(c,:); 
                        spks_lap = Restrict(spks_ave,lap);pos_lap = Restrict(pos_ave,lap);
                        [curveL , statsL] = FiringCurve(pos_lap , spks_lap , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                        %Spatial correlation
                        s = corrcoef(curveT.rate, curveL.rate);s = s(1,2);
                        stab_ave=[stab_ave,s];
%                         if mod(ii,2)% plot only odds
%                             subplot(nc,1,c+2);imagesc(curveL.rate), colormap 'jet'
%                         end    
                    end
                    stab_ave=stab_ave+ rescale;% rescaling to avoid 'rank deficient to within machine precision'
                    slope_ave = fitlm((1:size(stab_ave,2))',stab_ave');
                    
                    %Reward
                    %Construct a template with the last # template_laps laps 
                    template = in_lapR(end-template_laps+1:end,:); spks_t = Restrict(spks_rew,template);pos_t = Restrict(pos_rew,template);  
                    [curveT] = FiringCurve(pos_t , spks_t , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    %Compute the correlation between each lap and the template 
                    stab_rew=[];
                    for c=1:size(in_lapR,1)-template_laps
                        lap = in_lapR(c,:); 
                        spks_lap = Restrict(spks_rew,lap);pos_lap = Restrict(pos_rew,lap);
                        [curveL] = FiringCurve(pos_lap , spks_lap , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                        %Spatial correlation
                        s = corrcoef(curveT.rate, curveL.rate);s = s(1,2);
                        stab_rew=[stab_rew,s];
%                         if mod(ii,2)% plot only odds
%                             subplot(nc,1,c+2);imagesc(curveL.rate), colormap 'jet'
%                         end    
                    end
                    stab_rew=stab_rew+ rescale;% rescaling to avoid 'rank deficient to within machine precision'
                    slope_rew = fitlm((1:size(stab_rew,2))',stab_rew');
%                     figure, plot(slope)
                end
                
                if ~exist('stab_ave', 'var')
                    stab_ave = []; 
                end   
                if ~exist('stab_rew', 'var')
                    stab_rew = []; 
                end    
                
                %Save neuron info
                n.id = cluster; 
                n.nlap_ave=nlap_ave;
                n.nlap_rew=nlap_rew;
                n.stability.ave= stab_ave;
                n.stability.rew= stab_rew;
                n.slope.ave = slope_ave;
                n.slope.rew =slope_rew;
                dHPC_stab_subsampled{ii}= n;
        end
       
        %%%%%%%%%%%%%%%%%% VHPC%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('VHPC Stability calculation')
         try
           load('vHPC_pc_skaggs_pos.mat') 
            pc_vhpc = [];
            for pc=1:size(vHPC,2)
                pc_vhpc = [pc_vhpc;vHPC{pc}.id];
            end
        
        catch
           vHPC = []; disp('No vhpc pc')
        end 
        %Store pc identity
       
        for ii=1:size(pc_vhpc,1)

            cluster = pc_vhpc(ii);

                spks = spks_vHPC(spks_vHPC(:,1)==cluster,2); % select tspk from cluster

                % --- Aversive ---
                spks_ave = spks; 
 
                %Restrict spk and position to laps and movement periods
                spks_ave = Restrict(spks_ave,in_lapA);spks_ave = Restrict(spks_ave , movement.aversive); % Restrict spikes to movement periods and laps

                % --- Reward ---
                spks_rew = spks; 
               
                %Restrict spk and position to laps and movement periods
                spks_rew = Restrict(spks_rew,in_lapR); spks_rew = Restrict(spks_rew, movement.reward); % Restrict to movement periods
            
%                [curveA ] = FiringCurve(pos_ave , spks_ave , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                
               %%%%%%%%Control plot %%%%%%%%%%
%                figure(1);clf; hold on;subplot(1,2,1); plot(pos_ave(:,2),pos_ave(:,1),'color', [.3 .3 .3]);hold on; 
% % %                Find interpolated position of each spike:
%                xs = interp1(pos_ave(:,1),pos_ave(:,2),spks_ave);scatter(xs,spks_ave,'*'); 
%                title('Aversive');
%                    
%                subplot(1,2,2); plot(pos_rew(:,2),pos_rew(:,1),'color', [.3 .3 .3]);hold on; 
% % %            Find interpolated position of each spike:
%                xs = interp1(pos_rew(:,1),pos_rew(:,2),spks_rew); scatter(xs,spks_rew,'*');
%                title('Reward');
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            

               %STABILITY 
                if nlap_ave < nlap_rew % subsample rew laps 
                    %Aversive
                    %Construct a template with the last # template_laps laps 
                    template = in_lapA(end-template_laps+1:end,:); spks_t = Restrict(spks_ave,template);pos_t = Restrict(pos_ave,template);  
                    [curveT] = FiringCurve(pos_t , spks_t , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    %Compute the correlation between each lap and the template 
                    stab_ave=[];
                    for c=1:size(in_lapA,1)-template_laps
                    lap = in_lapA(c,:); 
                    spks_lap = Restrict(spks_ave,lap);pos_lap = Restrict(pos_ave,lap);
                    [curveL , statsL] = FiringCurve(pos_lap , spks_lap , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    %Spatial correlation
                    s = corrcoef(curveT.rate, curveL.rate);s = s(1,2);
                    stab_ave=[stab_ave,s];
%                         if mod(ii,2)% plot only odds
%                             subplot(nc,1,c+2);imagesc(curveL.rate), colormap 'jet'
%                         end    
                    end
                    slope_ave = fitlm((1:size(stab_ave,2))',stab_ave');
%                     figure, plot(slope_ave)

                    %Reward
                    %Subsampling 
                    f = round(nlap_rew/nlap_ave);
                    in_lapRs = in_lapR(1:f:end,:);
                    %Checking and fixing subsampling  
                    if size(in_lapRs,1)> nlap_ave
                        %Discard the last ones
                        in_lapRs = in_lapRs(1:nlap_ave,:);
                    elseif size(in_lapRs,1)< nlap_ave
                        la = nlap_ave - size(in_lapRs,1); 
                        %Randomly select la laps from the discarded ones 
                        idx =1:size(in_lapR,1); selected =idx(1:f:end); 
                        non_selected = idx(~ismember(idx,selected));
                        mix_idx=randperm(size(non_selected,2)); mix_idx=sort(mix_idx(1:la));% mixed not selected idx
                        extra_laps = in_lapR(non_selected(mix_idx),:);% picking randomly l laps from the non selected
                        in_lapRs = sort([in_lapRs;extra_laps]); 
                    end 
                    %Construct a template with the last # template_laps laps 
                    template = in_lapRs(end-template_laps+1:end,:); spks_t = Restrict(spks_rew,template);pos_t = Restrict(pos_rew,template);  
                    [curveT] = FiringCurve(pos_t , spks_t , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    %Compute the correlation between each lap and the template 
                    stab_rew=[];
                    for c=1:size(in_lapRs,1)-template_laps
                        lap = in_lapRs(c,:); 
                        spks_lap = Restrict(spks_rew,lap);pos_lap = Restrict(pos_rew,lap);
                        [curveL] = FiringCurve(pos_lap , spks_lap , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                        %Spatial correlation
                        s = corrcoef(curveT.rate, curveL.rate);s = s(1,2);
                        stab_rew=[stab_rew,s];
%                        
                    end
                    slope_rew = fitlm((1:size(stab_rew,2))',stab_rew');
%                     figure, plot(slope_rew)
                    
                    
                elseif nlap_ave > nlap_rew % subsample ave laps
                    %Reward
                    %Construct a template with the last # template_laps laps 
                    template = in_lapR(end-template_laps+1:end,:); spks_t = Restrict(spks_rew,template);pos_t = Restrict(pos_rew,template);  
                    [curveT] = FiringCurve(pos_t , spks_t , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    %Compute the correlation between each lap and the template 
                    stab_rew=[];
                    for c=1:size(in_lapR,1)-template_laps
                    lap = in_lapR(c,:); 
                    spks_lap = Restrict(spks_rew,lap);pos_lap = Restrict(pos_rew,lap);
                    [curveL] = FiringCurve(pos_lap , spks_lap , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    %Spatial correlation
                    s = corrcoef(curveT.rate, curveL.rate);s = s(1,2);
                    stab_rew=[stab_rew,s];
%                         if mod(ii,2)% plot only odds
%                             subplot(nc,1,c+2);imagesc(curveL.rate), colormap 'jet'
%                         end    
                    end
                    slope_rew = fitlm((1:size(stab_rew,2))',stab_rew');
%                     figure, plot(slope_rew)
                    
                    %Aversive
                    %Subsampling 
                    f = round(nlap_ave/nlap_rew);
                    in_lapAs = in_lapA(1:f:end,:);
                    %Checking and fixing subsampling  
                    if size(in_lapAs,1)> nlap_rew
                        %Discard the last ones
                        in_lapA = in_lapAs(1:nlap_rew,:);
                    elseif size(in_lapAs,1)< nlap_rew
                        la = nlap_rew - size(in_lapAs,1); 
                        %Randomly select la laps from the discarded ones 
                        idx =1:size(in_lapA,1); selected =idx(1:f:end); 
                        non_selected = idx(~ismember(idx,selected));
                        mix_idx=randperm(size(non_selected,2)); mix_idx=sort(mix_idx(1:la));% mixed not selected idx
                        extra_laps = in_lapA(non_selected(mix_idx),:);% picking randomly l laps from the non selected
                        in_lapAs = sort([in_lapAs;extra_laps]); 
                    end 
                    %Construct a template with the last # template_laps laps 
                    template = in_lapAs(end-template_laps+1:end,:); spks_t = Restrict(spks_ave,template);pos_t = Restrict(pos_ave,template);  
                    [curveT] = FiringCurve(pos_t , spks_t , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    %Compute the correlation between each lap and the template 
                    stab_ave=[];
                    for c=1:size(in_lapAs,1)-template_laps
                        lap = in_lapAs(c,:); 
                        spks_lap = Restrict(spks_ave,lap);pos_lap = Restrict(pos_ave,lap);
                        [curveL] = FiringCurve(pos_lap , spks_lap , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                        %Spatial correlation
                        s = corrcoef(curveT.rate, curveL.rate);s = s(1,2);
                        stab_ave=[stab_ave,s];
                    end
                    slope_ave = fitlm((1:size(stab_ave,2))',stab_ave');
%                     figure(), plot(slope_ave);

                elseif nlap_ave == nlap_rew % do not subsample
                    
                    %Aversive
                    %Construct a template with the last # template_laps laps 
                    template = in_lapA(end-template_laps+1:end,:); spks_t = Restrict(spks_ave,template);pos_t = Restrict(pos_ave,template);  
                    [curveT] = FiringCurve(pos_t , spks_t , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    %Compute the correlation between each lap and the template 
                    stab_ave=[];
                    for c=1:size(in_lapA,1)-template_laps
                    lap = in_lapA(c,:); 
                    spks_lap = Restrict(spks_ave,lap);pos_lap = Restrict(pos_ave,lap);
                    [curveL , statsL] = FiringCurve(pos_lap , spks_lap , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    %Spatial correlation
                    s = corrcoef(curveT.rate, curveL.rate);s = s(1,2);
                    stab_ave=[stab_ave,s];
%                         if mod(ii,2)% plot only odds
%                             subplot(nc,1,c+2);imagesc(curveL.rate), colormap 'jet'
%                         end    
                    end
                    slope_ave = fitlm((1:size(stab_ave,2))',stab_ave');
                    
                    %Reward
                    %Construct a template with the last # template_laps laps 
                    template = in_lapR(end-template_laps+1:end,:); spks_t = Restrict(spks_rew,template);pos_t = Restrict(pos_rew,template);  
                    [curveT] = FiringCurve(pos_t , spks_t , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    %Compute the correlation between each lap and the template 
                    stab_rew=[];
                    for c=1:size(in_lapR,1)-template_laps
                    lap = in_lapR(c,:); 
                    spks_lap = Restrict(spks_rew,lap);pos_lap = Restrict(pos_rew,lap);
                    [curveL] = FiringCurve(pos_lap , spks_lap , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    %Spatial correlation
                    s = corrcoef(curveT.rate, curveL.rate);s = s(1,2);
                    stab_rew=[stab_rew,s];
%                         if mod(ii,2)% plot only odds
%                             subplot(nc,1,c+2);imagesc(curveL.rate), colormap 'jet'
%                         end    
                    end
                    slope_rew = fitlm((1:size(stab_rew,2))',stab_rew');
%                     figure, plot(slope)
                end
                 
                if ~exist('stab_ave', 'var')
                    stab_ave = []; 
                end   
                if ~exist('stab_rew', 'var')
                    stab_rew = []; 
                end    
                
                %Save neuron info
                n.id = cluster; 
                n.pc =ismember(cluster, pc_vhpc);
                n.nlap_ave=nlap_ave;
                n.nlap_rew=nlap_rew;
                n.stability.ave= stab_ave;
                n.stability.rew= stab_rew;
                n.slope.ave = slope_ave;
                n.slope.rew =slope_rew;
                vHPC_stab_subsampled{ii}= n;
       
        end
     
        
         %% Saveing PC INFO 
        if ~isempty(dHPC_stab_subsampled)
            dHPC_stab_subsampled = dHPC_stab_subsampled(~cellfun('isempty',dHPC_stab_subsampled));
            save([cd,'\dHPC_stab_subsampled.mat'],'dHPC_stab_subsampled'); 
        end
        if ~isempty(vHPC_stab_subsampled)
            vHPC_stab_subsampled = vHPC_stab_subsampled(~cellfun('isempty',vHPC_stab_subsampled));
            save([cd,'\vHPC_stab_subsampled.mat'],'vHPC_stab_subsampled'); 
        end
        
        
        end 
    end
    disp(['-------- Finished rat#' , num2str(tt) , ' --------'])
    disp(' ')
end
%% Bar plots
clear
clc
close all
% Parameters
path = {'\\Maryjackson\e\Rat127\Ephys\pyr';'\\Maryjackson\e\Rat128\Ephys\in_pyr\ready';'\\Maryjackson\e\Rat103\usable';'\\Maryjackson\e\Rat132\recordings\in_pyr'; '\\Maryjackson\e\Rat165\in_pyr'; ...
    '\\Maryjackson\e\Rat126\Ephys\in_Pyr'};%List of folders from the path


%%%%%%% Plot spatial corr template vs lap %%%%%%%%%
%Output matrix 
 
dstab_= [];
vstab= [];

for tt = 1:length(path)
    %List of folders from the path
    files = dir(path{tt});
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
    clear files dirFlags

    for t = 1 : length(subFolders)-2
        disp([newline 'Initiating analysis of folder #' , num2str(t) , ' from rat #',num2str(tt) , ' --'])
        session = [subFolders(t+2).folder,'\',subFolders(t+2).name];
        cd([session,'\Spikesorting'])
        %Loading pc matrix of the sessions
        try
            load('dHPC_stab_subsampled.mat');
            %Save session parameters 
            for d=1:size(dHPC_stab_subsampled,2)
                clu = dHPC_stab_subsampled{d};
                dstab_ave= [dstab_ave;clu.stability.ave'];
                dstab_rew= [dstab_rew;clu.stability.rew'];

            end
        catch 
            dHPC = []; disp(['No dHPC_stab_subsampled.mat file in ',session]);
            cdhpc = cdhpc +1;
        end 
        
        try
            load('vHPC_stab_subsampled.mat');
            for d=1:size(vHPC_stab_subsampled,2)
                clu = vHPC_stab_subsampled{d}; 
                vstab_ave= [vstab_ave;clu.stability.ave'];
                vstab_rew= [vstab_rew;clu.stability.rew'];
               
            end
     
        catch
            vHPC= []; disp(['No vHPC_stab_subsampled.mat file in ',session]);
            cvhpc = cvhpc +1;
        end
        
    end 

end

%Ploting correlations 
[p,h] = ranksum(dstab_ave,dstab_rew)
figure(1);clf; hold on
boxplot([dstab_ave,dstab_rew]); 
title(['dHPC - temp vs lap p= ',num2str(p)])
ylabel('Spatial corr'); xlabel('Condition')
xticks([1 2])
xticklabels({'Ave', 'Rew'});

[p,h] = ranksum(vstab_ave,vstab_rew)
figure(2);clf; hold on
boxplot([vstab_ave,vstab_rew]); 
title(['vHPC - temp vs lap p= ',num2str(p)])
ylabel('Spatial corr'); xlabel('Condition')
xticks([1 2])
xticklabels({'Ave', 'Rew'});

%%%%%%%%%%%% Plot mean corr %%%%%%

%Plot corr
%Output matrix 
 
dstab_ave= [];dstab_rew= [];
vstab_ave= [];vstab_rew= [];

for tt = 1:length(path)
    %List of folders from the path
    files = dir(path{tt});
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
    clear files dirFlags

    for t = 1 : length(subFolders)-2
        disp([newline 'Initiating analysis of folder #' , num2str(t) , ' from rat #',num2str(tt) , ' --'])
        session = [subFolders(t+2).folder,'\',subFolders(t+2).name];
        cd([session,'\Spikesorting'])
        %Loading pc matrix of the sessions
        try
            load('dHPC_stab_subsampled.mat');
            %Save session parameters 
            for d=1:size(dHPC_stab_subsampled,2)
                clu = dHPC_stab_subsampled{d};
                dstab_ave= [dstab_ave;nanmean(clu.stability.ave)];
                dstab_rew= [dstab_rew;nanmean(clu.stability.rew)];

            end
        catch 
            dHPC = []; disp(['No dHPC_stab_subsampled.mat file in ',session]);
            cdhpc = cdhpc +1;
        end 
        
        try
            load('vHPC_stab_subsampled.mat');
            for d=1:size(vHPC_stab_subsampled,2)
                clu = vHPC_stab_subsampled{d}; 
                vstab_ave= [vstab_ave;nanmean(clu.stability.ave)];
                vstab_rew= [vstab_rew;nanmean(clu.stability.rew)];
               
            end
     
        catch
            vHPC= []; disp(['No vHPC_stab_subsampled.mat file in ',session]);
            cvhpc = cvhpc +1;
        end
        
    end 

end

%Ploting correlations 
[p,h] = ranksum(dstab_ave,dstab_rew)
figure(3);clf; hold on
boxplot([dstab_ave,dstab_rew]); 
title(['dHPC - mean temp vs lap p= ',num2str(p)])
ylabel('Spatial corr'); xlabel('Condition')
xticks([1 2])
xticklabels({'Ave', 'Rew'});

[p,h] = ranksum(vstab_ave,vstab_rew)
figure(4);clf; hold on
boxplot([vstab_ave,vstab_rew]); 
title(['vHPC - mean temp vs lap p= ',num2str(p)])
ylabel('Spatial corr'); xlabel('Condition')
xticks([1 2])
xticklabels({'Ave', 'Rew'});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cdhpc = 0; 
cvhpc = 0; 
%Preparing the data 
%Output matrix 
dstab_ave = []; 
dstab_rew = [];

vstab_ave = []; 
vstab_rew = [];

for tt = 1:length(path)
    %List of folders from the path
    files = dir(path{tt});
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
    clear files dirFlags

    for t = 1 : length(subFolders)-2
        disp([newline 'Initiating analysis of folder #' , num2str(t) , ' from rat #',num2str(tt) , ' --'])
        session = [subFolders(t+2).folder,'\',subFolders(t+2).name];
        cd([session,'\Spikesorting'])
        %Loading pc matrix of the sessions
        try
            load('dHPC_stab_subsampled.mat');
            %Save session parameters 
            for d=1:size(dHPC_stab_subsampled,2)
                clu = dHPC_stab_subsampled{d};
               
                    dstab_ave = [dstab_ave;clu.stability.ave']; 
                    dstab_rew = [dstab_rew;clu.stability.rew'];
               
            end
        catch 
            dHPC = []; disp(['No dHPC_stab_subsampled.mat file in ',session]);
            cdhpc = cdhpc +1;
        end 
        
        try
            load('vHPC_stab_subsampled.mat');
            for d=1:size(vHPC_stab_subsampled,2)
                clu = vHPC_stab_subsampled{d}; 
                    vstab_ave = [vstab_ave;clu.stability.ave']; 
                    vstab_rew = [vstab_rew;clu.stability.rew'];
               
            end
     
        catch
            vHPC= []; disp(['No vHPC_stab_subsampled.mat file in ',session]);
            cvhpc = cvhpc +1;
        end
        
    end 

end



%Ploting correlations 
dstab_ave
figure(1);hold on 
x = dstab_ave;
y = dstab_rew;
c = [.3,.3,.3];
scatter(x,y,30,c,"filled",'jitter','on', 'jitterAmount',0.3); xlim([0 2]);
ylabel('Spatial corr lap-template');
xticks([1 2])
xticklabels({'Ave', 'Rew'});


%Ploting 
%%%%%%%%%%%Dorsal 
figure(1);clf; hold on;
subplot(1,2,1);hold on;
title('Aversive stability dHPC')
for i=1:size(dstab_ave,2)
     c = [.3,.3,.3];
    scatter(ones(size(dstab_ave,1),1)*i,dstab_ave(:,i),[],c,"filled",'jitter','on', 'jitterAmount',0.1,'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3);
end
% bar((1:size(dstab_ave,2)),nanmean(dstab_ave),'r');
ylabel('Spatial correlation'); 
xlabel('Laps'); 

subplot(1,2,2);hold on;
title('Reward stability dHPC')
for i=1:size(dstab_rew,2)
     c = [.3,.3,.3];
    scatter(ones(size(dstab_rew,1),1)*i,dstab_rew(:,i),[],c,"filled",'jitter','on', 'jitterAmount',0.1,'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3);
end
% bar((1:size(dstab_rew,2)),nanmean(dstab_rew),'b');
ylabel('Spatial correlation'); 
xlabel('Laps'); 
sgtitle('Only PC dHPC - Subsampled');

%%%%%%%%%% Ventral 
figure(2);clf; hold on;
subplot(1,2,1);hold on;
title('Aversive stability vHPC')
for i=1:size(dstab_ave,2)
     c = [.3,.3,.3];
    scatter(ones(size(vstab_ave,1),1)*i,vstab_ave(:,i),[],c,"filled",'jitter','on', 'jitterAmount',0.1,'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3);
end
% bar((1:size(vstab_ave,2)),nanmean(vstab_ave),'r');
ylabel('Spatial correlation'); 
xlabel('Laps'); 
subplot(1,2,2);hold on;
title('Reward stability vHPC')
for i=1:size(vstab_rew,2)
     c = [.3,.3,.3];
    scatter(ones(size(vstab_rew,1),1)*i,vstab_rew(:,i),[],c,"filled",'jitter','on', 'jitterAmount',0.1,'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3);
end
% bar((1:size(vstab_rew,2)),nanmean(vstab_rew),'b');
ylabel('Spatial correlation'); 
xlabel('Laps'); 
sgtitle('Only PC vHPC - subsampled');


%% Distribution plots
clear
clc
close all
% Parameters
path = {'\\Maryjackson\e\Rat127\Ephys\pyr';'\\Maryjackson\e\Rat128\Ephys\in_pyr\ready';'\\Maryjackson\e\Rat103\usable';'\\Maryjackson\e\Rat132\recordings\in_pyr'; '\\Maryjackson\e\Rat165\in_pyr'; ...
    '\\Maryjackson\e\Rat126\Ephys\in_Pyr'};%List of folders from the path


cdhpc = 0; 
cvhpc = 0; 
%Preparing the data 
%Output matrix 
dstab_ave = []; 
dstab_rew = [];

vstab_ave = []; 
vstab_rew = [];

for tt = 1:length(path)
    %List of folders from the path
    files = dir(path{tt});
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
    clear files dirFlags

    for t = 1 : length(subFolders)-2
        disp([newline 'Initiating analysis of folder #' , num2str(t) , ' from rat #',num2str(tt) , ' --'])
        session = [subFolders(t+2).folder,'\',subFolders(t+2).name];
        cd([session,'\Spikesorting'])
        %Loading pc matrix of the sessions
        try
            load('dHPC_stab_subsampled.mat');
            %Save session parameters 
            for d=1:size(dHPC_stab_subsampled,2)
                clu = dHPC_stab_subsampled{d};
                if clu.pc ==1
                    dstab_ave = [dstab_ave;clu.slope.ave.Coefficients.Estimate(2)]; 
                    dstab_rew = [dstab_rew;clu.slope.rew.Coefficients.Estimate(2)];
                end
            end
        catch 
            dHPC = []; disp(['No dHPC_stab_subsampled.mat file in ',session]);
            cdhpc = cdhpc +1;
        end 
        
        try
            load('vHPC_stab_subsampled.mat');
            for d=1:size(vHPC_stab_subsampled,2)
                clu = vHPC_stab_subsampled{d}; 
                if clu.pc==1
                    vstab_ave = [vstab_ave;clu.slope.ave.Coefficients.Estimate(2)]; 
                    vstab_rew = [vstab_rew;clu.slope.rew.Coefficients.Estimate(2)];
                end
            end
     
        catch
            vHPC= []; disp(['No vHPC_stab_subsampled.mat file in ',session]);
            cvhpc = cvhpc +1;
        end
        
    end 

end

figure(1);clf; hold on
histogram(dstab_ave,50,'EdgeColor',[0 0 0],'FaceColor','r'); hold on
histogram(dstab_rew,50,'EdgeColor',[0 0 0],'FaceColor','b'); hold on
title('dHPC - slope stability')
ylabel('Counts'); xlabel('Slope')

figure(3);clf; hold on
boxplot([dstab_ave,dstab_rew]); 
title('dHPC - slope stability')
ylabel('R slope'); xlabel('Condition')
[p,h] = ranksum(dstab_ave,dstab_rew)

figure(2);clf; hold on
histogram(vstab_ave,50,'EdgeColor',[0 0 0],'FaceColor','r'); hold on
histogram(vstab_rew,50,'EdgeColor',[0 0 0],'FaceColor','b'); hold on
title('vHPC - slope stability')
ylabel('Counts'); xlabel('Slope')

[H p]= kstest2(vstab_ave,vstab_rew)


figure(4);clf; hold on
boxplot([vstab_ave,vstab_rew]); 
title('vHPC - slope stability')
ylabel('R slope'); xlabel('Condition')
[p,h] = ranksum(vstab_ave,vstab_rew)