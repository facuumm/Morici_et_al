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
Xedges = 60; %number of bins for RateMap construction - 3cm bin
sigma = 2;%round(15/(180/Xedges)); %defined for gauss kernel of 15cm
binSize = 0.001; % bin size for replay events detection
bin_size = 1; % to bin pos ans spks in between/within  
% Behavior
minimal_speed = 2.5;% minimal speed to detect quite periods
minimal_speed_time = 2; % minimal time to detect quite periods
%% SHOCK - PF DISTANCE  1st op
%Main loop

dhpc_ave = []; 
dhpc_rew = []; 

vhpc_ave = []; 
vhpc_rew = []; 

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
        
        %% Find interpolated positon of the shock 
        % Aversive 
        pos_ave = behavior.pos.aversive(:,1:2);
        pos_ave(:,2) =pos_ave(:,2)-min(pos_ave(:,2)); pos_ave(:,2) = pos_ave(:,2)/max(pos_ave(:,2)); %normalization of position
        in_maze = ToIntervals(pos_ave(:,1),and(pos_ave(:,2)>0.03 , pos_ave(:,2)<0.180));% eliminating the extrems of the maze
       
        x_shock = interp1(pos_ave(:,1),pos_ave(:,2),Shocks_filt);
        binx_shock = x_shock*60; % shock in bin scale 
        
%         figure(2);clf;hold on; 
%         plot(pos_ave(:,2),pos_ave(:,1),'color', [.3 .3 .3]);hold on; 
%         scatter(x_shock(1),Shocks_filt(1),'*','r'); figure(2)
        
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
        
        %%  Shock - center PF distance 
        disp('Distance calculation')
        cd([session,'\Spikesorting'])
        
        try
            load('dHPC_pc.mat');
            
            for d=1:size(dHPC,2)
                pc = dHPC{d}; 
                % Aversive
                pf_lim = pc.stats_ave.fieldX; 
                pc_frmap = pc.frMap_ave; 
                field = pc_frmap(pf_lim(1,1):pf_lim(1,2));
                
                %Center of mass of the field 
                idx = 1:size(field,2);                             
                com = sum(idx .* field) / sum(field);
                com = com + pf_lim(1,1); % back in general scale  
                
                %%%Control plot%%%
%                 figure(3);clf;hold on;
%                 imagesc(pc_frmap); colormap 'jet'
%                 axis tight 
%                 xline(pf_lim(1),'w');xline(pf_lim(2),'w') 
%                 scatter(find(pc_frmap==pc.stats_ave.peak),1,'*','w')
%                 scatter(com,1,'*','r')
                %%%%%%%%%%%%%%%%%%%
                
                %Distance shock-com
                dis = binx_shock-com; 
                dhpc_ave = [dhpc_ave;dis]; 
                
                % Reward
                pf_lim = pc.stats_rew.fieldX; 
                pc_frmap = pc.frMap_rew; 
                field = pc_frmap(pf_lim(1,1):pf_lim(1,2));
                
                %Center of mass of the field 
                idx = 1:size(field,2);                             
                com = sum(idx .* field) / sum(field);
                com = com + pf_lim(1,1); % back in general scale  
                
                %%%Control plot%%%
%                 figure(3);clf;hold on;
%                 imagesc(pc_frmap); colormap 'jet'
%                 axis tight 
%                 for p=1:size(pf_lim,1)
%                     xline(pf_lim(p,1),'w');xline(pf_lim(p,2),'w'); 
%                     scatter(find(pc_frmap==pc.stats_rew.peak(p)),1,'*','w')
%                 end
%                 scatter(com,1,'*','g')
                %%%%%%%%%%%%%%%%%%%
                
                %Distance shock-com
                dis = binx_shock-com; 
                dhpc_rew = [dhpc_rew;dis]; 
                
                
  
            end

        catch 
            dHPC = []; disp(['No dHPC_pc.mat file in ',session]);
        end
        
        try
            load('vHPC_pc.mat');
            
            for d=1:size(vHPC,2)
                pc = vHPC{d}; 
                % Aversive
                pf_lim = pc.stats_ave.fieldX; 
                pc_frmap = pc.frMap_ave; 
                field = pc_frmap(pf_lim(1,1):pf_lim(1,2));
                
                %Center of mass of the field 
                idx = 1:size(field,2);                             
                com = sum(idx .* field) / sum(field);
                com = com + pf_lim(1,1); % back in general scale  
                
                %%%Control plot%%%
%                 figure(3);clf;hold on;
%                 imagesc(pc_frmap); colormap 'jet'
%                 axis tight 
%                 xline(pf_lim(1),'w');xline(pf_lim(2),'w') 
%                 scatter(find(pc_frmap==pc.stats_ave.peak),1,'*','w')
%                 scatter(com,1,'*','r')
                %%%%%%%%%%%%%%%%%%%
                
                %Distance shock-com
                dis = binx_shock-com; 
                vhpc_ave = [dhpc_ave;dis]; 
                
                % Reward
                pf_lim = pc.stats_rew.fieldX; 
                pc_frmap = pc.frMap_rew; 
                field = pc_frmap(pf_lim(1,1):pf_lim(1,2));
                
                %Center of mass of the field 
                idx = 1:size(field,2);                             
                com = sum(idx .* field) / sum(field);
                com = com + pf_lim(1,1); % back in general scale  
                
                %%%Control plot%%%
%                 figure(3);clf;hold on;
%                 imagesc(pc_frmap); colormap 'jet'
%                 axis tight 
%                 for p=1:size(pf_lim,1)
%                     xline(pf_lim(p,1),'w');xline(pf_lim(p,2),'w'); 
%                     scatter(find(pc_frmap==pc.stats_rew.peak(p)),1,'*','w')
%                 end
%                 scatter(com,1,'*','g')
                %%%%%%%%%%%%%%%%%%%
                
                %Distance shock-com
                dis = binx_shock-com; 
                vhpc_rew = [dhpc_rew;dis]; 

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

 
 
figure(5); clf; hold on
h = histogram(dhpc_ave, 'Normalization', 'probability', 'EdgeColor', 'none', 'FaceColor',[1 0 0.1]);
histogram(dhpc_rew, 'Normalization', 'probability', 'EdgeColor', 'none', 'FaceColor',[0.1 0 1]);
xline(0); 
colororder([1 0 0;0 0 1]);
legend(["Aversive" "Reward"]);
title('dHPC')
 
figure(6); clf; hold on
histogram(vhpc_ave, 'Normalization', 'probability', 'EdgeColor', 'none', 'FaceColor',[1 0 0.1]);
histogram(vhpc_rew, 'Normalization', 'probability', 'EdgeColor', 'none', 'FaceColor',[0.1 0 1]);
xline(0);
colororder([1 0 0;0 0 1]);
legend(["Aversive" "Reward"]);
title('vHPC')

% pd = fitdist(dhpc_ave,'Kernel');
% x = min(dhpc_ave):0.01:max(dhpc_ave);
% p = cdf(pd,x);
% plot(x,p)
% hold on
 
%% SHOCK - PF DISTANCE  2nd op
%Output initialization
dhpc_den_ave=[]; 
dhpc_den_rew=[];

vhpc_den_ave=[]; 
vhpc_den_rew=[]; 

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
        
        %% Find interpolated positon of the shock 
        % Aversive 
        pos_ave = behavior.pos.aversive(:,1:2);
        pos_ave(:,2) =pos_ave(:,2)-min(pos_ave(:,2)); pos_ave(:,2) = pos_ave(:,2)/max(pos_ave(:,2)); %normalization of position
        in_maze = ToIntervals(pos_ave(:,1),and(pos_ave(:,2)>0.03 , pos_ave(:,2)<0.180));% eliminating the extrems of the maze
       
        x_shock = interp1(pos_ave(:,1),pos_ave(:,2),Shocks_filt);
        binx_shock = x_shock*60; % shock in bin scale 
        
%         figure(2);clf;hold on; 
%         plot(pos_ave(:,2),pos_ave(:,1),'color', [.3 .3 .3]);hold on; 
%         scatter(x_shock(1),Shocks_filt(1),'*','r'); figure(2)
        
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
        
        %%  Shock - center PF distance 
        disp('Distance calculation')
        cd([session,'\Spikesorting'])
        
        try
            load('dHPC_pc.mat');
            % Calculate session distances 
            dhpc_ave = [];
            dhpc_rew = [];
            for d=1:size(dHPC,2)
                pc = dHPC{d}; 
                % Aversive
                pf_lim = pc.stats_ave.fieldX; 
                pc_frmap = pc.frMap_ave; 
                field = pc_frmap(pf_lim(1,1):pf_lim(1,2));
                
                %Center of mass of the field 
                idx = 1:size(field,2);                             
                com = sum(idx .* field) / sum(field);
                com = com + pf_lim(1,1); % back in general scale  
                
                %%%Control plot%%%
%                 figure(3);clf;hold on;
%                 imagesc(pc_frmap); colormap 'jet'
%                 axis tight 
%                 xline(pf_lim(1),'w');xline(pf_lim(2),'w') 
%                 scatter(find(pc_frmap==pc.stats_ave.peak),1,'*','w')
%                 scatter(com,1,'*','r')
                %%%%%%%%%%%%%%%%%%%
                
                %Distance shock-com
                dis = binx_shock-com; 
                dhpc_ave = [dhpc_ave;dis]; 
                
                % Reward
                pf_lim = pc.stats_rew.fieldX; 
                pc_frmap = pc.frMap_rew; 
                field = pc_frmap(pf_lim(1,1):pf_lim(1,2));
                
                %Center of mass of the field 
                idx = 1:size(field,2);                             
                com = sum(idx .* field) / sum(field);
                com = com + pf_lim(1,1); % back in general scale  
                
                %%%Control plot%%%
%                 figure(3);clf;hold on;
%                 imagesc(pc_frmap); colormap 'jet'
%                 axis tight 
%                 for p=1:size(pf_lim,1)
%                     xline(pf_lim(p,1),'w');xline(pf_lim(p,2),'w'); 
%                     scatter(find(pc_frmap==pc.stats_rew.peak(p)),1,'*','w')
%                 end
%                 scatter(com,1,'*','g')
                %%%%%%%%%%%%%%%%%%%
                
                %Distance shock-com
                dis = binx_shock-com; 
                dhpc_rew = [dhpc_rew;dis]; 
                
                
  
            end
            
            % Calculate density 
            [N] = histcounts(dhpc_ave,120);
            density_ave = N/size(dhpc_ave,1);
%             figure;hold on; imagesc((Smooth(density, [1 1]))'); colormap 'gray';axis tight
%             xline(60,'r','LineWidth',1)
            
            [N] = histcounts(dhpc_rew,120);
            density_rew = N/size(dhpc_rew,1);
%             figure;hold on; imagesc((Smooth(density, [1 1]))'); colormap 'gray';axis tight
%             xline(60,'r','LineWidth',1)
            
            dhpc_den_ave=[dhpc_den_ave; density_ave]; 
            dhpc_den_rew=[dhpc_den_rew;density_rew];
            
        catch 
            dHPC = []; disp(['No dHPC_pc.mat file in ',session]);
        end
        
        try
            load('vHPC_pc.mat');
            % Calculate session distances 
            vhpc_ave = [];
            vhpc_rew = [];
            for d=1:size(vHPC,2)
                pc = vHPC{d}; 
                % Aversive
                pf_lim = pc.stats_ave.fieldX; 
                pc_frmap = pc.frMap_ave; 
                field = pc_frmap(pf_lim(1,1):pf_lim(1,2));
                
                %Center of mass of the field 
                idx = 1:size(field,2);                             
                com = sum(idx .* field) / sum(field);
                com = com + pf_lim(1,1); % back in general scale  
                
                %%%Control plot%%%
%                 figure(3);clf;hold on;
%                 imagesc(pc_frmap); colormap 'jet'
%                 axis tight 
%                 xline(pf_lim(1),'w');xline(pf_lim(2),'w') 
%                 scatter(find(pc_frmap==pc.stats_ave.peak),1,'*','w')
%                 scatter(com,1,'*','r')
                %%%%%%%%%%%%%%%%%%%
                
                %Distance shock-com
                dis = binx_shock-com; 
                vhpc_ave = [dhpc_ave;dis]; 
                
                % Reward
                pf_lim = pc.stats_rew.fieldX; 
                pc_frmap = pc.frMap_rew; 
                field = pc_frmap(pf_lim(1,1):pf_lim(1,2));
                
                %Center of mass of the field 
                idx = 1:size(field,2);                             
                com = sum(idx .* field) / sum(field);
                com = com + pf_lim(1,1); % back in general scale  
                
                %%%Control plot%%%
%                 figure(3);clf;hold on;
%                 imagesc(pc_frmap); colormap 'jet'
%                 axis tight 
%                 for p=1:size(pf_lim,1)
%                     xline(pf_lim(p,1),'w');xline(pf_lim(p,2),'w'); 
%                     scatter(find(pc_frmap==pc.stats_rew.peak(p)),1,'*','w')
%                 end
%                 scatter(com,1,'*','g')
                %%%%%%%%%%%%%%%%%%%
                
                %Distance shock-com
                dis = binx_shock-com; 
                vhpc_rew = [dhpc_rew;dis]; 

            end
            
             % Calculate density 
            [N] = histcounts(vhpc_ave,120);
            density_ave = N/size(vhpc_ave,1);
%             figure;hold on; imagesc((Smooth(density, [1 1]))'); colormap 'gray';axis tight
%             xline(60,'r','LineWidth',1)
            
            [N] = histcounts(vhpc_rew,120);
            density_rew = N/size(vhpc_rew,1);
%             figure;hold on; imagesc((Smooth(density, [1 1]))'); colormap 'gray';axis tight
%             xline(60,'r','LineWidth',1)
            
            vhpc_den_ave=[vhpc_den_ave;density_ave]; 
            vhpc_den_rew=[vhpc_den_rew;density_rew];
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



%Density heat maps plots per session
% Normalization of density plots
%%%dHPC%%%%
tempA=[];
tempR=[];
for i =1:size(dhpc_den_ave,1)
    A = dhpc_den_ave(i,:) - min(dhpc_den_ave(i,:));
    A = A ./ max(A);
    tempA=[tempA ; Smooth(A, [1 1])'];
    clear A
    
    R = dhpc_den_rew(i,:) - min(dhpc_den_rew(i,:));
    R = R ./ max(R);
    tempR=[tempR ;Smooth(R, [1 1])'];
    clear R    
end

[h idx] = max (tempA, [],2);
[m mm] = sort(idx); 

figure(1);clf;hold on; 
subplot(1,2,1);imagesc([0:1:120], [1:1:size(dhpc_den_ave,1)],tempA(mm,:)), colormap 'gray'; 
title('Aversive');xline(60,'r');xlabel(['Shock-COM distance',char(10),'Spatial bins']); 
ylabel('Sessions'); 

[h idx] = max (tempR, [],2);
[m mm] = sort(idx); 
subplot(1,2,2);imagesc([0:1:120], [1:1:size(dhpc_den_rew,1)],tempR(mm,:)), colormap 'gray'; 
title('Reward');xline(60,'r');xlabel(['Shock-COM distance',char(10),'Spatial bins']); 
ylabel('Sessions'); 
sgtitle('dHPC');

%window for stats
w_size = 2; % # of bins from the center. Each bin 3cm.
win_ave = nanmean(tempA(:, 60-w_size:60+w_size), 2); 
win_rew = nanmean(tempR(:, 60-w_size:60+w_size), 2); 
p = signrank(win_ave,win_rew) 

%%%vHPC%%%%
tempA=[];
tempR=[];
for i =1:size(vhpc_den_ave,1)
    A = vhpc_den_ave(i,:) - min(vhpc_den_ave(i,:));
    A = A ./ max(A);
    tempA=[tempA ; Smooth(A, [1 1])'];
    clear A
    
    R = vhpc_den_rew(i,:) - min(vhpc_den_rew(i,:));
    R = R ./ max(R);
    tempR=[tempR ;Smooth(R, [1 1])'];
    clear R    
end

[h idx] = max (tempA, [],2);
[m mm] = sort(idx); 

figure(2);clf;hold on; 
subplot(1,2,1);imagesc([0:1:120], [1:1:size(vhpc_den_ave,1)],tempA(mm,:)), colormap 'gray'; 
title('Aversive');xlabel(['Shock-COM distance',char(10),'Spatial bins']); 
ylabel('Sessions'); 
xline(60,'r',{'Shock'}); xline(56,'g',{'W4'}); xline(64,'g',{'W4'});
xline(50,'b',{'W10'}); xline(70,'b',{'W10'});

[h idx] = max (tempR, [],2);
[m mm] = sort(idx); 
subplot(1,2,2);imagesc([0:1:120], [1:1:size(vhpc_den_rew,1)],tempR(mm,:)), colormap 'gray'; 
title('Reward');xline(60,'r');xlabel(['Shock-COM distance',char(10),'Spatial bins']); 
xline(60,'r',{'Shock'}); xline(56,'g',{'W4'}); xline(64,'g',{'W4'});
xline(50,'b',{'W10'}); xline(70,'b',{'W10'});
ylabel('Sessions'); 
sgtitle('vHPC');


% Stats 
%window for stats
w_size = 3; % # of bins from the center. Each bin 3cm.
win_ave = nanmean(tempA(:, 60-w_size:60+w_size), 2); 
win_rew = nanmean(tempR(:, 60-w_size:60+w_size), 2); 
p = signrank(win_ave,win_rew)


%% SHOCK RESPONSIVE CELLS 

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
        % Definition of the emotional transition
        if aversiveTS_run(1)<rewardTS_run(1)
            cond = 1; % 1 if the order was Aversive -> Reward
        else
            cond = 2;% 1 if the order was Reward -> Aversive
        end
        
        %% Responsivness calculation
        %CCG parameters
         bin = 0.1;smooth=[1];th=2;window=4;
        
         dHPC_shock_norm={};
         vHPC_shock_norm={};
         
        %%%%%%% PYR 
%         group_dHPC = group_dHPC(group_dHPC(:,4)==1,:);
%         group_vHPC = group_vHPC(group_vHPC(:,4)==1,:);
%         
%         if ~isempty(group_dHPC)
%             %Aversive
%             curve=[]; 
%             limits=[behavior.speed.aversive(1,1),behavior.speed.aversive(end,1)];     
%             [curve , bins , responsive] = SU_responsivness(spks_dHPC,group_dHPC(:,1),Shocks_filt,...
%                                       limits,[0 1],window,bin,smooth,'gain',th); 
%              
%             %Save output
%             dhpc_shock=[dhpc_shock;curve']; 
%         end
%         
%         if ~isempty(group_vHPC)
%             disp('vHPC Responsivness')
%             %Aversive
%             curve=[]; 
%             limits=[behavior.speed.aversive(1,1),behavior.speed.aversive(end,1)];           
%             [curve , bins , responsive] = SU_responsivness(spks_vHPC,group_vHPC(:,1),Shocks_filt,...
%                                       limits,[0 1],4,bin,smooth,'gain',th);                                  
%             %Save output
%             vhpc_shock=[vhpc_shock;curve'];
%         end 
%         
        
%         %%%%%%% PC
        try
            load('dHPC_pc.mat');
            dhpc_pc = [];
            for d=1:size(dHPC,2)
                pc = dHPC{d}; 
                dhpc_pc = [dhpc_pc;pc.id]; 
            end
            
            group_dHPC = group_dHPC(ismember(group_dHPC(:,1),dhpc_pc),:); % keep only pc 
            disp('dHPC Responsivness')
            %Aversive
            curve=[]; 
            limits=[behavior.speed.aversive(1,1),behavior.speed.aversive(end,1)];     
            [curve , bins , responsive] = SU_responsivness(spks_dHPC,group_dHPC(:,1),Shocks_filt,...
                                      limits,[0 1],window,bin,smooth,'gain',th); 
             
            %Save output
            dHPC_shock_norm.id = group_dHPC(:,1); 
            dHPC_shock_norm.curve=curve'; 
            dHPC_shock_norm.resp=responsive'; 
            save([cd,'\dHPC_shock_norm.mat'],'dHPC_shock_norm'); 

        catch 
            dHPC = []; disp(['No dHPC_pc.mat file in ',session]);
        end 
        
        try
            load('vHPC_pc.mat');
            vhpc_pc = [];
            for d=1:size(vHPC,2)
                pc = vHPC{d}; 
                vhpc_pc = [vhpc_pc;pc.id]; 
            end
            group_vHPC = group_vHPC(ismember(group_vHPC(:,1),vhpc_pc),:);
            
            disp('vHPC Responsivness')
            %Aversive
            curve=[]; 
            limits=[behavior.speed.aversive(1,1),behavior.speed.aversive(end,1)];           
            [curve , bins , responsive] = SU_responsivness(spks_vHPC,group_vHPC(:,1),Shocks_filt,...
                                      limits,[0 1],4,bin,smooth,'gain',th);                                  
           
            %Save output
            vHPC_shock_norm.id = group_vHPC(:,1); 
            vHPC_shock_norm.curve=curve'; 
            vHPC_shock_norm.resp=responsive'; 
            save([cd,'\vHPC_shock_norm.mat'],'vHPC_shock_norm');
        
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


%%%%%%% Shock resp prop /pc  %%%%%%%
%Output matrix 
dHPC_pc = 0;dHPC_shock = 0;
vHPC_pc = 0;vHPC_shock = 0;
dHPC_id = [];vHPC_id = [];

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
            load('dHPC_shock_norm.mat');
            %Save session parameters 
            dHPC_pc = dHPC_pc +  size(dHPC_shock_norm.resp,1); 
            dHPC_shock = dHPC_shock +  sum(dHPC_shock_norm.resp);
            temp = [dHPC_shock_norm.id,dHPC_shock_norm.resp]; id_temp= temp(temp(:,2)==1,1); 
            dHPC_id = [dHPC_id;id_temp,ones(size(id_temp,1),1)*t];
        catch 
            dHPC = []; disp(['No dHPC_pc_lap.mat file in ',session]);
        end 
        
        try
            load('vHPC_shock_norm.mat');
            vHPC_pc = vHPC_pc +  size(vHPC_shock_norm.resp,1); 
            vHPC_shock = vHPC_shock +  sum(vHPC_shock_norm.resp);
            temp = [vHPC_shock_norm.id,vHPC_shock_norm.resp]; id_temp= temp(temp(:,2)==1,1); 
            vHPC_id = [vHPC_id;id_temp,ones(size(id_temp,1),1)*t];
        catch
            vHPC= []; disp(['No vHPC_shock_norm.mat file in ',session]);
        end
        
    end 

end 

%Proportion 
dhpc_prop = (dHPC_shock/dHPC_pc)*100;
vhpc_prop = (vHPC_shock/vHPC_pc)*100;

%%%%%%%%%%%%%%%%%%%% Density heat maps plots per %%%%%%%%%%%%%%%%%%%%%%%%
dhpc_shock = []; dhpc_resp = []; 
vhpc_shock = []; vhpc_resp = [];

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
        
        disp('Uploading dHPC session shock pc');
        try
            load('dHPC_shock_norm.mat');
            dhpc_shock = [dhpc_shock;dHPC_shock_norm.curve];
            dhpc_resp = [dhpc_resp; dHPC_shock_norm.resp]; 

        catch 
            dHPC = []; disp(['No dHPC file in ',session]);
        end 
        disp('Uploading vHPC session shock pc');
        try
            load('vHPC_shock_norm.mat');
            vhpc_shock = [vhpc_shock;vHPC_shock_norm.curve];
            vhpc_resp = [vhpc_resp; vHPC_shock_norm.resp]; 
     
        catch
            vHPC= []; disp(['No vHPC file in ',session]);
        end
        
    end 

end 


%%%dHPC%%%%

tempA=dhpc_shock; 

[h idx] = max (tempA, [],2);
[m mm] = sort(idx); 


figure(1);clf;hold on; 
subplot(2,1,1);imagesc([0:1:4], [1:1:size(dhpc_shock,1)],tempA(mm,:)), colormap 'gray'; axis tight
title('dHPC Aversive');xline(2,'r'); xline(3,'r'); xlabel(['Shock',char(10),'Time(sec)', char(10),num2str(dhpc_prop),'% resp']); 
ylabel('Nuerons');

subplot(2,1,2);hold on;plot([1:1:size(dhpc_shock,2)],nanmean(tempA(dhpc_resp==1,:))); 
ciplot(nanmean(tempA(dhpc_resp==1,:))-nansem(tempA(dhpc_resp==1,:)),nanmean(tempA(dhpc_resp==1,:))+nansem(tempA(dhpc_resp==1,:)),[1:1:size(dhpc_shock,2)]); alpha 0.1; hold on;
xline(20,'r');xline(30,'r');ylabel('Mean gain');title('Shock neurons');ylim([0.5 5.5]); xlim([0 40]); 

%window for stats
w_size = 2; % # of bins from the center. Each bin 3cm.
win_ave = nanmean(tempA(:, 20-w_size:60+w_size), 2); 
p = signrank(win_ave,win_rew) 

%%%vHPC%%%%
tempA=vhpc_shock;

[h idx] = max (tempA, [],2);
[m mm] = sort(idx); 

figure(2);clf;hold on  
subplot(2,1,1);imagesc([0:1:4], [1:1:size(vhpc_shock,1)],tempA(mm,:)), colormap 'gray'; axis tight
title('vHPC Aversive');xline(2,'r');xline(3,'r');xlabel(['Shock',char(10),'Time(sec)', char(10),num2str(vhpc_prop),'% resp']); 
ylabel('Neurons'); 
subplot(2,1,2);hold on; plot([1:1:size(vhpc_shock,2)],nanmean(tempA));
ciplot(nanmean(tempA)-nansem(tempA),nanmean(tempA)+nansem(tempA),[1:1:size(dhpc_shock,2)]); alpha 0.1; hold on;
xline(20,'r');xline(30,'r');ylim([0.5 5.5]);ylabel('Mean gain'); title('Shock neurons'); xlim([0 40]); 

% Stats 
%window for stats
w_size = 3; % # of bins from the center. Each bin 3cm.
win_ave = nanmean(tempA(:, 60-w_size:60+w_size), 2); 
p = signrank(win_ave,win_rew)
%%%%%%%%%%%%%%

%  Check the firng map of the shock responsivness cells 
%Output matrix 
dhpc_map_ave = []; dhpc_map_rew = [];
dhpc_no_ave = []; dhpc_no_rew = [];
vhpc_map_ave = []; vhpc_map_rew = [];
vhpc_no_ave = []; vhpc_no_rew = [];

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
            for id=1:size(dHPC_id,1)
                neuron = dHPC_id(id,1);
                ses= dHPC_id(id,2);
                if and(neuron==pc.id, ses==t)
                    A = pc.frMap_ave - min(pc.frMap_ave);A = A ./ max(A);
                    dhpc_map_ave = [dhpc_map_ave;A];
                    A = pc.frMap_rew - min(pc.frMap_rew);A = A ./ max(A);
                    dhpc_map_rew = [dhpc_map_rew;A];
                    
                else 
                    A = pc.frMap_ave - min(pc.frMap_ave);A = A ./ max(A);
                    dhpc_no_ave = [dhpc_no_ave;A];
                    A = pc.frMap_rew - min(pc.frMap_rew);A = A ./ max(A);
                    dhpc_no_rew = [dhpc_no_rew;A];
                end 
            end 
        end
        catch 
            dHPC = []; disp(['No dHPC_pc_lap.mat file in ',session]);
        end 
        
        try
            load('vHPC_pc.mat');
            for d=1:size(vHPC,2)
                pc = vHPC{d}; 
                for id=1:size(vHPC_id,1)
                    neuron = vHPC_id(id,1);
                    ses= vHPC_id(id,2);
                    if and(neuron==pc.id, ses==t)
                        A = pc.frMap_ave - min(pc.frMap_ave);A = A ./ max(A);
                        vhpc_map_ave = [vhpc_map_ave;A];
                        A = pc.frMap_rew - min(pc.frMap_rew);A = A ./ max(A);
                        vhpc_map_rew = [vhpc_map_rew;A];
                    else 
                        A = pc.frMap_ave - min(pc.frMap_ave);A = A ./ max(A);
                        vhpc_no_ave = [vhpc_no_ave;A];
                        A = pc.frMap_rew - min(pc.frMap_rew);A = A ./ max(A);
                        vhpc_no_rew = [vhpc_no_rew;A];
                    end 
                end 
            end 
     
        catch
            vHPC= []; disp(['No vHPC_pc_lap.mat file in ',session]);
        end
        
    end 

   
     
end 


% DHPC

[h idx] = max (dhpc_map_ave, [],2);
[m mm] = sort(idx); 

figure(1);clf;hold on; 
subplot(1,2,1);imagesc([0:50:150], [1:1:size(dhpc_map_ave,1)],dhpc_map_ave(mm,:)), colormap 'gray'; title('Aversive')

[h idx] = max (dhpc_map_rew, [],2);
[m mm] = sort(idx); 

subplot(1,2,2);imagesc([0:50:150], [1:1:size(dhpc_map_rew,1)],dhpc_map_rew(mm,:)), colormap 'gray'; title('Reward')
sgtitle('dHPC Firing maps SHOCK cells');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[h idx] = max (dhpc_no_ave, [],2);
[m mm] = sort(idx); 

figure(3);clf;hold on; 
subplot(1,2,1);imagesc([0:50:150], [1:1:size(dhpc_no_ave,1)],dhpc_no_ave(mm,:)), colormap 'gray'; title('Aversive')

[h idx] = max (dhpc_no_rew, [],2);
[m mm] = sort(idx); 

subplot(1,2,2);imagesc([0:50:150], [1:1:size(dhpc_no_rew,1)],dhpc_no_rew(mm,:)), colormap 'gray'; title('Reward')
sgtitle('dHPC Firing maps no SHOCK cells');


%%%%%%%%%%%%%%%VHPC

[h idx] = max (vhpc_map_ave, [],2);
[m mm] = sort(idx); 

figure(2);clf;hold on; 
subplot(1,2,1);imagesc([0:50:150], [1:1:size(vhpc_map_ave,1)],vhpc_map_ave(mm,:)), colormap 'gray'; title('Aversive')

[h idx] = max (vhpc_map_rew, [],2);
[m mm] = sort(idx); 

subplot(1,2,2);imagesc([0:50:150], [1:1:size(vhpc_map_rew,1)],vhpc_map_rew(mm,:)), colormap 'gray'; title('Reward')
sgtitle('VHPC Firing maps SHOCK cells');



%% Remapping parameters WITHOUT SHOCK resp 

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
        % Definition of the emotional transition
        if aversiveTS_run(1)<rewardTS_run(1)
            cond = 1; % 1 if the order was Aversive -> Reward
        else
            cond = 2;% 1 if the order was Reward -> Aversive
        end
        
        %% Firing curve calculation without shock resp
        %Session output
        
        dHPC = {};% one cell per pc
        vHPC = {};
    
        %%%%%% dHPC %%%%%%%
        try
            load('dHPC_shock_norm.mat');
            
            valid = dHPC_shock_norm.id(dHPC_shock_norm.resp==0);
            group_dHPC = group_dHPC(ismember(group_dHPC(:,1),valid),:); % keep only non responsive pc 
            
            disp('dHPC Remapping parameters computation')
             
            for ii=1:size(group_dHPC,1)

            cluster = group_dHPC(ii,1);
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
                    

            %Within aversive
            [withinA,withinA_tresh] = Within_pc(pos_ave,spks_ave,1,sigma,Xedges);
            %Within reward 
            [withinR,withinR_tresh] = Within_pc(pos_rew,spks_rew,1,sigma,Xedges);
            %Between
            [between] = Between_pc(pos_ave,spks_ave,pos_rew,spks_rew,bin_size,sigma,Xedges);

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
   
            dHPC{ii}= n;

           end
         
        catch 
            dHPC = []; disp(['No dHPC_pc.mat file in ',session]);
        end 
        
        %%%%%% vHPC %%%%%%%
        try
           load('vHPC_shock_norm.mat');
            
            valid = vHPC_shock_norm.id(vHPC_shock_norm.resp==0);
            group_vHPC = group_vHPC(ismember(group_vHPC(:,1),valid),:); % keep only non responsive pc 
            
            disp('vHPC Remapping parameters computation')
             
            for ii=1:size(group_vHPC,1)

            cluster = group_vHPC(ii,1);
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
                    

            %Within aversive
            [withinA,withinA_tresh] = Within_pc(pos_ave,spks_ave,1,sigma,Xedges);
            %Within reward 
            [withinR,withinR_tresh] = Within_pc(pos_rew,spks_rew,1,sigma,Xedges);
            %Between
            [between] = Between_pc(pos_ave,spks_ave,pos_rew,spks_rew,bin_size,sigma,Xedges);

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
   
            vHPC{ii}= n;

           end 
        
        catch
            vHPC= []; disp(['No vHPC_pc.mat file in ',session]);
        end

        %% Saveing PC INFO 
        if ~isempty(dHPC)
            dHPC = dHPC(~cellfun('isempty',dHPC));
            save([cd,'\dHPC_shock_remap.mat'],'dHPC'); 
        end
        if ~isempty(vHPC)
            vHPC = vHPC(~cellfun('isempty',vHPC));
            save([cd,'\vHPC_shock_remap.mat'],'vHPC'); 
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

%% PLOTS ans STATS: SHOCK PROP and REMAPPING WITHOUT SHOCK

%c1 = spatial c2= fr_change c3=  overlap c4: pf shift c5:  1 = between 2 = within aversive 3= within reward 
clear
clc
close all
% Parameters
path = {'\\Maryjackson\e\Rat127\Ephys\pyr';'\\Maryjackson\e\Rat128\Ephys\in_pyr\ready';'\\Maryjackson\e\Rat103\usable';'\\Maryjackson\e\Rat132\recordings\in_pyr'; '\\Maryjackson\e\Rat165\in_pyr'; ...
    '\\Maryjackson\e\Rat126\Ephys\in_Pyr'};%List of folders from the path

%Output matrix 
dhpc_sub = [];dhpc_n=0; 
vhpc_sub = [];vhpc_n=0;

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
            load('dHPC_shock_remap.mat');
             %Save session parameters 
        for d=1:size(dHPC,2)
            pc = dHPC{d}; 
            bet = [pc.between,1]; 
            wa = [pc.within_ave,2];
            wr = [pc.within_rew,3];
            dhpc_temp = [dhpc_temp;bet;wa;wr]; 
            dhpc_n = dhpc_n+1;
        end
        catch 
            dHPC = []; disp(['No dHPC_pc_lap.mat file in ',session]);
        end 
        
        try
            load('vHPC_shock_remap.mat');
            for d=1:size(vHPC,2)
                pc = vHPC{d}; 
                bet = [pc.between,1]; 
                wa = [pc.within_ave,2];
                wr = [pc.within_rew,3];
                vhpc_temp = [vhpc_temp;bet;wa;wr];
                vhpc_n = vhpc_n+1;
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
sgtitle('dHPC without shock neurons')
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
sgtitle('dHPC without shock cells')
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
sgtitle('vHPC without shock neurons')
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

%Only 3 plots
ylabels ={'Spatial correlation', 'Overlap', 'Pf shift'};
idx = [1,3,4]; 
xlabels = {'Between', 'Within Ave', 'Within Rew'}; 
figure(2);clf;hold on, 
sgtitle('vHPC without shock cells')
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


%% REMAPPING PARAMETERS without REWARD  
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
        % Definition of the emotional transition
        if aversiveTS_run(1)<rewardTS_run(1)
            cond = 1; % 1 if the order was Aversive -> Reward
        else
            cond = 2;% 1 if the order was Reward -> Aversive
        end
        
        %% Firing curve calculation without shock resp
        
        %Session output
        dHPC = {};% one cell per pc
        vHPC = {};
    
        %%%%%% dHPC %%%%%%%
        try
            load('dHPC_pc.mat');
            valid = []; 
            for d=1:size(dHPC,2)
                pc=dHPC{d};
                % Reward
                pf_lim = pc.stats_rew.fieldX; 
                pc_frmap = pc.frMap_rew; 
                field = pc_frmap(pf_lim(1,1):pf_lim(1,2));
                
                %Center of mass of the field 
                idx = 1:size(field,2);                             
                com = sum(idx .* field) / sum(field);
                com = com + pf_lim(1,1); % back in spatial bin scale
                
                th= 4;%# of bins to dicard 
                if or(com >1+th, com<60-th)   
                    valid = [valid;pc.id]; 
                end
            end 
            clear d pc pf_lim  pc_frmap field com idx
            
            group_dHPC = group_dHPC(ismember(group_dHPC(:,1),valid),:); % keep only non responsive pc 
            
            disp('dHPC Remapping parameters computation without reward neurons')
             
            for ii=1:size(group_dHPC,1)

            cluster = group_dHPC(ii,1);
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
                    

            %Within aversive
            [withinA,withinA_tresh] = Within_pc(pos_ave,spks_ave,1,sigma,Xedges);
            %Within reward 
            [withinR,withinR_tresh] = Within_pc(pos_rew,spks_rew,1,sigma,Xedges);
            %Between
            [between] = Between_pc(pos_ave,spks_ave,pos_rew,spks_rew,bin_size,sigma,Xedges);

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
   
            dHPC{ii}= n;

           end
         
        catch 
            dHPC = []; disp(['No dHPC_pc.mat file in ',session]);
        end 
        
        %%%%%% vHPC %%%%%%%
        try
           load('vHPC_pc.mat');
      
            valid = []; 
            for d=1:size(vHPC,2)
                pc=vHPC{d};
                % Reward
                pf_lim = pc.stats_rew.fieldX; 
                pc_frmap = pc.frMap_rew; 
                field = pc_frmap(pf_lim(1,1):pf_lim(1,2));
                
                %Center of mass of the field 
                idx = 1:size(field,2);                             
                com = sum(idx .* field) / sum(field);
                com = com + pf_lim(1,1); % back in spatial bin scale
                
                th= 4;%# of bins to dicard 
                if or(com >1+th, com<60-th)   
                    valid = [valid;pc.id]; 
                end
            end 
            clear d pc pf_lim  pc_frmap field com idx 

            group_vHPC = group_vHPC(ismember(group_vHPC(:,1),valid),:); % keep only non responsive pc 
            
            disp('vHPC Remapping parameters computation')
             
            for ii=1:size(group_vHPC,1)

            cluster = group_vHPC(ii,1);
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
                    

            %Within aversive
            [withinA,withinA_tresh] = Within_pc(pos_ave,spks_ave,1,sigma,Xedges);
            %Within reward 
            [withinR,withinR_tresh] = Within_pc(pos_rew,spks_rew,1,sigma,Xedges);
            %Between
            [between] = Between_pc(pos_ave,spks_ave,pos_rew,spks_rew,bin_size,sigma,Xedges);

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
   
            vHPC{ii}= n;

           end 
        
        catch
            vHPC= []; disp(['No vHPC_pc.mat file in ',session]);
        end

        %% Saveing PC INFO 
        if ~isempty(dHPC)
            dHPC = dHPC(~cellfun('isempty',dHPC));
            save([cd,'\dHPC_rew_remap.mat'],'dHPC'); 
        end
        if ~isempty(vHPC)
            vHPC = vHPC(~cellfun('isempty',vHPC));
            save([cd,'\vHPC_rew_remap.mat'],'vHPC'); 
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

%% Plots REMAPPING PARAMETERS without REWARD 
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
n_dhpc= 0; n_vhpc=0; 
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
            load('dHPC_rew_remap.mat');
             %Save session parameters 
        for d=1:size(dHPC,2)
            pc = dHPC{d}; 
            bet = [pc.between,1]; 
            wa = [pc.within_ave,2];
            wr = [pc.within_rew,3];
            dhpc_temp = [dhpc_temp;bet;wa;wr]; 
            n_dhpc= n_dhpc + 1;
            
        end
        catch 
            dHPC = []; disp(['No dHPC_pc_lap.mat file in ',session]);
        end 
        
        try
            load('vHPC_rew_remap.mat');
            for d=1:size(vHPC,2)
                pc = vHPC{d}; 
                bet = [pc.between,1]; 
                wa = [pc.within_ave,2];
                wr = [pc.within_rew,3];
                vhpc_temp = [vhpc_temp;bet;wa;wr];
                n_vhpc= n_dhpc + 1;
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
sgtitle('dHPC without reward cells')
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
sgtitle('vHPC without reward')
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
sgtitle('vHPC without reward cells')
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
[P,ANOVATAB,STATS] = kruskalwallis(vhpc_sub(:,5),vhpc_sub(:,6));
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

%%%%%%% Firing curves %%%%%%
%Output matrix 
dhpc_frmap_ave = [];dhpc_frmap_rew = [];
vhpc_frmap_ave = [];vhpc_frmap_rew = [];

for tt = 1:length(path)
    %List of folders from the path
    files = dir(path{tt});
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
    clear files dirFlags

    for t = 1 : length(subFolders)-2
        disp('');
        disp(['-- Initiating analysis of folder #' , num2str(t) , ' from rat #',num2str(tt) , ' --'])
        session = [subFolders(t+2).folder,'\',subFolders(t+2).name];
        cd([session,'\Spikesorting'])
        %Loading pc matrix of the sessions
        disp('Uploading session pc matrix');
        try
            load('dHPC_rew_remap.mat');
             %Save session parameters 
        for d=1:size(dHPC,2)
            pc = dHPC{d}; 
 
            A = pc.frMap_ave - min(pc.frMap_ave);
            A = A ./ max(A);
           
    
            R = pc.frMap_rew - min(pc.frMap_rew);
            R = R ./ max(R);
         
            dhpc_frmap_ave = [dhpc_frmap_ave;A];
            dhpc_frmap_rew = [dhpc_frmap_rew;R];
            
        end
        catch 
            dHPC = []; disp(['No dHPC file in ',session]);
        end 
        
        try
            load('vHPC_rew_remap.mat');
            for d=1:size(vHPC,2)
                pc = vHPC{d}; 
                A = pc.frMap_ave - min(pc.frMap_ave);
                A = A ./ max(A);
 
                R = pc.frMap_rew - min(pc.frMap_rew);
                R = R ./ max(R);
         
                vhpc_frmap_ave = [vhpc_frmap_ave;A];
                vhpc_frmap_rew = [vhpc_frmap_rew;R];
            end 
     
        catch
            vHPC= []; disp(['No vHPC_pc_lap.mat file in ',session]);
        end
        
    end 
     
end 


% DHPC

[h idx] = max (dhpc_frmap_rew, [],2);
[m mm] = sort(idx); 
figure(1);clf;hold on; 

subplot(1,2,1);imagesc([10:40:140], [1:1:size(dhpc_frmap_rew,1)],dhpc_frmap_rew(mm,:)), colormap 'gray'; title('Rewarded')
subplot(1,2,2);imagesc([10:40:140], [1:1:size(dhpc_frmap_ave,1)],dhpc_frmap_ave(mm,:)), colormap 'gray'; title('Aversive')
sgtitle('dHPC without "reward" cells');

%VHPC

[h idx] = max (vhpc_frmap_rew, [],2);
[m mm] = sort(idx); 
figure(2);clf;hold on; 

subplot(1,2,1);imagesc([10:40:140], [1:1:size(vhpc_frmap_rew,1)],vhpc_frmap_rew(mm,:)), colormap 'gray'; title('Rewarded')
subplot(1,2,2);imagesc([10:40:140], [1:1:size(vhpc_frmap_ave,1)],vhpc_frmap_ave(mm,:)), colormap 'gray'; title('Aversive')
sgtitle('vHPC without "reward" cells');

