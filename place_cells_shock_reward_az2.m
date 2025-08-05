clear
clc
close all
%% PARAMETERS - run this section before run any other
path = {'E:\Rat103\usable';'E:\Rat126\Ephys\in_Pyr';'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat132\recordings\in_pyr';'E:\Rat165\in_pyr\'};%List of folders from the path
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
         bin = 0.1;smooth=[1]; th=1 ;window=10;
        
         dHPC_shock_norm={};
         vHPC_shock_norm={};
       
         dHPC_valve_norm={};
         vHPC_valve_norm={};
         
        %%%%%%% PYR 
        group_dHPC = group_dHPC(group_dHPC(:,4)==1,:);
        group_vHPC = group_vHPC(group_vHPC(:,4)==1,:);
        
        if ~isempty(group_dHPC)
             disp('dHPC Responsivness')
            %Aversive
            curve=[]; 
            limits=[behavior.speed.aversive(1,1),behavior.speed.aversive(end,1)];     
            [curve , bins , responsive] = SU_responsivness_fixed(spks_dHPC,group_dHPC(:,1),Shocks_filt,...
                                      limits,[0 1],window,bin,smooth,'zscore',th); 

            %Save output
            dHPC_resp.id = group_dHPC(:,1); 
            dHPC_resp.curve_ave=curve'; 
            dHPC_resp.resp_ave=responsive'; 
            save([cd,'\dHPC_Shock_VF.mat'],'dHPC_resp'); 
        end
        
        if ~isempty(group_vHPC)
            disp('vHPC Responsivness')
            %Aversive
            curve=[]; 
            limits=[behavior.speed.aversive(1,1),behavior.speed.aversive(end,1)];           
            [curve , bins , responsive] = SU_responsivness_fixed(spks_vHPC,group_vHPC(:,1),Shocks_filt,...
                                      limits,[0 1],window,bin,smooth,'zscore',th);                                  
           
            %Save output
            vHPC_resp.id = group_vHPC(:,1); 
            vHPC_resp.curve_ave=curve'; 
            vHPC_resp.resp_ave=responsive'; 
            save([cd,'\vHPC_Shock_VF.mat'],'vHPC_resp');
                                  
        end 
        
        
%         %%%%%%% PC
%         try
%             load('dHPC_pc_skaggs_circular.mat');
%             dhpc_pc = [];
%             for d=1:size(dHPC,2)
%                 pc = dHPC{d}; 
%                 dhpc_pc = [dhpc_pc;pc.id]; 
%             end
%             
%             group_dHPC = group_dHPC(ismember(group_dHPC(:,1),dhpc_pc),:); % keep only pc 
%             
%             disp('dHPC Responsivness')
%             %Aversive
%             curve=[]; 
%             limits=[behavior.speed.aversive(1,1),behavior.speed.aversive(end,1)];     
%             [curve , bins , responsive] = SU_responsivness_fixed(spks_dHPC,group_dHPC(:,1),Shocks_filt,...
%                                       limits,[0 1],window,bin,smooth,'zscore',th); 
%              
%             %Reward-valve
%             curveR=[]; 
%             window = 7;
%             limits=[behavior.speed.reward(1,1),behavior.speed.reward(end,1)];     
%             [curveR , bins , responsiveR] = SU_responsivness_fixed(spks_dHPC,group_dHPC(:,1),sort(Rewards_filt(:,1)),...
%                                       limits,[0 5],window,bin,smooth,'zscore',th);                       
%             %Save output
%             dHPC_resp.id = group_dHPC(:,1); 
%             dHPC_resp.curve_ave=curve'; 
%             dHPC_resp.resp_ave=responsive'; 
%             dHPC_resp.curve_rew=curveR'; 
%             dHPC_resp.resp_rew=responsiveR';
%             save([cd,'\dHPC_responsivness.mat'],'dHPC_resp'); 
% 
%         catch 
%             dHPC = []; disp(['No dHPC_pc.mat file in ',session]);
%         end 
%         
%         try
%             load('vHPC_pc_skaggs_circular.mat');
%             vhpc_pc = [];
%             for d=1:size(vHPC,2)
%                 pc = vHPC{d}; 
%                 vhpc_pc = [vhpc_pc;pc.id]; 
%             end
%             group_vHPC = group_vHPC(ismember(group_vHPC(:,1),vhpc_pc),:);
%             
%             disp('vHPC Responsivness')
%             %Aversive
%             curve=[]; 
%             limits=[behavior.speed.aversive(1,1),behavior.speed.aversive(end,1)];           
%             [curve , bins , responsive] = SU_responsivness_fixed(spks_vHPC,group_vHPC(:,1),Shocks_filt,...
%                                       limits,[0 1],4,bin,smooth,'zscore',th);                                  
%            
%             curveR=[]; 
%             window = 7;
%             limits=[behavior.speed.reward(1,1),behavior.speed.reward(end,1)];     
%             [curveR , bins , responsiveR] = SU_responsivness_fixed(spks_vHPC,group_vHPC(:,1),sort(Rewards_filt(:,1)),...
%                                       limits,[0 5],window,bin,smooth,'zscore',th);                       
%             %Save output
%             vHPC_resp.id = group_vHPC(:,1); 
%             vHPC_resp.curve_ave=curve'; 
%             vHPC_resp.resp_ave=responsive'; 
%             vHPC_resp.curve_rew=curveR'; 
%             vHPC_resp.resp_rew=responsiveR';
%             save([cd,'\vHPC_responsivness.mat'],'vHPC_resp');
%         
%         catch
%             vHPC= []; disp(['No vHPC_pc.mat file in ',session]);
%         end


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
