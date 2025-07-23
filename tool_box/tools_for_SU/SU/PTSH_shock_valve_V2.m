clear
clc
close all

%Parameters

%%
% $$\$$

path = {'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat103\usable'};
% path = 'E:\Rat128\Ephys\in_pyr\ready';

%% Parameters
% ---------------
% --> SU
% ---------------
binsize = 0.05; % bin size for SU histogram in sec
ss = 2; %smooth level of CCG
d = 4; %duration of the CCG
n_SU_V = 0;
n_SU_D = 0;

celltype = 1; %to define which celltype to use (1: pyr / 2:int)

%for storing Peri-event Histogram
%Behavior
dHPC_Shock = []; vHPC_Shock = [];
dHPC_Valve = []; vHPC_Valve = [];
% Ripples
%dHPC
dHPC_dRipples = [];
dHPC_dRipplesA = []; dHPC_dRipplesR = []; dHPC_dRipplesB = [];
dHPC_dRipples_coordinatedA = []; dHPC_dRipples_coordinatedR = []; dHPC_dRipples_coordinatedB = [];
dHPC_dRipples_uncoordinatedA = []; dHPC_dRipples_uncoordinatedR = []; dHPC_dRipples_uncoordinatedB = [];
dHPC_vRipples = [];
dHPC_vRipples_coordinatedA = []; dHPC_vRipples_coordinatedR = []; dHPC_vRipples_coordinatedB = [];

%vHPC
vHPC_vRipples = [];
vHPC_vRipplesA = []; vHPC_vRipplesR = []; vHPC_vRipplesB = [];
vHPC_vRipples_coordinatedA = []; vHPC_vRipples_coordinatedR = []; vHPC_vRipples_coordinatedB = [];
vHPC_vRipples_uncoordinatedA = []; vHPC_vRipples_uncoordinatedR = []; vHPC_vRipples_uncoordinatedB = [];
vHPC_dRipples = [];
vHPC_dRipples_coordinatedA = []; vHPC_dRipples_coordinatedR = []; vHPC_dRipples_coordinatedB = [];


% Poisson for up or down modulated
dHPC_dRipples_poisson = [];
vHPC_vRipples_poisson = [];
dHPC_vRipples_poisson = [];
vHPC_dRipples_poisson = [];
cutting = 63; % final channel from vHPC

velocity_shock = [];%To store speed data locked to the events
velocity_valve = [];


% Parameters for CCG ripples-SU
b = 0.05; % binsize for ripples-SU modulation
dd = 1; % time window for ripples-SU modulation
sss = 1; % smooth level for ripples-SU modulation


% ---------------
% --> Behavior
% ---------------
pos_neck_valve = [];
pos_neck_shock = [];
pos_snout_valve = [];
pos_snout_shock = [];

% ---------------
% --> Behavior
% ---------------
CD = [];
CV = [];

%% -------------------------------------
%Main loop
for tt = 1:length(path)
    %List of folders from the path
    files = dir(path{tt});
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
    clear files dirFlags
    
    for t = 1 : length(subFolders)-2
        session = [subFolders(t+2).folder,'\',subFolders(t+2).name];
        cd(session)
%         load([cd,'\lfp.mat'])
%         Time = dHPC(:,1);
        
        %Loading TS of the sessions
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
        % Load digitalin:mat
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
        clear count deff shock
        
        %Rewards selection
        Rewards_filt = Restrict([leftvalve ; rightvalve],rewardTS_run./1000);
        
        if aversiveTS_run(1) > rewardTS_run
            load('laps1')
            posXY_reward = [posx,posy];
            posXY_snout_reward = [posx_snout,posy_snout];
            velocity_reward = velocity;
            clear posx posy velocity
            load('laps2')
            posXY_aversive = [posx,posy];
            posXY_snout_aversive = [posx_snout,posy_snout];
            velocity_aversive = velocity;
            
            clear posx posy velocity
        else
            load('laps2')
            posXY_reward = [posx,posy];
            posXY_snout_reward = [posx_snout,posy_snout];
            velocity_reward = velocity;
            clear posx posy velocity
            load('laps1')
            posXY_aversive = [posx,posy];
            posXY_snout_aversive = [posx_snout,posy_snout];
            velocity_aversive = velocity;
            clear posx posy velocity
        end
        clear ch
        
        %to restrict velocity sourrounding the events
        %Sync camare TTLs with speed and pos
        %aversive
        TTLs_camara_aversive = Restrict(camara,aversiveTS_run./1000);
        TTLs_camara_aversive = ((TTLs_camara_aversive(:,2) - TTLs_camara_aversive(:,1))/2)+TTLs_camara_aversive(:,1);
        velocity_aversive = [TTLs_camara_aversive(1:length(velocity_aversive)) , velocity_aversive];
        pos_aversive = [TTLs_camara_aversive(1:length(posXY_aversive)) , posXY_aversive(:,1)];
        pos_snout_aversive = [TTLs_camara_aversive(1:length(posXY_snout_aversive)) , posXY_snout_aversive(:,1)];
        %reward
        TTLs_camara_reward = Restrict(camara,rewardTS_run./1000);
        TTLs_camara_reward = ((TTLs_camara_reward(:,2) - TTLs_camara_reward(:,1))/2)+TTLs_camara_reward(:,1);
        velocity_reward = [TTLs_camara_reward(1:length(velocity_reward)) , velocity_reward];
        pos_reward = [TTLs_camara_reward(1:length(posXY_reward)) , posXY_reward(:,1)];
        pos_snout_reward = [TTLs_camara_reward(1:length(posXY_snout_reward)) , posXY_snout_reward(:,1)];
        
        
        %speed restriction
        %aversive
        tmp = [];
        for jj = 1:length(Shocks_filt)
            r = [Shocks_filt(jj)-d/2 , Shocks_filt(jj)+d/2];
            rr = Restrict(velocity_aversive,r);
            if (~isempty(rr) && length(rr)>=d*30 && sum(isnan(rr(:,2))) == 0)
                tmp = [tmp , rr(1:d*30,2)];
            end
            clear r rr
        end
        velocity_shock = [velocity_shock , mean(tmp,2)];
        pos_neck_shock = [pos_neck_shock ; Restrict(pos_aversive,[Shocks_filt(:,1) Shocks_filt(:,1)+1])];               
        pos_snout_shock = [pos_snout_shock ; Restrict(pos_snout_aversive,[Shocks_filt(:,1) Shocks_filt(:,1)+1])];               
        clear tmp velocity_aversive pos_aversive pos_snout_aversive
        
        %reward
        tmp = [];
        for jj = 1:length(Rewards_filt)
            r = [Rewards_filt(jj,1)-d/2 , Rewards_filt(jj,1)+d/2];
            rr = Restrict(velocity_reward,r);
            if (~isempty(rr) && length(rr)>=d*30 && sum(isnan(rr(:,2))) == 0)
                tmp = [tmp , rr(1:d*30,2)];
%                 figure,plot(rr(:,2))
            end
            clear r rr 
        end
        velocity_valve = [velocity_valve , mean(tmp,2)];
        pos_neck_valve = [pos_neck_valve ; Restrict(pos_reward,[Rewards_filt(:,1)+0.2 Rewards_filt(:,1)+1])];               
        pos_snout_valve = [pos_snout_valve ; Restrict(pos_snout_reward,[Rewards_filt(:,1)+0.2 Rewards_filt(:,1)+1])];               
        clear tmp velocity_reward pos_reward pos_snout_reward       
        
        %% Sleep
        x = dir([cd,'\*-states.mat']);    states = load([cd,'\',x.name]);    states = states.states;
        
        REM = ToIntervals(states==5);    NREM = ToIntervals(states==3);    WAKE = ToIntervals(states==1);
        clear x states
        
        % NREM events restriction according conditions
        NREM_B = Restrict(NREM,baselineTS./1000);
        NREM_R = Restrict(NREM,rewardTS./1000);
        NREM_A = Restrict(NREM,aversiveTS./1000);

        % REM events restriction according conditions
        REM_B = Restrict(REM,baselineTS./1000);
        REM_R = Restrict(REM,rewardTS./1000);
        REM_A = Restrict(REM,aversiveTS./1000);

                
%         load('detected_ripples.mat')
        ripplesD = table2array(readtable('ripplesD_customized2.csv'));
        ripplesV = table2array(readtable('ripplesV_customized2.csv'));
        
        % Coordinated dHPC ripples
        coordinated = [];
        coordinatedV = [];
        coordinatedV_refined = [];
        uncoordinated = [];
        for i = 1:length(ripplesD)
            r = ripplesD(i,:);
            tmp = sum(and(ripplesV(:,2)>= r(1,2)-0.2, ripplesV(:,2)<= r(1,2)+0.2));
            if tmp>0
                z = ripplesV(and(ripplesV(:,2)>= r(1,2)-0.2, ripplesV(:,2)<= r(1,2)+0.2),:);
                coordinatedV = [coordinatedV ; z];
                [p,indice] = min(abs(r(2)-z(:,2)));
                coordinatedV_refined = [coordinatedV_refined ; z(indice,:)];
                coordinated = [coordinated ; r];
                clear tmp2 tmp1
            end
            clear r
        end
        clear x tmp i
        
        coordinatedB = Restrict(coordinated,NREM_B);    coordinatedA = Restrict(coordinated,NREM_A);    coordinatedR = Restrict(coordinated,NREM_R);
        coordinatedB_V = Restrict(coordinatedV_refined,NREM_B);    coordinatedR_V = Restrict(coordinatedV_refined,NREM_R);    coordinatedA_V = Restrict(coordinatedV_refined,NREM_A);
        
        % Detection of uncoordinated ripples
        uncoordinated = ripplesD(~ismember(ripplesD(:,1),coordinated(:,1)),:);
        uncoordinatedV = ripplesV(~ismember(ripplesV(:,1),coordinatedV(:,1)),:);
        
        uncoordinatedB = Restrict(uncoordinated,NREM_B);    uncoordinatedA = Restrict(uncoordinated,NREM_A);    uncoordinatedR = Restrict(uncoordinated,NREM_R);
        uncoordinatedB_V = Restrict(uncoordinatedV,NREM_B);    uncoordinatedR_V = Restrict(uncoordinatedV,NREM_R);    uncoordinatedA_V = Restrict(uncoordinatedV,NREM_A);
    
        
    %% Spikes
    %Load Units
    cd 'Spikesorting'
    spks = double([readNPY('spike_times.npy') readNPY('spike_clusters.npy')]);
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
            n_SU_D = n_SU_D + length(group_dHPC);
            spks_dHPC = spks(ismember(spks(:,2),group_dHPC(:,1)),:); %keep spks from good clusters
        else
            n_SU_V = n_SU_V + length(group_vHPC);
            spks_vHPC = spks(ismember(spks(:,2),group_vHPC(:,1)),:); %keep spks from good clusters
        end
    end
    clear z
    spks_vHPC(:,1) = double(spks_vHPC(:,1))./20000;
    spks_dHPC(:,1) = double(spks_dHPC(:,1))./20000;
 
        %% Definition of cell type population to analize
        if celltype == 1
            tipoD = group_dHPC(:,3);
            tipoV = group_vHPC(:,3);
        elseif celltype == 2
            tipoD = group_dHPC(:,4);
            tipoV = group_vHPC(:,4);
        end
        
    %% PHIST - shock
    % -----------------------------------
    % Shock
    % -----------------------------------
    % dHPC
    base = InvertIntervals([Shocks_filt-4, Shocks_filt+4],aversiveTS_run(1)/1000 , aversiveTS_run(2)/1000);
    for i = 1 : length(group_dHPC)% ventral SU
        if logical(tipoD(i))
            y = Shocks_filt;
            x =  Restrict(spks_dHPC(spks_dHPC(:,2) == group_dHPC(i,1),1),aversiveTS_run./1000);
            if and(~isempty(x),length(x)>20)
                Base = Restrict(x,base);
                [times,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
                [ccg,TimeVector1] = CCG(times,groups,'binsize',binsize,'duration',d,'smooth',ss);
                FR = length(Base)/sum(base(:,2)-base(:,1));
                dHPC_Shock = [dHPC_Shock , ccg(:,1,2)./length(y)./binsize./FR];
                clear tmp FR Base ccg x y times ids groups
            else
                dHPC_Shock = [dHPC_Shock , nan(size(dHPC_Shock,1),1)];
                
            end
        end
    end
    clear base i    
    % vHPC
    base = InvertIntervals([Shocks_filt-4, Shocks_filt+4],aversiveTS_run(1)/1000 , aversiveTS_run(2)/1000);
    for i = 1 : length(group_vHPC)% ventral SU
        if logical(tipoV(i))
            y = Shocks_filt;
            x =  Restrict(spks_vHPC(spks_vHPC(:,2) == group_vHPC(i,1),1),aversiveTS_run./1000);
            if and(~isempty(x),length(x)>20)
                Base = Restrict(x,base);
                [times,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
                [ccg,TimeVector1] = CCG(times,groups,'binsize',binsize,'duration',d,'smooth',ss);
                FR = length(Base)/sum(base(:,2)-base(:,1));
                vHPC_Shock = [vHPC_Shock , ccg(:,1,2)./length(y)./binsize./FR];
                clear tmp FR Base ccg x y times ids groups
            else
                vHPC_Shock = [vHPC_Shock , nan(size(vHPC_Shock,1),1)];
            end
        end
    end
    clear base i
    
    % -----------------------------------
    % Valve
    % -----------------------------------
    % dHPC
    base = InvertIntervals([Rewards_filt-4, Rewards_filt+4],rewardTS_run(1)/1000 , rewardTS_run(2)/1000);
    for i = 1 : length(group_dHPC)% ventral SU
        if logical(tipoD(i))
            y = Rewards_filt(:,1);
            x =  Restrict(spks_dHPC(spks_dHPC(:,2) == group_dHPC(i,1),1),rewardTS_run./1000);
            if and(~isempty(x),length(x)>20)
                Base = Restrict(x,base);
                [times,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
                [ccg,TimeVector1] = CCG(times,groups,'binsize',binsize,'duration',d,'smooth',ss);
                FR = length(Base)/sum(base(:,2)-base(:,1));
                dHPC_Valve = [dHPC_Valve , ccg(:,1,2)./length(y)./binsize./FR];
                clear tmp FR Base ccg x y times ids groups
            else
                dHPC_Valve = [dHPC_Valve , nan(size(dHPC_Valve,1),1)];
            end
        end
    end
    clear base i
    
    % vHPC
    base = InvertIntervals([Rewards_filt-4, Rewards_filt+4],rewardTS_run(1)/1000 , rewardTS_run(2)/1000);
    for i = 1 : length(group_vHPC)% ventral SU
        if  logical(tipoV(i))
            y = Rewards_filt(:,1);
            x =  Restrict(spks_vHPC(spks_vHPC(:,2) == group_vHPC(i,1),1),rewardTS_run./1000);
            if and(~isempty(x),length(x)>20)
                Base = Restrict(x,base);
                [times,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
                [ccg,TimeVector1] = CCG(times,groups,'binsize',binsize,'duration',d,'smooth',ss);
                FR = length(Base)/sum(base(:,2)-base(:,1));
                vHPC_Valve = [vHPC_Valve , ccg(:,1,2)./length(y)./binsize./FR];
                clear tmp FR Base ccg x y times ids groups
            else
                vHPC_Valve = [vHPC_Valve , nan(size(vHPC_Valve,1),1)];
            end
        end
    end
    clear base i
    
    %% NOW PHIST surrounding ripples
    % -----------------------------------
    % Dorsal SU PHIST
    % -----------------------------------
    %dHPC SU - dRipples All
    base = InvertIntervals([ripplesD(:,1)-0.1, ripplesD(:,3)+0.1],NREM(:,1) , NREM(:,2));
    for i = 1 : length(group_dHPC)% ventral SU
        if logical(tipoD(i))
            y = ripplesD(:,2);
            x =  Restrict(spks_dHPC(spks_dHPC(:,2) == group_dHPC(i,1),1),NREM);
            Base = Restrict(x,base);
            [times,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,TimeVector2] = CCG(times,groups,'binsize',b,'duration',dd,'smooth',sss);
            FR = length(Base)/sum(base(:,2)-base(:,1));
            dHPC_dRipples = [dHPC_dRipples , ccg(:,1,2)./length(y)./b./FR];
            
            %Poisson
            totalrippletime = sum(ripplesD(:,3)-ripplesD(:,1));
            ripplespikes = Restrict(x,[ripplesD(:,1) ripplesD(:,3)]);
            nripplespikes = size(ripplespikes,1);
            
            ncellbaselinespikes = length(Base);
            ncellripplespikes = length(ripplespikes);
            totalbaselinetime = sum(base(:,2)-base(:,1));
            if ncellbaselinespikes~=0 & ncellripplespikes~=0
                [pInc pDec surp] = poissonTest(ncellbaselinespikes/totalbaselinetime,ncellripplespikes,totalrippletime);
                dHPC_dRipples_poisson = [dHPC_dRipples_poisson ; pInc pDec surp];
            else
                pInc = NaN;
                pDec = NaN;
                surp = NaN;
                dHPC_dRipples_poisson = [dHPC_dRipples_poisson ; pInc pDec surp];
            end
            clear tmp FR Base ccg x y times ids groups
            clear pInc pDec surp totalrippletime ripplespikes nripplespikes
            clear ncellbaselinespikes ncellripplespikes totalbaselinetime
        end
    end
    clear base i     
    
    %dHPC SU - coordinated dRipples
    base = InvertIntervals([ripplesD(:,1)-0.1, ripplesD(:,3)+0.1],NREM(:,1) , NREM(:,2));
    for i = 1 : length(group_dHPC)% ventral SU
        if logical(tipoD(i))
            %Baseline
            y = coordinatedB(:,2);
            x =  Restrict(spks_dHPC(spks_dHPC(:,2) == group_dHPC(i,1),1),NREM);
            Base = Restrict(x,base);
            [times,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,TimeVector2] = CCG(times,groups,'binsize',b,'duration',dd,'smooth',sss);
            FR = length(Base)/sum(base(:,2)-base(:,1));
            dHPC_dRipples_coordinatedB = [dHPC_dRipples_coordinatedB , ccg(:,1,2)./length(y)./b./FR];
            clear tmp FR Base ccg x y times ids groups
            
            %Reward
            y = coordinatedR(:,2);
            x =  Restrict(spks_dHPC(spks_dHPC(:,2) == group_dHPC(i,1),1),NREM);
            Base = Restrict(x,base);
            [times,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,TimeVector2] = CCG(times,groups,'binsize',b,'duration',dd,'smooth',sss);
            FR = length(Base)/sum(base(:,2)-base(:,1));
            dHPC_dRipples_coordinatedR = [dHPC_dRipples_coordinatedR , ccg(:,1,2)./length(y)./b./FR];
            clear tmp FR Base ccg x y times ids groups      
            
            %Aversive
            y = coordinatedA(:,2);
            x =  Restrict(spks_dHPC(spks_dHPC(:,2) == group_dHPC(i,1),1),NREM);
            Base = Restrict(x,base);
            [times,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,TimeVector2] = CCG(times,groups,'binsize',b,'duration',dd,'smooth',sss);
            FR = length(Base)/sum(base(:,2)-base(:,1));
            dHPC_dRipples_coordinatedA = [dHPC_dRipples_coordinatedA , ccg(:,1,2)./length(y)./b./FR];
            clear tmp FR Base ccg x y times ids groups            
        end
    end
    clear base i      
%     
%     %dHPC SU - uncoordinated dRipples
%     base = InvertIntervals([ripplesD(:,1)-0.1, ripplesD(:,3)+0.1],NREM(:,1) , NREM(:,2));
%     for i = 1 : length(group_dHPC)% ventral SU
%         if logical(tipoD(i))
%             %Baseline
%             y = uncoordinatedB(:,2);
%             x =  Restrict(spks_dHPC(spks_dHPC(:,2) == group_dHPC(i,1),1),NREM);
%             Base = Restrict(x,base);
%             [times,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
%             [ccg,TimeVector2] = CCG(times,groups,'binsize',b,'duration',dd,'smooth',sss);
%             FR = length(Base)/sum(base(:,2)-base(:,1));
%             dHPC_dRipples_uncoordinatedB = [dHPC_dRipples_uncoordinatedB , ccg(:,1,2)./length(y)./b./FR];
%             clear tmp FR Base ccg x y times ids groups
%             
%             %Reward
%             y = uncoordinatedR(:,2);
%             x =  Restrict(spks_dHPC(spks_dHPC(:,2) == group_dHPC(i,1),1),NREM);
%             Base = Restrict(x,base);
%             [times,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
%             [ccg,TimeVector2] = CCG(times,groups,'binsize',b,'duration',dd,'smooth',sss);
%             FR = length(Base)/sum(base(:,2)-base(:,1));
%             dHPC_dRipples_uncoordinatedR = [dHPC_dRipples_uncoordinatedR , ccg(:,1,2)./length(y)./b./FR];
%             clear tmp FR Base ccg x y times ids groups      
%             
%             %Aversive
%             y = uncoordinatedA(:,2);
%             x =  Restrict(spks_dHPC(spks_dHPC(:,2) == group_dHPC(i,1),1),NREM);
%             Base = Restrict(x,base);
%             [times,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
%             [ccg,TimeVector2] = CCG(times,groups,'binsize',b,'duration',dd,'smooth',sss);
%             FR = length(Base)/sum(base(:,2)-base(:,1));
%             dHPC_dRipples_uncoordinatedA = [dHPC_dRipples_uncoordinatedA , ccg(:,1,2)./length(y)./b./FR];
%             clear tmp FR Base ccg x y times ids groups            
%         end
%     end
%     clear base i      
%     
    %dHPC SU - coordinated vRipples
    base = InvertIntervals([ripplesD(:,1)-0.1, ripplesD(:,3)+0.1 ; ripplesV(:,1)-0.1, ripplesV(:,3)+0.1],NREM(:,1) , NREM(:,2));
    for i = 1 : length(group_dHPC)% ventral SU
        if logical(tipoD(i))
            %Baseline
            y = coordinatedB_V(:,2);
            x =  Restrict(spks_dHPC(spks_dHPC(:,2) == group_dHPC(i,1),1),NREM);
            Base = Restrict(x,base);
            [times,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,TimeVector2] = CCG(times,groups,'binsize',b,'duration',dd,'smooth',sss);
            FR = length(Base)/sum(base(:,2)-base(:,1));
            dHPC_vRipples_coordinatedB = [dHPC_vRipples_coordinatedB , ccg(:,1,2)./length(y)./b./FR];
            clear tmp FR Base ccg x y times ids groups
            
            %Reward
            y = coordinatedR_V(:,2);
            x =  Restrict(spks_dHPC(spks_dHPC(:,2) == group_dHPC(i,1),1),NREM);
            Base = Restrict(x,base);
            [times,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,TimeVector2] = CCG(times,groups,'binsize',b,'duration',dd,'smooth',sss);
            FR = length(Base)/sum(base(:,2)-base(:,1));
            dHPC_vRipples_coordinatedR = [dHPC_vRipples_coordinatedR , ccg(:,1,2)./length(y)./b./FR];
            clear tmp FR Base ccg x y times ids groups      
            
            %Aversive
            y = coordinatedA_V(:,2);
            x =  Restrict(spks_dHPC(spks_dHPC(:,2) == group_dHPC(i,1),1),NREM);
            Base = Restrict(x,base);
            [times,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,TimeVector2] = CCG(times,groups,'binsize',b,'duration',dd,'smooth',sss);
            FR = length(Base)/sum(base(:,2)-base(:,1));
            dHPC_vRipples_coordinatedA = [dHPC_vRipples_coordinatedA , ccg(:,1,2)./length(y)./b./FR];
            clear tmp FR Base ccg x y times ids groups            
        end
    end
    clear base i      
    
    %dHPC SU - dRipples
    base = InvertIntervals([ripplesD(:,1)-0.1, ripplesD(:,3)+0.1],NREM(:,1) , NREM(:,2));
    for i = 1 : length(group_dHPC)% ventral SU
        if logical(tipoD(i))
            %Baseline
            y = Restrict(ripplesD(:,2),NREM_B);
            x =  Restrict(spks_dHPC(spks_dHPC(:,2) == group_dHPC(i,1),1),NREM);
            Base = Restrict(x,base);
            [times,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,TimeVector2] = CCG(times,groups,'binsize',b,'duration',dd,'smooth',sss);
            FR = length(Base)/sum(base(:,2)-base(:,1));
            dHPC_dRipplesB = [dHPC_dRipplesB , ccg(:,1,2)./length(y)./b./FR];
            clear tmp FR Base ccg x y times ids groups
            
            %Reward
            y = Restrict(ripplesD(:,2),NREM_R);
            x =  Restrict(spks_dHPC(spks_dHPC(:,2) == group_dHPC(i,1),1),NREM);
            Base = Restrict(x,base);
            [times,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,TimeVector2] = CCG(times,groups,'binsize',b,'duration',dd,'smooth',sss);
            FR = length(Base)/sum(base(:,2)-base(:,1));
            dHPC_dRipplesR = [dHPC_dRipplesR , ccg(:,1,2)./length(y)./b./FR];
            clear tmp FR Base ccg x y times ids groups
            
            %Aversive
            y = Restrict(ripplesD(:,2),NREM_A);
            x =  Restrict(spks_dHPC(spks_dHPC(:,2) == group_dHPC(i,1),1),NREM);
            Base = Restrict(x,base);
            [times,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,TimeVector2] = CCG(times,groups,'binsize',b,'duration',dd,'smooth',sss);
            FR = length(Base)/sum(base(:,2)-base(:,1));
            dHPC_dRipplesA = [dHPC_dRipplesA , ccg(:,1,2)./length(y)./b./FR];
            clear tmp FR Base ccg x y times ids groups            
        end
    end
    clear base i      
 
%     %dHPC SU - vRipples
%     ev = [ripplesV(:,1)-0.1, ripplesV(:,3)+0.1 ; ripplesD(:,1)-0.1, ripplesD(:,3)+0.1];
%     base = InvertIntervals(ev,NREM(:,1) , NREM(:,2));
%     clear ev
%     ev = [ripplesV(:,1)-0.1, ripplesV(:,3)+0.1];
%     base1 = InvertIntervals(ev,NREM(:,1) , NREM(:,2));
%     clear ev
%     for i = 1 : length(group_dHPC)% ventral SU
%         if logical(tipoD(i))
%             y = ripplesV(:,2);
%             x =  Restrict(spks_dHPC(spks_dHPC(:,2) == group_dHPC(i,1),1),NREM);
%             Base = Restrict(x,base);
%             [times,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
%             [ccg,TimeVector2] = CCG(times,groups,'binsize',b,'duration',dd,'smooth',sss);
%             FR = length(Base)/sum(base(:,2)-base(:,1));
%             dHPC_vRipples = [dHPC_vRipples , ccg(:,1,2)./length(y)./b./FR];
%             
%             %Poisson
%             totalrippletime = sum(ripplesV(:,3)-ripplesV(:,1));
%             ripplespikes = Restrict(x,[ripplesV(:,1) ripplesV(:,3)]);
%             nripplespikes = size(ripplespikes,1);
%             Base = Restrict(x,base1);
%             ncellbaselinespikes = length(Base);
%             ncellripplespikes = length(ripplespikes);
%             totalbaselinetime = sum(base1(:,2)-base1(:,1));
%             if ncellbaselinespikes~=0 & ncellripplespikes~=0
%                 [pInc pDec surp] = poissonTest(ncellbaselinespikes/totalbaselinetime,ncellripplespikes,totalrippletime);
%                 dHPC_vRipples_poisson = [dHPC_vRipples_poisson ; pInc pDec surp];
%             else
%                 pInc = NaN;
%                 pDec = NaN;
%                 surp = NaN;
%                 dHPC_vRipples_poisson = [dHPC_vRipples_poisson ; pInc pDec surp];
%             end
%             clear tmp FR Base ccg x y times ids groups
%             clear pInc pDec surp totalrippletime ripplespikes nripplespikes
%             clear ncellbaselinespikes ncellripplespikes totalbaselinetime
%         end
%     end
%     clear base i base1   
%     
    % -----------------------------------
    % Ventral SU PHIST
    % -----------------------------------
    %vHPC SU - vRipples All
    base = InvertIntervals([ripplesV(:,1)-0.1, ripplesV(:,3)+0.1],NREM(:,1) , NREM(:,2));
    for i = 1 : length(group_vHPC)% ventral SU
        if logical(tipoV(i))
            y = ripplesV(:,2);
            x =  Restrict(spks_vHPC(spks_vHPC(:,2) == group_vHPC(i,1),1),NREM);
            Base = Restrict(x,base);
            [times,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,TimeVector2] = CCG(times,groups,'binsize',b,'duration',dd,'smooth',sss);
            FR = length(Base)/sum(base(:,2)-base(:,1));
            vHPC_vRipples = [vHPC_vRipples , ccg(:,1,2)./length(y)./b./FR];
            
            
            %Poisson
            totalrippletime = sum(ripplesV(:,3)-ripplesV(:,1));
            ripplespikes = Restrict(x,[ripplesV(:,1) ripplesV(:,3)]);
            nripplespikes = size(ripplespikes,1);
            
            ncellbaselinespikes = length(Base);
            ncellripplespikes = length(ripplespikes);
            totalbaselinetime = sum(base(:,2)-base(:,1));
            if ncellbaselinespikes~=0 & ncellripplespikes~=0
                [pInc pDec surp] = poissonTest(ncellbaselinespikes/totalbaselinetime,ncellripplespikes,totalrippletime);
                vHPC_vRipples_poisson = [vHPC_vRipples_poisson ; pInc pDec surp];
            else
                pInc = NaN;
                pDec = NaN;
                surp = NaN;
                vHPC_vRipples_poisson = [vHPC_vRipples_poisson ; pInc pDec surp];
            end
            clear tmp FR Base ccg x y times ids groups
            clear pInc pDec surp totalrippletime ripplespikes nripplespikes
            clear ncellbaselinespikes ncellripplespikes totalbaselinetime
        end
    end
    clear base i     
    
    %vHPC SU - coordinated vRipples
    base = InvertIntervals([ripplesV(:,1)-0.1, ripplesV(:,3)+0.1],NREM(:,1) , NREM(:,2));
    for i = 1 : length(group_vHPC)% ventral SU
        if logical(tipoV(i))
            %Baseline
            y = coordinatedB_V(:,2);
            x =  Restrict(spks_vHPC(spks_vHPC(:,2) == group_vHPC(i,1),1),NREM);
            Base = Restrict(x,base);
            [times,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,TimeVector2] = CCG(times,groups,'binsize',b,'duration',dd,'smooth',sss);
            FR = length(Base)/sum(base(:,2)-base(:,1));
            vHPC_vRipples_coordinatedB = [vHPC_vRipples_coordinatedB , ccg(:,1,2)./length(y)./b./FR];
            clear tmp FR Base ccg x y times ids groups
            
            %Reward
            y = coordinatedR_V(:,2);
            x =  Restrict(spks_vHPC(spks_vHPC(:,2) == group_vHPC(i,1),1),NREM);
            Base = Restrict(x,base);
            [times,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,TimeVector2] = CCG(times,groups,'binsize',b,'duration',dd,'smooth',sss);
            FR = length(Base)/sum(base(:,2)-base(:,1));
            vHPC_vRipples_coordinatedR = [vHPC_vRipples_coordinatedR , ccg(:,1,2)./length(y)./b./FR];
            clear tmp FR Base ccg x y times ids groups      
            
            %Aversive
            y = coordinatedA_V(:,2);
            x =  Restrict(spks_vHPC(spks_vHPC(:,2) == group_vHPC(i,1),1),NREM);
            Base = Restrict(x,base);
            [times,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,TimeVector2] = CCG(times,groups,'binsize',b,'duration',dd,'smooth',sss);
            FR = length(Base)/sum(base(:,2)-base(:,1));
            vHPC_vRipples_coordinatedA = [vHPC_vRipples_coordinatedA , ccg(:,1,2)./length(y)./b./FR];
            clear tmp FR Base ccg x y times ids groups            
        end
    end
    clear base i      
%     
%     %vHPC SU - uncoordinated vRipples
%     base = InvertIntervals([ripplesV(:,1)-0.1, ripplesV(:,3)+0.1],NREM(:,1) , NREM(:,2));
%     for i = 1 : length(group_vHPC)% ventral SU
%         if logical(tipoV(i))
%             %Baseline
%             y = uncoordinatedB_V(:,2);
%             x =  Restrict(spks_vHPC(spks_vHPC(:,2) == group_vHPC(i,1),1),NREM);
%             Base = Restrict(x,base);
%             [times,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
%             [ccg,TimeVector2] = CCG(times,groups,'binsize',b,'duration',dd,'smooth',sss);
%             FR = length(Base)/sum(base(:,2)-base(:,1));
%             vHPC_vRipples_uncoordinatedB = [vHPC_vRipples_uncoordinatedB , ccg(:,1,2)./length(y)./b./FR];
%             clear tmp FR Base ccg x y times ids groups
%             
%             %Reward
%             y = uncoordinatedR_V(:,2);
%             x =  Restrict(spks_vHPC(spks_vHPC(:,2) == group_vHPC(i,1),1),NREM);
%             Base = Restrict(x,base);
%             [times,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
%             [ccg,TimeVector2] = CCG(times,groups,'binsize',b,'duration',dd,'smooth',sss);
%             FR = length(Base)/sum(base(:,2)-base(:,1));
%             vHPC_vRipples_uncoordinatedR = [vHPC_vRipples_uncoordinatedR , ccg(:,1,2)./length(y)./b./FR];
%             clear tmp FR Base ccg x y times ids groups      
%             
%             %Aversive
%             y = uncoordinatedA_V(:,2);
%             x =  Restrict(spks_vHPC(spks_vHPC(:,2) == group_vHPC(i,1),1),NREM);
%             Base = Restrict(x,base);
%             [times,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
%             [ccg,TimeVector2] = CCG(times,groups,'binsize',b,'duration',dd,'smooth',sss);
%             FR = length(Base)/sum(base(:,2)-base(:,1));
%             vHPC_vRipples_uncoordinatedA = [vHPC_vRipples_uncoordinatedA , ccg(:,1,2)./length(y)./b./FR];
%             clear tmp FR Base ccg x y times ids groups            
%         end
%     end
%     clear base i      
%     
        %vHPC SU - coordinated dRipples
    base = InvertIntervals([ripplesV(:,1)-0.1, ripplesV(:,3)+0.1],NREM(:,1) , NREM(:,2));
    for i = 1 : length(group_vHPC)% ventral SU
        if logical(tipoV(i))
            %Baseline
            y = coordinatedB(:,2);
            x =  Restrict(spks_vHPC(spks_vHPC(:,2) == group_vHPC(i,1),1),NREM);
            Base = Restrict(x,base);
            [times,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,TimeVector2] = CCG(times,groups,'binsize',b,'duration',dd,'smooth',sss);
            FR = length(Base)/sum(base(:,2)-base(:,1));
            vHPC_dRipples_coordinatedB = [vHPC_dRipples_coordinatedB , ccg(:,1,2)./length(y)./b./FR];
            clear tmp FR Base ccg x y times ids groups
            
            %Reward
            y = coordinatedR(:,2);
            x =  Restrict(spks_vHPC(spks_vHPC(:,2) == group_vHPC(i,1),1),NREM);
            Base = Restrict(x,base);
            [times,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,TimeVector2] = CCG(times,groups,'binsize',b,'duration',dd,'smooth',sss);
            FR = length(Base)/sum(base(:,2)-base(:,1));
            vHPC_dRipples_coordinatedR = [vHPC_dRipples_coordinatedR , ccg(:,1,2)./length(y)./b./FR];
            clear tmp FR Base ccg x y times ids groups      
            
            %Aversive
            y = coordinatedA(:,2);
            x =  Restrict(spks_vHPC(spks_vHPC(:,2) == group_vHPC(i,1),1),NREM);
            Base = Restrict(x,base);
            [times,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,TimeVector2] = CCG(times,groups,'binsize',b,'duration',dd,'smooth',sss);
            FR = length(Base)/sum(base(:,2)-base(:,1));
            vHPC_dRipples_coordinatedA = [vHPC_dRipples_coordinatedA , ccg(:,1,2)./length(y)./b./FR];
            clear tmp FR Base ccg x y times ids groups            
        end
    end
    clear base i    
    
    %vHPC SU - vRipples
    base = InvertIntervals([ripplesD(:,1)-0.1, ripplesD(:,3)+0.1 ; ripplesV(:,1)-0.1, ripplesV(:,3)+0.1],NREM(:,1) , NREM(:,2));
    for i = 1 : length(group_vHPC)% ventral SU
        if logical(tipoV(i))
            %Baseline
            y = Restrict(ripplesV(:,2),NREM_B);
            x =  Restrict(spks_vHPC(spks_vHPC(:,2) == group_vHPC(i,1),1),NREM);
            Base = Restrict(x,base);
            [times,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,TimeVector2] = CCG(times,groups,'binsize',b,'duration',dd,'smooth',sss);
            FR = length(Base)/sum(base(:,2)-base(:,1));
            vHPC_vRipplesB = [vHPC_vRipplesB , ccg(:,1,2)./length(y)./b./FR];
            clear tmp FR Base ccg x y times ids groups
            
            %Reward
            y = Restrict(ripplesV(:,2),NREM_R);
            x =  Restrict(spks_vHPC(spks_vHPC(:,2) == group_vHPC(i,1),1),NREM);
            Base = Restrict(x,base);
            [times,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,TimeVector2] = CCG(times,groups,'binsize',b,'duration',dd,'smooth',sss);
            FR = length(Base)/sum(base(:,2)-base(:,1));
            vHPC_vRipplesR = [vHPC_vRipplesR , ccg(:,1,2)./length(y)./b./FR];
            clear tmp FR Base ccg x y times ids groups
            
            %Aversive
            y = Restrict(ripplesV(:,2),NREM_A);
            x =  Restrict(spks_vHPC(spks_vHPC(:,2) == group_vHPC(i,1),1),NREM);
            Base = Restrict(x,base);
            [times,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,TimeVector2] = CCG(times,groups,'binsize',b,'duration',dd,'smooth',sss);
            FR = length(Base)/sum(base(:,2)-base(:,1));
            vHPC_vRipplesA = [vHPC_vRipplesA , ccg(:,1,2)./length(y)./b./FR];
            clear tmp FR Base ccg x y times ids groups            
        end
    end
    clear base i    
%     
%     %vHPC SU - dRipples All
%     ev = [ripplesV(:,1)-0.1, ripplesV(:,3)+0.1 ; ripplesD(:,1)-0.1, ripplesD(:,3)+0.1];
%     base = InvertIntervals(ev,NREM(:,1) , NREM(:,2));
%     clear ev
%     ev = [ripplesD(:,1)-0.1, ripplesD(:,3)+0.1];
%     base1 = InvertIntervals(ev,NREM(:,1) , NREM(:,2));    
%     for i = 1 : length(group_vHPC)% ventral SU
%         if logical(tipoV(i))
%             y = ripplesD(:,2);
%             x =  Restrict(spks_vHPC(spks_vHPC(:,2) == group_vHPC(i,1),1),NREM);
%             Base = Restrict(x,base);
%             [times,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
%             [ccg,TimeVector2] = CCG(times,groups,'binsize',b,'duration',dd,'smooth',sss);
%             FR = length(Base)/sum(base(:,2)-base(:,1));
%             vHPC_dRipples = [vHPC_dRipples , ccg(:,1,2)./length(y)./b./FR];
%             
%             %Poisson
%             totalrippletime = sum(ripplesD(:,3)-ripplesD(:,1));
%             ripplespikes = Restrict(x,[ripplesD(:,1) ripplesD(:,3)]);
%             nripplespikes = size(ripplespikes,1);
%             Base = Restrict(x,base1);
%             ncellbaselinespikes = length(Base);
%             ncellripplespikes = length(ripplespikes);
%             totalbaselinetime = sum(base1(:,2)-base1(:,1));
%             if ncellbaselinespikes~=0 & ncellripplespikes~=0
%                 [pInc pDec surp] = poissonTest(ncellbaselinespikes/totalbaselinetime,ncellripplespikes,totalrippletime);
%                 vHPC_dRipples_poisson = [vHPC_dRipples_poisson ; pInc pDec surp];
%             else
%                 pInc = NaN;
%                 pDec = NaN;
%                 surp = NaN;
%                 vHPC_dRipples_poisson = [vHPC_dRipples_poisson ; pInc pDec surp];
%             end
%             clear tmp FR Base ccg x y times ids groups
%             clear pInc pDec surp totalrippletime ripplespikes nripplespikes
%             clear ncellbaselinespikes ncellripplespikes totalbaselinetime
%         end
%     end
%     clear base i base1   
    
    clear aversiveTS aversiveTS_run rewardTS rewardTS_run baselineTS
    clear camara coordinated coordinatedA coordinatedA_V coordinatedB
    clear coordinatedB_V coordinatedR coordinatedR_V coordinatedV
    clear uncoordinated uncoordinatedA uncoordinatedA_V uncoordinatedB
    clear uncoordinatedB_V uncoordinatedR uncoordinatedR_V uncoordinatedV
    clear Cell_type_classification leftvalve rightvalve Shocks_filt Rewards_filt
    clear NREM NREM_A NREM_B NREM_R REM REM_A REM_B REM_R laps K Kinfo numlaps
    clear ripplesD ripplesV spks spks_dHPC spks_vHPC
    clear TTLs_camara_aversive TTLs_camara_reward WAKE group_dHPC group_vHPC
    clear pos_aversive pos_reward
    
    t
    end
    tt
   
end
        

%% PHIST
% vHPC
ttt = [-2:binsize:2-binsize];
x = [];
xx = [];
criteriaV = [];
for iii = 1 : size(vHPC_Shock,2)
    if and(~isnan(vHPC_Shock(:,iii)) , sum(vHPC_Shock(:,iii))>0)
        x = [ x , zscore(vHPC_Shock(:,iii))];
        xx = [ xx , zscore(vHPC_Valve(:,iii))];
        
        [~ , i] = min(abs(ttt-0));
        [~ , ii] = min(abs(ttt-1));
        [~ , iiii] = min(abs(ttt-(-1)));
        m1 = mean(vHPC_Shock(i:ii,iii));
        m1b = mean(vHPC_Shock(1:iiii,iii));
        stdb1 = std(vHPC_Shock(1:iiii,iii));
        if m1>m1b+(stdb1*3)
            m1 = 1;
        elseif m1<m1b-(stdb1*3)
            m1 = 2;
        else
            m1 = 3;
        end
        
        [~ , i] = min(abs(ttt-(-1)));
        [~ , ii] = min(abs(ttt-1));
        [~ , iiii] = min(abs(ttt-(-1)));
        m2 = mean(vHPC_Valve(i:ii,iii));
        m2b = mean(vHPC_Valve(1:iiii,iii));
        stdb2 = std(vHPC_Valve(1:iiii,iii));
        
        if m2>m2b+(stdb2*3)
            m2 = 1;
        elseif m2<m2b-(stdb2*3)
            m2 = 2;
        else
            m2 = 3;
        end
        
        criteriaV = [criteriaV ; m1 , m2];
        clear m1 m2 m1b stdb1 m1  m2b stdb2 m2
    else
        criteriaV = [criteriaV ; nan , nan];
    end
end

y = mean(x(i:ii,:),1);
[~ , y] = sort(y','descend');

z = [];
zz = [];
for i = 1:size(x,2)
    z = [z , x(:,y(i))];
    zz = [zz , xx(:,y(i))];
end

figure
subplot(121), imagesc(ttt, [1:1:size(x,2)], z'), axis tight, xline(0,'LineWidth',1),xlim([-2 2]),caxis([0 4]),xlabel('Time(sec)'),ylabel('Neurons(id)')
xline(1,'--','LineWidth',1)
title('vHPC Shock')

subplot(122), imagesc(ttt, [1:1:size(xx,2)], zz'), axis tight, xline(0,'LineWidth',1),xlim([-2 2]),caxis([0 4]),xlabel('Time(sec)'),ylabel('Neurons(id)')
xline(1,'--','LineWidth',1)
title('vHPC Valve')

% dHPC
ttt = [-2:binsize:2-binsize];
[~ , i] = min(abs(ttt-0));
[~ , ii] = min(abs(ttt-1));

x = [];
xx = [];
criteriaD = [];
for iii = 1 : size(dHPC_Shock,2)
    if and(~isnan(dHPC_Shock(:,iii)) , sum(dHPC_Shock(:,iii))>0)
        x = [ x , dHPC_Shock(:,iii)];
        xx = [ xx , dHPC_Valve(:,iii)];

        [~ , i] = min(abs(ttt-0));
        [~ , ii] = min(abs(ttt-1));
        [~ , iiii] = min(abs(ttt-(-1)));
        m1 = mean(dHPC_Shock(i:ii,iii));
        m1b = mean(dHPC_Shock(1:iiii,iii));
        stdb1 = std(dHPC_Shock(1:iiii,iii));
        
        if m1>m1b+(stdb1*3)
            m1 = 1;
        elseif m1<m1b-(stdb1*3)
            m1 = 2;
        else
            m1 = 3;
        end
        
        [~ , i] = min(abs(ttt-(-1)));
        [~ , ii] = min(abs(ttt-1));
        [~ , iiii] = min(abs(ttt-(-1)));
        m2 = mean(dHPC_Valve(i:ii,iii));
        m2b = mean(dHPC_Valve(1:iiii,iii));
        stdb2 = std(dHPC_Valve(1:iiii,iii));
        
        if m2>m2b+(stdb2*3)
            m2 = 1;
        elseif m2<m2b-(stdb2*3)
            m2 = 2;
        else
            m2 = 3;
        end
        
        criteriaD = [criteriaD ; m1 , m2];
        clear m1 m2 m1b stdb1 m1  m2b stdb2 m2
    else
        criteriaD = [criteriaD ; nan , nan];
    end

end

y = mean(x(i:ii,:),1);
[~ , y] = sort(y','descend');

z = [];
zz = [];
for i = 1:size(x,2)
    z = [z , x(:,y(i))];
    zz = [zz , xx(:,y(i))];
end

figure
subplot(121), imagesc(ttt, [1:1:size(x,2)], z'), axis tight, xline(0,'LineWidth',1),xlim([-2 2]),caxis([-0 10]),xlabel('Time(sec)'),ylabel('Neurons(id)')
xline(1,'--','LineWidth',1)
title('dHPC Shock')

subplot(122), imagesc(ttt, [1:1:size(xx,2)], zz'), axis tight, xline(0,'LineWidth',1),xlim([-2 2]),caxis([-0 10]),xlabel('Time(sec)'),ylabel('Neurons(id)')
xline(1,'--','LineWidth',1)
title('dHPC Valve')


figure,
subplot(121),plot([-2:1/30:2-1/30],mean(velocity_shock,2)),ylim([4 30])
subplot(122),plot([-2:1/30:2-1/30],mean(velocity_valve,2)),ylim([4 30])

total = sum(not(isnan(criteriaD(:,1))));
percentage.dorsal.shock.up = sum(and(criteriaD(:,1)==1,criteriaD(:,2)==3))/total*100;
percentage.dorsal.shock.down = sum(and(criteriaD(:,1)==2,criteriaD(:,2)==3))/total*100;
percentage.dorsal.valve.up = sum(and(criteriaD(:,1)==3,criteriaD(:,2)==1))/total*100;
percentage.dorsal.valve.down = sum(and(criteriaD(:,1)==3,criteriaD(:,2)==2))/total*100;
percentage.dorsal.bidirectional.shock = sum(and(criteriaD(:,1)==1,criteriaD(:,2)==2))/total*100;
percentage.dorsal.bidirectional.valve = sum(and(criteriaD(:,1)==2,criteriaD(:,2)==1))/total*100;
percentage.dorsal.bidirectional.both.active = sum(and(criteriaD(:,1)==1,criteriaD(:,2)==1))/total*100;
percentage.dorsal.bidirectional.both.inactive = sum(and(criteriaD(:,1)==2,criteriaD(:,2)==2))/total*100;
percentage.dorsal.no_responsive = sum(and(criteriaD(:,1)==3,criteriaD(:,2)==3))/total*100;

a = sum([b c d])
b = sum([percentage.dorsal.shock.up percentage.dorsal.shock.down])
c = sum([percentage.dorsal.valve.up percentage.dorsal.valve.down])
d = sum([percentage.dorsal.bidirectional.shock percentage.dorsal.bidirectional.valve percentage.dorsal.bidirectional.both.active percentage.dorsal.bidirectional.both.inactive])
e = sum([percentage.dorsal.no_responsive])

total = sum(not(isnan(criteriaV(:,1))));
percentage.ventral.shock.up = sum(and(criteriaV(:,1)==1,criteriaV(:,2)==3))/total*100;
percentage.ventral.shock.down = sum(and(criteriaV(:,1)==2,criteriaV(:,2)==3))/total*100;
percentage.ventral.valve.up = sum(and(criteriaV(:,1)==3,criteriaV(:,2)==1))/total*100;
percentage.ventral.valve.down = sum(and(criteriaV(:,1)==3,criteriaV(:,2)==2))/total*100;
percentage.ventral.bidirectional.shock = sum(and(criteriaV(:,1)==1,criteriaV(:,2)==2))/total*100;
percentage.ventral.bidirectional.valve = sum(and(criteriaV(:,1)==2,criteriaV(:,2)==1))/total*100;
percentage.ventral.bidirectional.both.active = sum(and(criteriaV(:,1)==1,criteriaV(:,2)==1))/total*100;
percentage.ventral.bidirectional.both.inactive = sum(and(criteriaV(:,1)==2,criteriaV(:,2)==2))/total*100;
percentage.ventral.no_responsive = sum(and(criteriaV(:,1)==3,criteriaV(:,2)==3))/total*100;

a = sum([b c d])
b = sum([percentage.ventral.shock.up percentage.ventral.shock.down])
c = sum([percentage.ventral.valve.up percentage.ventral.valve.down])
d = sum([percentage.ventral.bidirectional.shock percentage.ventral.bidirectional.valve percentage.ventral.bidirectional.both.active percentage.ventral.bidirectional.both.inactive])
e = sum([percentage.ventral.no_responsive])

%% Plot dorsal and ventral responsive cells to the shock
% Shock
criterio = or(and(criteriaV(:,1)==1 , criteriaV(:,2)==3),and(criteriaV(:,1)==1 , criteriaV(:,2)==2));
y = vHPC_Shock(:,criterio');
x = vHPC_Shock(:,and(criteriaV(:,1)==3 , criteriaV(:,2)==3)');

m= mean(y(2:end,:),2);
figure,
subplot(222),plot(ttt,m),hold on
sem = std(y(2:end,:),0,2)/sqrt(size(y,2));
ciplot(m-sem , m+sem,ttt),alpha 0.05

m= mean(x(2:end,:),2);
plot(ttt,m)
sem = std(x(2:end,:),0,2)/sqrt(size(x,2));
ciplot(m-sem , m+sem,ttt),alpha 0.05

criterio = or(and(criteriaD(:,1)==1 , criteriaD(:,2)==3),and(criteriaD(:,1)==1 , criteriaD(:,2)==2));
y = dHPC_Shock(:,criterio');
x = dHPC_Shock(:,and(criteriaD(:,1)==3 , criteriaD(:,2)==3)');

m= mean(y(2:end,:),2);
subplot(221),plot(ttt,m),hold on
sem = std(y(2:end,:),0,2)/sqrt(size(y,2));
ciplot(m-sem , m+sem,ttt),alpha 0.05

m= mean(x(2:end,:),2);
plot(ttt,m)
sem = std(x(2:end,:),0,2)/sqrt(size(x,2));
ciplot(m-sem , m+sem,ttt),alpha 0.05


%Valve
criterio = or(and(criteriaD(:,1)==3 , criteriaD(:,2)==1),and(criteriaD(:,1)==2 , criteriaD(:,2)==1));
y = dHPC_Shock(:,criterio');
x = dHPC_Shock(:,and(criteriaD(:,1)==3 , criteriaD(:,2)==3)');

m= mean(y(2:end,:),2);
subplot(223),plot(ttt,m),hold on
sem = std(y(2:end,:),0,2)/sqrt(size(y,2));
ciplot(m-sem , m+sem,ttt),alpha 0.05

m= mean(x(2:end,:),2);
plot(ttt,m)
sem = std(x(2:end,:),0,2)/sqrt(size(x,2));
ciplot(m-sem , m+sem,ttt),alpha 0.05


criterio = or(and(criteriaV(:,1)==3 , criteriaV(:,2)==1),and(criteriaV(:,1)==2 , criteriaV(:,2)==1));
y = vHPC_Shock(:,criterio');
x = vHPC_Shock(:,and(criteriaV(:,1)==3 , criteriaV(:,2)==3)');

m= mean(y(2:end,:),2);
subplot(224),plot(ttt,m),hold on
sem = std(y(2:end,:),0,2)/sqrt(size(y,2));
ciplot(m-sem , m+sem,ttt),alpha 0.05

m= mean(x(2:end,:),2);
plot(ttt,m)
sem = std(x(2:end,:),0,2)/sqrt(size(x,2));
ciplot(m-sem , m+sem,ttt),alpha 0.05

%% PHIST coordinated ripples-SU
% vHPC
[~ , i] = min(abs(TimeVector2-(-0.2)));
[~ , ii] = min(abs(TimeVector2-0.2));

x = [];
xx = [];
xxx = [];
value = [];
value1 = [];
for iii = 1 : size(vHPC_vRipplesA,2)
    if or(and(criteriaV(iii,1)==1 , criteriaV(iii,2)==3), and(criteriaV(iii,1)==1 , criteriaV(iii,2)==2))
        [~ , i] = min(abs(TimeVector2-(-0.2)));
        [~ , ii] = min(abs(TimeVector2-0.2));
        x = [x , vHPC_vRipples_coordinatedB(:,iii)];
        xx = [xx , vHPC_vRipples_coordinatedR(:,iii)];
        xxx = [xxx , vHPC_vRipples_coordinatedA(:,iii)];
        
        value = [value ; mean( vHPC_vRipples_coordinatedB(i:ii,iii)) mean( vHPC_vRipples_coordinatedR(i:ii,iii)) mean( vHPC_vRipples_coordinatedA(i:ii,iii))];
        value1 = [value1 ; mean( vHPC_vRipples_coordinatedR(i:ii,iii))/mean( vHPC_vRipples_coordinatedB(i:ii,iii)) mean( vHPC_vRipples_coordinatedA(i:ii,iii))/mean( vHPC_vRipples_coordinatedB(i:ii,iii))];

    end
end

y = mean(xxx(i:ii,:),1);
[~ , y] = sort(y','descend');

z = [];
zz = [];
zzz = [];
for i = 1:size(x,2)
    z = [z , x(:,y(i))];
    zz = [zz , xx(:,y(i))];
    zzz = [zzz , xxx(:,y(i))];
end

figure
subplot(131), imagesc(TimeVector2, [1:1:size(x,2)], z'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default'),caxis([0 4])
title('Baseline')
subplot(132), imagesc(TimeVector2, [1:1:size(xx,2)], zz'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default'),caxis([0 4])
title('Reward')
subplot(133), imagesc(TimeVector2, [1:1:size(xxx,2)], zzz'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default'),caxis([0 4])
title('Aversive')

figure,
plot(TimeVector2,mean(z,2),'k'),hold on
plot(TimeVector2,mean(zz,2),'b'),hold on
plot(TimeVector2,mean(zzz,2),'r'),hold on

clear x xx xxx z zz zzz


% dHPC
[~ , i] = min(abs(TimeVector2-(-0.2)));
[~ , ii] = min(abs(TimeVector2-0.2));

x = [];
xx = [];
xxx = [];

for iii = 1 : size(dHPC_dRipplesA,2)
    if or(and(criteriaD(iii,1)==1 , criteriaD(iii,2)==3), and(criteriaD(iii,1)==1 , criteriaD(iii,2)==2))
        [~ , i] = min(abs(TimeVector2-(-0.05)));
        [~ , ii] = min(abs(TimeVector2-0.05));
        x = [x , dHPC_dRipples_coordinatedB(:,iii)];
        xx = [xx , dHPC_dRipples_coordinatedR(:,iii)];
        xxx = [xxx , dHPC_dRipples_coordinatedA(:,iii)];
        
        value = [value ; mean(dHPC_dRipples_coordinatedB(i:ii,iii)) mean(dHPC_dRipples_coordinatedR(i:ii,iii)) mean(dHPC_dRipples_coordinatedA(i:ii,iii))];
    end
end

y = mean(x(i:ii,:),1);
[~ , y] = sort(y','descend');

z = [];
zz = [];
zzz = [];
for i = 1:size(x,2)

    z = [z , x(:,y(i))];
    clear m mm mmm
    
    zz = [zz , xx(:,y(i))];
    clear m mm mmm
    
    zzz = [zzz , xxx(:,y(i))];
    clear m mm mmm
end


figure
subplot(131), imagesc(TimeVector2, [1:1:size(x,2)], z'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default'),caxis([0 10])
title('Baseline')
subplot(132), imagesc(TimeVector2, [1:1:size(xx,2)], zz'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default'),caxis([0 10])
title('Reward')
subplot(133), imagesc(TimeVector2, [1:1:size(xxx,2)], zzz'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default'),caxis([0 10])
title('Aversive')


figure,
plot(TimeVector2,mean(z,2),'k'),hold on
plot(TimeVector2,mean(zz,2),'b'),hold on
plot(TimeVector2,mean(zzz,2),'r'),hold on


%% PHIST ripples-SU
% vHPC
[~ , i] = min(abs(TimeVector2-(-0.05)));
[~ , ii] = min(abs(TimeVector2-0.05));

x = [];
xx = [];
xxx = [];
tmp1 = [];
for iii = 1 : size(vHPC_vRipples,2)
%     if and(criteriaV(iii,1)>1 , criteriaV(iii,2)<=1)
            x = [x , vHPC_dRipples(:,iii)];
            tmp1 = [tmp1 ; mean(vHPC_vRipples(i:ii,iii))>2];
%     end
end

y = mean(x(i:ii,:),1);
[~ , y] = sort(y','descend');

z = [];

for i = 1:size(x,2)
    z = [z , x(:,y(i))];
end

figure
subplot(121), imagesc(TimeVector2, [1:1:size(x,2)], z'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default'),caxis([0 6])
title('vHPC SU - vRipples')

% figure,
% plot(TimeVector2,mean(z,2),'k'),hold on

clear x xx xxx z zz zzz


% dHPC
[~ , i] = min(abs(TimeVector2-(-0.05)));
[~ , ii] = min(abs(TimeVector2-0.05));

x = [];

for iii = 1 : size(dHPC_dRipples,2)
%     if and(criteriaD(iii,1)>1 , criteriaD(iii,2)<=1)
            x = [x , dHPC_dRipples(:,iii)];
%     end
end

y = mean(x(i:ii,:),1);
[~ , y] = sort(y','descend');

z = [];

for i = 1:size(x,2)
    z = [z , x(:,y(i))];
end

subplot(122), imagesc(TimeVector2, [1:1:size(x,2)], z'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default'),caxis([0 6])
title('dHPC SU - dRipples')

%% Ripple Modulated neurons
up_modulated_vHPC = vHPC_vRipples_poisson(:,1)<0.01;
down_modulated_vHPC = vHPC_vRipples_poisson(:,2)<0.01;
no_modulated_vHPC = and(vHPC_vRipples_poisson(:,1)>0.01 , vHPC_vRipples_poisson(:,2)>0.01);
modulated = sum([sum(up_modulated_vHPC),sum(down_modulated_vHPC)]); 
total_vHPC = sum([sum(up_modulated_vHPC), sum(down_modulated_vHPC) , sum(no_modulated_vHPC)]);
percentage_vHPC =[sum(no_modulated_vHPC)/total_vHPC , sum([sum(up_modulated_vHPC),sum(down_modulated_vHPC)])/total_vHPC , sum(down_modulated_vHPC)/modulated , sum(up_modulated_vHPC)/modulated]*100;


up_modulated_dHPC = dHPC_dRipples_poisson(:,1)<0.01;
down_modulated_dHPC = dHPC_dRipples_poisson(:,2)<0.01;
no_modulated_dHPC = and(dHPC_dRipples_poisson(:,1)>0.01 , dHPC_dRipples_poisson(:,2)>0.01);
total_dHPC = sum([sum(up_modulated_dHPC), sum(down_modulated_dHPC) , sum(no_modulated_dHPC)]);
modulated = sum([sum(up_modulated_dHPC),sum(down_modulated_dHPC)]); 
percentage_dHPC =[sum(no_modulated_dHPC)/total_dHPC , sum([sum(up_modulated_dHPC),sum(down_modulated_dHPC)])/total_dHPC , sum(down_modulated_dHPC)/modulated , sum(up_modulated_dHPC)/modulated]*100;

[~ , i] = min(abs(TimeVector2-(-0.2)));
[~ , ii] = min(abs(TimeVector2-0.2));

%mean of up-modulated neurons
mean_dHPC = mean(dHPC_dRipples(i:ii,up_modulated_dHPC'));
mean_vHPC = mean(vHPC_vRipples(i:ii,up_modulated_vHPC'));
%quantile of upmodulated neurons
q_dHPC = quantile(mean_dHPC,3);
q_vHPC = quantile(mean_vHPC,3);

% vHPC
[~ , i] = min(abs(TimeVector2-(-0.2)));
[~ , ii] = min(abs(TimeVector2-0.2));

%mean of up-modulated neurons
mean_dHPC = mean(dHPC_dRipples(i:ii,:));
mean_vHPC = mean(vHPC_vRipples(i:ii,:));

x = [];
xx = [];
xxx = [];
value = [];
for iii = 1 : size(vHPC_vRipplesA,2)
%     if and(up_modulated_vHPC(iii),mean_vHPC(iii)>=q_vHPC(3))
    if up_modulated_vHPC(iii)
%         if and(criteriaD(iii,1)==1 , criteriaD(iii,2)==3)
            x = [x , vHPC_vRipples_coordinatedB(:,iii)];
            xx = [xx , vHPC_vRipples_coordinatedR(:,iii)];
            xxx = [xxx , vHPC_vRipples_coordinatedA(:,iii)];
            value = [value ; max(vHPC_vRipples_coordinatedB(i:ii,iii))  max(vHPC_vRipples_coordinatedR(i:ii,iii))  max(vHPC_vRipples_coordinatedA(i:ii,iii))];
%         end
    end
end

y = mean(xxx(i:ii,:),1);
[~ , y] = sort(y','descend');

z = [];
zz = [];
zzz = [];
for i = 1:size(x,2)
    z = [z , x(:,y(i))];
    zz = [zz , xx(:,y(i))];
    zzz = [zzz , xxx(:,y(i))];

end

figure
subplot(131), imagesc(TimeVector2, [1:1:size(x,2)], z'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default'),caxis([0 2])
title('Baseline')
subplot(132), imagesc(TimeVector2, [1:1:size(xx,2)], zz'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default'),caxis([0 2])
title('Reward')
subplot(133), imagesc(TimeVector2, [1:1:size(xxx,2)], zzz'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default'),caxis([0 2])
title('Aversive')

figure,
plot(TimeVector2,mean(z,2),'k'),hold on
SEM = std(z,0,2)/sqrt(size(z,2));
ciplot((mean(z,2)-SEM),(mean(z,2)+SEM),TimeVector2,'k')
alpha 0.3
plot(TimeVector2,mean(zz,2),'b'),hold on
SEM = std(zz,0,2)/sqrt(size(zz,2));
ciplot((mean(zz,2)-SEM),(mean(zz,2)+SEM),TimeVector2,'b')
alpha 0.3
plot(TimeVector2,mean(zzz,2),'r'),hold on
SEM = std(zzz,0,2)/sqrt(size(zzz,2));
ciplot((mean(zzz,2)-SEM),(mean(zzz,2)+SEM),TimeVector2,'r')
alpha 0.3
clear x xx xxx z zz zzz


% dHPC
[~ , i] = min(abs(TimeVector2-(-0.2)));
[~ , ii] = min(abs(TimeVector2-0.2));

x = [];
xx = [];
xxx = [];
value=[];
for iii = 1 : size(dHPC_dRipplesA,2)
%     if and(up_modulated_dHPC(iii),mean_dHPC(iii)>=q_dHPC(3))
    if up_modulated_dHPC(iii)
            x = [x , dHPC_dRipples_coordinatedB(:,iii)];
            xx = [xx , dHPC_dRipples_coordinatedR(:,iii)];
            xxx = [xxx , dHPC_dRipples_coordinatedA(:,iii)];
            
            value = [value ; max(dHPC_dRipples_coordinatedB(i:ii,iii))  max(dHPC_dRipples_coordinatedR(i:ii,iii))  max(dHPC_dRipples_coordinatedA(i:ii,iii))];
    end
end

y = mean(x(i:ii,:),1);
[~ , y] = sort(y','descend');

z = [];
zz = [];
zzz = [];
for i = 1:size(x,2)
    z = [z , x(:,y(i))];
    zz = [zz , xx(:,y(i))];
    zzz = [zzz , xxx(:,y(i))];
end


figure
subplot(131), imagesc(TimeVector2, [1:1:size(x,2)], z'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default'),caxis([0 4])
title('Baseline')
subplot(132), imagesc(TimeVector2, [1:1:size(xx,2)], zz'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default'),caxis([0 4])
title('Reward')
subplot(133), imagesc(TimeVector2, [1:1:size(xxx,2)], zzz'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default'),caxis([0 4])
title('Aversive')


figure,
plot(TimeVector2,mean(z,2),'k'),hold on
SEM = std(z,0,2)/sqrt(size(z,2));
ciplot((mean(z,2)-SEM),(mean(z,2)+SEM),TimeVector2,'k')
alpha 0.3
plot(TimeVector2,mean(zz,2),'b'),hold on
SEM = std(zz,0,2)/sqrt(size(zz,2));
ciplot((mean(zz,2)-SEM),(mean(zz,2)+SEM),TimeVector2,'b')
alpha 0.3
plot(TimeVector2,mean(zzz,2),'r'),hold on
SEM = std(zzz,0,2)/sqrt(size(zzz,2));
ciplot((mean(zzz,2)-SEM),(mean(zzz,2)+SEM),TimeVector2,'r')
alpha 0.3

%% Ripple cross-modulated neurons
up_modulated_vHPC = vHPC_dRipples_poisson(:,1)<0.001;
down_modulated_vHPC = vHPC_dRipples_poisson(:,2)<0.001;
no_modulated_vHPC = and(vHPC_dRipples_poisson(:,1)>0.001 , vHPC_dRipples_poisson(:,2)>0.01);
total_vHPC = sum([sum(up_modulated_vHPC), sum(down_modulated_vHPC) , sum(no_modulated_vHPC)]);
percentage_vHPC =[sum(no_modulated_vHPC)/total_vHPC , sum([sum(up_modulated_vHPC),sum(down_modulated_vHPC)])/total_vHPC , sum(down_modulated_vHPC)/total_vHPC , sum(up_modulated_vHPC)/total_vHPC]*100;


up_modulated_dHPC = dHPC_vRipples_poisson(:,1)<0.001;
down_modulated_dHPC = dHPC_vRipples_poisson(:,2)<0.001;
no_modulated_dHPC = and(dHPC_vRipples_poisson(:,1)>0.001 , dHPC_vRipples_poisson(:,2)>0.01);
total_dHPC = sum([sum(up_modulated_dHPC), sum(down_modulated_dHPC) , sum(no_modulated_dHPC)]);
percentage_dHPC =[sum(no_modulated_dHPC)/total_vHPC , sum([sum(up_modulated_dHPC),sum(down_modulated_dHPC)])/total_vHPC , sum(down_modulated_dHPC)/total_vHPC , sum(up_modulated_dHPC)/total_vHPC]*100;

[~ , i] = min(abs(TimeVector2-(-0.05)));
[~ , ii] = min(abs(TimeVector2-0.05));

%mean of up-modulated neurons
mean_dHPC = mean(dHPC_vRipples(i:ii,up_modulated_dHPC'));
mean_vHPC = mean(vHPC_dRipples(i:ii,up_modulated_vHPC'));
%quantile of upmodulated neurons
q_dHPC = quantile(mean_dHPC,3);
q_vHPC = quantile(mean_vHPC,3);

% vHPC
[~ , i] = min(abs(TimeVector2-(-0.05)));
[~ , ii] = min(abs(TimeVector2-0.05));

%mean of up-modulated neurons
mean_dHPC = mean(dHPC_vRipples(i:ii,:));
mean_vHPC = mean(vHPC_dRipples(i:ii,:));

x = [];
xx = [];
xxx = [];

for iii = 1 : size(vHPC_dRipples,2)
    if up_modulated_vHPC(iii)%,mean_vHPC(iii)>=q_vHPC(3))
            x = [x , vHPC_vRipples_coordinatedB(:,iii)];
            xx = [xx , vHPC_vRipples_coordinatedR(:,iii)];
            xxx = [xxx , vHPC_vRipples_coordinatedA(:,iii)];
    end
end

y = mean(x(i:ii,:),1);
[~ , y] = sort(y','descend');

z = [];
zz = [];
zzz = [];
for i = 1:size(x,2)
    z = [z , x(:,y(i))];
    zz = [zz , xx(:,y(i))];
    zzz = [zzz , xxx(:,y(i))];
end

figure
subplot(131), imagesc(TimeVector2, [1:1:size(x,2)], z'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default'),caxis([0 10])
title('Baseline')
subplot(132), imagesc(TimeVector2, [1:1:size(xx,2)], zz'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default'),caxis([0 10])
title('Reward')
subplot(133), imagesc(TimeVector2, [1:1:size(xxx,2)], zzz'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default'),caxis([0 10])
title('Aversive')

figure,
plot(TimeVector2,mean(z,2),'k'),hold on
SEM = std(z,0,2)/sqrt(size(z,2));
ciplot((mean(z,2)-SEM),(mean(z,2)+SEM),TimeVector2,'k')
alpha 0.3
plot(TimeVector2,mean(zz,2),'b'),hold on
SEM = std(zz,0,2)/sqrt(size(zz,2));
ciplot((mean(zz,2)-SEM),(mean(zz,2)+SEM),TimeVector2,'b')
alpha 0.3
plot(TimeVector2,mean(zzz,2),'r'),hold on
SEM = std(zzz,0,2)/sqrt(size(zzz,2));
ciplot((mean(zzz,2)-SEM),(mean(zzz,2)+SEM),TimeVector2,'r')
alpha 0.3
clear x xx xxx z zz zzz


% dHPC
[~ , i] = min(abs(TimeVector2-(-0.05)));
[~ , ii] = min(abs(TimeVector2-0.05));

x = [];
xx = [];
xxx = [];
for iii = 1 : size(dHPC_dRipplesA,2)
    if up_modulated_dHPC(iii)%,mean_dHPC(iii)>=q_dHPC(3))
            x = [x , dHPC_vRipples_coordinatedB(:,iii)];
            xx = [xx , dHPC_vRipples_coordinatedR(:,iii)];
            xxx = [xxx , dHPC_vRipples_coordinatedA(:,iii)];
    end
end

y = mean(xxx(i:ii,:),1);
[~ , y] = sort(y','descend');

z = [];
zz = [];
zzz = [];
for i = 1:size(x,2)
    z = [z , x(:,y(i))];
    zz = [zz , xx(:,y(i))];
    zzz = [zzz , xxx(:,y(i))];
end


figure
subplot(131), imagesc(TimeVector2, [1:1:size(x,2)], z'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default'),caxis([0 12])
title('Baseline')
subplot(132), imagesc(TimeVector2, [1:1:size(xx,2)], zz'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default'),caxis([0 12])
title('Reward')
subplot(133), imagesc(TimeVector2, [1:1:size(xxx,2)], zzz'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default'),caxis([0 12])
title('Aversive')


figure,
plot(TimeVector2,mean(z,2),'k'),hold on
SEM = std(z,0,2)/sqrt(size(z,2));
ciplot((mean(z,2)-SEM),(mean(z,2)+SEM),TimeVector2,'k')
alpha 0.3
plot(TimeVector2,mean(zz,2),'b'),hold on
SEM = std(zz,0,2)/sqrt(size(zz,2));
ciplot((mean(zz,2)-SEM),(mean(zz,2)+SEM),TimeVector2,'b')
alpha 0.3
plot(TimeVector2,mean(zzz,2),'r'),hold on
SEM = std(zzz,0,2)/sqrt(size(zzz,2));
ciplot((mean(zzz,2)-SEM),(mean(zzz,2)+SEM),TimeVector2,'r')
alpha 0.3


% Uncoordinated ripples
%vHPC
x = [];
xx = [];
xxx = [];

for iii = 1 : size(vHPC_vRipplesA,2)
    if and(up_modulated_vHPC(iii),mean_vHPC(iii)>=q_vHPC(3))
            x = [x , vHPC_vRipples_uncoordinatedB(:,iii)];
            xx = [xx , vHPC_vRipples_uncoordinatedR(:,iii)];
            xxx = [xxx , vHPC_vRipples_uncoordinatedA(:,iii)];
    end
end

y = mean(xxx(i:ii,:),1);
[~ , y] = sort(y','descend');

z = [];
zz = [];
zzz = [];
for i = 1:size(x,2)
    z = [z , x(:,y(i))];
    zz = [zz , xx(:,y(i))];
    zzz = [zzz , xxx(:,y(i))];
end

figure
subplot(131), imagesc(TimeVector2, [1:1:size(x,2)], z'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default'),caxis([0 10])
title('Baseline')
subplot(132), imagesc(TimeVector2, [1:1:size(xx,2)], zz'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default'),caxis([0 10])
title('Reward')
subplot(133), imagesc(TimeVector2, [1:1:size(xxx,2)], zzz'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default'),caxis([0 10])
title('Aversive')

figure,
plot(TimeVector2,mean(z,2),'k'),hold on
SEM = std(z,0,2)/sqrt(size(z,2));
ciplot((mean(z,2)-SEM),(mean(z,2)+SEM),TimeVector2,'k')
alpha 0.3
plot(TimeVector2,mean(zz,2),'b'),hold on
SEM = std(zz,0,2)/sqrt(size(zz,2));
ciplot((mean(zz,2)-SEM),(mean(zz,2)+SEM),TimeVector2,'b')
alpha 0.3
plot(TimeVector2,mean(zzz,2),'r'),hold on
SEM = std(zzz,0,2)/sqrt(size(zzz,2));
ciplot((mean(zzz,2)-SEM),(mean(zzz,2)+SEM),TimeVector2,'r')
alpha 0.3
clear x xx xxx z zz zzz


% dHPC
[~ , i] = min(abs(TimeVector2-(-0.05)));
[~ , ii] = min(abs(TimeVector2-0.05));

x = [];
xx = [];
xxx = [];
for iii = 1 : size(dHPC_dRipplesA,2)
    if and(up_modulated_dHPC(iii),mean_dHPC(iii)>=q_dHPC(3))
            x = [x , dHPC_dRipples_uncoordinatedB(:,iii)];
            xx = [xx , dHPC_dRipples_uncoordinatedR(:,iii)];
            xxx = [xxx , dHPC_dRipples_uncoordinatedA(:,iii)];
    end
end

y = mean(xxx(i:ii,:),1);
[~ , y] = sort(y','descend');

z = [];
zz = [];
zzz = [];
for i = 1:size(x,2)
    z = [z , x(:,y(i))];
    zz = [zz , xx(:,y(i))];
    zzz = [zzz , xxx(:,y(i))];
end


figure
subplot(131), imagesc(TimeVector2, [1:1:size(x,2)], z'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default'),caxis([0 12])
title('Baseline')
subplot(132), imagesc(TimeVector2, [1:1:size(xx,2)], zz'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default'),caxis([0 12])
title('Reward')
subplot(133), imagesc(TimeVector2, [1:1:size(xxx,2)], zzz'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default'),caxis([0 12])
title('Aversive')


figure,
plot(TimeVector2,mean(z,2),'k'),hold on
SEM = std(z,0,2)/sqrt(size(z,2));
ciplot((mean(z,2)-SEM),(mean(z,2)+SEM),TimeVector2,'k')
alpha 0.3
plot(TimeVector2,mean(zz,2),'b'),hold on
SEM = std(zz,0,2)/sqrt(size(zz,2));
ciplot((mean(zz,2)-SEM),(mean(zz,2)+SEM),TimeVector2,'b')
alpha 0.3
plot(TimeVector2,mean(zzz,2),'r'),hold on
SEM = std(zzz,0,2)/sqrt(size(zzz,2));
ciplot((mean(zzz,2)-SEM),(mean(zzz,2)+SEM),TimeVector2,'r')
alpha 0.3
