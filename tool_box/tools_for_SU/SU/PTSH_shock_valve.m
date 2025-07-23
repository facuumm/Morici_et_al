clear
clc
close all

%Parameters

%%
% $$\$$

path = {'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready'};
% path = 'E:\Rat128\Ephys\in_pyr\ready';

%% Parameters
% ---------------
% --> SU
% ---------------
cpval = 0.001/2; % p value to define if SU are ripple modulated
binsize = 0.05; % bin size for SU histogram in sec
ss = 1; %smooth level of CCG
d = 4; %duration of the CCG
n_SU_V = 0;
n_SU_D = 0;
FR_B_V = []; FR_R_V = [];  FR_A_V = []; % Firing Rate during NREM
FR_B_D = []; FR_R_D = [];  FR_A_D = []; % Firing Rate during NREM
poisson_dHPC_split = []; poisson_vHPC_split = []; %poisson results split by conditions

%for storing Peri-event Histogram
vHPC_neurons_shock_zscore = [];
vHPC_neurons_shock = [];
dHPC_neurons_shock_zscore = [];
dHPC_neurons_shock = [];

vHPC_neurons_valve_zscore = [];
vHPC_neurons_valve = [];
dHPC_neurons_valve_zscore = [];
dHPC_neurons_valve = [];

cutting = 63; % final channel from vHPC

velocity_shock = [];%To store speed data locked to the events
velocity_valve = [];

% ---------------
% --> Ripples
% ---------------
thresholdsD = [2 5]; % thresholds used for Ripple detection [1.5 5]
thresholdsV = [2 5]; % thresholds used for Ripple detection [1.5 5]

durations = [30 20 100]; % durations criteria for Ripple detection [Max , Min , Min_interval]
frequencies = [100 200]; % frequency band for ripple detection [100 250]

q = 0.25; %quantile to restrict above it ripples according to their peak amplitude

channels = [74 9 1 23 18 43];
channels = [channels ; 117 14 6 7 61 53];

ripples_coordinated_percentage = [];
rateV = []; rateD = []; % to store the ripples rate from dHPC and vHPC
deltaB = []; deltaR = []; deltaA = []; %to store all the time delta beween dorsal-ventral ripples

rate_in_time_V = [];
rate_in_time_D = [];

durationsB_dHPC = []; durationsR_dHPC = []; durationsA_dHPC = [];
amplitudesB_dHPC = []; amplitudesR_dHPC = []; amplitudesA_dHPC = [];


durationsB_vHPC = []; durationsR_vHPC = []; durationsA_vHPC = [];
amplitudesB_vHPC = []; amplitudesR_vHPC = []; amplitudesA_vHPC = [];

CCG_B = [];
CCG_R = [];
CCG_A = [];

% Parameters for CCG ripples-SU
b = 0.01; % binsize for ripples-SU modulation
dd = 1; % time window for ripples-SU modulation
sss = 1; % smooth level for ripples-SU modulation

dHPC_neurons_ripples_B = [];
dHPC_neurons_ripples_R = [];
dHPC_neurons_ripples_A = [];


vHPC_neurons_ripples_B = [];
vHPC_neurons_ripples_R = [];
vHPC_neurons_ripples_A = [];

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

%Behavior
dHPC_Shock = []; vHPC_Shock = [];
dHPC_Valve = []; vHPC_Valve = [];
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
        ripplesD = table2array(readtable('ripplesD_customized.csv'));
        ripplesV = table2array(readtable('ripplesV_customized.csv'));
        
        % Coordinated dHPC ripples
        coordinated = [];
        coordinatedV = [];
        for i = 1:length(ripplesD)
            r = ripplesD(i,:);
            tmp = sum(and(ripplesV(:,2)>= r(1,1)-0.05, ripplesV(:,2)<= r(1,3)+0.05));
            if tmp>0
                coordinatedV = [coordinatedV ; ripplesV(and(ripplesV(:,2)>= r(1,1)-0.05, ripplesV(:,2)<= r(1,3)+0.05),:)];
                coordinated = [coordinated ; r];
                clear tmp2 tmp1
            end
            clear r
        end
        clear x tmp i
        
        coordinatedB = Restrict(coordinated,NREM_B);    coordinatedA = Restrict(coordinated,NREM_A);    coordinatedR = Restrict(coordinated,NREM_R);
        coordinatedB_V = Restrict(coordinatedV,NREM_B);    coordinatedR_V = Restrict(coordinatedV,NREM_R);    coordinatedA_V = Restrict(coordinatedV,NREM_A);

        p_coordinatedB = (length(coordinatedB)*100)/length(Restrict(ripplesD,NREM_B));
        p_coordinatedA = (length(coordinatedA)*100)/length(Restrict(ripplesD,NREM_A));
        p_coordinatedR = (length(coordinatedR)*100)/length(Restrict(ripplesD,NREM_R));
        ripples_coordinated_percentage = [ripples_coordinated_percentage ; p_coordinatedB , p_coordinatedR , p_coordinatedA];
        clear p_coordinatedB p_coordinatedA p_coordinatedR
        
        rateV_B = length(Restrict(ripplesV,NREM_B)) / sum(NREM_B(:,2) - NREM_B(:,1));
        rateV_A = length(Restrict(ripplesV,NREM_A)) / sum(NREM_A(:,2) - NREM_A(:,1));
        rateV_R = length(Restrict(ripplesV,NREM_R)) / sum(NREM_R(:,2) - NREM_R(:,1));
        rateV = [rateV ; rateV_B , rateV_R , rateV_A];
        clear rateV_B rateV_A rateV_R
        
        rateD_B = length(Restrict(ripplesD,NREM_B)) / sum(NREM_B(:,2) - NREM_B(:,1));
        rateD_A = length(Restrict(ripplesD,NREM_A)) / sum(NREM_A(:,2) - NREM_A(:,1));
        rateD_R = length(Restrict(ripplesD,NREM_R)) / sum(NREM_R(:,2) - NREM_R(:,1));
        rateD = [rateD ; rateD_B , rateD_R , rateD_A];
        clear rateD_B rateD_A rateD_R
        
        x = Restrict(ripplesV(:,2),NREM_B);
        y = Restrict(ripplesD(:,2),NREM_B);
        duracion = sum(NREM_B(:,2) - NREM_B(:,2));
        [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
        [ccg,ttt] = CCG(s,ids,'binSize',binsize,'duration',d,'smooth',ss,'mode','ccg');
        CCG_B = [CCG_B , ccg(:,1,2)];
        clear x y s ids groups ccg duracion
        
        x = Restrict(ripplesV(:,2),NREM_R);
        y = Restrict(ripplesD(:,2),NREM_R);
        duracion = sum(NREM_R(:,2) - NREM_R(:,2));
        [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
        [ccg,ttt] = CCG(s,ids,'binSize',binsize,'duration',d,'smooth',ss,'mode','ccg');
        CCG_R = [CCG_R , ccg(:,1,2)];
        clear x y s ids groups ccg duracion 
        
        x = Restrict(ripplesV(:,2),NREM_A);
        y = Restrict(ripplesD(:,2),NREM_A);
        duracion = sum(NREM_A(:,2) - NREM_A(:,2));
        [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
        [ccg,ttt] = CCG(s,ids,'binSize',binsize,'duration',d,'smooth',ss,'mode','ccg');
        CCG_A = [CCG_A , ccg(:,1,2)];
        clear x y s ids groups ccg duracion        
    
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
    spks(:,1) = double(spks(:,1))./20000;
    

%     for i=1:length(group_vHPC)
%         [ST_V(:,i),bins]=binspikes(spks(spks(:,2)==group_vHPC(i),1),1/binsize,[0 segments.Var1(end)/1000]);
%     end
% 
%     for i=1:length(group_dHPC)
%         [ST_D(:,i),bins]=binspikes(spks(spks(:,2)==group_dHPC(i),1),1/binsize,[0 segments.Var1(end)/1000]);
%     end
%     
%     vHPC_Shock = [];
%     for i = 1 : length(group_vHPC)% ventral SU
%         if logical(group_vHPC(i,3))
%             tmp = [];
%             for ii = 1 : length(Shocks_filt)
%                 n = InIntervals(bins,[Shocks_filt(ii)-2, Shocks_filt(ii)+2]);
%                 tmp = [tmp , ST_V(n,i)];
%                 clear n
%             end
%             vHPC_Shock = [vHPC_Shock , Smooth(mean(tmp,2)./binsize,ss)];
%             clear tmp
%         end
%     end
%     
    
%     imagesc([-2:binsize:2-binsize], [1:1:size(vHPC_Shock,2)], vHPC_Shock'), axis tight, xline(0,'--','LineWidth',1),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default'),caxis([0 8])
    
        %% Peri-event firing rate histogram to awake stimuli
     % dHPC
    base = InvertIntervals([Shocks_filt-4, Shocks_filt+4],aversiveTS_run(1)/1000 , aversiveTS_run(2)/1000);
    criteriaD = [];
    for i = 1 : length(group_dHPC)% ventral SU
        if logical(group_dHPC(i,3))
            y = Shocks_filt;
            x =  Restrict(spks_dHPC(spks_dHPC(:,2) == group_dHPC(i,1),1),aversiveTS_run./1000);
            if and(~isempty(x),length(x)>20)
                Base = Restrict(x,base);
                [times,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
                [ccg,TimeVector1] = CCG(times,groups,'binsize',binsize,'duration',d,'smooth',ss);
                FR = length(Base)/sum(base(:,2)-base(:,1));
                dHPC_Shock = [dHPC_Shock , ccg(:,1,2)./length(y)./binsize./FR];
                
                [~,ind] = min(abs(TimeVector1-0));
                [~,ind2] = min(abs(TimeVector1-1));
                if mean(ccg(ind:ind2,1,2)./length(y)./binsize./FR)>2
                        criteriaD = [criteriaD ; true];
                else
                        criteriaD = [criteriaD ; false];
                end
                
                clear tmp FR Base ccg x y times ids groups
            else
                dHPC_Shock = [dHPC_Shock , nan(size(dHPC_Shock,1),1)];
                criteriaD = [criteriaD ; false];
                
            end
        else
                criteriaD = [criteriaD ; false];
        end
    end
    clear base i 
    
    % vHPC
    criteriaV = [];
    base = InvertIntervals([Shocks_filt-4, Shocks_filt+4],aversiveTS_run(1)/1000 , aversiveTS_run(2)/1000);
    for i = 1 : length(group_vHPC)% ventral SU
        if logical(group_vHPC(i,3))
            y = Shocks_filt;
            x =  Restrict(spks_vHPC(spks_vHPC(:,2) == group_vHPC(i,1),1),aversiveTS_run./1000);
            if and(~isempty(x),length(x)>20)
                Base = Restrict(x,base);
                [times,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
                [ccg,TimeVector1] = CCG(times,groups,'binsize',binsize,'duration',d,'smooth',ss);
                FR = length(Base)/sum(base(:,2)-base(:,1));
                vHPC_Shock = [vHPC_Shock , ccg(:,1,2)./length(y)./binsize./FR];
                
                [~,ind] = min(abs(TimeVector1-0));
                [~,ind2] = min(abs(TimeVector1-1));
                if mean(ccg(ind:ind2,1,2)./length(y)./binsize./FR)>2
                    criteriaV = [criteriaV ; true];
                else
                    criteriaV = [criteriaV ; false];
                end
                clear tmp FR Base ccg x y times ids groups
            else
                vHPC_Shock = [vHPC_Shock , nan(size(vHPC_Shock,1),1)];
                criteriaV = [criteriaV ; false];
                
            end
        else
                criteriaV = [criteriaV ; false];
        end
    end
    clear base i
    
    % -----------------------------------
    % Valve
    % -----------------------------------
    % dHPC
    base = InvertIntervals([Rewards_filt-4, Rewards_filt+4],rewardTS_run(1)/1000 , rewardTS_run(2)/1000);
    tmp = [];
    for i = 1 : length(group_dHPC)% ventral SU
        if logical(group_dHPC(i,3))
            y = Rewards_filt(:,1);
            x =  Restrict(spks_dHPC(spks_dHPC(:,2) == group_dHPC(i,1),1),rewardTS_run./1000);
            if and(~isempty(x),length(x)>20)
                Base = Restrict(x,base);
                [times,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
                [ccg,TimeVector1] = CCG(times,groups,'binsize',binsize,'duration',d,'smooth',ss);
                FR = length(Base)/sum(base(:,2)-base(:,1));
                dHPC_Valve = [dHPC_Valve , ccg(:,1,2)./length(y)./binsize./FR];
                
               [~,ind] = min(abs(TimeVector1-(-1)));
                [~,ind2] = min(abs(TimeVector1-1));
                if mean(ccg(ind:ind2,1,2)./length(y)./binsize./FR)>2
                        tmp = [tmp ; true];
                else
                        tmp = [tmp ; false];
                end                
                clear FR Base ccg x y times ids groups
            else
                dHPC_Valve = [dHPC_Valve , nan(size(dHPC_Valve,1),1)];
                        tmp = [tmp ; false];
            end
        else
            tmp = [tmp ; false];
        end
    end
    clear base i
    criteriaD = [criteriaD,tmp];

    % vHPC
    base = InvertIntervals([Rewards_filt-4, Rewards_filt+4],rewardTS_run(1)/1000 , rewardTS_run(2)/1000);
    tmp = [];
    for i = 1 : length(group_vHPC)% ventral SU
        if  logical(group_vHPC(i,3))
            y = Rewards_filt(:,1);
            x =  Restrict(spks_vHPC(spks_vHPC(:,2) == group_vHPC(i,1),1),rewardTS_run./1000);
            if and(~isempty(x),length(x)>20)
                Base = Restrict(x,base);
                [times,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
                [ccg,TimeVector1] = CCG(times,groups,'binsize',binsize,'duration',d,'smooth',ss);
                FR = length(Base)/sum(base(:,2)-base(:,1));
                vHPC_Valve = [vHPC_Valve , ccg(:,1,2)./length(y)./binsize./FR];
                
                [~,ind] = min(abs(TimeVector1-(-1)));
                [~,ind2] = min(abs(TimeVector1-1));
                if mean(ccg(ind:ind2,1,2)./length(y)./binsize./FR)>2
                        tmp = [tmp ; true];
                else
                        tmp = [tmp ; false];
                end
                clear FR Base ccg x y times ids groups
            else
                vHPC_Valve = [vHPC_Valve , nan(size(vHPC_Valve,1),1)];
                tmp = [tmp ; false];
            end
        else
            tmp = [tmp ; false];
        end
    end
    clear base i
    criteriaV = [criteriaV , tmp];
        save([cd,'\tag_shock_responsive_cells.mat'],'criteriaD','criteriaV')
        
                clear jj i numlaps R TTLs_camara_aversive TTLs_camara_reward
        clear aversiveTS aversiveTS_run rewardTS rewardTS_run baselineTS
        clear NREM NREM_A NREM_B NREM_R REM REM_A REM_B REM_R
        clear segments ripplesD ripplesV WAKE camara Cell_type_classification
        clear coordinated coordinatedA coordinatedA_V coordinatedB coordinatedB_V
        clear coordinatedR coordinatedR_V coordinatedV
        clear criteriaD criteriaV K Kinfo laps tmpD tmpV
        t
    end
    tt
end
        %% PHIST sourrouding dorsal ripples
        %Dorsal SU - Dorsal ripples
        event1 = Restrict(ripplesD,NREM_B);
        event1_NR = InvertIntervals([event1(:,1)-0.05 event1(:,3)+0.05] , baselineTS(1)/1000 , baselineTS(2)/1000);
        event1 = coordinatedB;
        event2 = Restrict(ripplesD,NREM_R);
        event2_NR = InvertIntervals([event2(:,1)-0.05 event2(:,3)+0.05] ,  rewardTS(1)/1000 , rewardTS(2)/1000);
        event2 = coordinatedR;
        event3 = Restrict(ripplesD,NREM_A);
        event3_NR = InvertIntervals([event3(:,1)-0.05 event3(:,3)+0.05] , aversiveTS(1)/1000 , aversiveTS(2)/1000);
        event3 = coordinatedA;
        
        tmpD = [];
        for i = 1 : length(group_dHPC)% dorsal SU
            %Check if is a Pyr
            if and(logical(group_dHPC(i,3)),logical(criteriaD(i,1)))
%             if and(and(logical(group_dHPC(i,3)),logical(criteriaD(i,1))),~logical(criteriaD(i,3)))
                spks = spks_dHPC(spks_dHPC(:,1)==group_dHPC(i,1),:);
                spks = spks(:,2);
                %Restrict Spikes within the ripples events
                % Baseline
                Spks1 = Restrict(spks,[min(event1(:,2))-10 max(event1(:,2)+1)+10]);
                m1 = (length(Restrict(spks,event1_NR)) / sum(event1_NR(:,2)-event1_NR(:,1)));
                cond1 = length(Spks1) > 20;%/(max(event1(:,2))-min(event1(:,2))) >= m1*0.2;
                % Reward
                Spks2 = Restrict(spks,[min(event2(:,2))-10 max(event2(:,2))+10]);
                m2 = (length(Restrict(spks,event2_NR)) / sum(event2_NR(:,2)-event2_NR(:,1)));
                cond2 = length(Spks2) > 20;%/(max(event2(:,2))-min(event2(:,2))) >= m2*0.2;
                % Aversive
                Spks3 = Restrict(spks,[min(event3(:,2))-10 max(event3(:,2))+10]);
                m3 = (length(Restrict(spks,event3_NR)) / sum(event3_NR(:,2)-event3_NR(:,1)));
                cond3 = length(Spks3) > 20;%/(max(event3(:,2))-min(event3(:,2))) >= m3*0.2;
                
                %check if FR within each interval is at least half of the mean FR
                if and(and(cond1 , cond2) , cond3)
                    %Baseline
                    [s,ids,groups] = CCGParameters(event1(:,2),ones(length(event1(:,2)),1),Spks1,ones(length(Spks1),1)*2);
                    [ccg,T] = CCG(s,ids,'binSize',b,'duration',dd,'smooth',sss,'mode','ccg');
                    %calculation of mean FRss
                    dHPC_neurons_ripples_B = [dHPC_neurons_ripples_B , ((ccg(:,1,2)./length(event1))./b)./m1];
                    clear s ids groups ccg T z m1 tmp ins Spks1
                    
                    %Reward
                    [s,ids,groups] = CCGParameters(event2(:,2),ones(length(event2(:,2)),1),Spks2,ones(length(Spks2),1)*2);
                    [ccg,T] = CCG(s,ids,'binSize',b,'duration',dd,'smooth',sss,'mode','ccg');
                    %calculation of mean FR
                    dHPC_neurons_ripples_R = [dHPC_neurons_ripples_R , ((ccg(:,1,2)./length(event2))./b)./m2];
                    clear s ids groups ccg T z m2 tmp ins Spks2
                    
                    %Aversive
                    [s,ids,groups] = CCGParameters(event3(:,2),ones(length(event3(:,2)),1),Spks3,ones(length(Spks3),1)*2);
                    [ccg,T] = CCG(s,ids,'binSize',b,'duration',dd,'smooth',sss,'mode','ccg');
                    %calculation of mean FR
                    dHPC_neurons_ripples_A = [dHPC_neurons_ripples_A , ((ccg(:,1,2)./length(event3))./b)./m3];
                    clear s ids groups ccg z m3 tmp ins Spks3
                    tmpD = [tmpD , criteriaD(i,2)];
                    
                end
            end
            clear spks rate cond1 cond2 cond3
        end
        
        clear i event1 event2 event3
        clear event1_NR event2_NR event3_NR
        
        event1 = Restrict(ripplesV,NREM_B);
        event1_NR = InvertIntervals([event1(:,1)-0.05 event1(:,3)+0.05] , baselineTS(1)/1000 , baselineTS(2)/1000);
        event1 = coordinatedB_V;
        event2 = Restrict(ripplesV,NREM_R);
        event2_NR = InvertIntervals([event2(:,1)-0.05 event2(:,3)+0.05] ,  rewardTS(1)/1000 , rewardTS(2)/1000);
        event2 = coordinatedR_V;
        event3 = Restrict(ripplesV,NREM_A);
        event3_NR = InvertIntervals([event3(:,1)-0.05 event3(:,3)+0.05] , aversiveTS(1)/1000 , aversiveTS(2)/1000);
        event3 = coordinatedA_V;
        tmpV = [];
        %Ventral SU - Dorsal ripples
        for i = 1 : length(group_vHPC)% ventral SU
            %Check if is a Pyr
            if and(logical(group_vHPC(i,3)),logical(criteriaV(i,1)))
%             if and(and(logical(group_vHPC(i,3)),logical(criteriaV(i,1))),~logical(criteriaV(i,3)))
                spks = spks_vHPC(spks_vHPC(:,1)==group_vHPC(i,1),:);
                spks = spks(:,2);
                %Restrict Spikes within the ripples events
                % Baseline
                Spks1 = Restrict(spks,[min(event1(:,2))-10 max(event1(:,2)+1)+10]);
                m1 = (length(Restrict(spks,event1_NR)) / sum(event1_NR(:,2)-event1_NR(:,1)));
                cond1 = length(Spks1) > 20;%/(max(event1(:,2))-min(event1(:,2))) >= m1*0.2;
                % Reward
                Spks2 = Restrict(spks,[min(event2(:,2))-10 max(event2(:,2))+10]);
                m2 = (length(Restrict(spks,event2_NR)) / sum(event2_NR(:,2)-event2_NR(:,1)));
                cond2 = length(Spks2) > 20;%/(max(event2(:,2))-min(event2(:,2))) >= m2*0.2;
                % Aversive
                Spks3 = Restrict(spks,[min(event3(:,2))-10 max(event3(:,2))+10]);
                m3 = (length(Restrict(spks,event3_NR)) / sum(event3_NR(:,2)-event3_NR(:,1)));
                cond3 = length(Spks3) > 20;%/(max(event3(:,2))-min(event3(:,2))) >= m3*0.2;
                
                %check if FR within each interval is equal or higer to mean FR
                if and(and(cond1 , cond2) , cond3)
                    %Baseline
                    [s,ids,groups] = CCGParameters(event1(:,2),ones(length(event1(:,2)),1),Spks1,ones(length(Spks1),1)*2);
                    [ccg,T] = CCG(s,ids,'binSize',b,'duration',dd,'smooth',sss,'mode','ccg');
                    %calculation of mean FRss
                    vHPC_neurons_ripples_B = [vHPC_neurons_ripples_B , ((ccg(:,1,2)./length(event1))./b)./m1];
                    clear s ids groups ccg T z m1 tmp ins Spks1
                    
                    %Reward
                    [s,ids,groups] = CCGParameters(event2(:,2),ones(length(event2(:,2)),1),Spks2,ones(length(Spks2),1)*2);
                    [ccg,T] = CCG(s,ids,'binSize',b,'duration',dd,'smooth',sss,'mode','ccg');
                    %calculation of mean FR
                    vHPC_neurons_ripples_R = [vHPC_neurons_ripples_R , ((ccg(:,1,2)./length(event2))./b)./m2];
                    clear s ids groups ccg T z m2 tmp ins Spks2
                    
                    %Aversive
                    [s,ids,groups] = CCGParameters(event3(:,2),ones(length(event3(:,2)),1),Spks3,ones(length(Spks3),1)*2);
                    [ccg,T] = CCG(s,ids,'binSize',b,'duration',dd,'smooth',sss,'mode','ccg');
                    %calculation of mean FR
                    vHPC_neurons_ripples_A = [vHPC_neurons_ripples_A , ((ccg(:,1,2)./length(event3))./b)./m3];
                    clear s ids groups ccg z m3 tmp ins Spks3
                    tmpV = [tmpV , criteriaV(i,2)];
                end
            end
            clear spks rate cond1 cond2 cond3
        end
        clear i event1 event2 event3
        clear event1_NR event2_NR event3_NR
        
        CD = [CD ,  tmpD];
        CV = [CV ,  tmpV];
        
        save([cd,'\tag_shock_responsive_cells.mat'],'criteriaD','criteriaV')
        t
        clear jj i numlaps R TTLs_camara_aversive TTLs_camara_reward
        clear aversiveTS aversiveTS_run rewardTS rewardTS_run baselineTS
        clear NREM NREM_A NREM_B NREM_R REM REM_A REM_B REM_R
        clear segments ripplesD ripplesV WAKE camara Cell_type_classification
        clear coordinated coordinatedA coordinatedA_V coordinatedB coordinatedB_V
        clear coordinatedR coordinatedR_V coordinatedV
        clear criteriaD criteriaV K Kinfo laps tmpD tmpV
    end
    tt
end

%% PHIST
% vHPC
[~ , i] = min(abs(ttt-0));
[~ , ii] = min(abs(ttt-0.5));

x = [];
xx = [];
vHPC_exc_S = [];
vHPC_exc_V = [];
vHPC_inh_S = [];
vHPC_inh_V = [];

for iii = 1 : size(vHPC_neurons_shock_zscore,2)
    if and(~isnan(sum(vHPC_neurons_shock_zscore(:,iii))) , ~isnan(sum(vHPC_neurons_valve_zscore(:,iii))))
        if and(not(mean(vHPC_neurons_shock_zscore(:,iii)) == 0) , not(mean(vHPC_neurons_valve_zscore(:,iii)) == 0))
            x = [ x , vHPC_neurons_shock_zscore(:,iii)];
            xx = [ xx , vHPC_neurons_valve_zscore(:,iii)];
            
            if mean(vHPC_neurons_shock_zscore(i:ii,iii))>1
                vHPC_exc_S = [vHPC_exc_S , vHPC_neurons_shock_zscore(:,iii)];
                vHPC_exc_V = [vHPC_exc_V , vHPC_neurons_valve_zscore(:,iii)];
            elseif mean(vHPC_neurons_shock_zscore(i:ii,iii))<-1
                vHPC_inh_S = [vHPC_inh_S , vHPC_neurons_shock_zscore(:,iii)];
                vHPC_inh_V = [vHPC_inh_V , vHPC_neurons_valve_zscore(:,iii)];    
            end
        end
    end
end


figure,
subplot(211)
m = mean(vHPC_exc_S,2); s = std(vHPC_exc_S,0,2)./sqrt(size(vHPC_exc_S,2));
time = [-2:binsize:2];
plot(time',m,'r','LineWidth',1),hold on, ciplot(m-s,m+s,time,'r'),alpha 0.05

m = mean(vHPC_exc_V,2); s = std(vHPC_exc_V,0,2)./sqrt(size(vHPC_exc_V,2));
time = [-2:binsize:2];
plot(time',m,'b','LineWidth',1),hold on, ciplot(m-s,m+s,time,'b'),alpha 0.05
ylabel('Neuronal Activity (z-score)'), xlabel('Time (sec)')

subplot(212)
m = mean(vHPC_inh_S,2); s = std(vHPC_inh_S,0,2)./sqrt(size(vHPC_exc_S,2));
time = [-2:binsize:2];
plot(time',m,'r','LineWidth',1),hold on, ciplot(m-s,m+s,time,'r'),alpha 0.05

m = mean(vHPC_inh_V,2); s = std(vHPC_inh_V,0,2)./sqrt(size(vHPC_exc_V,2));
time = [-2:binsize:2];
plot(time',m,'b','LineWidth',1),hold on, ciplot(m-s,m+s,time,'b'),alpha 0.05
ylabel('Neuronal Activity (z-score)'),xlabel('Time (sec)')



y = mean(x(i:ii,:),1);
[~ , y] = sort(y','descend');

z = [];
zz = [];
for i = 1:size(x,2)
    z = [z , x(:,y(i))];
    zz = [zz , xx(:,y(i))];
end

figure
subplot(121), imagesc(ttt, [1:1:size(x,2)], z'), axis tight, xline(0,'LineWidth',1),xlim([-2 2]),caxis([-2 2]),xlabel('Time(sec)'),ylabel('Neurons(id)')
xline(1,'--','LineWidth',1)
title('vHPC Shock')

subplot(122), imagesc(ttt, [1:1:size(x,2)], zz'), axis tight, xline(0,'LineWidth',1),xlim([-2 2]),caxis([-2 2]),xlabel('Time(sec)'),ylabel('Neurons(id)')
% xline(1,'--','LineWidth',1)
title('vHPC Valve')

%dHPC
[~ , i] = min(abs(ttt-0));
[~ , ii] = min(abs(ttt-1));

x = [];
xx = [];
dHPC_exc_S = [];
dHPC_exc_V = [];
dHPC_inh_S = [];
dHPC_inh_V = [];
for iii = 1 : size(dHPC_neurons_shock_zscore,2)
    if and(~isnan(sum(dHPC_neurons_shock_zscore(:,iii))) , ~isnan(sum(dHPC_neurons_valve_zscore(:,iii))))
        if and(not(mean(dHPC_neurons_shock_zscore(:,iii)) == 0) , not(mean(dHPC_neurons_valve_zscore(:,iii)) == 0))
            x = [ x , dHPC_neurons_shock_zscore(:,iii)];
            xx = [ xx , dHPC_neurons_valve_zscore(:,iii)];
            
            if mean(dHPC_neurons_shock_zscore(i:ii,iii))>1
                dHPC_exc_S = [dHPC_exc_S , dHPC_neurons_shock_zscore(:,iii)];
                dHPC_exc_V = [dHPC_exc_V , dHPC_neurons_valve_zscore(:,iii)];
            elseif mean(dHPC_neurons_shock_zscore(i:ii,iii))<-1
                dHPC_inh_S = [dHPC_inh_S , dHPC_neurons_shock_zscore(:,iii)];
                dHPC_inh_V = [dHPC_inh_V , dHPC_neurons_valve_zscore(:,iii)];    
            end        
        end
    end
end



figure,
subplot(211)
m = mean(dHPC_exc_S,2); s = std(dHPC_exc_S,0,2)./sqrt(size(dHPC_exc_S,2));
time = [-2:binsize:2];
plot(time',m,'r','LineWidth',1),hold on, ciplot(m-s,m+s,time,'r'),alpha 0.05

m = mean(dHPC_exc_V,2); s = std(dHPC_exc_V,0,2)./sqrt(size(dHPC_exc_V,2));
time = [-2:binsize:2];
plot(time',m,'b','LineWidth',1),hold on, ciplot(m-s,m+s,time,'b'),alpha 0.05
ylabel('Neuronal Activity (z-score)'), xlabel('Time (sec)')

subplot(212)
m = mean(dHPC_inh_S,2); s = std(vHPC_inh_S,0,2)./sqrt(size(vHPC_exc_S,2));
time = [-2:binsize:2];
plot(time',m,'r','LineWidth',1),hold on, ciplot(m-s,m+s,time,'r'),alpha 0.05

m = mean(dHPC_inh_V,2); s = std(vHPC_inh_V,0,2)./sqrt(size(vHPC_exc_V,2));
time = [-2:binsize:2];
plot(time',m,'b','LineWidth',1),hold on, ciplot(m-s,m+s,time,'b'),alpha 0.05
ylabel('Neuronal Activity (z-score)'),xlabel('Time (sec)')


y = mean(x(i:ii,:),1);
[~ , y] = sort(y','descend');

z = [];
zz = [];
for i = 1:size(x,2)
    z = [z , x(:,y(i))];
    zz = [zz , xx(:,y(i))];
end

figure
subplot(121), imagesc(ttt, [1:1:size(x,2)], z'), axis tight, xline(0,'LineWidth',1),xlim([-2 2]),caxis([-2 2]),xlabel('Time(sec)'),ylabel('Neurons(id)')
xline(1,'--','LineWidth',1)
title('dHPC Shock')

subplot(122), imagesc(ttt, [1:1:size(x,2)], zz'), axis tight, xline(0,'LineWidth',1),xlim([-2 2]),caxis([-2 2]),xlabel('Time(sec)'),ylabel('Neurons(id)')
% xline(1,'--','LineWidth',1)
title('dHPC Valve')



figure,
subplot(121),plot([-2:1/30:2-1/30],mean(velocity_shock,2)),ylim([4 30])
subplot(122),plot([-2:1/30:2-1/30],mean(velocity_valve,2)),ylim([4 30])


%% PHIST ripples-SU
% vHPC
[~ , i] = min(abs(T-(-0.05)));
[~ , ii] = min(abs(T-0.05));

x = [];
xx = [];
xxx = [];

for iii = 1 : size(vHPC_neurons_ripples_A,2)
    if and(and(~isnan(sum(vHPC_neurons_ripples_A(:,iii))) , ~isnan(sum(vHPC_neurons_ripples_R(:,iii)))),~isnan(sum(vHPC_neurons_ripples_B(:,iii))))
            x = [x , vHPC_neurons_ripples_B(:,iii)];
            xx = [xx , vHPC_neurons_ripples_R(:,iii)];
            xxx = [xxx , vHPC_neurons_ripples_A(:,iii)];
    end
end

y = mean(xxx(i:ii,:),1);
[~ , y] = sort(y','descend');

z = [];
zz = [];
zzz = [];
for i = 1:size(x,2)
%     m = min(x(:,y(i)));
%     mm = x(:,y(i)) - m;
%     mmm = mm./max(mm);
%     z = [z , mmm];
    z = [z , x(:,y(i))];
    clear m mm mmm
    
%     m = min(xx(:,y(i)));
%     mm = xx(:,y(i)) - m;
%     mmm = mm./max(mm);
%     zz = [zz , mmm];
    zz = [zz , xx(:,y(i))];
    clear m mm mmm
%     
%     m = min(xxx(:,y(i)));
%     mm = xxx(:,y(i)) - m;
%     mmm = mm./max(mm);
%     zzz = [zzz , mmm];
    zzz = [zzz , xxx(:,y(i))];
    clear m mm mmm
end

figure
subplot(131), imagesc(T, [1:1:size(x,2)], z'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default'),caxis([0 4])
title('Baseline')
subplot(132), imagesc(T, [1:1:size(xx,2)], zz'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default'),caxis([0 4])
title('Reward')
subplot(133), imagesc(T, [1:1:size(xxx,2)], zzz'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default'),caxis([0 4])
title('Aversive')

figure,
plot(T,mean(z,2),'k'),hold on
plot(T,mean(zz,2),'b'),hold on
plot(T,mean(zzz,2),'r'),hold on

clear x xx xxx z zz zzz


% dHPC
[~ , i] = min(abs(T-(-0.05)));
[~ , ii] = min(abs(T-0.05));

x = [];
xx = [];
xxx = [];

for iii = 1 : size(dHPC_neurons_ripples_A,2)
    if and(and(~isnan(sum(dHPC_neurons_ripples_A(:,iii))) , ~isnan(sum(dHPC_neurons_ripples_R(:,iii)))),~isnan(sum(dHPC_neurons_ripples_B(:,iii))))
            x = [x , dHPC_neurons_ripples_B(:,iii)];
            xx = [xx , dHPC_neurons_ripples_R(:,iii)];
            xxx = [xxx , dHPC_neurons_ripples_A(:,iii)];
    end
end

y = mean(xxx(i:ii,:),1);
[~ , y] = sort(y','descend');

z = [];
zz = [];
zzz = [];
for i = 1:size(x,2)
%     m = min(x(:,y(i)));
%     mm = x(:,y(i)) - m;
%     mmm = mm./max(mm);
%     z = [z , mmm];
    z = [z , x(:,y(i))];
    clear m mm mmm
    
%     m = min(xx(:,y(i)));
%     mm = xx(:,y(i)) - m;
%     mmm = mm./max(mm);
%     zz = [zz , mmm];
    zz = [zz , xx(:,y(i))];
    clear m mm mmm
    
%     m = min(xxx(:,y(i)));
%     mm = xxx(:,y(i)) - m;
%     mmm = mm./max(mm);
%     zzz = [zzz , mmm];
    zzz = [zzz , xxx(:,y(i))];
    clear m mm mmm
end


figure
subplot(131), imagesc(T, [1:1:size(x,2)], z'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default'),caxis([0 8])
title('Baseline')
subplot(132), imagesc(T, [1:1:size(xx,2)], zz'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default'),caxis([0 8])
title('Reward')
subplot(133), imagesc(T, [1:1:size(xxx,2)], zzz'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default'),caxis([0 8])
title('Aversive')


figure,
plot(T,mean(z,2),'k'),hold on
plot(T,mean(zz,2),'b'),hold on
plot(T,mean(zzz,2),'r'),hold on



%% Selection of SU according to their response during Shock
% vHPC
q = quantile(CV,4);
[~ , i] = min(abs(T-(-0.05)));
[~ , ii] = min(abs(T-0.05));

x = [];
xx = [];
xxx = [];

for iii = 1 : size(vHPC_neurons_ripples_A,2)
    if and(and(~isnan(sum(vHPC_neurons_ripples_A(:,iii))) , ~isnan(sum(vHPC_neurons_ripples_R(:,iii)))),~isnan(sum(vHPC_neurons_ripples_B(:,iii))))
        if CV(iii)>q(4)
            x = [x , vHPC_neurons_ripples_B(:,iii)];
            xx = [xx , vHPC_neurons_ripples_R(:,iii)];
            xxx = [xxx , vHPC_neurons_ripples_A(:,iii)];
        end
    end
end

y = mean(xxx(i:ii,:),1);
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
subplot(131), imagesc(T, [1:1:size(x,2)], z'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default'),caxis([0 4])
title('Baseline')
subplot(132), imagesc(T, [1:1:size(xx,2)], zz'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default'),caxis([0 4])
title('Reward')
subplot(133), imagesc(T, [1:1:size(xxx,2)], zzz'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default'),caxis([0 4])
title('Aversive')

figure,
plot(T,mean(z,2),'k'),hold on
plot(T,mean(zz,2),'b'),hold on
plot(T,mean(zzz,2),'r'),hold on

clear x xx xxx z zz zzz q


% dHPC
q = quantile(CD,4);
[~ , i] = min(abs(T-(-0.05)));
[~ , ii] = min(abs(T-0.05));

x = [];
xx = [];
xxx = [];

for iii = 1 : size(dHPC_neurons_ripples_A,2)
    if and(and(~isnan(sum(dHPC_neurons_ripples_A(:,iii))) , ~isnan(sum(dHPC_neurons_ripples_R(:,iii)))),~isnan(sum(dHPC_neurons_ripples_B(:,iii))))
        if CD(iii)>q(4)
            
            x = [x , dHPC_neurons_ripples_B(:,iii)];
            xx = [xx , dHPC_neurons_ripples_R(:,iii)];
            xxx = [xxx , dHPC_neurons_ripples_A(:,iii)];
        end
    end
end

y = mean(xxx(i:ii,:),1);
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
subplot(131), imagesc(T, [1:1:size(x,2)], z'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default'),caxis([0 8])
title('Baseline')
subplot(132), imagesc(T, [1:1:size(xx,2)], zz'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default'),caxis([0 8])
title('Reward')
subplot(133), imagesc(T, [1:1:size(xxx,2)], zzz'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default'),caxis([0 8])
title('Aversive')


figure,
plot(T,mean(z,2),'k'),hold on
plot(T,mean(zz,2),'b'),hold on
plot(T,mean(zzz,2),'r'),hold on