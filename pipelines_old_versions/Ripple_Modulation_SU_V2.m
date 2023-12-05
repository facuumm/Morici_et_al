clear
clc
close all

%Parameters

%%
% $$\$$

path = {'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat103\usable'};%;'E:\Rat128\Ephys\in_pyr\to_finish_spikesorting\usable'};
% path = 'E:\Rat128\Ephys\in_pyr\ready';

%% Parameters
% ---------------
% --> Ripples
% ---------------
ripples_coordinated_percentage = [];
ripples_coordinatedV_percentage = [];
rateV = []; rateD = []; % to store the ripples rate from dHPC and vHPC

CCG_B = []; CCG_BD = []; CCG_BV = [];
CCG_R = []; CCG_RD = []; CCG_RV = [];
CCG_A = []; CCG_AD = []; CCG_AV = [];

ACG_BD = []; ACG_RD = []; ACG_AD = [];
ACG_BV = []; ACG_RV = []; ACG_AV = [];

ACG_All_V = [];
ACG_All_D = [];

CCG_All = [];
totalTime = []; totalTimeB = []; totalTimeR = []; totalTimeA = [];

% Parameters for CCG ripples-SU
b = 0.01; % binsize for ripples-SU modulation
dd = 1; % time window for ripples-SU modulation
ss = 1; % smooth level for ripples-SU modulation

celltype = 1; %to define which celltype to use (1: pyr / 2:int)

% To store the CCG of each condition
dHPC_neurons_ripples_B = []; dHPC_neurons_ripples_R = []; dHPC_neurons_ripples_A = []; dHPC_neurons_ripples = [];
dHPC_neurons_bursts = []; dHPC_neurons_bursts_B = [];  dHPC_neurons_bursts_R = [];  dHPC_neurons_bursts_A = [];

vHPC_neurons_ripples_B = []; vHPC_neurons_ripples_R = []; vHPC_neurons_ripples_A = []; vHPC_neurons_ripples = [];
vHPC_neurons_bursts = []; vHPC_neurons_bursts_B = [];  vHPC_neurons_bursts_R = [];  vHPC_neurons_bursts_A = [];


% to store poisson result
dHPC_dRipples_poisson = []; vHPC_vRipples_poisson = [];

%% -------------------------------------
count = 0;
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
        count = count+1;
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
        
        ripplesD = table2array(readtable('ripplesD_customized2.csv'));
        ripplesV = table2array(readtable('ripplesV_customized2.csv'));
        
        coordinated = [];
        coordinatedV = [];
        for i = 1:length(ripplesD)
            r = ripplesD(i,:);
            tmp = sum(and(ripplesV(:,2)>= r(1,2)-0.1, ripplesV(:,2)<= r(1,2)+0.1));
            if tmp>0
                z = ripplesV(and(ripplesV(:,2)>= r(1,2)-0.2, ripplesV(:,2)<= r(1,2)+0.2),:);
                [p,indice] = min(abs(r(2)-z(:,2)));
                coordinatedV = [coordinatedV ; z(indice,:)];
                coordinated = [coordinated ; r];
                clear tmp2 tmp1 p indice z
            end
            clear r
        end
        clear x tmp i
        
        coordinatedB = Restrict(coordinated,NREM_B);    coordinatedA = Restrict(coordinated,NREM_A);    coordinatedR = Restrict(coordinated,NREM_R);
        coordinatedB_V = Restrict(coordinatedV,NREM_B);    coordinatedR_V = Restrict(coordinatedV,NREM_R);    coordinatedA_V = Restrict(coordinatedV,NREM_A);
        
        load('coordinated_ripple_bursts.mat')
        
        %% Spikes
        %Load Units
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
        
        %% Definition of cell type population to analize
        if celltype == 1
            tipoD = group_dHPC(:,3);
            tipoV = group_vHPC(:,3);
        elseif celltype == 2
            tipoD = group_dHPC(:,4);
            tipoV = group_vHPC(:,4);
        end
        
        %% Poisson test for Ripple modulation estimation
    %dHPC SU - dRipples All
    base = InvertIntervals([ripplesD(:,1)-0.1, ripplesD(:,3)+0.1],NREM(:,1) , NREM(:,2));
    for i = 1 : length(group_dHPC)% ventral SU
        if logical(tipoD(i))
            y = ripplesD(:,2);
            x =  Restrict(spks_dHPC(spks_dHPC(:,1) == group_dHPC(i,1),2),NREM);
            Base = Restrict(x,base);

            
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
    
    %vHPC SU - vRipples All
    base = InvertIntervals([ripplesV(:,1)-0.1, ripplesV(:,3)+0.1],NREM(:,1) , NREM(:,2));
    for i = 1 : length(group_vHPC)% ventral SU
        if logical(tipoV(i))
            y = ripplesV(:,2);
            x =  Restrict(spks_vHPC(spks_vHPC(:,1) == group_vHPC(i,1),2),NREM);
            Base = Restrict(x,base);
         
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
    
    %% PHIST sourrouding dorsal ripples
        %Dorsal SU - Dorsal ripples
        event1 = ripplesD;
        event1_NR = InvertIntervals([event1(:,1)-0.05 event1(:,3)+0.05] , baselineTS(1)/1000 , baselineTS(2)/1000);
        event1(isnan(event1(:,2)),:) = [];
        for i = 1 : length(group_dHPC)% ventral SU
            if logical(tipoD(i))
                spks = spks_dHPC(spks_dHPC(:,1)==group_dHPC(i,1),:);
                spks = spks(:,2);
                [s,ids,groups] = CCGParameters(event1(:,2),ones(length(event1(:,2)),1),spks,ones(length(spks),1)*2);
                [ccg,T] = CCG(s,ids,'binSize',b,'duration',dd,'smooth',0,'mode','ccg');
                %calculation of mean FRss
                m = (length(Restrict(spks,event1_NR)) / sum(event1_NR(:,2)-event1_NR(:,1)));
                dHPC_neurons_ripples = [dHPC_neurons_ripples , ((ccg(:,1,2)./length(event1))./b)./m];
                clear s ids groups ccg T z m tmp ins
            end
        end
        clear i event1 event1_NR
        
        event1 = ripplesV;
        event1_NR = InvertIntervals([event1(:,1)-0.05 event1(:,3)+0.05] , baselineTS(1)/1000 , baselineTS(2)/1000);
        event1(isnan(event1(:,2)),:) = [];
        %Ventral SU - Dorsal ripples
        for i = 1 : length(group_vHPC)% ventral SU
            if logical(tipoV(i))
                spks = spks_vHPC(spks_vHPC(:,1)==group_vHPC(i,1),:);
                spks = spks(:,2);
                [s,ids,groups] = CCGParameters(event1(:,2),ones(length(event1(:,2)),1),spks,ones(length(spks),1)*2);
                [ccg,T] = CCG(s,ids,'binSize',b,'duration',dd,'smooth',0,'mode','ccg');
                %calculation of mean FR
                m = (length(Restrict(spks,event1_NR)) / sum(event1_NR(:,2)-event1_NR(:,1)));
                vHPC_neurons_ripples = [vHPC_neurons_ripples , ((ccg(:,1,2)./length(event1))./b)./m];
                clear s ids groups ccg z m tmp ins T
            end
            clear spks rate
        end
        clear i event1 event1_NR
        
        %Dorsal SU - Dorsal ripples
        event1 = Restrict(ripplesD,NREM_B);
        event1_NR = InvertIntervals([event1(:,1)-0.05 event1(:,3)+0.05] , baselineTS(1)/1000 , baselineTS(2)/1000);
%         event1 = coordinatedB;
        event1(isnan(event1(:,2)),:) = [];
        
        event2 = Restrict(ripplesD,NREM_R);
        event2_NR = InvertIntervals([event2(:,1)-0.05 event2(:,3)+0.05] ,  rewardTS(1)/1000 , rewardTS(2)/1000);
%         event2 = coordinatedR;
        event2(isnan(event2(:,2)),:) = [];
        
        event3 = Restrict(ripplesD,NREM_A);
        event3_NR = InvertIntervals([event3(:,1)-0.05 event3(:,3)+0.05] , aversiveTS(1)/1000 , aversiveTS(2)/1000);
%         event3 = coordinatedA;
        event3(isnan(event3(:,2)),:) = [];
        for i = 1 : length(group_dHPC)% ventral SU
            if logical(tipoD(i))
                spks = spks_dHPC(spks_dHPC(:,1)==group_dHPC(i,1),:);
                spks = spks(:,2);
                % Baseline
                [s,ids,groups] = CCGParameters(event1(:,2),ones(length(event1(:,2)),1),spks,ones(length(spks),1)*2);
                [ccg,T] = CCG(s,ids,'binSize',b,'duration',dd,'smooth',ss,'mode','ccg');
                %calculation of mean FRss
                m = (length(Restrict(spks,event1_NR)) / sum(event1_NR(:,2)-event1_NR(:,1)));
                dHPC_neurons_ripples_B = [dHPC_neurons_ripples_B , ((ccg(:,1,2)./length(event1))./b)./m];
                clear s ids groups ccg T z m tmp ins
                
                % Reward
                [s,ids,groups] = CCGParameters(event2(:,2),ones(length(event2(:,2)),1),spks,ones(length(spks),1)*2);
                [ccg,T] = CCG(s,ids,'binSize',b,'duration',dd,'smooth',ss,'mode','ccg');
                %calculation of mean FR
                m = (length(Restrict(spks,event2_NR)) / sum(event2_NR(:,2)-event2_NR(:,1)));
                dHPC_neurons_ripples_R = [dHPC_neurons_ripples_R , ((ccg(:,1,2)./length(event2))./b)./m];
                clear s ids groups ccg T z m tmp ins
                
                % Aversive
                [s,ids,groups] = CCGParameters(event3(:,2),ones(length(event3(:,2)),1),spks,ones(length(spks),1)*2);
                [ccg,T] = CCG(s,ids,'binSize',b,'duration',dd,'smooth',ss,'mode','ccg');
                %calculation of mean FR
                m = (length(Restrict(spks,event3_NR)) / sum(event3_NR(:,2)-event3_NR(:,1)));
                dHPC_neurons_ripples_A = [dHPC_neurons_ripples_A , ((ccg(:,1,2)./length(event3))./b)./m];
                clear s ids groups ccg z m tmp ins spks
            end
        end
        clear i event1 event2 event3
        clear event1_NR event2_NR event3_NR
        
        event1 = Restrict(ripplesV,NREM_B);
        event1_NR = InvertIntervals([event1(:,1)-0.05 event1(:,3)+0.05] , baselineTS(1)/1000 , baselineTS(2)/1000);
%         event1 = coordinatedB_V;
        event1(isnan(event1(:,2)),:) = [];
        
        event2 = Restrict(ripplesV,NREM_R);
        event2_NR = InvertIntervals([event2(:,1)-0.05 event2(:,3)+0.05] ,  rewardTS(1)/1000 , rewardTS(2)/1000);
%         event2 = coordinatedR_V;
        event2(isnan(event2(:,2)),:) = [];
        
        event3 = Restrict(ripplesV,NREM_A);
        event3_NR = InvertIntervals([event3(:,1)-0.05 event3(:,3)+0.05] , aversiveTS(1)/1000 , aversiveTS(2)/1000);
%         event3 = coordinatedA_V;
        event3(isnan(event3(:,2)),:) = [];
        
        %Ventral SU - Dorsal ripples
        for i = 1 : length(group_vHPC)% ventral SU
            if logical(tipoV(i))
                spks = spks_vHPC(spks_vHPC(:,1)==group_vHPC(i,1),:);
                spks = spks(:,2);
                % Baseline
                [s,ids,groups] = CCGParameters(event1(:,2),ones(length(event1(:,2)),1),spks,ones(length(spks),1)*2);
                [ccg,T] = CCG(s,ids,'binSize',b,'duration',dd,'smooth',ss,'mode','ccg');
                %calculation of mean FR
                m = (length(Restrict(spks,event1_NR)) / sum(event1_NR(:,2)-event1_NR(:,1)));
                vHPC_neurons_ripples_B = [vHPC_neurons_ripples_B , ((ccg(:,1,2)./length(event1))./b)./m];
                clear s ids groups ccg z m tmp ins T
                
                % Reward
                [s,ids,groups] = CCGParameters(event2(:,2),ones(length(event2(:,2)),1),spks,ones(length(spks),1)*2);
                [ccg,T] = CCG(s,ids,'binSize',b,'duration',dd,'smooth',ss,'mode','ccg');
                %calculation of mean FR
                m = (length(Restrict(spks,event2_NR)) / sum(event2_NR(:,2)-event2_NR(:,1)));
                vHPC_neurons_ripples_R = [vHPC_neurons_ripples_R , ((ccg(:,1,2)./length(event2))./b)./m];
                clear s ids groups ccg z m tmp ins T
                
                % Aversive
                [s,ids,groups] = CCGParameters(event3(:,2),ones(length(event3(:,2)),1),spks,ones(length(spks),1)*2);
                [ccg,T] = CCG(s,ids,'binSize',b,'duration',dd,'smooth',ss,'mode','ccg');
                %calculation of mean FR
                m = (length(Restrict(spks,event3_NR)) / sum(event3_NR(:,2)-event3_NR(:,1)));
                vHPC_neurons_ripples_A = [vHPC_neurons_ripples_A , ((ccg(:,1,2)./length(event3))./b)./m];
                clear s ids groups ccg z m tmp ins Spks
            end
            clear spks rate
        end
        clear i event1 event2 event3
        clear event1_NR event2_NR event3_NR    
        
        %% For coordinated burst of activity
        %All
        %Dorsal
        event1 = coordinated_ripple_bursts;
        event1_NR = InvertIntervals([event1(:,1)-0.05 event1(:,3)+0.05] , baselineTS(1)/1000 , baselineTS(2)/1000);
        event1(isnan(event1(:,2)),:) = [];
        for i = 1 : length(group_dHPC)% ventral SU
            if logical(tipoD(i))
                spks = spks_dHPC(spks_dHPC(:,1)==group_dHPC(i,1),:);
                spks = spks(:,2);
                [s,ids,groups] = CCGParameters(event1(:,2),ones(length(event1(:,2)),1),spks,ones(length(spks),1)*2);
                [ccg,T] = CCG(s,ids,'binSize',b,'duration',dd,'smooth',0,'mode','ccg');
                %calculation of mean FRss
                m = (length(Restrict(spks,event1_NR)) / sum(event1_NR(:,2)-event1_NR(:,1)));
                dHPC_neurons_bursts = [dHPC_neurons_bursts , ((ccg(:,1,2)./length(event1))./b)./m];
                clear s ids groups ccg T z m tmp ins
            end
        end
        clear i event1 event1_NR    
        
        event1 = coordinated_ripple_bursts;
        event1_NR = InvertIntervals([event1(:,1)-0.05 event1(:,3)+0.05] , baselineTS(1)/1000 , baselineTS(2)/1000);
        event1(isnan(event1(:,2)),:) = [];
        %Ventral SU - Dorsal ripples
        for i = 1 : length(group_vHPC)% ventral SU
            if logical(tipoV(i))
                spks = spks_vHPC(spks_vHPC(:,1)==group_vHPC(i,1),:);
                spks = spks(:,2);
                [s,ids,groups] = CCGParameters(event1(:,2),ones(length(event1(:,2)),1),spks,ones(length(spks),1)*2);
                [ccg,T] = CCG(s,ids,'binSize',b,'duration',dd,'smooth',0,'mode','ccg');
                %calculation of mean FR
                m = (length(Restrict(spks,event1_NR)) / sum(event1_NR(:,2)-event1_NR(:,1)));
                vHPC_neurons_bursts = [vHPC_neurons_bursts , ((ccg(:,1,2)./length(event1))./b)./m];
                clear s ids groups ccg z m tmp ins T
            end
            clear spks rate
        end
        clear i event1 event1_NR                
        
        % per conditions
        event1 = Restrict(coordinated_ripple_bursts,NREM_B);
        event1_NR = InvertIntervals([event1(:,1)-0.05 event1(:,3)+0.05] , baselineTS(1)/1000 , baselineTS(2)/1000);
        event1(isnan(event1(:,2)),:) = [];
        
        event2 = Restrict(coordinated_ripple_bursts,NREM_R);
        event2_NR = InvertIntervals([event2(:,1)-0.05 event2(:,3)+0.05] ,  rewardTS(1)/1000 , rewardTS(2)/1000);
        event2(isnan(event2(:,2)),:) = [];
        
        event3 = Restrict(coordinated_ripple_bursts,NREM_A);
        event3_NR = InvertIntervals([event3(:,1)-0.05 event3(:,3)+0.05] , aversiveTS(1)/1000 , aversiveTS(2)/1000);
        event3(isnan(event3(:,2)),:) = [];
        
        %Dorsal SU
        for i = 1 : length(group_dHPC)
            if logical(tipoD(i))
                spks = spks_dHPC(spks_dHPC(:,1)==group_dHPC(i,1),:);
                spks = spks(:,2);
                % Baseline
                [s,ids,groups] = CCGParameters(event1(:,2),ones(length(event1(:,2)),1),spks,ones(length(spks),1)*2);
                [ccg,T] = CCG(s,ids,'binSize',b,'duration',dd,'smooth',ss,'mode','ccg');
                %calculation of mean FRss
                m = (length(Restrict(spks,event1_NR)) / sum(event1_NR(:,2)-event1_NR(:,1)));
                dHPC_neurons_bursts_B = [dHPC_neurons_bursts_B , ((ccg(:,1,2)./length(event1))./b)./m];
                clear s ids groups ccg T z m tmp ins
                
                % Reward
                [s,ids,groups] = CCGParameters(event2(:,2),ones(length(event2(:,2)),1),spks,ones(length(spks),1)*2);
                [ccg,T] = CCG(s,ids,'binSize',b,'duration',dd,'smooth',ss,'mode','ccg');
                %calculation of mean FR
                m = (length(Restrict(spks,event2_NR)) / sum(event2_NR(:,2)-event2_NR(:,1)));
                dHPC_neurons_bursts_R = [dHPC_neurons_bursts_R , ((ccg(:,1,2)./length(event2))./b)./m];
                clear s ids groups ccg T z m tmp ins
                
                % Aversive
                [s,ids,groups] = CCGParameters(event3(:,2),ones(length(event3(:,2)),1),spks,ones(length(spks),1)*2);
                [ccg,T] = CCG(s,ids,'binSize',b,'duration',dd,'smooth',ss,'mode','ccg');
                %calculation of mean FR
                m = (length(Restrict(spks,event3_NR)) / sum(event3_NR(:,2)-event3_NR(:,1)));
                dHPC_neurons_bursts_A = [dHPC_neurons_bursts_A , ((ccg(:,1,2)./length(event3))./b)./m];
                clear s ids groups ccg z m tmp ins spks
            end
        end        
        
        %Ventral SU
        for i = 1 : length(group_vHPC)% ventral SU
            if logical(tipoV(i))
                spks = spks_vHPC(spks_vHPC(:,1)==group_vHPC(i,1),:);
                spks = spks(:,2);
                % Baseline
                [s,ids,groups] = CCGParameters(event1(:,2),ones(length(event1(:,2)),1),spks,ones(length(spks),1)*2);
                [ccg,T] = CCG(s,ids,'binSize',b,'duration',dd,'smooth',ss,'mode','ccg');
                %calculation of mean FR
                m = (length(Restrict(spks,event1_NR)) / sum(event1_NR(:,2)-event1_NR(:,1)));
                vHPC_neurons_bursts_B = [vHPC_neurons_bursts_B , ((ccg(:,1,2)./length(event1))./b)./m];
                clear s ids groups ccg z m tmp ins T
                
                % Reward
                [s,ids,groups] = CCGParameters(event2(:,2),ones(length(event2(:,2)),1),spks,ones(length(spks),1)*2);
                [ccg,T] = CCG(s,ids,'binSize',b,'duration',dd,'smooth',ss,'mode','ccg');
                %calculation of mean FR
                m = (length(Restrict(spks,event2_NR)) / sum(event2_NR(:,2)-event2_NR(:,1)));
                vHPC_neurons_bursts_R = [vHPC_neurons_bursts_R , ((ccg(:,1,2)./length(event2))./b)./m];
                clear s ids groups ccg z m tmp ins T
                
                % Aversive
                [s,ids,groups] = CCGParameters(event3(:,2),ones(length(event3(:,2)),1),spks,ones(length(spks),1)*2);
                [ccg,T] = CCG(s,ids,'binSize',b,'duration',dd,'smooth',ss,'mode','ccg');
                %calculation of mean FR
                m = (length(Restrict(spks,event3_NR)) / sum(event3_NR(:,2)-event3_NR(:,1)));
                vHPC_neurons_bursts_A = [vHPC_neurons_bursts_A , ((ccg(:,1,2)./length(event3))./b)./m];
                clear s ids groups ccg z m tmp ins Spks
            end
            clear spks rate
        end
        clear i event1 event2 event3
        clear event1_NR event2_NR event3_NR        
        t
    end
    tt
end


%% PHIST ripples-SU
%All ripples pooled
% vHPC
[~ , i] = min(abs(T-(-0.05)));
[~ , ii] = min(abs(T-0.05));

x = [];
xx = [];
xxx = [];

for iii = 1 : size(vHPC_neurons_ripples_A,2)
    if ~isnan(sum(vHPC_neurons_ripples(:,iii)))
            x = [x , vHPC_neurons_ripples(:,iii)];
    end
end

y = mean(x(i:ii,:),1);
[~ , y] = sort(y','descend');

z = [];
for i = 1:size(x,2)
    z = [z , x(:,y(i))];
    clear m mm mmm
end

figure
subplot(122), imagesc(T, [1:1:size(x,2)], z'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),caxis([0 5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default')
title('vHPC SU')


% dHPC
[~ , i] = min(abs(T-(-0.05)));
[~ , ii] = min(abs(T-0.05));

xx = [];

for iii = 1 : size(dHPC_neurons_ripples_A,2)
    if ~isnan(sum(dHPC_neurons_ripples(:,iii)))
            xx = [xx , dHPC_neurons_ripples(:,iii)];
    end
end

y = mean(xx(i:ii,:),1);
[~ , y] = sort(y','descend');

z = [];
for i = 1:size(xx,2)
    z = [z , xx(:,y(i))];
    clear m mm mmm
end

subplot(121), imagesc(T, [1:1:size(xx,2)], z'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),caxis([0 5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default')
title('dHPC SU')

p = mean(x,2);
pp = mean(xx,2);
figure,plot(p),hold on, plot(pp)

%% PHIST burst-SU
%All burst pooled
% vHPC
[~ , i] = min(abs(T-(-0.05)));
[~ , ii] = min(abs(T-0.05));

x = [];
xx = [];
xxx = [];

for iii = 1 : size(vHPC_neurons_bursts,2)
    if ~isnan(sum(vHPC_neurons_bursts(:,iii)))
            x = [x , vHPC_neurons_bursts(:,iii)];
    end
end

y = mean(x(i:ii,:),1);
[~ , y] = sort(y','descend');

z = [];
for i = 1:size(x,2)
    z = [z , x(:,y(i))];
    clear m mm mmm
end

figure
subplot(122), imagesc(T, [1:1:size(x,2)], z'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),caxis([0 5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default')
title('vHPC SU')


% dHPC
[~ , i] = min(abs(T-(-0.05)));
[~ , ii] = min(abs(T-0.05));

xx = [];

for iii = 1 : size(dHPC_neurons_bursts,2)
    if ~isnan(sum(dHPC_neurons_bursts(:,iii)))
            xx = [xx , dHPC_neurons_bursts(:,iii)];
    end
end

y = mean(xx(i:ii,:),1);
[~ , y] = sort(y','descend');

z = [];
for i = 1:size(xx,2)
    z = [z , xx(:,y(i))];
    clear m mm mmm
end

subplot(121), imagesc(T, [1:1:size(xx,2)], z'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),caxis([0 5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default')
title('dHPC SU')

p = mean(x,2);
pp = mean(xx,2);
figure,plot(p),hold on, plot(pp)

%% Split into conditions
% vHPC
[~ , i] = min(abs(T-(-0.2)));
[~ , ii] = min(abs(T-0.2));

x = [];
xx = [];
xxx = [];

for iii = 1 : size(vHPC_neurons_bursts_A,2)
    if and(and(~isnan(sum(vHPC_neurons_bursts_A(:,iii))) , ~isnan(sum(vHPC_neurons_bursts_R(:,iii)))),~isnan(sum(vHPC_neurons_bursts_B(:,iii))))
            x = [x , vHPC_neurons_bursts_B(:,iii)];
            xx = [xx , vHPC_neurons_bursts_R(:,iii)];
            xxx = [xxx , vHPC_neurons_bursts_A(:,iii)];
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
subplot(131), imagesc(T, [1:1:size(x,2)], z'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),caxis([0 5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default')
title('Baseline')
subplot(132), imagesc(T, [1:1:size(xx,2)], zz'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),caxis([0 5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default')
title('Reward')
subplot(133), imagesc(T, [1:1:size(xxx,2)], zzz'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),caxis([0 5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default')
title('Aversive')

figure,
plot(T,mean(z,2),'k'),hold on
plot(T,mean(zz,2),'b'),hold on
plot(T,mean(zzz,2),'r'),hold on

clear x xx xxx z zz zzz


% dHPC
[~ , i] = min(abs(T-(-0.2)));
[~ , ii] = min(abs(T-0.2));

x = [];
xx = [];
xxx = [];

for iii = 1 : size(dHPC_neurons_bursts_A,2)
    if and(and(~isnan(sum(dHPC_neurons_bursts_A(:,iii))) , ~isnan(sum(dHPC_neurons_bursts_R(:,iii)))),~isnan(sum(dHPC_neurons_bursts_B(:,iii))))
            x = [x , dHPC_neurons_bursts_B(:,iii)];
            xx = [xx , dHPC_neurons_bursts_R(:,iii)];
            xxx = [xxx , dHPC_neurons_bursts_A(:,iii)];
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
subplot(131), imagesc(T, [1:1:size(z,2)], x'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),caxis([0 10]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default')
title('Baseline')
subplot(132), imagesc(T, [1:1:size(zz,2)], xx'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),caxis([0 10]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default')
title('Reward')
subplot(133), imagesc(T, [1:1:size(zzz,2)], xxx'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),caxis([0 10]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default')
title('Aversive')


figure,
plot(T,mean(z,2),'k'),hold on
plot(T,mean(zz,2),'b'),hold on
plot(T,mean(zzz,2),'r'),hold on



%% Split into conditions
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

    z = [z , x(:,y(i))];
    clear m mm mmm
    
    zz = [zz , xx(:,y(i))];
    clear m mm mmm
    
    zzz = [zzz , xxx(:,y(i))];
    clear m mm mmm
end

figure
subplot(131), imagesc(T, [1:1:size(x,2)], z'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),caxis([0 5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default')
title('Baseline')
subplot(132), imagesc(T, [1:1:size(xx,2)], zz'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),caxis([0 5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default')
title('Reward')
subplot(133), imagesc(T, [1:1:size(xxx,2)], zzz'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),caxis([0 5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default')
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
    z = [z , x(:,y(i))];
    clear m mm mmm
    
    zz = [zz , xx(:,y(i))];
    clear m mm mmm
    
    zzz = [zzz , xxx(:,y(i))];
    clear m mm mmm
end

figure
subplot(131), imagesc(T, [1:1:size(z,2)], x'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),caxis([0 5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default')
title('Baseline')
subplot(132), imagesc(T, [1:1:size(zz,2)], xx'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),caxis([0 5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default')
title('Reward')
subplot(133), imagesc(T, [1:1:size(zzz,2)], xxx'), axis tight, xline(0,'--','LineWidth',1),xlim([-0.5 0.5]),caxis([0 5]),xlabel('Time(sec)'),ylabel('Neurons(id)'),colormap('default')
title('Aversive')


figure,
plot(T,mean(z,2),'k'),hold on
plot(T,mean(zz,2),'b'),hold on
plot(T,mean(zzz,2),'r'),hold on



tmp = [];
tmpB = [];
tmpR = [];
tmpA = [];

for i = 1:size(CCG_All,2)
    tmp = [tmp,CCG_All(:,i)./sum(CCG_All(:,i))];
    tmpB = [tmpB,CCG_B(:,i)./sum(CCG_B(:,i))];
    tmpR = [tmpR,CCG_R(:,i)./sum(CCG_R(:,i))];
    tmpA = [tmpB,CCG_A(:,i)./sum(CCG_A(:,i))];
end





tmpD = [];
tmpB = [];
tmpR = [];
tmpA = [];
burst_index_D=[];
[o ,ii] = min(abs(ttt-0));
[o ,iii] = min(abs(ttt-0.2));

for i = 1:size(ACG_All_D,2)
    tmpD = [tmpD,ACG_All_D(:,i)./max(ACG_All_D(:,i))];
    tmpB = [tmpB,ACG_BD(:,i)./max(ACG_BD(:,i))];
    tmpR = [tmpR,ACG_RD(:,i)./max(ACG_RD(:,i))];
    tmpA = [tmpA,ACG_AD(:,i)./max(ACG_AD(:,i))];
    
    burst_index_D=[burst_index_D ; max(ACG_BD(ii:iii,i))/mean(ACG_BD(ii:end,i)) , max(ACG_RD(ii:iii,i))/mean(ACG_RD(ii:end,i)) , max(ACG_AD(ii:iii,i))/mean(ACG_AD(ii:end,i))];
end
figure,
subplot(311)
plot(ttt,tmpB,'K')
subplot(312)
plot(ttt,tmpR,'b')
subplot(313)
plot(ttt,tmpA,'r')



tmpV = [];
tmpB = [];
tmpR = [];
tmpA = [];
burst_index_V = [];

for i = 1:size(ACG_All_V,2)
    tmpV = [tmpV,ACG_All_V(:,i)./max(ACG_All_V(:,i))];
    tmpB = [tmpB,ACG_BV(:,i)./max(ACG_BV(:,i))];
    tmpR = [tmpR,ACG_RV(:,i)./max(ACG_RV(:,i))];
    tmpA = [tmpA,ACG_AV(:,i)./max(ACG_AV(:,i))];
    
    burst_index_V=[burst_index_V ; max(ACG_BV(ii:iii,i))/mean(ACG_BV(ii:end,i)) , max(ACG_RV(ii:iii,i))/mean(ACG_RV(ii:end,i)) , max(ACG_AV(ii:iii,i))/mean(ACG_AV(ii:end,i))];
end
figure,
subplot(311)
plot(ttt,tmpB,'K')
subplot(312)
plot(ttt,tmpR,'b')
subplot(313)
plot(ttt,tmpA,'r')


figure,subplot(211),plot(ttt,tmpD,'k'),subplot(212),plot(ttt,tmpV,'k')


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

