clear
clc
close all

%Parameters

%%
% $$\$$

path = 'E:\Rat127\Ephys\pyr';
% path = 'E:\Rat128\Ephys\in_pyr\ready';

%%
fs = 1250; %Sampling frequency

dt = 1/fs;

%List of folders from the path
files = dir(path);
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags);

final = cell(length(subFolders)-2,3);
clear files dirFlags

% Ripples
thresholdsD = [1.5 3]; % thresholds used for Ripple detection [1.5 5]
thresholdsV = [1.5 3]; % thresholds used for Ripple detection [1.5 5]

durations = [30 20 100]; % durations criteria for Ripple detection [Max , Min , Min_interval]
frequencies = [100 200]; % frequency band for ripple detection [100 250]

q = 0.25; %quantile to restrict above it ripples according to their peak amplitude

ch = [74 9 1 23 18 43];

Ripples= cell(length(subFolders)-2,1);

ripples_coordinated_percentage = [];
rateV = []; rateD = []; % to store the ripples rate from dHPC and vHPC
deltaB = []; deltaR = []; deltaA = []; %to store all the time delta beween dorsal-ventral ripples

rate_in_time_V = [];
rate_in_time_D = [];

durationsB_dHPC = []; durationsR_dHPC = []; durationsA_dHPC = [];
amplitudesB_dHPC = []; amplitudesR_dHPC = []; amplitudesA_dHPC = [];


durationsB_vHPC = []; durationsR_vHPC = []; durationsA_vHPC = [];
amplitudesB_vHPC = []; amplitudesR_vHPC = []; amplitudesA_vHPC = [];

TS_ventral_Ripples_B = [];
TS_ventral_Ripples_R = [];
TS_ventral_Ripples_A = [];
coordinated_ventral_ripples = cell(length(subFolders)-2,3);

% for SU
pval = 0.001/2; % p value to define if SU are ripple modulated
binsize = 0.005; % bin size for MU histogram in sec
ss = 2; %smooth level of CCG
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


corrB = []; corrB1 = []; corrB2 = [];
corrA = []; corrA1 = []; corrA2 = [];
corrR = []; corrR1 = []; corrR2 = [];

% cross-correlation
corrB_MU_dHPC = []; corrR_MU_dHPC = []; corrA_MU_dHPC = [];
corrB_MU_vHPC = []; corrR_MU_vHPC = []; corrA_MU_vHPC = [];

% cross-correlograms Rippls-SU
corrB_SU_dHPC = []; corrR_SU_dHPC = []; corrA_SU_dHPC = [];
corrB_SU_vHPC = []; corrR_SU_vHPC = []; corrA_SU_vHPC = [];

% phase
Phase_dHPC = {};
Phase_vHPC = {};

%Loading LFP
for t = 1 : length(subFolders)-2
    session = [subFolders(t+2).folder,'\',subFolders(t+2).name];
    cd(session)
    load([cd,'\lfp.mat'])
    Time = dHPC(:,1);
    
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
    Rewards_filt = Restrict([leftvalve ; rightvalve],rewardTS_run ./1000);
    
    %% Sleep
    x = dir([cd,'\*-states.mat']);    states = load([cd,'\',x.name]);    states = states.states;
    
    REM = ToIntervals(states==5);    NREM = ToIntervals(states==3);    WAKE = ToIntervals(states==1);
    clear x states
    
    % NREM events restriction according conditions
    NREM_B = NREM(NREM(:,2)<baselineTS(1,2)/1000,:);
    NREM_A = NREM(NREM(:,2)>aversiveTS(1,1)/1000 & NREM(:,2)<aversiveTS(1,2)/1000,:);
    NREM_R = NREM(NREM(:,2)>rewardTS(1,1)/1000 & NREM(:,2)<rewardTS(1,2)/1000,:);
    % REM events restriction according conditions
    REM_B = REM(REM(:,2)<baselineTS(1,2)/1000,:);
    REM_A = REM(REM(:,2)>aversiveTS(1,1)/1000 & REM(:,2)<aversiveTS(1,2)/1000,:);
    REM_R = REM(REM(:,2)>rewardTS(1,1)/1000 & REM(:,2)<rewardTS(1,2)/1000,:);
    
    store_ripples = cell(1,6); %variable where I will store the ripples
    [ripplesV1,~,~] = FindRipples(FilterLFP(Restrict(vHPC1,NREM),'passband',frequencies),'thresholds',thresholdsV,'durations', durations);
    lengths = length(ripplesV1);    store_ripples{1} = ripplesV1;
    
    [ripplesV2,~,~] = FindRipples(FilterLFP(Restrict(vHPC2,NREM),'passband',frequencies),'thresholds',thresholdsV,'durations', durations);
    lengths = [lengths , length(ripplesV2)];    store_ripples{2} = ripplesV2;
    
    [ripplesV3,~,~] = FindRipples(FilterLFP(Restrict(vHPC3,NREM),'passband',frequencies),'thresholds',thresholdsV,'durations', durations);
    lengths = [lengths , length(ripplesV3)];    store_ripples{3} = ripplesV3;
    
    [ripplesV4,~,~] = FindRipples(FilterLFP(Restrict(vHPC4,NREM),'passband',frequencies),'thresholds',thresholdsV,'durations', durations);
    lengths = [lengths , length(ripplesV4)];    store_ripples{4} = ripplesV4;
    
    [ripplesV5,~,~] = FindRipples(FilterLFP(Restrict(vHPC5,NREM),'passband',frequencies),'thresholds',thresholdsV,'durations', durations);
    lengths = [lengths , length(ripplesV5)];    store_ripples{5} = ripplesV5;
    
    [ripplesD,~,~] = FindRipples(FilterLFP(Restrict(dHPC,NREM),'passband',frequencies),'thresholds',thresholdsD,'durations', durations);
    store_ripples{6} = ripplesD;
    
    [~,ind] = max(lengths);
    
    ripplesV = store_ripples{ind};    ripplesD = store_ripples{6};    Ripples{t} = store_ripples;
    
%     % Saveing the Ripples
%     SaveRippleEvents(['Rat127-20221109.VRP.evt'],ripplesV,ch(ind),'overwrite','on')
%     SaveRippleEvents(['Rat127-20221109.DRP.evt'],ripplesD,ch(1),'overwrite','on')
%     
   ripplesV = ripplesV(ripplesV(:,4) > quantile(ripplesV(:,4),q) , :);
   ripplesD = ripplesD(ripplesD(:,4) > quantile(ripplesD(:,4),q) , :);
    
    % Coordinated dHPC ripples
    coordinated = [];
    coordinatedV = [];
    for i = 1:length(ripplesD)
        r = ripplesD(i,:);
        tmp = sum(and(ripplesV(:,2)>= r(1,1)-0.05, ripplesV(:,2)<= r(1,3)+0.05));
        if tmp>0
            coordinatedV = [coordinatedV ; ripplesV(and(ripplesV(:,2)>= r(1,1)-0.05, ripplesV(:,2)<= r(1,3)+0.05),2)];
            coordinated = [coordinated ; r(1,2)];
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
    
    %% Spikes
    %Load Units
    cd 'Spikesorting'
    spks = double([readNPY('spike_times.npy') readNPY('spike_clusters.npy')]);
    K = tdfread('cluster_group.tsv'); % import clusters ID and groups (MUA,Noise,Good)
    Kinfo = tdfread('cluster_info.tsv'); % import information of clusters
    K = [K.cluster_id(K.group(:,1) == 'g') , Kinfo.ch(K.group(:,1) == 'g')]; % defining good clusters
    
    % Load neuronal classification
    load([cd,'\SU_classification'])
    
    %Loop to select dorsal or ventral LFP and SU
    % z=1 --> dorsal
    % z=2 --> ventral
    for z = 1:2
        if z == 1
            %             group_dHPC = K(K(:,2) > cutting , 1); % keep only dorsal hippocampal SU.
            n_SU_D = n_SU_D + length(group_dHPC);
            spks_dHPC = spks(ismember(spks(:,2),group_dHPC),:); %keep spks from good clusters
        else
            %             group_vHPC = K(K(:,2) <= cutting , 1); % keep only dorsal hippocampal SU.
            n_SU_V = n_SU_V + length(group_vHPC);
            spks_vHPC = spks(ismember(spks(:,2),group_vHPC),:); %keep spks from good clusters
        end
    end
    clear z
    spks_vHPC(:,1) = double(spks_vHPC(:,1))./20000;
    spks_dHPC(:,1) = double(spks_dHPC(:,1))./20000;
    
    %% Peri-event firing rate histogram to awake stimuli
    %Shock
    criteriaV = [];
    for i = 1 : length(group_vHPC)% ventral SU
        spks = spks_vHPC(spks_vHPC(:,2)==group_vHPC(i),:);
        spks = Restrict(spks(:,1),aversiveTS_run./1000);
        rate = (length(spks)/40 ) > 1;
        if rate
            [s,ids,groups] = CCGParameters(Shocks_filt,ones(length(Shocks_filt),1),spks,ones(length(spks),1)*2);
            [ccg,ttt] = CCG(s,ids,'binSize',0.1,'duration',40,'smooth',3,'mode','ccg');
            
            vHPC_neurons_shock = [vHPC_neurons_shock , ccg(:,1,2)./0.1];
            z = zscore(ccg(:,1,2)./0.1);
            vHPC_neurons_shock_zscore = [vHPC_neurons_shock_zscore , z];
            
            [~,ii] = min(abs(ttt-0));
            [~,iii] = min(abs(ttt-1));
            if or(mean(z(ii:iii)) > 3 , mean(z(ii:iii)) < -3)
                criteriaV = [criteriaV ; true];
            else
                criteriaV = [criteriaV ; false];
            end
            clear z ii iii
        else
            vHPC_neurons_shock = [vHPC_neurons_shock , nan(length(vHPC_neurons_shock),1)];
            vHPC_neurons_shock_zscore = [vHPC_neurons_shock_zscore,nan(length(vHPC_neurons_shock),1)];
            criteriaV = [criteriaV ; false];
        end
        clear spks rate
    end
    
    criteriaD = [];
    for i = 1 : length(group_dHPC)% ventral SU
        %         if classificationV(i,4)==1
        spks = spks_dHPC(spks_dHPC(:,2)==group_dHPC(i),:);
        spks = Restrict(spks(:,1),aversiveTS_run./1000);
        rate = (length(spks)/1 ) > 10;
        if rate
            [s,ids,groups] = CCGParameters(Shocks_filt,ones(length(Shocks_filt),1),spks,ones(length(spks),1)*2);
            [ccg,ttt] = CCG(s,ids,'binSize',0.1,'duration',40,'smooth',3,'mode','ccg');
            vHPC_neurons_shock = [vHPC_neurons_shock , ccg(:,1,2)./0.1];
            z = zscore(ccg(:,1,2)./0.1);
            
            dHPC_neurons_shock_zscore = [dHPC_neurons_shock_zscore , z];
            
            [~,ii] = min(abs(ttt - 0));
            [~,iii] = min(abs(ttt - 1));
            if or(mean(z(ii:iii)) > 3 , mean(z(ii:iii)) < -3)
                criteriaD = [criteriaD ; true];
            else
                criteriaD = [criteriaD ; false];
            end
            clear z ii iii
        else
            dHPC_neurons_shock = [dHPC_neurons_shock , nan(length(dHPC_neurons_shock),1)];
            dHPC_neurons_shock_zscore = [dHPC_neurons_shock_zscore , nan(length(dHPC_neurons_shock),1)];
            criteriaD = [criteriaD ; false];
        end
        clear spks rate
    end
    
    %Valve
    tmp = [];
    for i = 1 : length(group_vHPC)% ventral SU
        spks = spks_vHPC(spks_vHPC(:,2)==group_vHPC(i),:);
        spks = Restrict(spks(:,1),rewardTS_run./1000);
        rate = (length(spks)/40 ) > 1;
        if rate
            [s,ids,groups] = CCGParameters(Rewards_filt(:,1),ones(length(Rewards_filt(:,1)),1),spks,ones(length(spks),1)*2);
            [ccg,ttt] = CCG(s,ids,'binSize',0.1,'duration',40,'smooth',3,'mode','ccg');
%             figure,bar(ttt,ccg(:,1,2))
            vHPC_neurons_valve = [vHPC_neurons_valve , ccg(:,1,2)./0.1];
            z = zscore(ccg(:,1,2)./0.1);
            vHPC_neurons_valve_zscore = [vHPC_neurons_valve_zscore , z];
            
            [~,ii] = min(abs(ttt-(-1)));
            [~,iii] = min(abs(ttt-1));
            if or(mean(z(ii:iii)) > 3 , mean(z(ii:iii)) < -3)
                tmp = [tmp ; true];
            else
                tmp = [tmp ; false];
            end
            clear z ii iii
        else
            vHPC_neurons_valve = [vHPC_neurons_valve , nan(length(vHPC_neurons_valve),1)];
            vHPC_neurons_valve_zscore = [vHPC_neurons_valve_zscore,nan(length(vHPC_neurons_valve),1)];
            tmp = [tmp ; false];
        end
        clear spks rate
    end
    criteriaV = [criteriaV , tmp];
    clear tmp
    
    tmp = [];
    for i = 1 : length(group_dHPC)% ventral SU
        spks = spks_dHPC(spks_dHPC(:,2)==group_dHPC(i),:);
        spks = Restrict(spks(:,1),rewardTS_run./1000);
        rate = (length(spks)/40 ) > 1;
        if rate
            [s,ids,groups] = CCGParameters(Rewards_filt(:,1),ones(length(Rewards_filt(:,1)),1),spks,ones(length(spks),1)*2);
            [ccg,ttt] = CCG(s,ids,'binSize',0.1,'duration',40,'smooth',3,'mode','ccg');
%             figure,bar(ttt,ccg(:,1,2))
            dHPC_neurons_valve = [dHPC_neurons_valve , ccg(:,1,2)./0.1];
            z = zscore(ccg(:,1,2)./0.1);
            dHPC_neurons_valve_zscore = [dHPC_neurons_valve_zscore , z];
            
            [~,ii] = min(abs(ttt-(-1)));
            [~,iii] = min(abs(ttt-1));
            if or(mean(z(ii:iii)) > 3 , mean(z(ii:iii)) < -3)
                tmp = [tmp ; true];
            else
                tmp = [tmp ; false];
            end
            clear z ii iii
        else
            dHPC_neurons_valve = [dHPC_neurons_valve , nan(length(dHPC_neurons_valve),1)];
            dHPC_neurons_valve_zscore = [dHPC_neurons_valve_zscore,nan(length(dHPC_neurons_valve),1)];
            tmp = [tmp ; false];
        end
        clear spks rate
    end
    criteriaD = [criteriaD , tmp];
    clear tmp
    
    %     %% Poisson test
    %     % For selection of SU for further analyses
    %     pIncV = [];
    %     pDecV = [];
    %     surpV = [];
    %     for i = 1 : length(group_vHPC)% ventral Ripples-SU
    %         spks = spks_vHPC(spks_vHPC(:,2) == group_vHPC(i),:);    % Restrict to sws
    %         spks=Restrict(spks,NREM);
    %
    %         % baseline epochs
    %         bufferedripples=[ripplesV(:,1)-0.1 ripplesV(:,3)+0.1];
    %         [baseline,ind]=SubtractIntervals(NREM,bufferedripples,'strict','on');
    %         totalbaselinetime=sum(baseline(:,2)-baseline(:,1));
    %
    %         baselinespikes=Restrict(spks(:,1),baseline);
    %         % Restrict further to the third of ripples with the highest amplitude
    %
    %         totalrippletime=sum(ripplesV(:,3)-ripplesV(:,1));
    %         ripplespikes=Restrict(spks(:,1),[ripplesV(:,1) ripplesV(:,3)]);
    %
    %         ncellbaselinespikes=length(baselinespikes);
    %         ncellripplespikes=length(ripplespikes);
    %         if ncellbaselinespikes~=0 & ncellripplespikes~=0
    %             [pIncV(i) pDecV(i) surpV(i)] = poissonTest(ncellbaselinespikes/totalbaselinetime,ncellripplespikes,totalrippletime);
    %         else
    %             pIncV(i)=NaN;
    %             pDecV(i)=NaN;
    %             surpv(i)=NaN;
    %         end
    %         clear spks ncellbaselinespikes ncellripplespikes totalrippletime ripplespikes
    %         clear baselinespikes totalbaselinetime baseline bufferedripples
    %     end
    %     pIncD = [];
    %     pDecD = [];
    %     surpD = [];
    %     for i = 1 : length(group_dHPC)% ventral Ripples-SU
    %         spks = spks_dHPC(spks_dHPC(:,2) == group_dHPC(i),:);    % Restrict to sws
    %         spks=Restrict(spks,NREM);
    %
    %         % baseline epochs
    %         bufferedripples=[ripplesD(:,1)-0.1 ripplesD(:,3)+0.1];
    %         [baseline,ind]=SubtractIntervals(NREM,bufferedripples,'strict','on');
    %         totalbaselinetime=sum(baseline(:,2)-baseline(:,1));
    %
    %         baselinespikes=Restrict(spks(:,1),baseline);
    %         % Restrict further to the third of ripples with the highest amplitude
    %
    %         totalrippletime=sum(ripplesD(:,3)-ripplesD(:,1));
    %         ripplespikes=Restrict(spks(:,1),[ripplesD(:,1) ripplesD(:,3)]);
    %
    %         ncellbaselinespikes=length(baselinespikes);
    %         ncellripplespikes=length(ripplespikes);
    %         if ncellbaselinespikes~=0 & ncellripplespikes~=0
    %             [pIncD(i) pDecD(i) surpD(i)] = poissonTest(ncellbaselinespikes/totalbaselinetime,ncellripplespikes,totalrippletime);
    %         else
    %             pIncD(i)=NaN;
    %             pDecD(i)=NaN;
    %             surpD(i)=NaN;
    %         end
    %         clear spks ncellbaselinespikes ncellripplespikes totalrippletime ripplespikes
    %         clear baselinespikes totalbaselinetime baseline bufferedripples
    %     end
    %
    %     % To quantify percentage across conditions
    %     pInc = [];
    %     pDec = [];
    %     surp = [];
    %     for i = 1 : length(group_vHPC)% ventral SU
    %         spks1 = spks_vHPC(spks_vHPC(:,2) == group_vHPC(i),:);    % Restrict to sws
    %
    %         for ii = 1:3
    %             if ii ==1 %Baseline sleep
    %                 spks=Restrict(spks1,NREM_B);
    %                 iii = Restrict(ripplesV,NREM_B);
    %                 % baseline epochs
    %                 bufferedripples=[iii(:,1)-0.1 iii(:,3)+0.1];
    %                 [baseline,ind]=SubtractIntervals(NREM_B,bufferedripples,'strict','on');
    %                 totalbaselinetime=sum(baseline(:,2)-baseline(:,1));
    %                 baselinespikes=Restrict(spks(:,1),baseline);
    %                 % Restrict further to the third of ripples with the highest amplitude
    %                 totalrippletime=sum(iii(:,3)-iii(:,1));
    %                 ripplespikes=Restrict(spks(:,1),[iii(:,1) iii(:,3)]);
    %                 ncellbaselinespikes=length(baselinespikes);
    %                 ncellripplespikes=length(ripplespikes);
    %             elseif ii == 2 %Reward sleep
    %                 spks=Restrict(spks1,NREM_R);
    %                 iii = Restrict(ripplesV,NREM_R);
    %                 % baseline epochs
    %                 bufferedripples=[iii(:,1)-0.1 iii(:,3)+0.1];
    %                 [baseline,ind]=SubtractIntervals(NREM_R,bufferedripples,'strict','on');
    %                 totalbaselinetime=sum(baseline(:,2)-baseline(:,1));
    %                 baselinespikes=Restrict(spks(:,1),baseline);
    %                 % Restrict further to the third of ripples with the highest amplitude
    %                 totalrippletime=sum(iii(:,3)-iii(:,1));
    %                 ripplespikes=Restrict(spks(:,1),[iii(:,1) iii(:,3)]);
    %                 ncellbaselinespikes=length(baselinespikes);
    %                 ncellripplespikes=length(ripplespikes);
    %             else %Aversive sleep
    %                 spks=Restrict(spks1,NREM_A);
    %                 iii = Restrict(ripplesV,NREM_A);
    %                 % baseline epochs
    %                 bufferedripples=[iii(:,1)-0.1 iii(:,3)+0.1];
    %                 [baseline,ind]=SubtractIntervals(NREM_A,bufferedripples,'strict','on');
    %                 totalbaselinetime=sum(baseline(:,2)-baseline(:,1));
    %                 baselinespikes=Restrict(spks(:,1),baseline);
    %                 % Restrict further to the third of ripples with the highest amplitude
    %                 totalrippletime=sum(iii(:,3)-iii(:,1));
    %                 ripplespikes=Restrict(spks(:,1),[iii(:,1) iii(:,3)]);
    %                 ncellbaselinespikes=length(baselinespikes);
    %                 ncellripplespikes=length(ripplespikes);
    %             end
    %
    %             if ncellbaselinespikes~=0 & ncellripplespikes~=0
    %                 [pI(ii) pD(ii) sur(ii)] = poissonTest(ncellbaselinespikes/totalbaselinetime,ncellripplespikes,totalrippletime);
    %             else
    %                 pI(ii)=NaN;
    %                 pD(ii)=NaN;
    %                 sur(ii)=NaN;
    %             end
    %             clear spks ncellbaselinespikes ncellripplespikes totalrippletime ripplespikes
    %             clear baselinespikes totalbaselinetime baseline bufferedripples
    %         end
    %         pInc = [pInc ; pI];
    %         pDec = [pDec ; pD];
    %         surp = [surp ;sur];
    %     end
    %     poisson_vHPC_split = [poisson_vHPC_split ; pInc , pDec];
    %     clear pInc pDec surp
    %
    %     pInc = [];
    %     pDec = [];
    %     surp = [];
    %     for i = 1 : length(group_dHPC)% dorsal SU
    %         spks1 = spks_dHPC(spks_dHPC(:,2) == group_dHPC(i),:);    % Restrict to sws
    %
    %         for ii = 1:3
    %             if ii ==1 %Baseline sleep
    %                 spks=Restrict(spks1,NREM_B);
    %                 iii = Restrict(ripplesD,NREM_B);
    %                 % baseline epochs
    %                 bufferedripples=[iii(:,1)-0.1 iii(:,3)+0.1];
    %                 [baseline,ind]=SubtractIntervals(NREM_B,bufferedripples,'strict','on');
    %                 totalbaselinetime=sum(baseline(:,2)-baseline(:,1));
    %                 baselinespikes=Restrict(spks(:,1),baseline);
    %                 % Restrict further to the third of ripples with the highest amplitude
    %                 totalrippletime=sum(iii(:,3)-iii(:,1));
    %                 ripplespikes=Restrict(spks(:,1),[iii(:,1) iii(:,3)]);
    %                 ncellbaselinespikes=length(baselinespikes);
    %                 ncellripplespikes=length(ripplespikes);
    %             elseif ii == 2 %Reward sleep
    %                 spks=Restrict(spks1,NREM_R);
    %                 iii = Restrict(ripplesD,NREM_R);
    %                 % baseline epochs
    %                 bufferedripples=[iii(:,1)-0.1 iii(:,3)+0.1];
    %                 [baseline,ind]=SubtractIntervals(NREM_R,bufferedripples,'strict','on');
    %                 totalbaselinetime=sum(baseline(:,2)-baseline(:,1));
    %                 baselinespikes=Restrict(spks(:,1),baseline);
    %                 % Restrict further to the third of ripples with the highest amplitude
    %                 totalrippletime=sum(iii(:,3)-iii(:,1));
    %                 ripplespikes=Restrict(spks(:,1),[iii(:,1) iii(:,3)]);
    %                 ncellbaselinespikes=length(baselinespikes);
    %                 ncellripplespikes=length(ripplespikes);
    %             else %Aversive sleep
    %                 spks=Restrict(spks1,NREM_A);
    %                 iii = Restrict(ripplesD,NREM_A);
    %                 % baseline epochs
    %                 bufferedripples=[iii(:,1)-0.1 iii(:,3)+0.1];
    %                 [baseline,ind]=SubtractIntervals(NREM_A,bufferedripples,'strict','on');
    %                 totalbaselinetime=sum(baseline(:,2)-baseline(:,1));
    %                 baselinespikes=Restrict(spks(:,1),baseline);
    %                 % Restrict further to the third of ripples with the highest amplitude
    %                 totalrippletime=sum(iii(:,3)-iii(:,1));
    %                 ripplespikes=Restrict(spks(:,1),[iii(:,1) iii(:,3)]);
    %                 ncellbaselinespikes=length(baselinespikes);
    %                 ncellripplespikes=length(ripplespikes);
    %             end
    %
    %             if ncellbaselinespikes~=0 & ncellripplespikes~=0
    %                 [pI(ii) pD(ii) sur(ii)] = poissonTest(ncellbaselinespikes/totalbaselinetime,ncellripplespikes,totalrippletime);
    %             else
    %                 pI(ii)=NaN;
    %                 pD(ii)=NaN;
    %                 sur(ii)=NaN;
    %             end
    %             clear spks ncellbaselinespikes ncellripplespikes totalrippletime ripplespikes
    %             clear baselinespikes totalbaselinetime baseline bufferedripples
    %         end
    %         pInc = [pInc ; pI];
    %         pDec = [pDec ; pD];
    %         surp = [surp ;sur];
    %     end
    %     poisson_dHPC_split = [poisson_dHPC_split ; pInc , pDec];
    %     clear pInc pDec surp
    %
    %     % Determination of dHPC selectivity to one specific session
    %     x = poisson_dHPC_split<pval;
    %     countD_exc = [];
    %     for i = 1:length(x)
    %         if x(i,1)==1 && x(i,2)==1 && x(i,3)==1
    %             countD_exc = [countD_exc ; 1];%all
    %         elseif x(i,1)==1 && x(i,2)==0 && x(i,3)==0
    %             countD_exc = [countD_exc ; 2];%baseline
    %         elseif x(i,1)==0 && x(i,2)==1 && x(i,3)==0
    %             countD_exc = [countD_exc ; 3];%reward
    %         elseif x(i,1)==1 && x(i,2)==0 && x(i,3)==1
    %             countD_exc = [countD_exc ; 4];%no-reward
    %         elseif x(i,1)==0 && x(i,2)==0 && x(i,3)==1
    %             countD_exc = [countD_exc ; 5];%aversive
    %         elseif x(i,1)==1 && x(i,2)==1 && x(i,3)==0
    %             countD_exc = [countD_exc ; 6];%no-aversive
    %         elseif x(i,1)==0 && x(i,2)==1 && x(i,3)==1
    %             countD_exc = [countD_exc ; 7];%just emotional
    %         elseif x(i,1)==0 && x(i,2)==0 && x(i,3)==0
    %             countD_exc = [countD_exc ; 8];%non
    %         end
    %     end
    %     clear x
    %
    %     x = poisson_vHPC_split<pval;
    %     countV_exc = [];
    %     for i = 1:length(x)
    %         if x(i,1)==1 && x(i,2)==1 && x(i,3)==1
    %             countV_exc = [countV_exc ; 1];%all
    %         elseif x(i,1)==1 && x(i,2)==0 && x(i,3)==0
    %             countV_exc = [countV_exc ; 2];%baseline
    %         elseif x(i,1)==0 && x(i,2)==1 && x(i,3)==0
    %             countV_exc = [countV_exc ; 3];%reward
    %         elseif x(i,1)==1 && x(i,2)==0 && x(i,3)==1
    %             countV_exc = [countV_exc ; 4];%no-reward
    %         elseif x(i,1)==0 && x(i,2)==0 && x(i,3)==1
    %             countV_exc = [countV_exc ; 5];%aversive
    %         elseif x(i,1)==1 && x(i,2)==1 && x(i,3)==0
    %             countV_exc = [countV_exc ; 6];%no-aversive
    %         elseif x(i,1)==0 && x(i,2)==1 && x(i,3)==1
    %             countV_exc = [countV_exc ; 7];%just emotional
    %         elseif x(i,1)==0 && x(i,2)==0 && x(i,3)==0
    %             countV_exc = [countV_exc ; 8];%non
    %         end
    %     end
    %     clear x
    %
    %% Ripples - SU cross-correlations
    criteria = criteriaV(:,1)-criteriaV(:,2);
    group_vHPC_filt = group_vHPC(criteria==1);
    xB = FilterLFP(Restrict(vHPC1,baselineTS./1000),'passband',[6 10]);
    xR = FilterLFP(Restrict(vHPC1,rewardTS./1000),'passband',[6 10]);
    xA = FilterLFP(Restrict(vHPC1,aversiveTS./1000),'passband',[6 10]);
    if ~isempty(group_vHPC_filt)
        for i = 1 : length(group_vHPC_filt)% ventral Ripples-SU
            %         if countV_exc(i) == 1
            %         if and(or(pIncV(i)<pval,pDecV(i)<pval),~isnan(and(pIncV(i),pDecV(i))))
            %         if classificationV(i,4)==1
            spks = spks_vHPC(spks_vHPC(:,2)==group_vHPC_filt(i),:);
            z = Restrict(spks(:,1),baselineTS./1000);
            zz = Restrict(spks(:,1),rewardTS./1000);
            zzz = Restrict(spks(:,1),aversiveTS./1000);
            if ~isempty(z) && ~isempty(zz) && ~isempty(zzz)
                % ---- Baseline ----
                                x = Restrict(ripplesV(:,2),baselineTS./1000);
%                 x = Restrict(coordinatedV,baselineTS./1000);
                y = Restrict(spks(:,1),baselineTS./1000);
                duration = sum(NREM_B(:,2)-NREM_B(:,1));
                [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                [ccg,ti] = CCG(s,ids,'binSize',binsize,'duration',2,'smooth',ss,'mode','ccv');%,'totalTime',duration);
                corrB_SU_vHPC = [corrB_SU_vHPC , ccg(:,1,2)];
                FR_B_V = [FR_B_V ; length(y)/duration];
                
                clear s ids groups x ccg y duration p
                
                % ---- Reward ----
                                x = Restrict(ripplesV(:,2),rewardTS./1000);
%                 x = Restrict(coordinatedV,rewardTS./1000);
                y = Restrict(spks(:,1),rewardTS./1000);
                duration = sum(NREM_R(:,2)-NREM_R(:,1));
                [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                [ccg,ti] = CCG(s,ids,'binSize',binsize,'duration',2,'smooth',ss,'mode','ccv');%,'totalTime',duration);
                corrR_SU_vHPC = [corrR_SU_vHPC,ccg(:,1,2)];
                FR_R_V = [FR_R_V ; length(y)/duration];
                clear s ids groups x ccg y duration
                
                % ---- Aversive ----
                x = Restrict(ripplesV(:,2),aversiveTS./1000);
                x = Restrict(coordinatedV,aversiveTS./1000);
                y = Restrict(spks(:,1),aversiveTS./1000);
                duration = sum(NREM_A(:,2)-NREM_A(:,1));
                [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                [ccg,ti] = CCG(s,ids,'binSize',binsize,'duration',2,'smooth',ss,'mode','ccv');%,'totalTime',duration);
                corrA_SU_vHPC = [corrA_SU_vHPC,ccg(:,1,2)];
                FR_A_V = [FR_A_V ; length(y)/duration];
                clear s ids groups x ccg y duration
            end
            clear z zz zzz
            %         end
            clear spks
        end
    end
    
    criteria = criteriaD(:,1)-criteriaD(:,2);
    group_dHPC_filt = group_dHPC(criteria==1);
    if ~isempty(group_dHPC_filt)
        %         group_dHPC_filt = group_dHPC(logical(criteriaD));
        for i = 1 : length(group_dHPC_filt)% dorsal Ripples-SU
            %         if countD_exc(i) == 1
            %         if and(or(pIncD(i)<pval,pDecD(i)<pval),~isnan(and(pIncD(i),pDecD(i))))
            %         if classificationD(i,4)==1
            spks = spks_dHPC(spks_dHPC(:,2) == group_dHPC_filt(i),:);
            z = Restrict(spks(:,1),baselineTS./1000);
            zz = Restrict(spks(:,1),rewardTS./1000);
            zzz = Restrict(spks(:,1),aversiveTS./1000);
            if ~isempty(z) && ~isempty(zz) && ~isempty(zzz)
                % ---- Baseline ----
                x = Restrict(ripplesD(:,2),baselineTS./1000);
                %                 x = Restrict(coordinated,baselineTS./1000);
                y = Restrict(spks(:,1),baselineTS./1000);
                duration = sum(NREM_B(:,2)-NREM_B(:,1));
                [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                [ccg,ti] = CCG(s,ids,'binSize',binsize,'duration',2,'smooth',ss,'mode','ccv');%,'totalTime',duration);
                corrB_SU_dHPC = [corrB_SU_dHPC , ccg(:,1,2)];
                FR_B_D = [FR_B_D ; length(y)/duration];
                %             figure,
                %             subplot(3,1,1),bar(ti,ccg(:,1,2)./sum(ccg(:,1,2)),'FaceColor','k'),xlim([-0.06 0.06])
                clear s ids groups x ccg y duration
                
                % ---- Reward ----
                x = Restrict(ripplesD(:,2),rewardTS./1000);
                %                 x = Restrict(coordinated,rewardTS./1000);
                y = Restrict(spks(:,1),rewardTS./1000);
                duration = sum(NREM_R(:,2)-NREM_R(:,1));
                [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                [ccg,ti] = CCG(s,ids,'binSize',binsize,'duration',2,'smooth',ss,'mode','ccv');%,'totalTime',duration);
                corrR_SU_dHPC = [corrR_SU_dHPC,ccg(:,1,2)];
                FR_R_D = [FR_R_D ; length(y)/duration];
                %             subplot(3,1,2),bar(ti,ccg(:,1,2)./sum(ccg(:,1,2)),'FaceColor','b'),xlim([-0.06 0.06])
                clear s ids groups x ccg y duration
                
                % ---- Aversive ----
                x = Restrict(ripplesD(:,2),aversiveTS./1000);
                %                 x = Restrict(coordinated,aversiveTS./1000);
                y = Restrict(spks(:,1),aversiveTS./1000);
                duration = sum(NREM_A(:,2)-NREM_A(:,1));
                [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                [ccg,ti] = CCG(s,ids,'binSize',binsize,'duration',2,'smooth',ss,'mode','ccv');%,'totalTime',duration);
                corrA_SU_dHPC = [corrA_SU_dHPC,ccg(:,1,2)];
                FR_A_D = [FR_A_D ; length(y)/duration];
                %             subplot(3,1,3),bar(ti,ccg(:,1,2)./sum(ccg(:,1,2)),'FaceColor','r'),xlim([-0.06 0.06])
                clear s ids groups x ccg duration
            end
            clear z zz zzz
            %         end
            clear spks
        end
    end
    t
end
    %% REM theta phase locking
    criteria = criteriaV(:,1)-criteriaV(:,2);
    group_vHPC_filt = group_vHPC(criteria==1);
    xB = FilterLFP(Restrict(vHPC1,baselineTS./1000),'passband',[6 10]);
    xR = FilterLFP(Restrict(vHPC1,rewardTS./1000),'passband',[6 10]);
    xA = FilterLFP(Restrict(vHPC1,aversiveTS./1000),'passband',[6 10]);
    if ~isempty(group_vHPC_filt)
        for i = 1 : length(group_vHPC_filt)% ventral Ripples-SU
            spks = spks_vHPC(spks_vHPC(:,2)==group_vHPC_filt(i),:);
            z = Restrict(spks(:,1),baselineTS./1000);
            zz = Restrict(spks(:,1),rewardTS./1000);
            zzz = Restrict(spks(:,1),aversiveTS./1000);
            if ~isempty(z) && ~isempty(zz) && ~isempty(zzz)
                % ---- Baseline ----
                y = Restrict(spks(:,1),baselineTS./1000);
                [phaseB,~,~]=Phase(xB,y);
                clear y
                
                % ---- Reward ----
                y = Restrict(spks(:,1),rewardTS./1000);
                [phaseR,~,~]=Phase(xR,y);
                clear y
                
                % ---- Aversive ----
                y = Restrict(spks(:,1),aversiveTS./1000);
                [phaseA,~,~]=Phase(xA,y);
                clear y
                Phase_vHPC{end+1,1} = phaseB(:,2);
                Phase_vHPC{end,2} = phaseR(:,2);
                Phase_vHPC{end,3} = phaseA(:,2);
                figure,subplot(131),rose(phaseB(:,2))
                subplot(132),rose(phaseR(:,2))
                subplot(133),rose(phaseA(:,2))
                clear phaseB phaseR phaseA
            end
            clear z zz zzz spks
        end
    end    
    
    criteria = criteriaD(:,1)-criteriaD(:,2);
    group_dHPC_filt = group_dHPC(criteria==1);
    xB = FilterLFP(Restrict(dHPC,baselineTS./1000),'passband',[6 10]);
    xR = FilterLFP(Restrict(dHPC,rewardTS./1000),'passband',[6 10]);
    xA = FilterLFP(Restrict(dHPC,aversiveTS./1000),'passband',[6 10]);
    if ~isempty(group_dHPC_filt)
        for i = 1 : length(group_dHPC_filt)% ventral Ripples-SU
            spks = spks_dHPC(spks_dHPC(:,2)==group_dHPC_filt(i),:);
            z = Restrict(spks(:,1),baselineTS./1000);
            zz = Restrict(spks(:,1),rewardTS./1000);
            zzz = Restrict(spks(:,1),aversiveTS./1000);
            if ~isempty(z) && ~isempty(zz) && ~isempty(zzz)
                % ---- Baseline ----
                y = Restrict(spks(:,1),baselineTS./1000);
                [phaseB,~,~]=Phase(xB,y);
                clear y
                
                % ---- Reward ----
                y = Restrict(spks(:,1),rewardTS./1000);
                [phaseR,~,~]=Phase(xR,y);
                clear y
                
                % ---- Aversive ----
                y = Restrict(spks(:,1),aversiveTS./1000);
                [phaseA,~,~]=Phase(xA,y);
                clear y
                Phase_dHPC{end+1,1} = phaseB(:,2);
                Phase_dHPC{end,2} = phaseR(:,2);
                Phase_dHPC{end,3} = phaseA(:,2);
                figure,subplot(131),rose(phaseB(:,2))
                subplot(132),rose(phaseR(:,2))
                subplot(133),rose(phaseA(:,2))
                clear phaseB phaseR phaseA
            end
            clear z zz zzz spks
        end
    end       
    
    t
end
    %% Ripples MU cross-correlograms
    spks = spks_vHPC;
    % ---- Baseline ----
    %     x = Restrict(ripplesV(:,2),baselineTS./1000);
    x = Restrict(coordinatedV,baselineTS./1000);
    y = Restrict(spks(:,1),baselineTS./1000);
    duration = sum(NREM_B(:,2)-NREM_B(:,1));
    [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
    [ccg,ti] = CCG(s,ids,'binSize',binsize,'duration',1,'smooth',2,'mode','ccg','totalTime',duration);
    corrB_MU_vHPC = [corrB_MU_vHPC , ccg(:,1,2)];
    %             figure,
    %             subplot(3,1,1),bar(ti,ccg(:,1,2)./sum(ccg(:,1,2)),'FaceColor','k'),xlim([-0.06 0.06])
    FR_B_V = [FR_B_V ; length(y)/duration];
    clear s ids groups x ccg y duration
    
    % ---- Reward ----
    %     x = Restrict(ripplesV(:,2),rewardTS./1000);
    x = Restrict(coordinatedV,rewardTS./1000);
    y = Restrict(spks(:,1),rewardTS./1000);
    duration = sum(NREM_R(:,2)-NREM_R(:,1));
    [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
    [ccg,ti] = CCG(s,ids,'binSize',binsize,'duration',1,'smooth',2,'mode','ccg','totalTime',duration);
    corrR_MU_vHPC = [corrR_MU_vHPC,ccg(:,1,2)];
    %             subplot(3,1,2),bar(ti,ccg(:,1,2)./sum(ccg(:,1,2)),'FaceColor','b'),xlim([-0.06 0.06])
    FR_R_V = [FR_R_V ; length(y)/duration];
    clear s ids groups x ccg y duration
    
    % ---- Aversive ----
    %     x = Restrict(ripplesV(:,2),aversiveTS./1000);
    x = Restrict(coordinatedV,aversiveTS./1000);
    y = Restrict(spks(:,1),aversiveTS./1000);
    duration = sum(NREM_A(:,2)-NREM_A(:,1));
    [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
    [ccg,ti] = CCG(s,ids,'binSize',binsize,'duration',1,'smooth',2,'mode','ccg','totalTime',duration);
    corrA_MU_vHPC = [corrA_MU_vHPC,ccg(:,1,2)];
    %             subplot(3,1,3),bar(ti,ccg(:,1,2)./sum(ccg(:,1,2)),'FaceColor','r'),xlim([-0.06 0.06])
    FR_A_V = [FR_A_V ; length(y)/duration];
    clear s ids groups x ccg y duration
    clear spks
    
    
    spks = spks_dHPC;
    % ---- Baseline ----
    %     x = Restrict(ripplesD(:,2),baselineTS./1000);
    x = Restrict(coordinated,baselineTS./1000);
    y = Restrict(spks(:,1),baselineTS./1000);
    duration = sum(NREM_B(:,2)-NREM_B(:,1));
    [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
    [ccg,ti] = CCG(s,ids,'binSize',binsize,'duration',1,'smooth',2,'mode','ccg','totalTime',duration);
    corrB_MU_dHPC = [corrB_MU_dHPC , ccg(:,1,2)];
    FR_B_D = [FR_B_D ; length(y)/duration];
    %             figure,
    %             subplot(3,1,1),bar(ti,ccg(:,1,2)./sum(ccg(:,1,2)),'FaceColor','k'),xlim([-0.06 0.06])
    clear s ids groups x ccg y duration
    
    % ---- Reward ----
    %     x = Restrict(ripplesD(:,2),rewardTS./1000);
    x = Restrict(coordinated,rewardTS./1000);
    y = Restrict(spks(:,1),rewardTS./1000);
    duration = sum(NREM_R(:,2)-NREM_R(:,1));
    [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
    [ccg,ti] = CCG(s,ids,'binSize',binsize,'duration',1,'smooth',2,'mode','ccg','totalTime',duration);
    corrR_MU_dHPC = [corrR_MU_dHPC,ccg(:,1,2)];
    FR_R_D = [FR_R_D ; length(y)/duration];
    %             subplot(3,1,2),bar(ti,ccg(:,1,2)./sum(ccg(:,1,2)),'FaceColor','b'),xlim([-0.06 0.06])
    clear s ids groups x ccg y duration
    
    % ---- Aversive ----
    %     x = Restrict(ripplesD(:,2),aversiveTS./1000);
    x = Restrict(coordinated,aversiveTS./1000);
    y = Restrict(spks(:,1),aversiveTS./1000);
    duration = sum(NREM_A(:,2)-NREM_A(:,1));
    [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
    [ccg,ti] = CCG(s,ids,'binSize',binsize,'duration',1,'smooth',2,'mode','ccg','totalTime',duration);
    corrA_MU_dHPC = [corrA_MU_dHPC,ccg(:,1,2)];
    FR_A_D = [FR_A_D ; length(y)/duration];
    %             subplot(3,1,3),bar(ti,ccg(:,1,2)./sum(ccg(:,1,2)),'FaceColor','r'),xlim([-0.06 0.06])
    clear s ids groups x ccg duration
    clear spks
    t
end

%% Plots Ripples-SU coordination
%vHPC neurons
tmp1 = []; tmp2 = []; tmp3 = [];
lagB = []; lagR = []; lagA = [];
CB_SU_vHPC = []; CR_SU_vHPC = []; CA_SU_vHPC = [];
for i = 1:size(corrB_SU_vHPC,2)
    [~,ii] = min(abs(ti - 0.1));
    [~,iii] = min(abs(ti - (-0.1)));
%     if sum(corrB_SU_vHPC(iii:ii,i)) >0 && sum(corrR_SU_vHPC(iii:ii,i))>0 & sum(corrA_SU_vHPC(iii:ii,i))>0
        B = corrB_SU_vHPC(:,i);% ./ sum(corrB_SU_vHPC(:,i));
        R = corrR_SU_vHPC(:,i);% ./ sum(corrR_SU_vHPC(:,i));
        A = corrA_SU_vHPC(:,i);% ./ sum(corrA_SU_vHPC(:,i));
        
        CB_SU_vHPC = [CB_SU_vHPC , B];
        CR_SU_vHPC = [CR_SU_vHPC , R];
        CA_SU_vHPC = [CA_SU_vHPC , A];
        
        [mB,iB] = max(B(iii:ii,1));
        lagB = [lagB ; ti(iB+iii)];
        
        [mR,iR] = max(R(iii:ii,1));
        lagR = [lagR ; ti(iR+iii)];
        
        [mA,iA] = max(A(iii:ii,1));
        lagA = [lagA ; ti(iA+iii)];
        
        tmp1 = [tmp1 ; mB];
        tmp2 = [tmp2 ; mR];
        tmp3 = [tmp3 ; mA];
        m = max([mB , mR , mA]);
        mm = min([min(B(:,1)) , min(R(:,1)) , min(A(:,1))]);
        
        figure,
        subplot(1,3,1),bar(ti,B,'FaceColor','k','FaceAlpha',0.7),xlim([-0.1 0.1])%,ylim([0 m])
        subplot(1,3,2),bar(ti,R,'FaceColor','b','FaceAlpha',0.7),xlim([-0.1 0.1])%,ylim([0 m])
        subplot(1,3,3),bar(ti,A,'FaceColor','r','FaceAlpha',0.7),xlim([-0.1 0.1])%,ylim([0 m])
        clear mB mR mA iB iR iA A R B
%     else
% %         figure
%     end
end
figure,
subplot(1,3,1),bar(ti,mean(CB_SU_vHPC,2),1,'k','LineStyle','none'),xline(0,'--'),xlim([-0.05 0.05])
subplot(1,3,2),bar(ti,mean(CR_SU_vHPC,2),1,'b','LineStyle','none'),xline(0,'--'),xlim([-0.05 0.05])
subplot(1,3,3),bar(ti,mean(CA_SU_vHPC,2),1,'r','LineStyle','none'),xline(0,'--'),xlim([-0.05 0.05])

figure,
boxplot([lagB lagR lagA]),ylim([-0.05 0.05])
[h , p] = ttest(lagB)
[h , p] = ttest(lagR)
[h , p] = ttest(lagA)

figure,
boxplot([tmp1 tmp2 tmp3]),ylim([0 0.04])
[p,table,stats] = anova1([tmp1 tmp2 tmp3]);
[c,m,h,nms] = multcompare(stats,'alpha',0.05,'ctype','bonferroni');



%dHPC neurons
tmp1 = []; tmp2 = []; tmp3 = [];
lagB = []; lagR = []; lagA = [];
CB_SU_dHPC = []; CR_SU_dHPC = []; CA_SU_dHPC = [];
for i = 1:size(corrB_SU_dHPC,2)
    [~,ii] = min(abs(ti - 0.05));
    [~,iii] = min(abs(ti - (-0.05)));
    if sum(corrB_SU_dHPC(iii:ii,i)) >0 && sum(corrR_SU_dHPC(iii:ii,i))>0 & sum(corrA_SU_dHPC(iii:ii,i))>0
        B = corrB_SU_dHPC(:,i);% ./ sum(corrB_SU_dHPC(:,i));
        R = corrR_SU_dHPC(:,i);% ./ sum(corrR_SU_dHPC(:,i));
        A = corrA_SU_dHPC(:,i);% ./ sum(corrA_SU_dHPC(:,i));
        
        CB_SU_dHPC = [CB_SU_dHPC , B];
        CR_SU_dHPC = [CR_SU_dHPC , R];
        CA_SU_dHPC = [CA_SU_dHPC , A];
        
        [mB,iB] = max(B(iii:ii,1));
        lagB = [lagB ; ti(iB+iii)];
        
        [mR,iR] = max(R(iii:ii,1));
        lagR = [lagR ; ti(iR+iii)];
        
        [mA,iA] = max(A(iii:ii,1));
        lagA = [lagA ; ti(iA+iii)];
        
        tmp1 = [tmp1 ; mB];
        tmp2 = [tmp2 ; mR];
        tmp3 = [tmp3 ; mA];
        m = max([mB , mR , mA]);
        mm = min([min(B(:,1)) , min(R(:,1)) , min(A(:,1))]);
        
%         figure,
%         subplot(1,3,1),bar(ti,B,'FaceColor','k','FaceAlpha',0.7),xlim([-0.1 0.1]),ylim([mm-0.001 m])
%         subplot(1,3,2),bar(ti,R,'FaceColor','b','FaceAlpha',0.7),xlim([-0.1 0.1]),ylim([mm-0.001 m])
%         subplot(1,3,3),bar(ti,A,'FaceColor','r','FaceAlpha',0.7),xlim([-0.1 0.1]),ylim([mm-0.001 m])
        clear mB mR mA iB iR iA A R B
    else
%         figure
    end
end
figure,
subplot(1,3,1),bar(ti,mean(CB_SU_dHPC,2),1,'k','LineStyle','none'),xline(0,'--'),ylim([0 0.06]),xlim([-0.1 0.1])
subplot(1,3,2),bar(ti,mean(CR_SU_dHPC,2),1,'b','LineStyle','none'),xline(0,'--'),ylim([0 0.06]),xlim([-0.1 0.1])
subplot(1,3,3),bar(ti,mean(CA_SU_dHPC,2),1,'r','LineStyle','none'),xline(0,'--'),ylim([0 0.06]),xlim([-0.1 0.1])

figure,
boxplot([lagB lagR lagA]),ylim([-0.05 0.05])
[h , p] = ttest(lagB)
[h , p] = ttest(lagR)
[h , p] = ttest(lagA)

figure,
boxplot([tmp1 tmp2 tmp3]),ylim([0 0.015])
[p,table,stats] = anova1([tmp1 tmp2 tmp3]);
c
%% Plots Ripples-MU coordination
%vHPC neurons
tmp1 = []; tmp2 = []; tmp3 = [];
lagB = []; lagR = []; lagA = [];
CB_MU_vHPC = []; CR_MU_vHPC = []; CA_MU_vHPC = [];
for i = 1:size(corrB_MU_vHPC,2)
    [~,ii] = min(abs(ti - 0.05));
    [~,iii] = min(abs(ti - (-0.05)));
    if sum(corrB_MU_vHPC(iii:ii,i)) >0 && sum(corrR_MU_vHPC(iii:ii,i))>0 & sum(corrA_MU_vHPC(iii:ii,i))>0
        B = corrB_MU_vHPC(:,i) ./ sum(corrB_MU_vHPC(:,i));
        R = corrR_MU_vHPC(:,i) ./ sum(corrR_MU_vHPC(:,i));
        A = corrA_MU_vHPC(:,i) ./ sum(corrA_MU_vHPC(:,i));
        
        CB_MU_vHPC = [CB_MU_vHPC , B];
        CR_MU_vHPC = [CR_MU_vHPC , R];
        CA_MU_vHPC = [CA_MU_vHPC , A];
        
        [mB,iB] = max(B(iii:ii,1));
        lagB = [lagB ; ti(iB+iii)];
        
        [mR,iR] = max(R(iii:ii,1));
        lagR = [lagR ; ti(iR+iii)];
        
        [mA,iA] = max(A(iii:ii,1));
        lagA = [lagA ; ti(iA+iii)];
        
        tmp1 = [tmp1 ; mB];
        tmp2 = [tmp2 ; mR];
        tmp3 = [tmp3 ; mA];
        m = max([mB , mR , mA]);
        mm = min([min(B(:,1)) , min(R(:,1)) , min(A(:,1))]);
        
        figure,
        subplot(1,3,1),bar(ti,B,'FaceColor','k','FaceAlpha',0.7),xlim([-0.05 0.05]),ylim([mm-0.001 m])
        subplot(1,3,2),bar(ti,R,'FaceColor','b','FaceAlpha',0.7),xlim([-0.05 0.05]),ylim([mm-0.001 m])
        subplot(1,3,3),bar(ti,A,'FaceColor','r','FaceAlpha',0.7),xlim([-0.05 0.05]),ylim([mm-0.001 m])
        clear mB mR mA iB iR iA A R B
    else
%         figure
    end
end
figure,
subplot(1,3,1),bar(ti,mean(CB_MU_vHPC,2),1,'k','LineStyle','none'),xline(0,'--'),ylim([0.002 0.04]),xlim([-0.05 0.05])
subplot(1,3,2),bar(ti,mean(CR_MU_vHPC,2),1,'b','LineStyle','none'),xline(0,'--'),ylim([0.002 0.04]),xlim([-0.05 0.05])
subplot(1,3,3),bar(ti,mean(CA_MU_vHPC,2),1,'r','LineStyle','none'),xline(0,'--'),ylim([0.002 0.04]),xlim([-0.05 0.05])

figure,
boxplot([lagB lagR lagA]),ylim([-0.05 0.05])
[h , p] = ttest(lagB)
[h , p] = ttest(lagR)
[h , p] = ttest(lagA)

figure,
boxplot([tmp1 tmp2 tmp3]),ylim([0.005 0.025])
[p,table,stats] = anova1([tmp1 tmp2 tmp3]);
[c,m,h,nms] = multcompare(stats,'alpha',0.05,'ctype','bonferroni');



%dHPC neurons
tmp1 = []; tmp2 = []; tmp3 = [];
lagB = []; lagR = []; lagA = [];
CB_MU_dHPC = []; CR_MU_dHPC = []; CA_MU_dHPC = [];
for i = 1:size(corrB_MU_dHPC,2)
    [~,ii] = min(abs(ti - 0.1));
    [~,iii] = min(abs(ti - (-0.1)));
    if sum(corrB_MU_dHPC(iii:ii,i)) >0 && sum(corrR_MU_dHPC(iii:ii,i))>0 & sum(corrA_MU_dHPC(iii:ii,i))>0
        B = corrB_MU_dHPC(:,i) ./ sum(corrB_MU_dHPC(:,i));
        R = corrR_MU_dHPC(:,i) ./ sum(corrR_MU_dHPC(:,i));
        A = corrA_MU_dHPC(:,i) ./ sum(corrA_MU_dHPC(:,i));
        
        CB_MU_dHPC = [CB_MU_dHPC , B];
        CR_MU_dHPC = [CR_MU_dHPC , R];
        CA_MU_dHPC = [CA_MU_dHPC , A];
        
        [mB,iB] = max(B(iii:ii,1));
        lagB = [lagB ; ti(iB+iii)];
        
        [mR,iR] = max(R(iii:ii,1));
        lagR = [lagR ; ti(iR+iii)];
        
        [mA,iA] = max(A(iii:ii,1));
        lagA = [lagA ; ti(iA+iii)];
        
        tmp1 = [tmp1 ; mB];
        tmp2 = [tmp2 ; mR];
        tmp3 = [tmp3 ; mA];
        m = max([mB , mR , mA]);
        mm = min([min(B(:,1)) , min(R(:,1)) , min(A(:,1))]);
        
%         figure,
%         subplot(1,3,1),bar(ti,B,'FaceColor','k','FaceAlpha',0.7),xlim([-0.05 0.05]),ylim([mm-0.001 m])
%         subplot(1,3,2),bar(ti,R,'FaceColor','b','FaceAlpha',0.7),xlim([-0.05 0.05]),ylim([mm-0.001 m])
%         subplot(1,3,3),bar(ti,A,'FaceColor','r','FaceAlpha',0.7),xlim([-0.05 0.05]),ylim([mm-0.001 m])
        clear mB mR mA iB iR iA A R B
    else
%         figure
    end
end
figure,
subplot(1,3,1),bar(ti,mean(CB_MU_dHPC,2),1,'k','LineStyle','none'),xline(0,'--'),ylim([0 0.06]),xlim([-0.05 0.05])
subplot(1,3,2),bar(ti,mean(CR_MU_dHPC,2),1,'b','LineStyle','none'),xline(0,'--'),ylim([0 0.06]),xlim([-0.05 0.05])
subplot(1,3,3),bar(ti,mean(CA_MU_dHPC,2),1,'r','LineStyle','none'),xline(0,'--'),ylim([0 0.06]),xlim([-0.05 0.05])

figure,
boxplot([lagB lagR lagA]),ylim([-0.05 0.05])
[h , p] = ttest(lagB)
[h , p] = ttest(lagR)
[h , p] = ttest(lagA)

figure,
boxplot([tmp1 tmp2 tmp3]),ylim([0 0.08])
[p,table,stats] = anova1([tmp1 tmp2 tmp3]);
[c,m,h,nms] = multcompare(stats,'alpha',0.05,'ctype','bonferroni');

x = poisson_dHPC_split<0.001;
countD_exc = zeros(1,8);
for i = 1:length(x)
    if x(i,1)==1 && x(i,2)==1 && x(i,3)==1
        countD_exc(1) = countD_exc(1)+1;
        %all
    elseif x(i,1)==1 && x(i,2)==0 && x(i,3)==0 
        countD_exc(2) = countD_exc(2)+1;
        %baseline
    elseif x(i,1)==0 && x(i,2)==1 && x(i,3)==0 
        countD_exc(3) = countD_exc(3)+1;
        %reward
    elseif x(i,1)==1 && x(i,2)==0 && x(i,3)==1          
        countD_exc(4) = countD_exc(4)+1;
        %no-reward        
    elseif x(i,1)==0 && x(i,2)==0 && x(i,3)==1 
        countD_exc(5) = countD_exc(5)+1;
        %aversive
    elseif x(i,1)==1 && x(i,2)==1 && x(i,3)==0 
        countD_exc(6) = countD_exc(6)+1;
        %no-aversive        
    elseif x(i,1)==0 && x(i,2)==1 && x(i,3)==1          
        countD_exc(7) = countD_exc(7)+1;
        %just emotional
    elseif x(i,1)==0 && x(i,2)==0 && x(i,3)==0 
        countD_exc(8) = countD_exc(8)+1;
        %non
    end
end
percentage_dHPC_exc = (countD_exc./sum(countD_exc))*100;

x = poisson_vHPC_split<0.001;
countV_exc = zeros(1,8);
for i = 1:length(x)
    if x(i,1)==1 && x(i,2)==1 && x(i,3)==1
        countV_exc(1) = countV_exc(1)+1;
        %all
    elseif x(i,1)==1 && x(i,2)==0 && x(i,3)==0 
        countV_exc(2) = countV_exc(2)+1;
        %baseline
    elseif x(i,1)==0 && x(i,2)==1 && x(i,3)==0 
        countV_exc(3) = countV_exc(3)+1;
        %reward
    elseif x(i,1)==1 && x(i,2)==0 && x(i,3)==1          
        countV_exc(4) = countV_exc(4)+1;
        %no-reward        
    elseif x(i,1)==0 && x(i,2)==0 && x(i,3)==1 
        countV_exc(5) = countV_exc(5)+1;
        %aversive
    elseif x(i,1)==1 && x(i,2)==1 && x(i,3)==0 
        countV_exc(6) = countV_exc(6)+1;
        %no-aversive
    elseif x(i,1)==0 && x(i,2)==1 && x(i,3)==1          
        countV_exc(7) = countV_exc(7)+1;
        %just emotional
    elseif x(i,1)==0 && x(i,2)==0 && x(i,3)==0 
        countV_exc(8) = countV_exc(8)+1;
        %non
    end
end
percentage_vHPC_exc = (countV_exc./sum(countV_exc))*100;


%% PHIST
[~ , i] = min(abs(ttt-0));
[~ , ii] = min(abs(ttt-2));

x = [];
for iii = 1 : size(vHPC_neurons_shock_zscore,2)
    if ~isnan(sum(vHPC_neurons_shock_zscore(:,iii)))
        x = [ x , vHPC_neurons_shock_zscore(:,iii)];
    end
end
y = mean(x(i:ii,:),1);
[~ , y] = sort(y','descend');

z = [];
for i = 1:size(x,2)
    z = [z , x(:,y(i))];
end

figure
subplot(121), imagesc(ttt, [1:1:size(x,2)], z'), axis tight, xline(0,'LineWidth',1),xlim([-5 5]),caxis([-2 10]),xlabel('Time(sec)'),ylabel('Neurons(id)')
xline(1,'--','LineWidth',1)
title('vHPC SU')

[~ , i] = min(abs(ttt-0));
[~ , ii] = min(abs(ttt-2));

x = [];
for iii = 1 : size(dHPC_neurons_shock_zscore,2)
    if ~isnan(sum(dHPC_neurons_shock_zscore(:,iii)))
        x = [ x , dHPC_neurons_shock_zscore(:,iii)];
    end
end
y = mean(x(i:ii,:),1);
[~ , y] = sort(y','descend');

z = [];
for i = 1:size(x,2)
    z = [z , x(:,y(i))];
end

subplot(122), imagesc(ttt, [1:1:size(x,2)], z'), axis tight, xline(0,'LineWidth',1),xlim([-5 5]),caxis([-2 10]),xlabel('Time(sec)'),ylabel('Neurons(id)')
xline(1,'--','LineWidth',1)
title('dHPC SU')

%% PHIST
[~ , i] = min(abs(ttt-(-2)));
[~ , ii] = min(abs(ttt-2));

x = [];
for iii = 1 : size(vHPC_neurons_valve_zscore,2)
    if ~isnan(sum(vHPC_neurons_valve_zscore(:,iii)))
        x = [ x , vHPC_neurons_valve_zscore(:,iii)];
    end
end
y = mean(x(i:ii,:),1);
[~ , y] = sort(y','descend');

z = [];
for i = 1:size(x,2)
    z = [z , x(:,y(i))];
end

figure
subplot(121), imagesc(ttt, [1:1:size(x,2)], z'), axis tight, xline(0,'LineWidth',1),xlim([-5 5]),caxis([-2 10]),xlabel('Time(sec)'),ylabel('Neurons(id)')
xline(1,'--','LineWidth',1)
title('vHPC SU')

[~ , i] = min(abs(ttt-());
[~ , ii] = min(abs(ttt-2));

x = [];
for iii = 1 : size(dHPC_neurons_valve_zscore,2)
    if ~isnan(sum(dHPC_neurons_valve_zscore(:,iii)))
        x = [ x , dHPC_neurons_valve_zscore(:,iii)];
    end
end
y = mean(x(i:ii,:),1);
[~ , y] = sort(y','descend');

z = [];
for i = 1:size(x,2)
    z = [z , x(:,y(i))];
end

subplot(122), imagesc(ttt, [1:1:size(x,2)], z'), axis tight, xline(0,'LineWidth',1),xlim([-5 5]),caxis([-2 10]),xlabel('Time(sec)'),ylabel('Neurons(id)')
xline(1,'--','LineWidth',1)
title('dHPC SU')
