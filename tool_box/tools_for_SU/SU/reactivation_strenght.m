clear
clc
close all

%% Parameters
path = {'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat103\usable'};%List of folders from the path

%Sleep
time_criteria = 600; %time criteria to define the maximal time of sleep to include

% Ripples
thresholdsD = [1.5 5]; % thresholds used for Ripple detection [1.5 5]
thresholdsV = [1.5 5]; % thresholds used for Ripple detection [1.5 5]
durations = [30 20 100]; % durations criteria for Ripple detection [Max , Min , Min_interval]
frequencies = [100 200]; % frequency band for ripple detection [100 250]
q = 0.25; %quantile to restrict above it ripples according to their peak amplitude
ch = [74 9 1 23 18 43];
Ripples= cell(14,1);
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

% for SU
criteria_fr = 0; %criteria to include or not a SU into the analysis
criteria_n = 3; % minimal number of neurons from each structure
pval = 0.001/2; % p value to define if SU are ripple modulated
ss = 2; %smooth level of CCG
n_SU_V = 0;
n_SU_D = 0;
FR_B_V = []; FR_R_V = [];  FR_A_V = []; % Firing Rate during NREM
FR_B_D = []; FR_R_D = [];  FR_A_D = []; % Firing Rate during NREM
poisson_dHPC_split = []; poisson_vHPC_split = []; %poisson results split by conditions

%For EV and REV
binSize = 0.05;
EV_A = []; REV_A = [];
EV_R = []; REV_R = [];

window = [-2 2];
%% Main lo op, to iterate across sessions
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
    
%     load([cd,'\lfp.mat'])
%     Time = dHPC(:,1);
    
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
    clear count deff shock i
    
    %Rewards selection
    Rewards_filt = Restrict([leftvalve ; rightvalve],rewardTS_run ./1000);
    
    % Load laps from DLC data **remember that the SF is 1/30**
    load([cd,'\laps1.mat'],'laps','posx','posy','velocity')
    laps1 = laps;
    pos1 = [posx,posy];
    velocity1 = velocity;
    clear laps posx posy velocity
    
    load([cd,'\laps2.mat'],'laps','posx','posy','velocity')
    laps2 = laps;
    pos2 = [posx,posy];
    velocity2 = velocity;
    clear laps posx posy velocity
    
    %Organizing TS deppending on their condition
    if aversiveTS(1,1) > rewardTS(1,1)
        lapsA = laps2;
        posA = pos2;
        velA = velocity2;
        
        lapsR = laps1;
        posR = pos1;
        velR = velocity1;
    else
        lapsA = laps1;
        posA = pos1;
        velA = velocity1;
        
        lapsR = laps2;
        posR = pos2;
        velR = velocity2;
    end
    
    clear laps1 laps2 pos1 pos2 velocity1 velocity2    
    
    %find the begining of the camara
    [camaraA1,~] = find((camara(:,1)-aversiveTS_run(1)/1000)>0,1,'first'); %TimeStamp of the begining of aversive
    [camaraA2,~] = find((camara(:,1)-aversiveTS_run(2)/1000)<0,1,'last'); %TimeStamp of the ending of aversive
    timeA = ((camara(:,2) - camara(:,1))/2) + camara(:,1); timeA = timeA(camaraA1:camaraA2);
    if length(timeA) == length(posA),   posAA = [timeA posA]; clear timeA
    elseif length(timeA) < length(posA),   posAA = [timeA posA(1:length(timeA),:)]; clear timeA
    else,   posAA = [timeA(1:length(posA),:) posA]; clear timeA
    end %same length for camara TTLs velocity
    
    timeA = camara(camaraA1:camaraA2,:);
    if length(timeA) < length(posA),  posA = posA(1:length(timeA,:));  %clear timeA
    else, timeA = timeA(1:length(posA),:);
    end %same length for camara TTLs velocity
    
    V = LinearVelocity(posAA,3);
    speedA = ToIntervals(V(:,1),V(:,2)>2);
    speedA = speedA((speedA(:,2)-speedA(:,1))>1,:);
    clear posA velA posAA timeA V
    
    %reward
    [camaraR1,~] = find((camara(:,1)-rewardTS(1)/1000)>0,1,'first'); %TimeStamp of the begining of reward
    [camaraR2,~] = find((camara(:,1)-rewardTS(2)/1000)<0,1,'last'); %TimeStamp of the ending of reward
    timeR = ((camara(:,2) - camara(:,1))/2) + camara(:,1); timeR = timeR(camaraR1:camaraR2);
    if length(timeR) == length(posR),   posRR = [timeR posR]; clear timeR
    elseif length(timeR) < length(posR),   posRR = [timeR posR(1:length(timeR),:)]; clear timeR
    else,   posRR = [timeR(1:length(posR),:) posR]; clear timeR
    end %same length for camara TTLs velocity
        
    timeR = camara(camaraR1:camaraR2,:);
    if length(timeR) < length(posR),  posR = posR(1:length(timeR),:);  %clear timeA
    else, timeR = timeR(1:length(posR),:);
    end %same length for camara TTLs velocit
    
    V = LinearVelocity(posRR,3);
    speedR = ToIntervals(V(:,1),V(:,2)>2);
    speedR = speedR((speedR(:,2)-speedR(:,1))>1,:);
    clear posR velR posR timeR V
    
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
    %         load('detected_ripples.mat')
    ripplesD = table2array(readtable('ripplesD_customized2.csv'));
    ripplesV = table2array(readtable('ripplesV_customized2.csv'));
    
    % Coordinated dHPC ripples
    coordinated = [];
    coordinatedV = [];
    uncoordinated = [];
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
    
    % Detection of uncoordinated ripples
    uncoordinated = ripplesD(~ismember(ripplesD(:,1),coordinated(:,1)),:);
    uncoordinatedV = ripplesV(~ismember(ripplesV(:,1),coordinatedV(:,1)),:);
    
    uncoordinatedB = Restrict(uncoordinated,NREM_B);    uncoordinatedA = Restrict(uncoordinated,NREM_A);    uncoordinatedR = Restrict(uncoordinated,NREM_R);
    uncoordinatedB_V = Restrict(uncoordinatedV,NREM_B);    uncoordinatedR_V = Restrict(uncoordinatedV,NREM_R);    uncoordinatedA_V = Restrict(uncoordinatedV,NREM_A);
    
    load('coordinated_ripple_bursts.mat')
    
    ripple_burst_B = Restrict(coordinated_ripple_bursts,NREM_B);
    ripple_burst_R = Restrict(coordinated_ripple_bursts,NREM_R);
    ripple_burst_A = Restrict(coordinated_ripple_bursts,NREM_A);

    
    %% Spikes
    %Load Units
    cd 'Spikesorting'
    spks = double([readNPY('spike_clusters.npy') readNPY('spike_times.npy')]);
    K = tdfread('cluster_group.tsv'); % import clusters ID and groups (MUA,Noise,Good)
    Kinfo = tdfread('cluster_info.tsv'); % import information of clusters
    K = [K.cluster_id(K.group(:,1) == 'g') , Kinfo.ch(K.group(:,1) == 'g'),]; % defining good clusters
    % Load neuronal classification
    load('Cell_type_classification')
    load('tag_shock_responsive_cells.mat')
    
    K = [K , Cell_type_classification(:,6:7)];
    group_dHPC = K(K(:,2) > 63,:);
    group_vHPC = K(K(:,2) <= 63,:);
    %Loop to select dorsal or ventral LFP and SU
    % z=1 --> dorsal
    % z=2 --> ventral
    for z = 1:2
        if z == 1
            n_SU_D = n_SU_D + length(group_dHPC);
            spks_dHPC = spks(ismember(spks(:,1),group_dHPC(:,1)),:); %keep spks from good clusters
        else
            n_SU_V = n_SU_V + length(group_vHPC);
            spks_vHPC = spks(ismember(spks(:,1),group_vHPC(:,1)),:); %keep spks from good clusters
        end
    end
    clear z
    spks_vHPC(:,2) = double(spks_vHPC(:,2))./20000;
    spks_dHPC(:,2) = double(spks_dHPC(:,2))./20000;
    
        
%     load('tag_shock_responsive_cells.mat')
    
    % constructing Spiketrains
    freq = 1/binSize;
    limits = [0 segments.Var1(end)/1000];
    spiketrains_dHPC = [];
    spiketrains_vHPC = [];
    
    for ii=1:length(group_dHPC)
        cluster = group_dHPC(ii,1);
        cellulartype = logical(Cell_type_classification(Cell_type_classification(:,1) == cluster,6));
        if and(cellulartype,criteriaD(ii,1))
            spks = spks_dHPC(spks_dHPC(:,1)==cluster,2);
            [tmp,bins]=binspikes(spks,freq,limits);
            m1 = mean(movmean(tmp(InIntervals(bins,NREM_B)),freq));
            m2 = mean(movmean(tmp(InIntervals(bins,NREM_R)),freq));
            m3 = mean(movmean(tmp(InIntervals(bins,NREM_A)),freq));
            m4 = mean([m1 m2 m3]);
            if m4 > criteria_fr
                if logical(group_dHPC(ii,3)) %check if is pyr
                    spiketrains_dHPC = [spiketrains_dHPC , tmp];
                end
            end
            clear spks tmp m1 m2 m3 m4
        end
    end
    
    for ii=1:length(group_vHPC)
        cluster = group_vHPC(ii,1);
        cellulartype = logical(Cell_type_classification(Cell_type_classification(:,1) == cluster,6));
        if and(cellulartype,criteriaV(ii,1))
            spks = spks_vHPC(spks_vHPC(:,1)==cluster,2);
            [tmp,bins]=binspikes(spks,freq,limits);
            m1 = mean(movmean(tmp(InIntervals(bins,NREM_B)),freq));
            m2 = mean(movmean(tmp(InIntervals(bins,NREM_R)),freq));
            m3 = mean(movmean(tmp(InIntervals(bins,NREM_A)),freq));
            m4 = mean([m1 m2 m3]);
            if m4 > criteria_fr
                if logical(group_vHPC(ii,3)) %check if is pyr
                    spiketrains_vHPC = [spiketrains_vHPC , tmp];
                end
            end
            clear spks tmp m1 m2 m3 m4
        end
    end
    clear freq limits
    
    if and(size(spiketrains_vHPC,2) >= criteria_n,size(spiketrains_dHPC,2) >= criteria_n)
        %Restricting bins inside each condition
        is.baseline.sws = InIntervals(bins,NREM_B);
        is.aversive.sws = InIntervals(bins,NREM_A);
        is.reward.sws = InIntervals(bins,NREM_R);
        
        is.baseline.rem = InIntervals(bins,REM_B);
        is.aversive.rem = InIntervals(bins,REM_A);
        is.reward.rem = InIntervals(bins,REM_R);
        
        is.aversive.run = InIntervals(bins,aversiveTS_run ./ 1000);
        is.reward.run = InIntervals(bins,rewardTS_run ./ 1000);
        
%         is.aversive.run = InIntervals(bins,speedA);
%         is.reward.run = InIntervals(bins,speedR);        
        
        is.baseline.ripplesD = InIntervals(bins,Restrict([ripplesD(:,2)-2 ripplesD(:,2)+2],NREM_B));
        is.aversive.ripplesD = InIntervals(bins,Restrict([ripplesD(:,2)-2 ripplesD(:,2)+2],NREM_A));
        is.reward.ripplesD = InIntervals(bins,Restrict([ripplesD(:,2)-2 ripplesD(:,2)+2],NREM_R));
        
        is.baseline.ripplesV = InIntervals(bins,Restrict([ripplesV(:,2)-2 ripplesV(:,2)+2],NREM_B));
        is.aversive.ripplesV = InIntervals(bins,Restrict([ripplesV(:,2)-2 ripplesV(:,2)+2],NREM_A));
        is.reward.ripplesV = InIntervals(bins,Restrict([ripplesV(:,2)-2 ripplesV(:,2)+2],NREM_R));
        
        is.baseline.burst = InIntervals(bins,Restrict([ripple_burst_B(:,2)-2 ripple_burst_B(:,2)+2],NREM_B));
        is.aversive.ripplesD = InIntervals(bins,Restrict([ripple_burst_A(:,2)-2 ripple_burst_A(:,2)+2],NREM_A));
        is.reward.ripplesD = InIntervals(bins,Restrict([ripple_burst_R(:,2)-2 ripple_burst_R(:,2)+2],NREM_R));
        
        
        %Zscoring and Restrincing Spike trians
        dHPCtrains.run.reward = zscore(spiketrains_dHPC(is.reward.run,:),0,1);
        dHPCtrains.run.aversive = zscore(spiketrains_dHPC(is.aversive.run,:),0,1);
        dHPCtrains.sleep.baseline = zscore(spiketrains_dHPC(is.baseline.sws,:),0,1);
        dHPCtrains.sleep.reward = zscore(spiketrains_dHPC(is.reward.sws,:),0,1);
        dHPCtrains.sleep.aversive = zscore(spiketrains_dHPC(is.aversive.sws,:),0,1);
        
        vHPCtrains.run.reward = zscore(spiketrains_vHPC(is.reward.run,:),0,1);
        vHPCtrains.run.aversive = zscore(spiketrains_vHPC(is.aversive.run,:),0,1);
        vHPCtrains.sleep.baseline = zscore(spiketrains_vHPC(is.baseline.sws,:),0,1);
        vHPCtrains.sleep.reward = zscore(spiketrains_vHPC(is.reward.sws,:),0,1);
        vHPCtrains.sleep.aversive = zscore(spiketrains_vHPC(is.aversive.sws,:),0,1);
        
        %Correlation Matrix
        %RUN z-scored correlation matrix
        corrM.dHPC.reward = dHPCtrains.run.reward' * dHPCtrains.run.reward / sum(is.reward.run);
        corrM.dHPC.aversive = dHPCtrains.run.aversive' * dHPCtrains.run.aversive / sum(is.aversive.run);
        corrM.vHPC.reward = vHPCtrains.run.reward' * vHPCtrains.run.reward / sum(is.reward.run);
        corrM.vHPC.aversive = vHPCtrains.run.aversive' * vHPCtrains.run.aversive / sum(is.aversive.run);
        corrM.DV.reward = (1/sum(is.reward.run)) * dHPCtrains.run.reward' * vHPCtrains.run.reward;
        corrM.DV.aversive = (1/sum(is.aversive.run)) * dHPCtrains.run.aversive' * vHPCtrains.run.aversive;
        
        %Calculation of Reactivation Strenght(R) across sleep sessions
        clear R
        for i = 1:sum(is.baseline.sws)
            R.dHPC.baseline.sws.reward(i) = dHPCtrains.sleep.baseline(i,:) * corrM.dHPC.reward * dHPCtrains.sleep.baseline(i,:)' - dHPCtrains.sleep.baseline(i,:) * dHPCtrains.sleep.baseline(i,:)';
            R.dHPC.baseline.sws.aversive(i) = dHPCtrains.sleep.baseline(i,:) * corrM.dHPC.aversive * dHPCtrains.sleep.baseline(i,:)' - dHPCtrains.sleep.baseline(i,:) * dHPCtrains.sleep.baseline(i,:)';
            R.vHPC.baseline.sws.reward(i) = vHPCtrains.sleep.baseline(i,:) * corrM.vHPC.reward * vHPCtrains.sleep.baseline(i,:)' - vHPCtrains.sleep.baseline(i,:) * vHPCtrains.sleep.baseline(i,:)';
            R.vHPC.baseline.sws.aversive(i) = vHPCtrains.sleep.baseline(i,:) * corrM.vHPC.aversive * vHPCtrains.sleep.baseline(i,:)' - vHPCtrains.sleep.baseline(i,:) * vHPCtrains.sleep.baseline(i,:)';
            R.DV.baseline.sws.reward(i) = dHPCtrains.sleep.baseline(i,:) * corrM.DV.reward * vHPCtrains.sleep.baseline(i,:)';
            R.DV.baseline.sws.aversive(i) = dHPCtrains.sleep.baseline(i,:) * corrM.DV.aversive * vHPCtrains.sleep.baseline(i,:)';
        end
        
        for i = 1:sum(is.reward.sws)
            R.dHPC.post.sws.reward(i) = dHPCtrains.sleep.reward(i,:) * corrM.dHPC.reward * dHPCtrains.sleep.reward(i,:)' - dHPCtrains.sleep.reward(i,:) * dHPCtrains.sleep.reward(i,:)';
            R.vHPC.post.sws.reward(i) = vHPCtrains.sleep.reward(i,:) * corrM.vHPC.reward * vHPCtrains.sleep.reward(i,:)' - vHPCtrains.sleep.reward(i,:) * vHPCtrains.sleep.reward(i,:)';
            R.DV.post.sws.reward(i) = dHPCtrains.sleep.reward(i,:) * corrM.DV.reward * vHPCtrains.sleep.reward(i,:)';
        end
        
        for i = 1:sum(is.aversive.sws)
            R.dHPC.post.sws.aversive(i) = dHPCtrains.sleep.aversive(i,:) * corrM.dHPC.aversive * dHPCtrains.sleep.aversive(i,:)' - dHPCtrains.sleep.aversive(i,:) * dHPCtrains.sleep.aversive(i,:)';
            R.vHPC.post.sws.aversive(i) = vHPCtrains.sleep.aversive(i,:) * corrM.vHPC.aversive * vHPCtrains.sleep.aversive(i,:)' - vHPCtrains.sleep.aversive(i,:) * vHPCtrains.sleep.aversive(i,:)';
            R.DV.post.sws.aversive(i) = dHPCtrains.sleep.aversive(i,:) * corrM.DV.aversive * vHPCtrains.sleep.aversive(i,:)';
        end
        
        
        event = ripple_burst_B(:,2);
        [syncR.dHPC.baseline.reward , idc.dHPC.baseline.reward] = Sync([bins(is.baseline.sws)'  R.DV.baseline.sws.reward'],event,'durations',window);
        figure,PlotSync(syncR.dHPC.baseline.reward,idc.dHPC.baseline.reward)
        
        event = ripple_burst_R(:,2);
        [syncR.dHPC.baseline.reward , idc.dHPC.baseline.reward] = Sync([bins(is.reward.sws)'  R.DV.post.sws.reward'],event,'durations',window);
        figure,PlotSync(syncR.dHPC.baseline.reward,idc.dHPC.baseline.reward)
        
        event = ripple_burst_A(:,2);
        [syncR.dHPC.baseline.reward , idc.dHPC.baseline.reward] = Sync([bins(is.aversive.sws)'  R.DV.post.sws.aversive'],event,'durations',window);
        figure,PlotSync(syncR.dHPC.baseline.reward,idc.dHPC.baseline.reward)
        
        
        %Baseline
        %dHPC
        event = Restrict(ripplesD(:,2),NREM_B);
        [syncR.dHPC.baseline.reward , idc.dHPC.baseline.reward] = Sync([bins(is.baseline.sws)'  R.dHPC.baseline.sws.reward'],event,'durations',window);
%         figure,PlotSync(syncR.dHPC.baseline.reward,idc.dHPC.baseline.reward)
        
        [syncR.dHPC.baseline.aversive , idc.dHPC.baseline.aversive] = Sync([bins(is.baseline.sws)'  R.dHPC.baseline.sws.aversive'],event,'durations',window);
%         figure,PlotSync(syncR.dHPC.baseline.aversive , idc.dHPC.baseline.aversive)
        
        %vHPC
        event = Restrict(ripplesV(:,2),NREM_B);
        [syncR.vHPC.baseline.reward , idc.vHPC.baseline.reward] = Sync([bins(is.baseline.sws)'  R.vHPC.baseline.sws.reward'],event,'durations',window);
%         figure,PlotSync(syncR.dHPC.baseline.reward,idc.dHPC.baseline.reward)
        
        [syncR.vHPC.baseline.aversive , idc.vHPC.baseline.aversive] = Sync([bins(is.baseline.sws)'  R.vHPC.baseline.sws.aversive'],event,'durations',window);
%         figure,PlotSync(syncR.vHPC.baseline.aversive , idc.vHPC.baseline.aversive)
        
        % dHPC - vHPC
        event = Restrict(ripplesD(:,2),NREM_B);
        [syncR.DV.baseline.reward.dRipples , idc.DV.baseline.reward.dRipples] = Sync([bins(is.baseline.sws)'  R.DV.baseline.sws.reward'],event,'durations',window);
%         figure,PlotSync(syncR.DV.baseline.reward , idc.DV.baseline.reward)
        
        [syncR.DV.baseline.aversive.dRipples , idc.DV.baseline.aversive.dRipples] = Sync([bins(is.baseline.sws)'  R.DV.baseline.sws.aversive'],event,'durations',window);
%         figure,PlotSync(syncR.DV.baseline.aversive , idc.DV.baseline.aversive)
        
        event = Restrict(ripplesV(:,2),NREM_B);
        [syncR.DV.baseline.reward.vRipples , idc.DV.baseline.reward.vRipples] = Sync([bins(is.baseline.sws)'  R.DV.baseline.sws.reward'],event,'durations',window);
%         figure,PlotSync(syncR.DV.baseline.reward , idc.DV.baseline.reward)
        
        [syncR.DV.baseline.aversive.vRipples , idc.DV.baseline.aversive.vRipples] = Sync([bins(is.baseline.sws)'  R.DV.baseline.sws.aversive'],event,'durations',window);
%         figure,PlotSync(syncR.DV.baseline.aversive , idc.DV.baseline.aversive)

        %Reward
        event = Restrict(ripplesD(:,2),NREM_R);
        [syncR.dHPC.post.reward , idc.dHPC.post.reward] = Sync([bins(is.reward.sws)'  R.dHPC.post.sws.reward'],event,'durations',window);
%         figure,PlotSync(syncR.dHPC.post.reward , idc.dPC.post.reward)
        
        %vHPC
        event = Restrict(ripplesV(:,2),NREM_R);
        [syncR.vHPC.post.reward , idc.vHPC.post.reward] = Sync([bins(is.reward.sws)'  R.vHPC.post.sws.reward'],event,'durations',window);
%         figure,PlotSync(syncR.vHPC.post.reward , idc.vHPC.post.reward)
        
        %dHPC - vHPC
        event = Restrict(ripplesD(:,2),NREM_R);
        [syncR.DV.post.reward.dRipples , idc.DV.post.reward.dRipples] = Sync([bins(is.reward.sws)'  R.DV.post.sws.reward'],event,'durations',window);
%         figure,PlotSync(syncR.DV.post.reward.dRipples , idc.DV.post.reward.dRipples)
        
        event = Restrict(ripplesV(:,2),NREM_R);
        [syncR.DV.post.reward.vRipples , idc.DV.post.reward.vRipples] = Sync([bins(is.reward.sws)'  R.DV.post.sws.reward'],event,'durations',window);
%         figure,PlotSync(syncR.DV.post.reward.vRipples , idc.DV.post.reward.vRipples)
        
        
        %Aversive
        event = Restrict(ripplesD(:,2),NREM_A);
        [syncR.dHPC.post.aversive , idc.dHPC.post.aversive] = Sync([bins(is.aversive.sws)'  R.dHPC.post.sws.aversive'],event,'durations',window);
%         figure,PlotSync(syncR.dHPC.post.aversive , idc.dPC.post.aversive)
        
        %vHPC
        event = Restrict(ripplesV(:,2),NREM_A);
        [syncR.vHPC.post.aversive , idc.vHPC.post.aversive] = Sync([bins(is.aversive.sws)'  R.vHPC.post.sws.aversive'],event,'durations',window);
%         figure,PlotSync(syncR.vHPC.post.aversive , idc.vHPC.post.aversive)
        
        %dHPC - vHPC
        event = Restrict(ripplesD(:,2),NREM_A);
        [syncR.DV.post.aversive.dRipples , idc.DV.post.aversive.dRipples] = Sync([bins(is.aversive.sws)'  R.DV.post.sws.aversive'],event,'durations',window);
%         figure,PlotSync(syncR.DV.post.aversive.dRipples , idc.DV.post.aversive.dRipples)
        
        event = Restrict(ripplesV(:,2),NREM_A);
        [syncR.DV.post.aversive.vRipples , idc.DV.post.aversive.vRipples] = Sync([bins(is.aversive.sws)'  R.DV.post.sws.aversive'],event,'durations',window);
%         figure,PlotSync(syncR.DV.post.aversive.vRipples , idc.DV.post.aversive.vRipples)
        
        
        [meanR.dHPC.pre.aversive,err.dHPC.pre.aversive,tb]=SyncHist(syncR.dHPC.baseline.aversive , idc.dHPC.baseline.aversive,'mode','mean','durations',window,'smooth',2);
        [meanR.dHPC.pre.reward,err.dHPC.pre.reward,tb]=SyncHist(syncR.dHPC.baseline.reward , idc.dHPC.baseline.reward,'mode','mean','durations',window,'smooth',2);
       
        [meanR.dHPC.post.aversive,err.dHPC.post.aversive,tb]=SyncHist(syncR.dHPC.post.aversive , idc.dHPC.post.aversive,'mode','mean','durations',window,'smooth',2);
        [meanR.dHPC.post.reward,err.dHPC.post.reward,tb]=SyncHist(syncR.dHPC.post.reward , idc.dHPC.post.reward,'mode','mean','durations',window,'smooth',2);
        
        figure,subplot(121),plot(meanR.dHPC.pre.aversive,'k'),hold on
        plot(meanR.dHPC.post.aversive,'r')
        
        subplot(122),plot(meanR.dHPC.pre.reward,'k'),hold on
        plot(meanR.dHPC.post.reward,'b')
        
        
        [meanR.vHPC.pre.aversive,err.vHPC.pre.aversive,tb]=SyncHist(syncR.vHPC.baseline.aversive , idc.vHPC.baseline.aversive,'mode','mean','durations',window,'smooth',2);
        [meanR.vHPC.pre.reward,err.vHPC.pre.reward,tb]=SyncHist(syncR.vHPC.baseline.reward , idc.vHPC.baseline.reward,'mode','mean','durations',window,'smooth',2);
       
        [meanR.vHPC.post.aversive,err.vHPC.post.aversive,tb]=SyncHist(syncR.vHPC.post.aversive , idc.vHPC.post.aversive,'mode','mean','durations',window,'smooth',2);
        [meanR.vHPC.post.reward,err.vHPC.post.reward,tb]=SyncHist(syncR.vHPC.post.reward , idc.vHPC.post.reward,'mode','mean','durations',window,'smooth',2);
        
        figure,subplot(121),plot(meanR.vHPC.pre.aversive,'k'),hold on
        plot(meanR.vHPC.post.aversive,'r')
        
        subplot(122),plot(meanR.vHPC.pre.reward,'k'),hold on
        plot(meanR.vHPC.post.reward,'b')
        
        
        [meanR.DV.pre.aversive,err.DV.pre.aversive,tb]=SyncHist(syncR.DV.baseline.aversive.dRipples , idc.DV.baseline.aversive.dRipples,'mode','mean','durations',window,'smooth',2);
        [meanR.DV.pre.reward,err.DV.pre.reward,tb]=SyncHist(syncR.DV.baseline.reward.dRipples , idc.DV.baseline.reward.dRipples,'mode','mean','durations',window,'smooth',2);
       
        [meanR.DV.post.aversive,err.DV.post.aversive,tb]=SyncHist(syncR.DV.post.aversive.dRipples , idc.DV.post.aversive.dRipples,'mode','mean','durations',window,'smooth',2);
        [meanR.DV.post.reward,err.vHPC.post.reward,tb]=SyncHist(syncR.DV.post.reward.dRipples , idc.DV.post.reward.dRipples,'mode','mean','durations',window,'smooth',2);
        
        figure,subplot(121),plot(meanR.DV.pre.aversive,'k'),hold on
        plot(meanR.DV.post.aversive,'r')
        
        subplot(122),plot(meanR.DV.pre.reward,'k'),hold on
        plot(meanR.DV.post.reward,'b')
        
        
        [meanR.DV.pre.aversive,err.DV.pre.aversive,tb]=SyncHist(syncR.DV.baseline.aversive.vRipples , idc.DV.baseline.aversive.vRipples,'mode','mean','durations',window,'smooth',2);
        [meanR.DV.pre.reward,err.DV.pre.reward,tb]=SyncHist(syncR.DV.baseline.reward.vRipples , idc.DV.baseline.reward.vRipples,'mode','mean','durations',window,'smooth',2);
       
        [meanR.DV.post.aversive,err.DV.post.aversive,tb]=SyncHist(syncR.DV.post.aversive.dRipples , idc.DV.post.aversive.dRipples,'mode','mean','durations',window,'smooth',2);
        [meanR.DV.post.reward,err.vHPC.post.reward,tb]=SyncHist(syncR.DV.post.reward.vRipples , idc.DV.post.reward.vRipples,'mode','mean','durations',window,'smooth',2);
        
        figure,subplot(121),plot(meanR.DV.pre.aversive,'k'),hold on
        plot(meanR.DV.post.aversive,'r')
        
        subplot(122),plot(meanR.DV.pre.reward,'k'),hold on
        plot(meanR.DV.post.reward,'b')
        
        
        %% Coordinated ripples
        %Baseline
        %dHPC
        event = Restrict(coordinated(:,2),NREM_B);
        [syncR.dHPC.baseline.reward , idc.dHPC.baseline.reward] = Sync([bins(is.baseline.sws)'  R.dHPC.baseline.sws.reward'],event,'durations',window);
%         figure,PlotSync(syncR.dHPC.baseline.reward,idc.dHPC.baseline.reward)
        
        [syncR.dHPC.baseline.aversive , idc.dHPC.baseline.aversive] = Sync([bins(is.baseline.sws)'  R.dHPC.baseline.sws.aversive'],event,'durations',window);
%         figure,PlotSync(syncR.dHPC.baseline.aversive , idc.dHPC.baseline.aversive)
        
        %vHPC
        event = Restrict(coordinatedV(:,2),NREM_B);
        [syncR.vHPC.baseline.reward , idc.vHPC.baseline.reward] = Sync([bins(is.baseline.sws)'  R.vHPC.baseline.sws.reward'],event,'durations',window);
%         figure,PlotSync(syncR.dHPC.baseline.reward,idc.dHPC.baseline.reward)
        
        [syncR.vHPC.baseline.aversive , idc.vHPC.baseline.aversive] = Sync([bins(is.baseline.sws)'  R.vHPC.baseline.sws.aversive'],event,'durations',window);
%         figure,PlotSync(syncR.vHPC.baseline.aversive , idc.vHPC.baseline.aversive)
        
        % dHPC - vHPC
        event = Restrict(coordinated(:,2),NREM_B);
        [syncR.DV.baseline.reward.dRipples , idc.DV.baseline.reward.dRipples] = Sync([bins(is.baseline.sws)'  R.DV.baseline.sws.reward'],event,'durations',window);
%         figure,PlotSync(syncR.DV.baseline.reward , idc.DV.baseline.reward)
        
        [syncR.DV.baseline.aversive.dRipples , idc.DV.baseline.aversive.dRipples] = Sync([bins(is.baseline.sws)'  R.DV.baseline.sws.aversive'],event,'durations',window);
%         figure,PlotSync(syncR.DV.baseline.aversive , idc.DV.baseline.aversive)

        event = Restrict(coordinatedV(:,2),NREM_B);
        [syncR.DV.baseline.reward.vRipples , idc.DV.baseline.reward.vRipples] = Sync([bins(is.baseline.sws)'  R.DV.baseline.sws.reward'],event,'durations',window);
%         figure,PlotSync(syncR.DV.baseline.reward , idc.DV.baseline.reward)
        
        [syncR.DV.baseline.aversive.vRipples , idc.DV.baseline.aversive.vRipples] = Sync([bins(is.baseline.sws)'  R.DV.baseline.sws.aversive'],event,'durations',window);
%         figure,PlotSync(syncR.DV.baseline.aversive , idc.DV.baseline.aversive)        
        
        %Reward
        event = Restrict(coordinated(:,2),NREM_R);
        [syncR.dHPC.post.reward , idc.dHPC.post.reward] = Sync([bins(is.reward.sws)'  R.dHPC.post.sws.reward'],event,'durations',window);
%         figure,PlotSync(syncR.dHPC.post.reward , idc.dPC.post.reward)
        
        %vHPC
        event = Restrict(coordinatedV(:,2),NREM_R);
        [syncR.vHPC.post.reward , idc.vHPC.post.reward] = Sync([bins(is.reward.sws)'  R.vHPC.post.sws.reward'],event,'durations',window);
%         figure,PlotSync(syncR.vHPC.post.reward , idc.vHPC.post.reward)
        
        %dHPC - vHPC
        event = Restrict(coordinated(:,2),NREM_R);
        [syncR.DV.post.reward.dRipples , idc.DV.post.reward.dRipples] = Sync([bins(is.reward.sws)'  R.DV.post.sws.reward'],event,'durations',window);
%         figure,PlotSync(syncR.DV.post.reward.dRipples , idc.DV.post.reward.dRipples)
        
        event = Restrict(coordinatedV(:,2),NREM_R);
        [syncR.DV.post.reward.vRipples , idc.DV.post.reward.vRipples] = Sync([bins(is.reward.sws)'  R.DV.post.sws.reward'],event,'durations',window);
%         figure,PlotSync(syncR.DV.post.reward.vRipples , idc.DV.post.reward.vRipples)
        
        
        %Aversive
        event = Restrict(coordinated(:,2),NREM_A);
        [syncR.dHPC.post.aversive , idc.dHPC.post.aversive] = Sync([bins(is.aversive.sws)'  R.dHPC.post.sws.aversive'],event,'durations',window);
%         figure,PlotSync(syncR.dHPC.post.aversive , idc.dPC.post.aversive)
        
        %vHPC
        event = Restrict(coordinatedV(:,2),NREM_A);
        [syncR.vHPC.post.aversive , idc.vHPC.post.aversive] = Sync([bins(is.aversive.sws)'  R.vHPC.post.sws.aversive'],event,'durations',window);
%         figure,PlotSync(syncR.vHPC.post.aversive , idc.vHPC.post.aversive)
        
        %dHPC - vHPC
        event = Restrict(coordinated(:,2),NREM_A);
        [syncR.DV.post.aversive.dRipples , idc.DV.post.aversive.dRipples] = Sync([bins(is.aversive.sws)'  R.DV.post.sws.aversive'],event,'durations',window);
%         figure,PlotSync(syncR.DV.post.aversive.dRipples , idc.DV.post.aversive.dRipples)
        
        event = Restrict(coordinatedV(:,2),NREM_A);
        [syncR.DV.post.aversive.vRipples , idc.DV.post.aversive.vRipples] = Sync([bins(is.aversive.sws)'  R.DV.post.sws.aversive'],event,'durations',window);
%         figure,PlotSync(syncR.DV.post.aversive.vRipples , idc.DV.post.aversive.vRipples)        
        
    
        [meanR.dHPC.pre.aversive,err.dHPC.pre.aversive,tb]=SyncHist(syncR.dHPC.baseline.aversive , idc.dHPC.baseline.aversive,'mode','mean','durations',window,'smooth',2);
        [meanR.dHPC.pre.reward,err.dHPC.pre.reward,tb]=SyncHist(syncR.dHPC.baseline.reward , idc.dHPC.baseline.reward,'mode','mean','durations',window,'smooth',2);
       
        [meanR.dHPC.post.aversive,err.dHPC.post.aversive,tb]=SyncHist(syncR.dHPC.post.aversive , idc.dHPC.post.aversive,'mode','mean','durations',window,'smooth',2);
        [meanR.dHPC.post.reward,err.dHPC.post.reward,tb]=SyncHist(syncR.dHPC.post.reward , idc.dHPC.post.reward,'mode','mean','durations',window,'smooth',2);
        
        figure,subplot(121),plot(meanR.dHPC.pre.aversive,'k'),hold on
        plot(meanR.dHPC.post.aversive,'r')
        
        subplot(122),plot(meanR.dHPC.pre.reward,'k'),hold on
        plot(meanR.dHPC.post.reward,'b')
        
        [meanR.vHPC.pre.aversive,err.vHPC.pre.aversive,tb]=SyncHist(syncR.vHPC.baseline.aversive , idc.vHPC.baseline.aversive,'mode','mean','durations',window,'smooth',2);
        [meanR.vHPC.pre.reward,err.vHPC.pre.reward,tb]=SyncHist(syncR.vHPC.baseline.reward , idc.vHPC.baseline.reward,'mode','mean','durations',window,'smooth',2);
       
        [meanR.vHPC.post.aversive,err.vHPC.post.aversive,tb]=SyncHist(syncR.vHPC.post.aversive , idc.vHPC.post.aversive,'mode','mean','durations',window,'smooth',2);
        [meanR.vHPC.post.reward,err.vHPC.post.reward,tb]=SyncHist(syncR.vHPC.post.reward , idc.vHPC.post.reward,'mode','mean','durations',window,'smooth',2);
        
        figure,subplot(121),plot(meanR.vHPC.pre.aversive,'k'),hold on
        plot(meanR.vHPC.post.aversive,'r')
        
        subplot(122),plot(meanR.vHPC.pre.reward,'k'),hold on
        plot(meanR.vHPC.post.reward,'b')
        
        
        [meanR.DV.pre.aversive,err.DV.pre.aversive,tb]=SyncHist(syncR.DV.baseline.aversive.dRipples , idc.DV.baseline.aversive.dRipples,'mode','mean','durations',window,'smooth',2);
        [meanR.DV.pre.reward,err.DV.pre.reward,tb]=SyncHist(syncR.DV.baseline.reward.dRipples , idc.DV.baseline.reward.dRipples,'mode','mean','durations',window,'smooth',2);
       
        [meanR.DV.post.aversive,err.DV.post.aversive,tb]=SyncHist(syncR.DV.post.aversive.dRipples , idc.DV.post.aversive.dRipples,'mode','mean','durations',window,'smooth',2);
        [meanR.DV.post.reward,err.vHPC.post.reward,tb]=SyncHist(syncR.DV.post.reward.dRipples , idc.DV.post.reward.dRipples,'mode','mean','durations',window,'smooth',2);
        
        figure,subplot(121),plot(meanR.DV.pre.aversive,'k'),hold on
        plot(meanR.DV.post.aversive,'r')
        
        subplot(122),plot(meanR.DV.pre.reward,'k'),hold on
        plot(meanR.DV.post.reward,'b')
        
        
        [meanR.DV.pre.aversive,err.DV.pre.aversive,tb]=SyncHist(syncR.DV.baseline.aversive.vRipples , idc.DV.baseline.aversive.vRipples,'mode','mean','durations',window,'smooth',2);
        [meanR.DV.pre.reward,err.DV.pre.reward,tb]=SyncHist(syncR.DV.baseline.reward.vRipples , idc.DV.baseline.reward.vRipples,'mode','mean','durations',window,'smooth',2);
       
        [meanR.DV.post.aversive,err.DV.post.aversive,tb]=SyncHist(syncR.DV.post.aversive.dRipples , idc.DV.post.aversive.dRipples,'mode','mean','durations',window,'smooth',2);
        [meanR.DV.post.reward,err.vHPC.post.reward,tb]=SyncHist(syncR.DV.post.reward.vRipples , idc.DV.post.reward.vRipples,'mode','mean','durations',window,'smooth',2);
        
        figure,subplot(121),plot(meanR.DV.pre.aversive,'k'),hold on
        plot(meanR.DV.post.aversive,'r')
        
        subplot(122),plot(meanR.DV.pre.reward,'k'),hold on
        plot(meanR.DV.post.reward,'b')
        
    end
    
    clear spiketrains_dHPC spiketrains_vHPC is
    clear Sx Sy Sz T S1 S2 CorrMatrixR_RUN CorrMatrixA_RUN
    t
    end
% end
tt

end



% subplot(1,3,1),imagesc([1:size(S1,2)], [1:size(S1,1)], CorrMatrixR_RUN),caxis([-0.04 0.07]), axis tight,xlabel('dorsal SU'),ylabel('ventral SU'),title('Baseline Sleep')
% subplot(1,3,2),imagesc([1:size(T,2)], [1:size(T,1)], CorrMatrixR_RUN), axis tight,xlabel('dorsal SU'),ylabel('ventral SU'),title('Reward Sleep')
% subplot(1,3,3),imagesc([1:size(S2,2)], [1:size(S2,1)], CorrMatrixA_RUN), axis tight,xlabel('dorsal SU'),ylabel('ventral SU'),title('Aversive Sleep')
% 
% boxplot([REV_R , EV_R , REV_A , EV_A])



subplot(121)
boxplot([REV_A*100 , EV_A*100]),ylim([0 10])
subplot(122)
boxplot([REV_R*100 , EV_R*100]),ylim([0 10])

[h p] = ttest(REV_A , EV_A)