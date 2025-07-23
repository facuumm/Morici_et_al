clear
clc
close all

%% Parameters
path = {'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready'};%List of folders from the path

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
%% Main loop, to iterate across sessions
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
    
    %% Sleep
    x = dir([cd,'\*-states.mat']);    states = load([cd,'\',x.name]);    states = states.states;
    
    REM = ToIntervals(states==5);    NREM = ToIntervals(states==3);    WAKE = ToIntervals(states==1);
    clear x states
    
    % NREM events restriction according conditions
    NREM_B = NREM(NREM(:,2)<baselineTS(1,2)/1000,:);
    x = (NREM_B(:,2)-NREM_B(:,1));
    NREM_B = NREM_B(x>180,:);
    x = (NREM_B(:,2)-NREM_B(:,1));
    d = 0;
    tmp = [];
    for i = 1:length(NREM_B)
        d = d + x(i);
        if d < time_criteria
            tmp = [tmp ; NREM_B(i,:)];
        else
            y = d - time_criteria;
            tmp = [tmp ; NREM_B(i,1) , NREM_B(i,2)-y];
            clear y
            break
        end
    end
    NREM_B = tmp;
    clear x d tmp
    
    NREM_A = NREM(NREM(:,2)>aversiveTS(1,1)/1000 & NREM(:,2)<aversiveTS(1,2)/1000,:);
    x = (NREM_A(:,2)-NREM_A(:,1));
    NREM_A = NREM_A(x>180,:);
    x = (NREM_A(:,2)-NREM_A(:,1));
    d = 0;
    tmp = [];
    for i = 1:length(NREM_A)
        d = d + x(i);
        if d < time_criteria
            tmp = [tmp ; NREM_A(i,:)];
        else
            y = d - time_criteria;
            tmp = [tmp ; NREM_A(i,1) , NREM_A(i,2)-y];
            clear y
            break
        end
    end
    NREM_A = tmp;
    clear x d tmp
        
    NREM_R = NREM(NREM(:,2)>rewardTS(1,1)/1000 & NREM(:,2)<rewardTS(1,2)/1000,:);
    x = (NREM_R(:,2)-NREM_R(:,1));
    NREM_R = NREM_R(x>180,:);
    x = (NREM_R(:,2)-NREM_R(:,1));
    d = 0;
    tmp = [];
    for i = 1:length(NREM_R)
        d = d + x(i);
        if d < time_criteria
            tmp = [tmp ; NREM_R(i,:)];
        else
            y = d - time_criteria;
            tmp = [tmp ; NREM_R(i,1) , NREM_R(i,2)-y];
            clear y
            break
        end
    end
    NREM_R = tmp;
    clear x d tmp
    
    
    % REM events restriction according conditions
    REM_B = REM(REM(:,2)<baselineTS(1,2)/1000,:);
    REM_A = REM(REM(:,2)>aversiveTS(1,1)/1000 & REM(:,2)<aversiveTS(1,2)/1000,:);
    REM_R = REM(REM(:,2)>rewardTS(1,1)/1000 & REM(:,2)<rewardTS(1,2)/1000,:);
%     % Coordinated dHPC ripples
%     coordinated = [];
%     coordinatedV = [];
%     middle = [];
%     for i = 1:length(ripplesD)
%         r = ripplesD(i,:);
%         tmp = sum(and(ripplesV(:,2)>= r(1,2)-0.05, ripplesV(:,2)<= r(1,2)+0.05));
%         if tmp>0
%             v = ripplesV(and(ripplesV(:,2)>= r(1,2)-0.05, ripplesV(:,2)<= r(1,2)+0.05),2);
%             coordinatedV = [coordinatedV ; v];
%             coordinated = [coordinated ; r(1,2)];
%             clear tmp2 tmp1 v m
%         end
%         clear r
%     end
%     clear x tmp i
%     
%     coordinatedB = Restrict(coordinated,NREM_B);    coordinatedA = Restrict(coordinated,NREM_A);    coordinatedR = Restrict(coordinated,NREM_R);
%     coordinatedB_V = Restrict(coordinatedV,NREM_B);    coordinatedR_V = Restrict(coordinatedV,NREM_R);    coordinatedA_V = Restrict(coordinatedV,NREM_A);
%     
%     p_coordinatedB = (length(coordinatedB)*100)/length(Restrict(ripplesD,NREM_B));
%     p_coordinatedA = (length(coordinatedA)*100)/length(Restrict(ripplesD,NREM_A));
%     p_coordinatedR = (length(coordinatedR)*100)/length(Restrict(ripplesD,NREM_R));
%     ripples_coordinated_percentage = [ripples_coordinated_percentage ; p_coordinatedB , p_coordinatedR , p_coordinatedA];
%     %     clear p_coordinatedB p_coordinatedA p_coordinatedR
%     
%     rateV_B = length(Restrict(ripplesV,NREM_B)) / sum(NREM_B(:,2) - NREM_B(:,1));
%     rateV_A = length(Restrict(ripplesV,NREM_A)) / sum(NREM_A(:,2) - NREM_A(:,1));
%     rateV_R = length(Restrict(ripplesV,NREM_R)) / sum(NREM_R(:,2) - NREM_R(:,1));
%     rateV = [rateV ; rateV_B , rateV_R , rateV_A];
%     %     clear rateV_B rateV_A rateV_R
%     
%     rateD_B = length(Restrict(ripplesD,NREM_B)) / sum(NREM_B(:,2) - NREM_B(:,1));
%     rateD_A = length(Restrict(ripplesD,NREM_A)) / sum(NREM_A(:,2) - NREM_A(:,1));
%     rateD_R = length(Restrict(ripplesD,NREM_R)) / sum(NREM_R(:,2) - NREM_R(:,1));
%     rateD = [rateD ; rateD_B , rateD_R , rateD_A];
%     %     clear rateD_B rateD_A rateD_R
%     
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
    
    for ii=1:length(group_vHPC)
        cluster = group_vHPC(ii,1);
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
        
        %Correlation Matrix
        [CorrMatrixR_RUN , PvalMatrixR_RUN]=corr(spiketrains_dHPC(is.reward.run,:),spiketrains_vHPC(is.reward.run,:));
        [CorrMatrixA_RUN , PvalMatrixA_RUN]=corr(spiketrains_dHPC(is.aversive.run,:),spiketrains_vHPC(is.aversive.run,:));
        
        %EV Aversive
        if aversiveTS_run(1) < rewardTS_run(1)
%             x = [spiketrains_dHPC(is.baseline.sws,:) ; spiketrains_dHPC(is.reward.sws,:)];
%             y = [spiketrains_vHPC(is.baseline.sws,:); spiketrains_vHPC(is.reward.sws,:)];
%             [S1 , PvalMatrixB_NREM]=corr(x,y);
%             [S2 , PvalMatrixA_NREM]=corr(spiketrains_dHPC(is.aversive.sws,:),spiketrains_vHPC(is.aversive.sws,:));
%             T = CorrMatrixA_RUN;
%         else
            x = [spiketrains_dHPC(is.baseline.sws,:)];
            y = [spiketrains_vHPC(is.baseline.sws,:)];
            [S1 , PvalMatrixB_NREM]=corr(x,y);
            [S2 , PvalMatrixA_NREM]=corr(spiketrains_dHPC(is.aversive.sws,:),spiketrains_vHPC(is.aversive.sws,:));
            T = CorrMatrixA_RUN;
            
            Sx = corrcoef(S2,T,'rows','complete');     Sx = Sx(1,2);
            Sy = corrcoef(S1,T,'rows','complete');     Sy = Sy(1,2);
            Sz = corrcoef(S1,S2,'rows','complete');     Sz = Sz(1,2);
            
            EV = (Sx-Sy*Sz);
            EV = (EV/sqrt((1-Sy^2)*(1-Sz^2)))^2;
            EV_A = [EV_A ; EV]; clear EV
            
            REV = (Sy-Sx*Sz);
            REV = (REV/sqrt((1-Sy^2)*(1-Sz^2)))^2;
            REV_A = [REV_A ; REV];
            clear Sx Sy Sz REV
        end
        
        %EV Reward
        if   rewardTS_run(1) < aversiveTS_run(1)
            x = [spiketrains_dHPC(is.baseline.sws,:)];
            y = [spiketrains_vHPC(is.baseline.sws,:)];
            [S1 , PvalMatrixB_NREM]=corr(x,y);
            [S2 , PvalMatrixA_NREM]=corr(spiketrains_dHPC(is.reward.sws,:),spiketrains_vHPC(is.reward.sws,:));
            T = CorrMatrixR_RUN;
%         else
%             x = [spiketrains_dHPC(is.baseline.sws,:) ; spiketrains_dHPC(is.aversive.sws,:)];
%             y = [spiketrains_vHPC(is.baseline.sws,:) ; spiketrains_vHPC(is.aversive.sws,:)];
%             [S1 , PvalMatrixB_NREM]=corr(x,y);
%             [S2 , PvalMatrixA_NREM]=corr(spiketrains_dHPC(is.reward.sws,:),spiketrains_vHPC(is.reward.sws,:));
%             T = CorrMatrixR_RUN;
            
            Sx = corrcoef(S2,T,'rows','complete');     Sx = Sx(1,2);
            Sy = corrcoef(S1,T,'rows','complete');     Sy = Sy(1,2);
            Sz = corrcoef(S1,S2,'rows','complete');     Sz = Sz(1,2);
            
            EV = (Sx-Sy*Sz);
            EV = (EV/sqrt((1-Sy^2)*(1-Sz^2)))^2;
            EV_R = [EV_R ; EV]; clear EV
            

            REV = (Sy-Sx*Sz);
            REV = (REV/sqrt((1-Sy^2)*(1-Sz^2)))^2;
            REV_R = [REV_R ; REV];
            clear Sx Sy Sz REV
        end
        
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