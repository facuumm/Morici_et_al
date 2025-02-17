clear
clc
close all

%% Parameters
path = {'E:\Rat103\usable';'E:\Rat126\Ephys\in_Pyr';'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat132\recordings\in_pyr';'E:\Rat165\in_pyr\'};%List of folders from the path

%Sleep
time_criteria = 600; %time criteria to define the maximal time of sleep to include

% Ripples
q = 0.25; %quantile to restrict above it ripples according to their peak amplitude
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
n_SU_V = [];
n_SU_D = [];
FR_B_V = []; FR_R_V = [];  FR_A_V = []; % Firing Rate during NREM
FR_B_D = []; FR_R_D = [];  FR_A_D = []; % Firing Rate during NREM
poisson_dHPC_split = []; poisson_vHPC_split = []; %poisson results split by conditions

%For EV and REV
binSize = 0.1;
EV_A = []; REV_A = [];
EV_R = []; REV_R = [];

window = [-2 2];

FR_dHPC_pyr = [];
FR_dHPC_int = [];
FR_vHPC_pyr = [];
FR_vHPC_int = [];


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
        
        %% Sleep
        x = dir([cd,'\*-states.mat']);    states = load([cd,'\',x.name]);    states = states.states;
        
        REM = ToIntervals(states==5);    NREM = ToIntervals(states==3);    WAKE = ToIntervals(states==1);
        clear x states
        
        %keep only WAKE in HomeCage
        %         WAKE = Restrict(WAKE, [aversiveTS ; rewardTS ; baselineTS] ./1000);
        
        % NREM events restriction according conditions
        NREM_B = NREM(NREM(:,2)<baselineTS(1,2)/1000,:);
        NREM_A = NREM(NREM(:,2)>aversiveTS(1,1)/1000 & NREM(:,2)<aversiveTS(1,2)/1000,:);
        NREM_R = NREM(NREM(:,2)>rewardTS(1,1)/1000 & NREM(:,2)<rewardTS(1,2)/1000,:);
        
        
        % REM events restriction according conditions
        REM_B = REM(REM(:,2)<baselineTS(1,2)/1000,:);
        REM_A = REM(REM(:,2)>aversiveTS(1,1)/1000 & REM(:,2)<aversiveTS(1,2)/1000,:);
        REM_R = REM(REM(:,2)>rewardTS(1,1)/1000 & REM(:,2)<rewardTS(1,2)/1000,:);
        
        %% Spikes
        %Load Units
        cd 'Spikesorting'
        spks = double([readNPY('spike_clusters.npy') readNPY('spike_times.npy')]);
        K = tdfread('cluster_group.tsv'); % import clusters ID and groups (MUA,Noise,Good)
        Kinfo = tdfread('cluster_info.tsv'); % import information of clusters
        K = [K.cluster_id(K.group(:,1) == 'g') , Kinfo.ch(K.group(:,1) == 'g'),]; % defining good clusters
        % Load neuronal classification
        load('Cell_type_classification')
        %         load('tag_shock_responsive_cells.mat')
        
        K = [K , Cell_type_classification(:,5:7)];
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
%                 n_SU_V = [n_SU_V ; length(group_vHPC)];
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
        spiketrains_dHPC_pyr = [];        spiketrains_dHPC_int = [];
        spiketrains_vHPC_pyr = [];        spiketrains_vHPC_int = [];
        
        for ii=1:size(group_dHPC,1)
            cluster = group_dHPC(ii,1);
            cellulartype = logical(Cell_type_classification(Cell_type_classification(:,1) == cluster,7));
            iterator = Cell_type_classification(Cell_type_classification(:,1) == cluster,:);
            %             if and(cellulartype,criteriaD(ii,1))
            spks = spks_dHPC(spks_dHPC(:,1)==cluster,2);
            [tmp,bins]=binspikes(spks,freq,limits);
            m1 = length(Restrict(spks,NREM_B))/sum(NREM_B(:,2)-NREM_B(:,1));
            m2 = length(Restrict(spks,NREM_R))/sum(NREM_R(:,2)-NREM_R(:,1));
            m3 = length(Restrict(spks,NREM_A))/sum(NREM_A(:,2)-NREM_A(:,1));
            
            %                 m1 = mean(movmean(tmp(InIntervals(bins,NREM_B))./binSize,freq));
            %                 m2 = mean(movmean(tmp(InIntervals(bins,NREM_R))./binSize,freq));
            %                 m3 = mean(movmean(tmp(InIntervals(bins,NREM_A))./binSize,freq));
            m4 = mean([m1 m2 m3]);
            
            m5 = length(Restrict(spks,REM_B))/sum(REM_B(:,2)-REM_B(:,1));
            m6 = length(Restrict(spks,REM_R))/sum(REM_R(:,2)-REM_R(:,1));
            m7 = length(Restrict(spks,REM_A))/sum(REM_A(:,2)-REM_A(:,1));
            
            %                 m5 = mean(movmean(tmp(InIntervals(bins,REM_B))./binSize,freq));
            %                 m6 = mean(movmean(tmp(InIntervals(bins,REM_B))./binSize,freq));
            %                 m7 = mean(movmean(tmp(InIntervals(bins,REM_A))./binSize,freq));
            m8 = mean([m5 m6 m7]);
            
            m9 = length(Restrict(spks,rewardTS_run./1000))/sum(rewardTS_run(2)./1000 - rewardTS_run(1)./1000);
            m10 = length(Restrict(spks,aversiveTS_run./1000))/sum(aversiveTS_run(2)./1000 - aversiveTS_run(1)./1000);
            
            %                 m9 = mean(movmean(tmp(InIntervals(bins,rewardTS_run./1000))./binSize,freq));
            %                 m10 = mean(movmean(tmp(InIntervals(bins,aversiveTS_run./1000))./binSize,freq));
            m11 = mean([m9 m10]);
            %                 m12 = mean(movmean(tmp(InIntervals(bins,WAKE))./binSize,freq));
            m12 = length(Restrict(spks,WAKE))/sum(WAKE(:,2)-WAKE(:,1));
            m13 = length(spks)/(segments.Var1(end)/1000);
            
            %                 if m4 > criteria_fr
            if logical(group_dHPC(ii,5)) %check if is pyr
                spiketrains_dHPC_pyr = [spiketrains_dHPC_pyr , zscore(tmp./binSize,0,1)];
                FR_dHPC_pyr = [FR_dHPC_pyr ; m4 m1 m2 m3 m8 m5 m6 m7 m11 m9 m10 m12 m13 iterator(6) iterator(5)];
            else
                spiketrains_dHPC_int = [spiketrains_dHPC_int , zscore(tmp./binSize,0,1)];
                FR_dHPC_int = [FR_dHPC_int ; m4 m1 m2 m3 m8 m5 m6 m7 m11 m9 m10 m12 m13 iterator(6) iterator(5)];
            end
            %                 end
            clear spks tmp m1 m2 m3 m4 m5 m6 m7 m8 m9 m10 m11 m12 m13
            %             end
        end
        countP = 0;
        countI = 0;
        for ii=1:size(group_vHPC,1)
            cluster = group_vHPC(ii,1);
            cellulartype = logical(Cell_type_classification(Cell_type_classification(:,1) == cluster,6));
            iterator = Cell_type_classification(Cell_type_classification(:,1) == cluster,:);
            %             if and(cellulartype,criteriaV(ii,1))
            spks = spks_vHPC(spks_vHPC(:,1)==cluster,2);
            [tmp,bins]=binspikes(spks,freq,limits);
            m1 = length(Restrict(spks,NREM_B))/sum(NREM_B(:,2)-NREM_B(:,1));
            m2 = length(Restrict(spks,NREM_R))/sum(NREM_R(:,2)-NREM_R(:,1));
            m3 = length(Restrict(spks,NREM_A))/sum(NREM_A(:,2)-NREM_A(:,1));
            
            %                 m1 = mean(movmean(tmp(InIntervals(bins,NREM_B))./binSize,freq));
            %                 m2 = mean(movmean(tmp(InIntervals(bins,NREM_R))./binSize,freq));
            %                 m3 = mean(movmean(tmp(InIntervals(bins,NREM_A))./binSize,freq));
            m4 = mean([m1 m2 m3]);
            
            m5 = length(Restrict(spks,REM_B))/sum(REM_B(:,2)-REM_B(:,1));
            m6 = length(Restrict(spks,REM_R))/sum(REM_R(:,2)-REM_R(:,1));
            m7 = length(Restrict(spks,REM_A))/sum(REM_A(:,2)-REM_A(:,1));
            
            %                 m5 = mean(movmean(tmp(InIntervals(bins,REM_B))./binSize,freq));
            %                 m6 = mean(movmean(tmp(InIntervals(bins,REM_B))./binSize,freq));
            %                 m7 = mean(movmean(tmp(InIntervals(bins,REM_A))./binSize,freq));
            m8 = mean([m5 m6 m7]);
            
            m9 = length(Restrict(spks,rewardTS_run./1000))/sum(rewardTS_run(2)./1000 - rewardTS_run(1)./1000);
            m10 = length(Restrict(spks,aversiveTS_run./1000))/sum(aversiveTS_run(2)./1000 - aversiveTS_run(1)./1000);
            
            %                 m9 = mean(movmean(tmp(InIntervals(bins,rewardTS_run./1000))./binSize,freq));
            %                 m10 = mean(movmean(tmp(InIntervals(bins,aversiveTS_run./1000))./binSize,freq));
            m11 = mean([m9 m10]);
            %                 m12 = mean(movmean(tmp(InIntervals(bins,WAKE))./binSize,freq));
            m12 = length(Restrict(spks,WAKE))/sum(WAKE(:,2)-WAKE(:,1));
            m13 = length(spks)/(segments.Var1(end)/1000);
            %                 if m4 > criteria_fr
            if logical(group_vHPC(ii,5)) %check if is pyr
                spiketrains_vHPC_pyr = [spiketrains_vHPC_pyr , zscore(tmp./binSize,0,1)];
                FR_vHPC_pyr = [FR_vHPC_pyr ; m4 m1 m2 m3 m8 m5 m6 m7 m11 m9 m10 m12 m13 iterator(6) iterator(5)];
                countP = countP + 1;
            else
                spiketrains_vHPC_int = [spiketrains_vHPC_int , zscore(tmp./binSize,0,1)];
                FR_vHPC_int = [FR_vHPC_int ; m4 m1 m2 m3 m8 m5 m6 m7 m11 m9 m10 m12 m13 iterator(6) iterator(5)];
                countI = countI + 1;
            end
            %                 end
            clear spks tmp m1 m2 m3 m4 m5 m6 m7 m8 m9 m10 m11 m12 m13
            %             end
        end
        clear freq limits
        n_SU_V = [n_SU_V ; countP countI]; clear countI countP
        %         if and(size(spiketrains_vHPC,2) >= criteria_n,size(spiketrains_dHPC,2) >= criteria_n)
        %Restricting bins inside each condition
        is.baseline.sws = InIntervals(bins,NREM_B);
        is.aversive.sws = InIntervals(bins,NREM_A);
        is.reward.sws = InIntervals(bins,NREM_R);
        
        is.baseline.rem = InIntervals(bins,REM_B);
        is.aversive.rem = InIntervals(bins,REM_A);
        is.reward.rem = InIntervals(bins,REM_R);
        
        is.aversive.run = InIntervals(bins,aversiveTS_run ./ 1000);
        is.reward.run = InIntervals(bins,rewardTS_run ./ 1000);
        
        
        clear spiketrains_dHPC_int spiketrains_dHPC_pyr spiketrains_vHPC_int spiketrains_vHPC_pyr
        t
    end
    tt
end


figure,
subplot(231),
scatter(FR_dHPC_pyr(:,12),FR_dHPC_pyr(:,1),5,'filled','MarkerFaceColor',[0.1986 0.7214 0.6310]),hold on
scatter(FR_dHPC_int(:,12),FR_dHPC_int(:,1),5,'filled','MarkerFaceColor',[0.6 0.0 0.6])
set(gca,'xscale','log','yscale','log'),xlim([0.001 100]),ylim([0.001 100]),
x = linspace(0.001,100);
y = linspace(0.001,100);
plot(x,y),ylabel('NREM rate (HZ)') , xlabel('WAKE rate (Hz)')
title('dHPC pyr + int')

subplot(232),
scatter(FR_dHPC_pyr(:,12),FR_dHPC_pyr(:,5),5,'filled','MarkerFaceColor',[0.1986 0.7214 0.6310]),hold on
scatter(FR_dHPC_int(:,12),FR_dHPC_int(:,5),5,'filled','MarkerFaceColor',[0.6 0.0 0.6])
set(gca,'xscale','log','yscale','log'),xlim([0.001 100]),ylim([0.001 100]),
plot(x,y),ylabel('REM rate (HZ)') , xlabel('WAKE rate (Hz)')
title('dHPC pyr + int')

subplot(233),
scatter(FR_dHPC_pyr(:,1),FR_dHPC_pyr(:,5),5,'filled','MarkerFaceColor',[0.1986 0.7214 0.6310]),hold on
scatter(FR_dHPC_int(:,1),FR_dHPC_int(:,5),5,'filled','MarkerFaceColor',[0.6 0.0 0.6])
set(gca,'xscale','log','yscale','log'),xlim([0.001 100]),ylim([0.001 100]),
plot(x,y),ylabel('REM rate (HZ)') , xlabel('NREM rate (Hz)')
title('dHPC pyr + int')

subplot(234),
scatter(FR_vHPC_pyr(:,12),FR_vHPC_pyr(:,1),5,'filled','MarkerFaceColor',[0.1986 0.7214 0.6310]),hold on
scatter(FR_vHPC_int(:,12),FR_vHPC_int(:,1),5,'filled','MarkerFaceColor',[0.6 0.0 0.6])
set(gca,'xscale','log','yscale','log'),xlim([0.001 100]),ylim([0.001 100]),
plot(x,y),ylabel('NREM rate (HZ)') , xlabel('WAKE rate (Hz)')
title('vHPC pyr + int')

subplot(235),
scatter(FR_vHPC_pyr(:,12),FR_vHPC_pyr(:,5),5,'filled','MarkerFaceColor',[0.1986 0.7214 0.6310]),hold on
scatter(FR_vHPC_int(:,12),FR_vHPC_int(:,5),5,'filled','MarkerFaceColor',[0.6 0.0 0.6])
set(gca,'xscale','log','yscale','log'),xlim([0.001 100]),ylim([0.001 100]),
plot(x,y),ylabel('REM rate (HZ)') , xlabel('WAKE rate (Hz)')
title('vHPC pyr + int')

subplot(236),
scatter(FR_vHPC_pyr(:,1),FR_vHPC_pyr(:,5),5,'filled','MarkerFaceColor',[0.1986 0.7214 0.6310]),hold on
scatter(FR_vHPC_int(:,1),FR_vHPC_int(:,5),5,'filled','MarkerFaceColor',[0.6 0.0 0.6])
set(gca,'xscale','log','yscale','log'),xlim([0.001 100]),ylim([0.001 100]),
plot(x,y),ylabel('REM rate (HZ)') , xlabel('NREM rate (Hz)')
title('vHPC pyr + int')

%% Ratio calculations and plot
ratio.NREM.dHPC.pyr = ((FR_dHPC_pyr(:,1)-FR_dHPC_pyr(:,12))./(FR_dHPC_pyr(:,12)+FR_dHPC_pyr(:,1)));
ratio.NREM.dHPC.int = ((FR_dHPC_int(:,1)-FR_dHPC_int(:,12))./(FR_dHPC_int(:,12)+FR_dHPC_int(:,1)));
ratio.NREM.vHPC.pyr = ((FR_vHPC_pyr(:,1)-FR_vHPC_pyr(:,12))./(FR_vHPC_pyr(:,12)+FR_vHPC_pyr(:,1)));
ratio.NREM.vHPC.int = ((FR_vHPC_int(:,1)-FR_vHPC_int(:,12))./(FR_vHPC_int(:,12)+FR_vHPC_int(:,1)));

ratio.REM.dHPC.pyr = ((FR_dHPC_pyr(:,5)-FR_dHPC_pyr(:,12))./(FR_dHPC_pyr(:,12)+FR_dHPC_pyr(:,5)));
ratio.REM.dHPC.int = ((FR_dHPC_int(:,5)-FR_dHPC_int(:,12))./(FR_dHPC_int(:,12)+FR_dHPC_int(:,5)));
ratio.REM.vHPC.pyr = ((FR_vHPC_pyr(:,5)-FR_vHPC_pyr(:,12))./(FR_vHPC_pyr(:,12)+FR_vHPC_pyr(:,5)));
ratio.REM.vHPC.int = ((FR_vHPC_int(:,5)-FR_vHPC_int(:,12))./(FR_vHPC_int(:,12)+FR_vHPC_int(:,5)));

ratio.sleep.dHPC.pyr = ((FR_dHPC_pyr(:,5)-FR_dHPC_pyr(:,1))./(FR_dHPC_pyr(:,1)+FR_dHPC_pyr(:,5)));
ratio.sleep.dHPC.int = ((FR_dHPC_int(:,5)-FR_dHPC_int(:,1))./(FR_dHPC_int(:,1)+FR_dHPC_int(:,5)));
ratio.sleep.vHPC.pyr = ((FR_vHPC_pyr(:,5)-FR_vHPC_pyr(:,1))./(FR_vHPC_pyr(:,1)+FR_vHPC_pyr(:,5)));
ratio.sleep.vHPC.int = ((FR_vHPC_int(:,5)-FR_vHPC_int(:,1))./(FR_vHPC_int(:,1)+FR_vHPC_int(:,5)));

figure,
subplot(231)
[f,x] = ksdensity(ratio.NREM.dHPC.pyr,'Bandwidth',0.2);
plot(x,f./max(f),'Color',[0.1986 0.7214 0.6310],'LineWidth',2);,hold on
xline(nanmedian(ratio.NREM.dHPC.pyr),'Color',[0.1986 0.7214 0.6310],'LineWidth',1)
clear f x
[f,x] = ksdensity(ratio.NREM.dHPC.int,'Bandwidth',0.2);
plot(x,f./max(f),'Color',[0.6 0.0 0.6],'LineWidth',2);xlim([-1 1])
xline(nanmedian(ratio.NREM.dHPC.int),'Color',[0.6 0.0 0.6],'LineWidth',1)
clear f x
xline(0,'--')
ylabel('Normalized count') , xlabel('Firing rate ratio')
title('dHPC pyr + int NREM')
p = signrank(ratio.NREM.dHPC.pyr)
p = signrank(ratio.NREM.dHPC.int)


subplot(232)
[f,x] = ksdensity(ratio.REM.dHPC.pyr,'Bandwidth',0.2);
plot(x,f./max(f),'Color',[0.1986 0.7214 0.6310],'LineWidth',2);,hold on
xline(nanmedian(ratio.REM.dHPC.pyr),'Color',[0.1986 0.7214 0.6310],'LineWidth',1)
clear f x
[f,x] = ksdensity(ratio.REM.dHPC.int,'Bandwidth',0.2);
plot(x,f./max(f),'Color',[0.6 0.0 0.6],'LineWidth',2);xlim([-1 1])
xline(nanmedian(ratio.REM.dHPC.int),'Color',[0.6 0.0 0.6],'LineWidth',1)
clear f x
xline(0,'--')
ylabel('Normalized count') , xlabel('Firing rate ratio')
title('dHPC pyr + int REM')
p = signrank(ratio.REM.dHPC.pyr)
p = signrank(ratio.REM.dHPC.int)

subplot(233)
[f,x] = ksdensity(ratio.sleep.dHPC.pyr,'Bandwidth',0.2);
plot(x,f./max(f),'Color',[0.1986 0.7214 0.6310],'LineWidth',2);,hold on
xline(nanmedian(ratio.sleep.dHPC.pyr),'Color',[0.1986 0.7214 0.6310],'LineWidth',1)
clear f x
[f,x] = ksdensity(ratio.sleep.dHPC.int,'Bandwidth',0.2);
plot(x,f./max(f),'Color',[0.6 0.0 0.6],'LineWidth',2);xlim([-1 1])
xline(nanmedian(ratio.sleep.dHPC.int),'Color',[0.6 0.0 0.6],'LineWidth',1)
clear f x
xline(0,'--')
ylabel('Normalized count') , xlabel('Firing rate ratio')
title('dHPC pyr + int REM-NREM')
p = signrank(ratio.sleep.dHPC.pyr)
p = signrank(ratio.sleep.dHPC.int)

subplot(234)
[f,x] = ksdensity(ratio.NREM.vHPC.pyr,'Bandwidth',0.2);
plot(x,f./max(f),'Color',[0.1986 0.7214 0.6310],'LineWidth',2);,hold on
xline(nanmedian(ratio.NREM.vHPC.pyr),'Color',[0.1986 0.7214 0.6310],'LineWidth',1)
clear f x
[f,x] = ksdensity(ratio.NREM.vHPC.int,'Bandwidth',0.2);
plot(x,f./max(f),'Color',[0.6 0.0 0.6],'LineWidth',2);xlim([-1 1])
xline(nanmedian(ratio.NREM.vHPC.int),'Color',[0.6 0.0 0.6],'LineWidth',1)
clear f x
xline(0,'--')
ylabel('Normalized count') , xlabel('Firing rate ratio')
title('vHPC pyr + int NREM')
p = signrank(ratio.NREM.vHPC.pyr)
p = signrank(ratio.NREM.vHPC.int)

subplot(235)
[f,x] = ksdensity(ratio.REM.vHPC.pyr,'Bandwidth',0.2);
plot(x,f./max(f),'Color',[0.1986 0.7214 0.6310],'LineWidth',2);,hold on
xline(nanmedian(ratio.REM.vHPC.pyr),'Color',[0.1986 0.7214 0.6310],'LineWidth',1)
clear f x
[f,x] = ksdensity(ratio.REM.vHPC.int,'Bandwidth',0.2);
plot(x,f./max(f),'Color',[0.6 0.0 0.6],'LineWidth',2);xlim([-1 1])
xline(nanmedian(ratio.REM.vHPC.int),'Color',[0.6 0.0 0.6],'LineWidth',1)
clear f x
xline(0,'--')
ylabel('Normalized count') , xlabel('Firing rate ratio')
title('vHPC pyr + int REM')
p = signrank(ratio.REM.vHPC.pyr)
p = signrank(ratio.REM.vHPC.int)

subplot(236)
[f,x] = ksdensity(ratio.sleep.vHPC.pyr,'Bandwidth',0.2);
plot(x,f./max(f),'Color',[0.1986 0.7214 0.6310],'LineWidth',2);,hold on
xline(nanmedian(ratio.sleep.vHPC.pyr),'Color',[0.1986 0.7214 0.6310],'LineWidth',1)
clear f x
[f,x] = ksdensity(ratio.sleep.vHPC.int,'Bandwidth',0.2);
plot(x,f./max(f),'Color',[0.6 0.0 0.6],'LineWidth',2);xlim([-1 1])
xline(nanmedian(ratio.sleep.vHPC.int),'Color',[0.6 0.0 0.6],'LineWidth',1)
clear f x
xline(0,'--')
ylabel('Normalized count') , xlabel('Firing rate ratio')
title('vHPC pyr + int REM-NREM')
p = signrank(ratio.sleep.vHPC.pyr)
p = signrank(ratio.sleep.vHPC.int)

%% NREM across conditions
figure,
subplot(231),
scatter(FR_dHPC_pyr(:,12),FR_dHPC_pyr(:,2),5,'filled','MarkerFaceColor',[0.1986 0.7214 0.6310]),hold on
scatter(FR_dHPC_int(:,12),FR_dHPC_int(:,2),5,'filled','MarkerFaceColor',[0.6 0.0 0.6])
set(gca,'xscale','log','yscale','log'),xlim([0.001 100]),ylim([0.001 100]),
x = linspace(0.001,100);
y = linspace(0.001,100);
plot(x,y)

subplot(232),
scatter(FR_dHPC_pyr(:,12),FR_dHPC_pyr(:,3),5,'filled','MarkerFaceColor',[0.1986 0.7214 0.6310]),hold on
scatter(FR_dHPC_int(:,12),FR_dHPC_int(:,3),5,'filled','MarkerFaceColor',[0.6 0.0 0.6])
set(gca,'xscale','log','yscale','log'),xlim([0.001 100]),ylim([0.001 100]),
plot(x,y)

subplot(233),
scatter(FR_dHPC_pyr(:,12),FR_dHPC_pyr(:,4),5,'filled','MarkerFaceColor',[0.1986 0.7214 0.6310]),hold on
scatter(FR_dHPC_int(:,12),FR_dHPC_int(:,4),5,'filled','MarkerFaceColor',[0.6 0.0 0.6])
set(gca,'xscale','log','yscale','log'),xlim([0.001 100]),ylim([0.001 100]),
plot(x,y)

subplot(234),
scatter(FR_vHPC_pyr(:,12),FR_vHPC_pyr(:,2),5,'filled','MarkerFaceColor',[0.1986 0.7214 0.6310]),hold on
scatter(FR_vHPC_int(:,12),FR_vHPC_int(:,2),5,'filled','MarkerFaceColor',[0.6 0.0 0.6])
set(gca,'xscale','log','yscale','log'),xlim([0.001 100]),ylim([0.001 100]),
x = linspace(0.001,100);
y = linspace(0.001,100);
plot(x,y)

subplot(235),
scatter(FR_vHPC_pyr(:,12),FR_vHPC_pyr(:,3),5,'filled','MarkerFaceColor',[0.1986 0.7214 0.6310]),hold on
scatter(FR_vHPC_int(:,12),FR_vHPC_int(:,3),5,'filled','MarkerFaceColor',[0.6 0.0 0.6])
set(gca,'xscale','log','yscale','log'),xlim([0.001 100]),ylim([0.001 100]),
plot(x,y)

subplot(236),
scatter(FR_vHPC_pyr(:,12),FR_vHPC_pyr(:,4),5,'filled','MarkerFaceColor',[0.1986 0.7214 0.6310]),hold on
scatter(FR_vHPC_int(:,12),FR_vHPC_int(:,4),5,'filled','MarkerFaceColor',[0.6 0.0 0.6])
set(gca,'xscale','log','yscale','log'),xlim([0.001 100]),ylim([0.001 100]),
plot(x,y)

ratio.NREM.dHPC.pyrB = ((FR_dHPC_pyr(:,2)-FR_dHPC_pyr(:,12))./(FR_dHPC_pyr(:,12)+FR_dHPC_pyr(:,2)));
ratio.NREM.dHPC.intB = ((FR_dHPC_int(:,2)-FR_dHPC_int(:,12))./(FR_dHPC_int(:,12)+FR_dHPC_int(:,2)));
ratio.NREM.vHPC.pyrB = ((FR_vHPC_pyr(:,2)-FR_vHPC_pyr(:,12))./(FR_vHPC_pyr(:,12)+FR_vHPC_pyr(:,2)));
ratio.NREM.vHPC.intB = ((FR_vHPC_int(:,2)-FR_vHPC_int(:,12))./(FR_vHPC_int(:,12)+FR_vHPC_int(:,2)));

ratio.NREM.dHPC.pyrR = ((FR_dHPC_pyr(:,3)-FR_dHPC_pyr(:,12))./(FR_dHPC_pyr(:,12)+FR_dHPC_pyr(:,3)));
ratio.NREM.dHPC.intR = ((FR_dHPC_int(:,3)-FR_dHPC_int(:,12))./(FR_dHPC_int(:,12)+FR_dHPC_int(:,3)));
ratio.NREM.vHPC.pyrR = ((FR_vHPC_pyr(:,3)-FR_vHPC_pyr(:,12))./(FR_vHPC_pyr(:,12)+FR_vHPC_pyr(:,3)));
ratio.NREM.vHPC.intR = ((FR_vHPC_int(:,3)-FR_vHPC_int(:,12))./(FR_vHPC_int(:,12)+FR_vHPC_int(:,3)));

ratio.NREM.dHPC.pyrA = ((FR_dHPC_pyr(:,4)-FR_dHPC_pyr(:,12))./(FR_dHPC_pyr(:,12)+FR_dHPC_pyr(:,4)));
ratio.NREM.dHPC.intA = ((FR_dHPC_int(:,4)-FR_dHPC_int(:,12))./(FR_dHPC_int(:,12)+FR_dHPC_int(:,4)));
ratio.NREM.vHPC.pyrA = ((FR_vHPC_pyr(:,4)-FR_vHPC_pyr(:,12))./(FR_vHPC_pyr(:,12)+FR_vHPC_pyr(:,4)));
ratio.NREM.vHPC.intA = ((FR_vHPC_int(:,4)-FR_vHPC_int(:,12))./(FR_vHPC_int(:,12)+FR_vHPC_int(:,4)));

figure,
subplot(121),boxplot([ratio.NREM.dHPC.pyrB ratio.NREM.dHPC.pyrR ratio.NREM.dHPC.pyrA ])
subplot(122),boxplot([ratio.NREM.vHPC.pyrB ratio.NREM.vHPC.pyrR ratio.NREM.vHPC.pyrA ])
figure,
subplot(121),boxplot([ratio.NREM.dHPC.intB ratio.NREM.dHPC.intR ratio.NREM.dHPC.intA ])
subplot(122),boxplot([ratio.NREM.vHPC.intB ratio.NREM.vHPC.intR ratio.NREM.vHPC.intA ])

% REM across conditions
figure,
subplot(231),
scatter(FR_dHPC_pyr(:,12),FR_dHPC_pyr(:,6),5,'filled','MarkerFaceColor',[0.1986 0.7214 0.6310]),hold on
scatter(FR_dHPC_int(:,12),FR_dHPC_int(:,6),5,'filled','MarkerFaceColor',[0.6 0.0 0.6])
set(gca,'xscale','log','yscale','log'),xlim([0.001 100]),ylim([0.001 100]),
x = linspace(0.001,100);
y = linspace(0.001,100);
plot(x,y)

subplot(232),
scatter(FR_dHPC_pyr(:,12),FR_dHPC_pyr(:,7),5,'filled','MarkerFaceColor',[0.1986 0.7214 0.6310]),hold on
scatter(FR_dHPC_int(:,12),FR_dHPC_int(:,7),5,'filled','MarkerFaceColor',[0.6 0.0 0.6])
set(gca,'xscale','log','yscale','log'),xlim([0.001 100]),ylim([0.001 100]),
plot(x,y)

subplot(233),
scatter(FR_dHPC_pyr(:,12),FR_dHPC_pyr(:,8),5,'filled','MarkerFaceColor',[0.1986 0.7214 0.6310]),hold on
scatter(FR_dHPC_int(:,12),FR_dHPC_int(:,8),5,'filled','MarkerFaceColor',[0.6 0.0 0.6])
set(gca,'xscale','log','yscale','log'),xlim([0.001 100]),ylim([0.001 100]),
plot(x,y)

subplot(234),
scatter(FR_vHPC_pyr(:,12),FR_vHPC_pyr(:,6),5,'filled','MarkerFaceColor',[0.1986 0.7214 0.6310]),hold on
scatter(FR_vHPC_int(:,12),FR_vHPC_int(:,6),5,'filled','MarkerFaceColor',[0.6 0.0 0.6])
set(gca,'xscale','log','yscale','log'),xlim([0.001 100]),ylim([0.001 100]),
x = linspace(0.001,100);
y = linspace(0.001,100);
plot(x,y)

subplot(235),
scatter(FR_vHPC_pyr(:,12),FR_vHPC_pyr(:,7),5,'filled','MarkerFaceColor',[0.1986 0.7214 0.6310]),hold on
scatter(FR_vHPC_int(:,12),FR_vHPC_int(:,7),5,'filled','MarkerFaceColor',[0.6 0.0 0.6])
set(gca,'xscale','log','yscale','log'),xlim([0.001 100]),ylim([0.001 100]),
plot(x,y)

subplot(236),
scatter(FR_vHPC_pyr(:,12),FR_vHPC_pyr(:,8),5,'filled','MarkerFaceColor',[0.1986 0.7214 0.6310]),hold on
scatter(FR_vHPC_int(:,12),FR_vHPC_int(:,8),5,'filled','MarkerFaceColor',[0.6 0.0 0.6])
set(gca,'xscale','log','yscale','log'),xlim([0.001 100]),ylim([0.001 100]),
plot(x,y)

ratio.REM.dHPC.pyrB = ((FR_dHPC_pyr(:,6)-FR_dHPC_pyr(:,12))./(FR_dHPC_pyr(:,12)+FR_dHPC_pyr(:,6)));
ratio.REM.dHPC.intB = ((FR_dHPC_int(:,6)-FR_dHPC_int(:,12))./(FR_dHPC_int(:,12)+FR_dHPC_int(:,6)));
ratio.REM.vHPC.pyrB = ((FR_vHPC_pyr(:,6)-FR_vHPC_pyr(:,12))./(FR_vHPC_pyr(:,12)+FR_vHPC_pyr(:,6)));
ratio.REM.vHPC.intB = ((FR_vHPC_int(:,6)-FR_vHPC_int(:,12))./(FR_vHPC_int(:,12)+FR_vHPC_int(:,6)));

ratio.REM.dHPC.pyrR = ((FR_dHPC_pyr(:,7)-FR_dHPC_pyr(:,12))./(FR_dHPC_pyr(:,12)+FR_dHPC_pyr(:,7)));
ratio.REM.dHPC.intR = ((FR_dHPC_int(:,7)-FR_dHPC_int(:,12))./(FR_dHPC_int(:,12)+FR_dHPC_int(:,7)));
ratio.REM.vHPC.pyrR = ((FR_vHPC_pyr(:,7)-FR_vHPC_pyr(:,12))./(FR_vHPC_pyr(:,12)+FR_vHPC_pyr(:,7)));
ratio.REM.vHPC.intR = ((FR_vHPC_int(:,7)-FR_vHPC_int(:,12))./(FR_vHPC_int(:,12)+FR_vHPC_int(:,7)));

ratio.REM.dHPC.pyrA = ((FR_dHPC_pyr(:,8)-FR_dHPC_pyr(:,12))./(FR_dHPC_pyr(:,12)+FR_dHPC_pyr(:,8)));
ratio.REM.dHPC.intA = ((FR_dHPC_int(:,8)-FR_dHPC_int(:,12))./(FR_dHPC_int(:,12)+FR_dHPC_int(:,8)));
ratio.REM.vHPC.pyrA = ((FR_vHPC_pyr(:,8)-FR_vHPC_pyr(:,12))./(FR_vHPC_pyr(:,12)+FR_vHPC_pyr(:,8)));
ratio.REM.vHPC.intA = ((FR_vHPC_int(:,8)-FR_vHPC_int(:,12))./(FR_vHPC_int(:,12)+FR_vHPC_int(:,8)));

figure,
subplot(121),boxplot([ratio.REM.dHPC.pyrB ratio.REM.dHPC.pyrR ratio.REM.dHPC.pyrA ])
subplot(122),boxplot([ratio.REM.vHPC.pyrB ratio.REM.vHPC.pyrR ratio.REM.vHPC.pyrA ])
figure,
subplot(121),boxplot([ratio.REM.dHPC.intB ratio.REM.dHPC.intR ratio.REM.dHPC.intA ])
subplot(122),boxplot([ratio.REM.vHPC.intB ratio.REM.vHPC.intR ratio.REM.vHPC.intA ])



%% scatter of mean FR
figure,
grps = [ones(size(FR_vHPC_pyr,1),1) ; ones(size(FR_vHPC_int,1),1)*2 ; ones(size(FR_dHPC_pyr,1),1)*3 ; ones(size(FR_dHPC_int,1),1)*4];
x = [FR_vHPC_pyr(:,13) ; FR_vHPC_int(:,13) ; FR_dHPC_pyr(:,13) ; FR_dHPC_int(:,13)];
scatter(grps,x,"filled",'jitter','on', 'jitterAmount',0.1),xlim([0 5]),hold on
scatter([1 2 3 4] , [nanmean(FR_vHPC_pyr(:,13)) ; nanmean(FR_vHPC_int(:,13)) ; nanmean(FR_dHPC_pyr(:,13)) ; nanmean(FR_dHPC_int(:,13))],"filled")
% set(gca, 'YScRale', 'log')
boxplot(x,grps)
ylabel('Firing Rate (Hz)')



x = [FR_dHPC_pyr(:,13) ; FR_dHPC_int(:,13) ; FR_vHPC_pyr(:,13) ; FR_vHPC_int(:,13)];

group = [ones(sum([length(FR_dHPC_pyr) , length(FR_dHPC_int)]),1) ; ones(sum([length(FR_vHPC_pyr) , length(FR_vHPC_int)]),1)*2];
group = [group , [ones(length(FR_dHPC_pyr),1) ; ones(length(FR_dHPC_int),1)*2 ; ones(length(FR_vHPC_pyr),1) ; ones(length(FR_vHPC_int),1)*2]];


[p,tbl,stats] = anovan(x,group,'model','full','varnames',{'structure','type'})
[results,~,~,gnames] = multcompare(stats,"Dimension",[1 2]);

figure,
grps = [ones(size(FR_vHPC_pyr,1),1) ; ones(size(FR_vHPC_int,1),1)*2 ; ones(size(FR_dHPC_pyr,1),1)*3 ; ones(size(FR_dHPC_int,1),1)*4];
x = [FR_vHPC_pyr(:,15) ; FR_vHPC_int(:,15) ; FR_dHPC_pyr(:,15) ; FR_dHPC_int(:,15)];
scatter(grps,x,"filled",'jitter','on', 'jitterAmount',0.1),xlim([0 5]),hold on
scatter([1 2 3 4] , [nanmean(FR_vHPC_pyr(:,15)) ; nanmean(FR_vHPC_int(:,15)) ; nanmean(FR_dHPC_pyr(:,15)) ; nanmean(FR_dHPC_int(:,15))],"filled")
% set(gca, 'YScRale', 'log')
ylabel('Burst Index')


x = [FR_dHPC_pyr(:,15) ; FR_dHPC_int(:,15) ; FR_vHPC_pyr(:,15) ; FR_vHPC_int(:,15)];

group = [ones(sum([length(FR_dHPC_pyr) , length(FR_dHPC_int)]),1) ; ones(sum([length(FR_vHPC_pyr) , length(FR_vHPC_int)]),1)*2];
group = [group , [ones(length(FR_dHPC_pyr),1) ; ones(length(FR_dHPC_int),1)*2 ; ones(length(FR_vHPC_pyr),1) ; ones(length(FR_vHPC_int),1)*2]];


[p,tbl,stats] = anovan(x,group,'model','full','varnames',{'structure','type'})
[results,~,~,gnames] = multcompare(stats,"Dimension",[1 2]);






%% plot 100% stacked chart
n = sum([length(FR_vHPC_pyr) , length(FR_vHPC_int)]);
p = (length(FR_vHPC_pyr)/n)*100;
bar(1,100),hold on
bar(1,p)

n = sum([length(FR_dHPC_pyr) , length(FR_dHPC_int)]);
p = (length(FR_dHPC_pyr)/n)*100;
bar(2,100),hold on
bar(2,p)
