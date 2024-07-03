clear
clc
close all

%% Parameters
path = {'\\Maryjackson\e\Rat127\Ephys\pyr';'\\Maryjackson\e\Rat128\Ephys\in_pyr\ready';'\\Maryjackson\e\Rat103\usable';'\\Maryjackson\e\Rat132\recordings\in_pyr'; '\\Maryjackson\e\Rat165\in_pyr'; ...
    '\\Maryjackson\e\Rat126\Ephys\in_Pyr'};%List of folders from the path
%Sleep
time_criteria = 1; % minimal time to include a NREM epoch (in min)

% for SU
criteria_fr = 0; %criteria to include or not a SU into the analysis
criteria_n = [3 3]; % minimal number of neurons from each structure [vHPC dHPC]
criteria_type = 0; %criteria for celltype (0:pyr, 1:int, 2:all)
binSize = 0.025; %for qssemblie detection qnd qxctivity strength
n_SU_V = 0;
n_SU_D = 0;

win = 300; % time window for bin construction

% Behavior
minimal_speed = 7; % minimal speed to detect quite periods
minimal_speed_time = 2; % minimal time to detect quite periods

th = 3; % threshold for peak detection
s = true; % true if I wanna save maps at each folder

Number_of_assemblies.aversive = [];
Number_of_assemblies.reward = [];

map.bothA.aversive = []; map.bothA.reward = [];
map.bothR.aversive = []; map.bothR.reward = [];
map.dHPCA.aversive = []; map.dHPCA.reward = [];
map.dHPCR.aversive = []; map.dHPCR.reward = [];
map.vHPCA.aversive = []; map.vHPCA.reward = [];
map.vHPCR.aversive = []; map.vHPCR.reward = [];


Within.bothA.aversive = []; Within.bothA.reward = [];
Within.bothR.aversive = []; Within.bothR.reward = [];
Within.dHPCA.aversive = []; Within.dHPCA.reward = [];
Within.dHPCR.aversive = []; Within.dHPCR.reward = [];
Within.vHPCA.aversive = []; Within.vHPCA.reward = [];
Within.vHPCR.aversive = []; Within.vHPCR.reward = [];

Between.dHPCA = []; Between.dHPCR = [];
Between.vHPCA = []; Between.vHPCR = [];

% Sacar el filtro que puse del FR en el counts de neuronas
%% Main loop, to iterate across sessions
for tt = 1:length(path)
    %List of folders from the path
    files = dir(path{tt});
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
    clear files dirFlags
    num_assembliesR = [];
    num_assembliesA = [];
    for t = 1 : length(subFolders)-2
        disp(['-- Initiating analysis of folder #' , num2str(t) , ' from rat #',num2str(tt) , ' --'])
        session = [subFolders(t+2).folder,'\',subFolders(t+2).name];
        cd(session)
        
        %% Loading TS of the sessions
        disp('Uploading session time stamps')
        load('behavioral_data.mat')
        load('session_organization.mat')

        %% Load sleep states
        disp('Uploading sleep scoring')
        x = dir([cd,'\*-states.mat']);    states = load([cd,'\',x.name]);    states = states.states;
        REM.all = ToIntervals(states==5);    NREM.all = ToIntervals(states==3);    WAKE.all = ToIntervals(states==1);
        %         REM.all(REM.all(:,2)-REM.all(:,1)<60,:) = [];
        clear x states
        
        NREM.all(NREM.all(:,2)-NREM.all(:,1)<60*time_criteria,:)=[];
        NREM.baseline = Restrict(NREM.all,baselineTS./1000);
        NREM.aversive = Restrict(NREM.all,aversiveTS./1000);
        NREM.reward = Restrict(NREM.all,rewardTS./1000);
        
        REM.baseline = Restrict(REM.all,baselineTS./1000);
        REM.aversive = Restrict(REM.all,aversiveTS./1000);
        REM.reward = Restrict(REM.all,rewardTS./1000);
        
        WAKE.baseline = Restrict(WAKE.all,baselineTS./1000);
        WAKE.aversive = Restrict(WAKE.all,aversiveTS./1000);
        WAKE.reward = Restrict(WAKE.all,rewardTS./1000);
        
        %% Spikes
        % Load Units
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
                spks_dHPC = spks(ismember(spks(:,1),group_dHPC(:,1)),:); %keep spks from good clusters
            else
                spks_vHPC = spks(ismember(spks(:,1),group_vHPC(:,1)),:); %keep spks from good clusters
            end
        end
        clear z
        spks_vHPC(:,2) = double(spks_vHPC(:,2))./20000;
        spks_dHPC(:,2) = double(spks_dHPC(:,2))./20000;
        spks(:,2) = double(spks(:,2))./20000;
        
        % Selection of celltype to analyze
        if criteria_type == 0 %pyr
            cellulartype = [K(:,1) , K(:,4)];
        elseif criteria_type == 1 % int
            cellulartype = [K(:,1) , not(K(:,4))];
        elseif criteria_type == 2 % all
            cellulartype = [K(:,1) , ones(length(K),1)];
        end
        
        %% Counting the Number f SU
        numberD = 0;
        clusters.all = [];
        clusters.dHPC = [];
        for ii=1:size(group_dHPC,1)
            cluster = group_dHPC(ii,1);
            Cellulartype = logical(cellulartype(cellulartype(:,1) == cluster,2));
            if Cellulartype
                a = length(Restrict(spks(spks(:,1)==cluster,2),aversiveTS_run./1000)) / ((aversiveTS_run(2)-aversiveTS_run(1))/1000);
                r = length(Restrict(spks(spks(:,1)==cluster,2),rewardTS_run./1000)) / ((rewardTS_run(2)-rewardTS_run(1))/1000);
                if or(a > criteria_fr ,  r > criteria_fr)
                    numberD = numberD+1;
                    clusters.all = [clusters.all ; cluster];
                    clusters.dHPC = [clusters.dHPC ; cluster];
                end
            end
            clear tmp cluster Cellulartype fr1 fr2 fr3 fr4 fr5 r a
        end
        
        numberV = 0;
        clusters.vHPC = [];
        for ii=1:size(group_vHPC,1)
            cluster = group_vHPC(ii,1);
            Cellulartype = logical(cellulartype(cellulartype(:,1) == cluster,2));
            if Cellulartype
                a = length(Restrict(spks(spks(:,1)==cluster,2),aversiveTS_run./1000)) / ((aversiveTS_run(2)-aversiveTS_run(1))/1000);
                r = length(Restrict(spks(spks(:,1)==cluster,2),rewardTS_run./1000)) / ((rewardTS_run(2)-rewardTS_run(1))/1000);
                if or(a > criteria_fr ,  r > criteria_fr)
                    numberV = numberV+1;
                    clusters.all = [clusters.all ; cluster];
                    clusters.vHPC = [clusters.vHPC ; cluster];
                end
            end
            clear tmp cluster Cellulartype fr1 fr2 fr3 fr4 fr5 r a
        end
        clear freq limits
        clear camara shock rightvalve leftvalve
        clear ejeX ejeY dX dY dX_int dY_int
        
        criteria_n = [3 3];
        if or(numberD >= 3 , numberV >= 3)
            disp('Lets go for the assemblies')
            if isfile('separated_assemblies.mat')
                load('separated_assemblies.mat')
                
                % SpikeTrains construction
                limits = [0 segments.Var1(end)/1000];
                events = [];
                if numberD >= 3
                    [Spikes.dHPC , bins , Clusters.dHPC] = spike_train_construction([spks_dHPC], clusters.dHPC, cellulartype, binSize, limits, events, false, true);
                end 
                if numberV >= 3
                    [Spikes.vHPC , bins , Clusters.vHPC] = spike_train_construction([spks_vHPC], clusters.vHPC, cellulartype, binSize, limits, events, false, true);
                    
                end 
                clear limits events
                
                if isfield(patterns,'aversive')
                    % --- AVERSIVE ---
                    %dHPC
                    if isfield(patterns.aversive,'dHPC')
                        if and(not(isempty(patterns.aversive.dHPC)), size(patterns.aversive.dHPC,1)>2)
                        pos = [behavior.pos.aversive(:,1:2) ; behavior.pos.reward(:,1:2)];
                        [~, xx] = sort(pos(:,1));
                        pos = pos(xx,:);
                        events = cell(2,1);
                        events{1} = movement.aversive;
                        events{2}  = movement.reward;
                        
                        [Maps pc between within] = FiringMap_Assemblies(patterns.aversive.dHPC , logical(ones(1,size(patterns.aversive.dHPC,2))) , [bins' , Spikes.dHPC] , th , pos , events , 60 , true , true);
                        clear pos x xx
                        
                        
                        pc = or(pc.cond1 , pc.cond2); %logical to select maps to save
                        
                        map.dHPCA.aversive = [map.dHPCA.aversive ; Maps.cond1(pc,:)];
                        map.dHPCA.reward = [map.dHPCA.reward ; Maps.cond2(pc,:)];
                        
                        Within.dHPCA.aversive = [Within.dHPCA.aversive ; within.cond1(pc,:)];
                        Within.dHPCA.reward = [Within.dHPCA.reward ; within.cond2(pc,:)];
                        
                        Between.dHPCA = [Between.dHPCA ; between(pc,:)];
                        
                        clear within between pc Maps pos x xx events
                        end
                    end
                    
                    % vHPC
                    if isfield(patterns.aversive,'vHPC')
                        if and(not(isempty(patterns.aversive.vHPC)), size(patterns.aversive.vHPC,1)>2)
                        pos = [behavior.pos.aversive(:,1:2) ; behavior.pos.reward(:,1:2)];
                        [x xx] = sort(pos(:,1));
                        pos = pos(xx,:);
                        events = cell(2,1);
                        events{1} = movement.aversive;
                        events{2}  = movement.reward;
                        
                        [Maps pc between within] = FiringMap_Assemblies(patterns.aversive.vHPC , logical(ones(1,size(patterns.aversive.vHPC,2))) , [bins' , Spikes.vHPC] , th , pos , events , 60 , true , true);
                        clear pos x xx
                        
                        
                        
                        pc = or(pc.cond1 , pc.cond2); %logical to select maps to save
                        
                        map.vHPCA.aversive = [map.vHPCA.aversive ; Maps.cond1(pc,:)];
                        map.vHPCA.reward = [map.vHPCA.reward ; Maps.cond2(pc,:)];
                        
                        Within.vHPCA.aversive = [Within.vHPCA.aversive ; within.cond1(pc,:)];
                        Within.vHPCA.reward = [Within.vHPCA.reward ; within.cond2(pc,:)];
                        
                        Between.vHPCA = [Between.vHPCA ; between(pc,:)];
                        
                        clear within between pc Maps pos x xx events
                        end 
                    end
                end
                
                
                
                if isfield(patterns,'reward')
                    % --- Reward ---
                    %dHPC
                    if isfield(patterns.reward,'dHPC')
                        if and(not(isempty(patterns.reward.dHPC)), size(patterns.reward.dHPC,1)>2)
                        pos = [behavior.pos.reward(:,1:2) ; behavior.pos.aversive(:,1:2)];
                        [~, xx] = sort(pos(:,1));
                        pos = pos(xx,:);
                        events = cell(2,1);
                        events{1} = movement.reward;
                        events{2}  = movement.aversive;
                        
                        [Maps pc between within] = FiringMap_Assemblies(patterns.reward.dHPC , logical(ones(1,size(patterns.reward.dHPC,2))) , [bins' , Spikes.dHPC] , th , pos , events , 60 , true , true);
                        clear pos x xx
                        
                        
                        pc = or(pc.cond1 , pc.cond2); %logical to select maps to save
                        
                        map.dHPCR.reward = [map.dHPCR.reward ; Maps.cond1(pc,:)];
                        map.dHPCR.aversive = [map.dHPCR.aversive ; Maps.cond2(pc,:)];
                        
                        Within.dHPCR.reward = [Within.dHPCR.reward ; within.cond1(pc,:)];
                        Within.dHPCR.aversive = [Within.dHPCR.aversive ; within.cond2(pc,:)];
                        
                        Between.dHPCR = [Between.dHPCR ; between(pc,:)];
                        
                        clear within between pc Maps pos x xx events
                        end
                    end
                    
                    % vHPC
                    if isfield(patterns.reward,'vHPC')
                        if and(not(isempty(patterns.reward.vHPC)), size(patterns.reward.vHPC,1)>2)
                        pos = [behavior.pos.reward(:,1:2) ; behavior.pos.aversive(:,1:2)];
                        [x xx] = sort(pos(:,1));
                        pos = pos(xx,:);
                        events = cell(2,1);
                        events{1} = movement.reward;
                        events{2}  = movement.aversive;
                        
                        [Maps pc between within] = FiringMap_Assemblies(patterns.reward.vHPC , logical(ones(1,size(patterns.reward.vHPC,2))) , [bins' , Spikes.vHPC] , th , pos , events , 60 , true , true);
                        clear pos x xx
                        
                        
                        pc = or(pc.cond1 , pc.cond2); %logical to select maps to save
                        
                        map.vHPCR.reward = [map.vHPCR.reward ; Maps.cond1(pc,:)];
                        map.vHPCR.aversive = [map.vHPCR.aversive ; Maps.cond2(pc,:)];
                        
                        Within.vHPCR.reward = [Within.vHPCR.reward ; within.cond1(pc,:)];
                        Within.vHPCR.aversive = [Within.vHPCR.aversive ; within.cond2(pc,:)];
                        
                        Between.vHPCR = [Between.vHPCR ; between(pc,:)];
                        
                        clear within between pc Maps pos x xx events
                         end 
                    end
                end                
                
            end
 
        end
        disp(' ')
        
        %% Clear Work space
        clear A aversiveTS aversiveTS_run baselineTS rewardTS rewardTS_run
        clear behavior bins Cell_type_classification cellulartype cond
        clear is K Kinfo group_dHPC group_vHPC matrixC matrixD matrixV
        clear NREM REM WAKE segmentation tmp cond config
        clear spiketrains_dHPC spiketrains_vHPC opts MUA
        clear patterns Thresholded i  ii numberD numberV movement cross crossN
        clear Spikes bins Clusters Shocks_filt Rewards_filt config n_SU_D n_SU_V
        clear clusters coordinated coordinated_ripple_bursts coordinatedV
        clear cooridnated_event coordinatedV_refined coordinatedV_refined
    end
    
    Number_of_assemblies.aversive = [Number_of_assemblies.aversive ; sum(num_assembliesA)];
    Number_of_assemblies.reward = [Number_of_assemblies.reward ; sum(num_assembliesR)];
    clear num_assembliesA num_assembliesR
    
end
%% Plotting section
% The following section is used to plot the outputs of the script
% Please, be sure that you undersand the structure of the matrix you are
% intending to plot.

%% Both
% figure
% x = map.bothA.aversive - min(map.bothA.aversive,[],2);
% x = x./max(x,[],2);
% 
% y = map.bothA.reward - min(map.bothA.reward,[],2);
% y = y./max(y,[],2);
% 
% [c i] = max(x,[],2);
% [c i] = sort(i);
% 
% subplot(121),imagesc(x(i,:))
% subplot(122),imagesc(y(i,:))
% 
% 
% figure
% x = map.bothR.aversive - min(map.bothR.aversive,[],2);
% x = x./max(x,[],2);
% 
% y = map.bothR.reward - min(map.bothR.reward,[],2);
% y = y./max(y,[],2);
% 
% [c i] = max(y,[],2);
% [c i] = sort(i);
% 
% subplot(121),imagesc(x(i,:))
% subplot(122),imagesc(y(i,:))

%% dHPC
figure
x = map.dHPCA.aversive - min(map.dHPCA.aversive,[],2);
x = x./max(x,[],2);

y = map.dHPCA.reward - min(map.dHPCA.reward,[],2);
y = y./max(y,[],2);

[c i] = max(x,[],2);
[c i] = sort(i);

subplot(121),imagesc(x(i,:)); colormap 'gray'; title('Aversive')
subplot(122),imagesc(y(i,:));colormap 'gray';title('Reward')
sgtitle('dHPC-Aversive');

figure
x = map.dHPCR.aversive - min(map.dHPCR.aversive,[],2);
x = x./max(x,[],2);

y = map.dHPCR.reward - min(map.dHPCR.reward,[],2);
y = y./max(y,[],2);

[c i] = max(y,[],2);
[c i] = sort(i);
sgtitle('dHPC-Reward');
subplot(121),imagesc(x(i,:));colormap 'gray'; title('Aversive')
subplot(122),imagesc(y(i,:));colormap 'gray';title('Reward')

%% vHPC
figure
x = map.vHPCA.aversive - min(map.vHPCA.aversive,[],2);
x = x./max(x,[],2);

y = map.vHPCA.reward - min(map.vHPCA.reward,[],2);
y = y./max(y,[],2);

[c i] = max(x,[],2);
[c i] = sort(i);

subplot(121),imagesc(x(i,:)); colormap 'gray'; title('Aversive')
subplot(122),imagesc(y(i,:));colormap 'gray';title('Reward')
sgtitle('vHPC-Aversive');


figure
x = map.vHPCR.aversive - min(map.vHPCR.aversive,[],2);
x = x./max(x,[],2);

y = map.vHPCR.reward - min(map.vHPCR.reward,[],2);
y = y./max(y,[],2);

[c i] = max(y,[],2);
[c i] = sort(i);

subplot(121),imagesc(x(i,:)); colormap 'gray'; title('Aversive')
subplot(122),imagesc(y(i,:));colormap 'gray';title('Reward')
sgtitle('vHPC-Reward');


%% Plot parameters
% Spatial Correlation
% figure
% subplot(3,2,1)
% y = [Between.bothA(:,1) ; Within.bothA.reward(:,1) ; Within.bothA.aversive(:,1)];
% x = [ones(length(Between.bothA(:,1)),1) ; ones(length(Within.bothA.reward(:,1)),1)*2 ; ones(length(Within.bothA.aversive(:,1)),1)*3];
% scatter(x,y,'filled','jitter','on', 'jitterAmount',0.1),hold on,xlim([0 4]),ylim([-0.55 1.05])
% scatter([1 2 3],[nanmedian(Between.bothA(:,1)) , nanmedian( Within.bothA.reward(:,1)), nanmedian( Within.bothA.aversive(:,1))],'filled')
% % boxplot(y,x),ylim([-0.55 1.05])
% 
% [p,tbl,stats] = kruskalwallis(y,x);
% c = multcompare(stats)
% 
% 
% subplot(3,2,2)
% y = [Between.bothR(:,1) ; Within.bothR.reward(:,1) ; Within.bothR.aversive(:,1)];
% x = [ones(length(Between.bothR(:,1)),1) ; ones(length(Within.bothR.reward(:,1)),1)*2 ; ones(length(Within.bothR.aversive(:,1)),1)*3];
% scatter(x,y,'filled','jitter','on', 'jitterAmount',0.1),hold on,xlim([0 4]),ylim([-0.55 1.05])
% scatter([1 2 3],[nanmedian(Between.bothR(:,1)) , nanmedian( Within.bothR.reward(:,1)), nanmedian( Within.bothR.aversive(:,1))],'filled')
% % boxplot(y,x),ylim([-0.55 1.05])
% 
% [p,tbl,stats] = kruskalwallis(y,x);
% c = multcompare(stats)

color = [.3,.3,.3];

figure(1);clf; hold on
subplot(1,2,1)

y = [Between.dHPCA(:,1);Within.dHPCA.aversive(:,1); Within.dHPCA.reward(:,1)];
x = [ones(length(Between.dHPCA(:,1)),1) ; ones(length(Within.dHPCA.aversive(:,1)),1)*2 ; ones(length(Within.dHPCA.reward(:,1)),1)*3];
scatter(x,y,[],color,"filled",'jitter','on', 'jitterAmount',0.1),hold on,xlim([0 4]),ylim([-0.55 1.05])
ylabel('Spatial correlation');
xticks([1 2 3])
xticklabels({'Bet', 'WA', 'WR'});
scatter([1 2 3],[nanmedian(Between.dHPCA(:,1)) , nanmedian( Within.dHPCA.aversive(:,1)), nanmedian( Within.dHPCA.reward(:,1))],'filled')
% boxplot(y,x),ylim([-0.55 1.05])


% Paired non-parametric anova - Friedman test
%Prepare data to test
data = [Between.dHPCA(:,1), Within.dHPCA.aversive(:,1),Within.dHPCA.reward(:,1)];
data(any(isnan(data), 2), :) = [];

[p,~,stats] =  friedman(data,1); 
c = multcompare(stats);
tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])

subplot(1,2,2)
y = [Between.dHPCR(:,1) ; Within.dHPCR.aversive(:,1);Within.dHPCR.reward(:,1)];
x = [ones(length(Between.dHPCR(:,1)),1) ; ones(length(Within.dHPCR.aversive(:,1)),1)*2 ; ones(length(Within.dHPCR.reward(:,1)),1)*3];
scatter(x,y,[],color,"filled",'jitter','on', 'jitterAmount',0.1),hold on,xlim([0 4]),ylim([-0.55 1.05])
scatter([1 2 3],[nanmedian(Between.dHPCR(:,1)) , nanmedian( Within.dHPCR.aversive(:,1)), nanmedian( Within.dHPCR.reward(:,1))],'filled')
% boxplot(y,x),ylim([-0.55 1.05])
ylabel('Spatial correlation');
xticks([1 2 3])
xticklabels({'Bet', 'WA', 'WR'});

% Paired non-parametric anova - Friedman test
%Prepare data to test
data = [Between.dHPCR(:,1), Within.dHPCR.reward(:,1), Within.dHPCR.aversive(:,1)];
data(any(isnan(data), 2), :) = [];

[p,~,stats] =  friedman(data,1); 
c = multcompare(stats);
tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])

figure(2);clf; hold on
subplot(1,2,1)
y = [Between.vHPCA(:,1) ; Within.vHPCA.aversive(:,1) ; Within.vHPCA.reward(:,1)];
x = [ones(length(Between.vHPCA(:,1)),1) ; ones(length(Within.vHPCA.aversive(:,1)),1)*2 ; ones(length(Within.vHPCA.reward(:,1)),1)*3];
scatter(x,y,[],color,"filled",'jitter','on', 'jitterAmount',0.1),hold on,xlim([0 4]),ylim([-0.55 1.05])
scatter([1 2 3],[nanmedian(Between.vHPCA(:,1)) , nanmedian( Within.vHPCA.aversive(:,1)), nanmedian( Within.vHPCA.reward(:,1))],'filled')
% boxplot(y,x),ylim([-0.55 1.05])
ylabel('Spatial correlation');
xticks([1 2 3])
xticklabels({'Bet', 'WA', 'WR'});

% stats
data = [Between.vHPCA(:,1), Within.vHPCA.reward(:,1), Within.vHPCA.aversive(:,1)];
data(any(isnan(data), 2), :) = [];

[p,~,stats] =  friedman(data,1); 
c = multcompare(stats);
tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])

subplot(1,2,2)
y = [Between.vHPCR(:,1) ; Within.vHPCR.aversive(:,1) ; Within.vHPCR.reward(:,1)];
x = [ones(length(Between.vHPCR(:,1)),1) ; ones(length(Within.vHPCR.aversive(:,1)),1)*2 ; ones(length(Within.vHPCR.reward(:,1)),1)*3];
scatter(x,y,[],color,"filled",'jitter','on', 'jitterAmount',0.1),hold on,xlim([0 4]),ylim([-0.55 1.05])
scatter([1 2 3],[nanmedian(Between.vHPCR(:,1)) , nanmedian( Within.vHPCR.aversive(:,1)), nanmedian( Within.vHPCR.reward(:,1))],'filled')
% boxplot(y,x),ylim([-0.55 1.05])
ylabel('Spatial correlation');
xticks([1 2 3])
xticklabels({'Bet', 'WA', 'WR'});
sgtitle('vHPC Assemblies')


data = [Between.vHPCR(:,1), Within.vHPCR.reward(:,1), Within.vHPCR.aversive(:,1)];
data(any(isnan(data), 2), :) = [];

[p,~,stats] =  friedman(data,1); 
c = multcompare(stats);
tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])


%% Comparing Between groups across type of assemblies
figure
subplot(121)
y = [Between.bothA(:,1) ; Between.dHPCA(:,1) ;  Between.vHPCA(:,1)];
x = [ones(length(Between.bothA(:,1)),1) ; ones(length(Between.dHPCA(:,1)),1)*2 ; ones(length(Between.vHPCA(:,1)),1)*3];
scatter(x,y,'filled','jitter','on', 'jitterAmount',0.1),hold on,xlim([0 4]),ylim([-0.55 1.05])
scatter([1 2 3],[nanmedian(Between.bothA(:,1)) , nanmedian(Between.dHPCA(:,1)), nanmedian(Between.vHPCA(:,1))],'filled')

[p,tbl,stats] = kruskalwallis(y,x);
c = multcompare(stats)


subplot(122)
y = [Between.bothR(:,1) ; Between.dHPCR(:,1) ;  Between.vHPCR(:,1)];
x = [ones(length(Between.bothR(:,1)),1) ; ones(length(Between.dHPCR(:,1)),1)*2 ; ones(length(Between.vHPCR(:,1)),1)*3];
scatter(x,y,'filled','jitter','on', 'jitterAmount',0.1),hold on,xlim([0 4]),ylim([-0.55 1.05])
scatter([1 2 3],[nanmedian(Between.bothR(:,1)) , nanmedian(Between.dHPCR(:,1)), nanmedian(Between.vHPCR(:,1))],'filled')

[p,tbl,stats] = kruskalwallis(y,x);
c = multcompare(stats)