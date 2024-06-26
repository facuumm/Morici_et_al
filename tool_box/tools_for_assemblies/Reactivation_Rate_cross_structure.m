function [Reactivation , Cumulative] = Reactivation_Rate_cross_structure(path)
% This function a CCG between dorsal and ventral SUs. It does it during RUN
% periods, and if they are coordinated (>95th percentile of a surrogate)
% they store their CCG during PRE and POST sleep.
%
% INPUTS
% path: cell, inside each cell you should put the path of each session you
%       want to analyze
%
% OUTPUT
% cross: structure, it contains the CCG for aversive and reward
%                .pre.aversive
%                    .reward
%
%                .post.aversive
%                     .reward
%
%
%                .run.reward
%                .run.aversive.up
%                             .down
%
% time: column vector, it contains the time axis for the CCG plot.
%
% other functions: CCG from FMA toolbox
% Morci Juan Facundo 01/2024

% Variabbles
minimal_speed = 7;
minimal_speed_time = 2;

% for peak detection
tempo = 45*60; %time to restrict the sleep epochs
th = 1; % SD

% Initialization of structures that wil contain the outputs
Reactivation.aversive = [];
Reactivation.reward = [];

Cumulative.aversive = [];
Cumulative.reward= [];

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
        load('session_organization.mat')
        
        %% Awake
        disp('Uploading behavioral data')
        % Load digitalin.mat
        load('behavioral_data.mat')

        %% load sleep states
        disp('Uploading sleep scoring')
        x = dir([cd,'\*-states.mat']);    states = load([cd,'\',x.name]);    states = states.states;
        REM.all = ToIntervals(states==5);    NREM.all = ToIntervals(states==3);    WAKE.all = ToIntervals(states==1);
        %         REM.all(REM.all(:,2)-REM.all(:,1)<60,:) = [];
        clear x states
        
        %         NREM.all(NREM.all(:,2)-NREM.all(:,1)<60*time_criteria,:)=[];
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
        cellulartype = [K(:,1) , K(:,4)];
        
        
        %% Counting the Number f SU
        numberD = 0;
        clusters.all = [];
        clusters.dHPC = [];
        clusters.int.dHPC = [];
        for ii=1:size(group_dHPC,1)
            cluster = group_dHPC(ii,1);
            Cellulartype = logical(cellulartype(cellulartype(:,1) == cluster,2));
            if Cellulartype
                a = length(Restrict(spks(spks(:,1)==cluster,2),aversiveTS_run./1000)) / ((aversiveTS_run(2)-aversiveTS_run(1))/1000);
                r = length(Restrict(spks(spks(:,1)==cluster,2),rewardTS_run./1000)) / ((rewardTS_run(2)-rewardTS_run(1))/1000);
                if or(a > 0 ,  r > 0)
                    numberD = numberD+1;
                    clusters.all = [clusters.all ; cluster];
                    clusters.dHPC = [clusters.dHPC ; cluster];
                end
            else
                clusters.int.dHPC = [clusters.int.dHPC ; cluster];
            end
            clear tmp cluster Cellulartype fr1 fr2 fr3 fr4 fr5 r a
        end
        
        numberV = 0;
        clusters.vHPC = [];
        clusters.int.vHPC = [];
        for ii=1:size(group_vHPC,1)
            cluster = group_vHPC(ii,1);
            Cellulartype = logical(cellulartype(cellulartype(:,1) == cluster,2));
            if Cellulartype
                a = length(Restrict(spks(spks(:,1)==cluster,2),aversiveTS_run./1000)) / ((aversiveTS_run(2)-aversiveTS_run(1))/1000);
                r = length(Restrict(spks(spks(:,1)==cluster,2),rewardTS_run./1000)) / ((rewardTS_run(2)-rewardTS_run(1))/1000);
                if or(a > 0 ,  r > 0)
                    numberV = numberV+1;
                    clusters.all = [clusters.all ; cluster];
                    clusters.vHPC = [clusters.vHPC ; cluster];
                end
            else
                clusters.int.vHPC = [clusters.int.vHPC ; cluster];
            end
            clear tmp cluster Cellulartype fr1 fr2 fr3 fr4 fr5 r a
        end
        clear freq limits
        clear camara shock rightvalve leftvalve
        clear ejeX ejeY dX dY dX_int dY_int
        
        if or(numberD >3 , numberV > 3)
            % --- Aversive ---
            disp('Lets go for the assemblies')
            if isfile('dorsalventral_assemblies_aversiveVF.mat')
                disp('Loading Aversive template')
                load('dorsalventral_assemblies_aversiveVF.mat')
            end
            
            Thresholded.aversive.all = Th;
            patterns.all.aversive = pat;
            clear cond Th pat
            
            % Detection of members
            if not(isempty(Thresholded.aversive.all))
                if numberD>0
                    cond1 =  sum(Thresholded.aversive.all(1:size(clusters.dHPC,1),:),1)>0; %checking of dHPC SU
                    cond2 =  sum(Thresholded.aversive.all(size(clusters.dHPC,1)+1:end,:),1)>0; %checking of vHPC SU
                    cond.dHPC.aversive = and(cond1 , not(cond2));
                    cond.vHPC.aversive = and(cond2 , not(cond1));
                    cond.both.aversive = and(cond1 , cond2); clear cond1 cond2
                else
                    cond1 =  logical(zeros(1,size(Thresholded.aversive.all,2))); %checking of dHPC SU
                    cond2 =  logical(ones(1,size(Thresholded.aversive.all,2))); %checking of vHPC SU
                    cond.dHPC.aversive = and(cond1 , not(cond2));
                    cond.vHPC.aversive = and(cond2 , not(cond1));
                    cond.both.aversive = and(cond1 , cond2); clear cond1 cond2
                end
            else
                cond1 =  false; %checking of dHPC SU
                cond2 =  logical(0); %checking of vHPC SU
                cond.dHPC.aversive = and(cond1 , not(cond2));
                cond.vHPC.aversive = and(cond2 , not(cond1));
                cond.both.aversive = and(cond1 , cond2); clear cond1 cond2
            end
            
            % --- Reward ---
            disp('Loading Reward template')
            if isfile('dorsalventral_assemblies_rewardVF.mat')
                load('dorsalventral_assemblies_rewardVF.mat')
            end
            
            Thresholded.reward.all = Th;
            patterns.all.reward = pat;
            clear Th pat
            
            % Detection of members using
            if not(isempty(Thresholded.reward.all))
                if numberD>0
                    cond1 =  sum(Thresholded.reward.all(1:size(clusters.dHPC,1),:),1)>0; %checking of dHPC SU
                    cond2 =  sum(Thresholded.reward.all(size(clusters.dHPC,1)+1:end,:),1)>0; %checking of vHPC SU
                    cond.dHPC.reward = and(cond1 , not(cond2));
                    cond.vHPC.reward = and(cond2 , not(cond1));
                    cond.both.reward = and(cond1 , cond2); clear cond1 cond2
                else
                    cond1 =  logical(zeros(1,size(Thresholded.reward.all,2))); %checking of dHPC SU
                    cond2 =  logical(ones(1,size(Thresholded.reward.all,2))); %checking of vHPC SU
                    cond.dHPC.reward = and(cond1 , not(cond2));
                    cond.vHPC.reward = and(cond2 , not(cond1));
                    cond.both.reward = and(cond1 , cond2); clear cond1 cond2
                end
            else
                cond1 =  logical(0); %checking of dHPC SU
                cond2 =  logical(0); %checking of vHPC SU
                cond.dHPC.reward = and(cond1 , not(cond2));
                cond.vHPC.reward = and(cond2 , not(cond1));
                cond.both.reward = and(cond1 , cond2); clear cond1 cond2
            end
            
            %% SpikeTrains construction
            limits = [0 segments.Var1(end)/1000];
            events = [];
            [Spikes , bins , Clusters] = spike_train_construction([spks_dHPC;spks_vHPC], clusters.all, cellulartype, 0.025, limits, events, false, true);
            clear limits events
            
            if sum(cond.both.aversive)>=1
                if aversiveTS(1)< rewardTS(1)
                    timestamps.PreSleep = Restrict(NREM.baseline,[NREM.baseline(end,2)-tempo NREM.baseline(end,2)]);
                    timestamps.PostSleep = Restrict(NREM.aversive,[NREM.aversive(1,1) NREM.aversive(1,1)+tempo]);
                else
                    timestamps.PreSleep = Restrict(NREM.reward,[NREM.reward(end,2)-tempo NREM.reward(end,2)]);
                    timestamps.PostSleep = Restrict(NREM.aversive,[NREM.aversive(1,1) NREM.aversive(1,1)+tempo]);                  
                end
                p = patterns.all.aversive;
                [P R] = reactivation_strength_members(p , cond.both.aversive , [bins' , Spikes] , th , Thresholded.aversive.all , [numberD numberV] , timestamps); clear templates
                Reactivation.aversive = [Reactivation.aversive ; R]; 
                Cumulative.aversive = [Cumulative.aversive , P]; clear R p P
            end
            
            if sum(cond.both.reward)>=1
                if aversiveTS(1)< rewardTS(1)
                    timestamps.PreSleep = Restrict(NREM.aversive,[NREM.aversive(end,2)-tempo NREM.aversive(end,2)]);
                    timestamps.PostSleep = Restrict(NREM.reward,[NREM.reward(1,1) NREM.reward(1,1)+tempo]);                
                else
                    timestamps.PreSleep = Restrict(NREM.baseline,[NREM.baseline(end,2)-tempo NREM.baseline(end,2)]);
                    timestamps.PostSleep = Restrict(NREM.reward,[NREM.reward(1,1) NREM.reward(1,1)+tempo]);                    
                end
                p = patterns.all.aversive;
                [P R] = reactivation_strength_members(p , cond.both.aversive , [bins' , Spikes] , th , Thresholded.aversive.all , [numberD numberV] , timestamps); clear templates
                Reactivation.reward = [Reactivation.reward ; R];
                Cumulative.reward = [Cumulative.reward , P]; clear R p P
            end
            
            clear spks NREM REM WAKE group_dHPC group_vHPC K Kinfo
            clear rewardTS rewardTS_run aversiveTS aversiveTS_run baselineTS
            clear numberD numberV patterns Thresholded cond Shocks_filt Spikes
            clear segments spks_dHPC spks_vHPC timestamps behavior Cell_type_classification
            clear config clusters movement Rewards_filt 
            disp(' ')
        end
    end
end

%% Plor reactivation
figure
x = Reactivation.reward;
y = Reactivation.aversive;
    
kstest(x)
kstest(y)
[h, p] = ranksum(x,y)
[h, p] = signrank(y,1)
[h, p] = signrank(x,1)
    
grps = [ones(size(x,1),1) ; ones(size(y,1),1)*2];
Y = [x;y];
scatter(grps,Y,"filled",'jitter','on', 'jitterAmount',0.1),hold on
scatter([1 2] , [nanmean(x) nanmean(y)],'filled'),xlim([0 3])%,ylim([-0.5 0.5])

%% plot FR in function of time
figure
c = [];
coefficients.aversive = [];
RA = [];
for i = 1 : size(Cumulative.aversive,2)
    t = Cumulative.aversive(:,i);
%     s1 = corrcoef(t',[1:1:size(t,1)]);
    f = fitlm(t',[1:1:size(t,1)],'linear');
    coefficients.aversive = [coefficients.aversive, f.Coefficients{:, 'Estimate'}];
    c = [c , (Smooth(Cumulative.aversive(:,i),1))]; clear f 
%     RA = [RA ; nanmean(cumulative.aversive.dvHPC(1:5,i))  nanmean(cumulative.aversive.dvHPC(11:20,i))  nanmean(cumulative.aversive.dvHPC(21:30,i))  nanmean(cumulative.aversive.dvHPC(31:40,i))];
end
plot(nanmean(c,2),'r'),hold on
ciplot(nanmean(c',1)-nansem(c'),nanmean(c',1)+nansem(c'),[1:45],'r'),alpha 0.5

c = [];
coefficients.reward = [];
RR = [];
for i = 1 : size(Cumulative.reward,2)
    t = Cumulative.reward(:,i);
%     s1 = corrcoef(t',[1:1:size(t,1)]);
    f = fitlm(t',[1:1:size(t,1)],'linear');
    coefficients.reward = [coefficients.reward, f.Coefficients{:, 'Estimate'}];
    c = [c , (Smooth(Cumulative.reward(:,i),1))]; clear f 
%     RR = [RR ; nanmean(Cumulative.reward(1:5,i))  nanmean(Cumulative.reward(11:20,i))  nanmean(Cumulative.reward(21:30,i))  nanmean(Cumulative.reward(31:40,i))];
end
plot(nanmean(c,2),'b'),hold on
ciplot(nanmean(c',1)-nansem(c'),nanmean(c',1)+nansem(c'),[1:45],'b'),alpha 0.5


end