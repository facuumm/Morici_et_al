function [Pre Post Run] = Reactivation_Peaks(path)
% This function calculates the Pre and Post sleep assemblies rate.
%
% INPUTS
% path: cell, inside each cell you should put the path of each session you
%       want to analyze
%
% OUTPUT
% Pre and Post: Structure, it store the Assemblies Rate for each type of
%               assemblies.
%
%               Architecture of each output:
%                   Pre.dHPC
%                      .vHPC.aversive
%                           .reward
%                   Post.dHPC
%                       .vHPC.aversive
%                            .reward
%
% Morci Juan Facundo 04/2024

%variables to use in the script
criteria_fr = 0;
th = 5; % threshold for detecting peak assemblies
% 3 SD funciona

% storage variables
Pre.dHPC.aversive = [];    Pre.vHPC.aversive = [];
Post.dHPC.aversive = [];   Post.vHPC.aversive = [];
Pre.dHPC.reward = [];      Pre.vHPC.reward = [];
Post.dHPC.reward = [];     Post.vHPC.reward = [];
Run.dHPC.reward = [];      Run.dHPC.aversive = [];
Run.vHPC.reward = [];      Run.vHPC.aversive = [];

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
        disp(['-- Initiating analysis of folder #' , num2str(t) , ' from rat #',num2str(tt) , ' --'])
        session = [subFolders(t+2).folder,'\',subFolders(t+2).name];
        cd(session)
        
        %Loading TS of the sessions
        disp('Uploading session time stamps')
        load('session_organization.mat')
        load('behavioral_data.mat')
        
        %% load sleep states
        disp('Uploading sleep scoring')
        x = dir([cd,'\*-states.mat']);    states = load([cd,'\',x.name]);    states = states.states;
        REM.all = ToIntervals(states==5);    NREM.all = ToIntervals(states==3);
        clear x states
        
        NREM.baseline = Restrict(NREM.all,baselineTS./1000);
        NREM.aversive = Restrict(NREM.all,aversiveTS./1000);
        NREM.reward = Restrict(NREM.all,rewardTS./1000);
        
        REM.baseline = Restrict(REM.all,baselineTS./1000);
        REM.aversive = Restrict(REM.all,aversiveTS./1000);
        REM.reward = Restrict(REM.all,rewardTS./1000);

        
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
        K = [K , Cell_type_classification(:,6:8)];
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
        % pyr
        cellulartype = [K(:,1) , K(:,4)];
        
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
        
        %% Assemblies detection
        if or(numberD > 3 , numberV > 3)
            % --- Aversive ---
            disp('Lets go for the assemblies')
            if isfile('separated_assemblies.mat')
                load('separated_assemblies.mat')
            end
            
            % Definition of limits for assemblie activity strength
            limits = [0 segments.Var1(end)/1000];
            events = [];
            
            % Spiketrains construction
            if numberD>0 %for dHPC SUs
                [SpksTrains.dHPC , Bins , Cluster] = spike_train_construction(spks_dHPC, clusters.dHPC, cellulartype, 0.005, limits, events, false,true);
            end
            if numberV>0 %for vHPC SUs
                [SpksTrains.vHPC , Bins , Cluster] = spike_train_construction(spks_vHPC, clusters.vHPC, cellulartype, 0.005, limits, events, false,true);
            end
            
            if isfield(patterns,'aversive') %checking if there are aversive assemblies
                if aversiveTS_run(1)<rewardTS_run(1)
%                     TS.pre = REM.baseline;
%                     TS.post = REM.aversive;
                    TS.pre = Restrict(NREM.baseline,[NREM.baseline(end,2)-1800 NREM.baseline(end,2)]);
                    TS.post = Restrict(NREM.aversive,[NREM.aversive(1,1) NREM.aversive(1,1)+1800]);
                    TS.run.run1 = movement.aversive;
                    TS.run.run2 = movement.reward;
                else
%                     TS.pre = REM.reward;
%                     TS.post = REM.aversive;
                    TS.pre = Restrict(NREM.reward,[NREM.reward(end,2)-1800 NREM.reward(end,2)]);
                    TS.post = Restrict(NREM.aversive,[NREM.aversive(1,1) NREM.aversive(1,1)+1800]);   
                    TS.run.run1 = movement.aversive;
                    TS.run.run2 = movement.reward;
                end
                
                if numberV>0
                    pks.vHPC.aversive = assemblies_peaks([Bins' SpksTrains.vHPC] , patterns.aversive.vHPC , th);
                    [pre1 , post1] = Reactivation_Rate(pks.vHPC.aversive,TS);
                    [pre2 , post2] = Reactivation_Mean([Bins' SpksTrains.vHPC] , patterns.aversive.vHPC , TS);
                    
                    [run1 run2] = Activation_Run_Rate(pks.vHPC.aversive,TS.run);
                    [run3 run4] = Activation_Run_Mean([Bins' SpksTrains.vHPC] , patterns.aversive.vHPC , TS.run);                    
                    
                    Pre.vHPC.aversive = [Pre.vHPC.aversive ; pre1  pre2]; clear pre1 pre2
                    Post.vHPC.aversive = [Post.vHPC.aversive ; post1 post2]; clear post1 post2
                    Run.vHPC.aversive = [Run.vHPC.aversive ; run1 run2 run3 run4]; clear run1 run2 run3 run4
                end
                
                if numberD>0
                    pks.dHPC.aversive = assemblies_peaks([Bins' SpksTrains.dHPC] , patterns.aversive.dHPC , th);
                    [pre1 , post1] = Reactivation_Rate(pks.dHPC.aversive,TS);
                    [pre2 , post2] = Reactivation_Mean([Bins' SpksTrains.dHPC] , patterns.aversive.dHPC , TS);
                    
                    [run1 run2] = Activation_Run_Rate(pks.dHPC.aversive,TS.run);
                    [run3 run4] = Activation_Run_Mean([Bins' SpksTrains.dHPC] , patterns.aversive.dHPC , TS.run); 
                    
                    Pre.dHPC.aversive = [Pre.dHPC.aversive ; pre1 pre2]; clear pre1 pre2
                    Post.dHPC.aversive = [Post.dHPC.aversive ; post1 post2]; clear post1 post2
                    Run.dHPC.aversive = [Run.dHPC.aversive ; run1 run2 run3 run4]; clear run1 run2 run3 run4
                end
            end
            
            
            if isfield(patterns,'reward') %checking if there are reward assemblies
                if aversiveTS_run(1)<rewardTS_run(1)
%                     TS.pre = REM.aversive;
%                     TS.post =REM.reward;
                    TS.pre = Restrict(NREM.aversive,[NREM.aversive(end,2)-1800 NREM.aversive(end,2)]);
                    TS.post = Restrict(NREM.reward,[NREM.reward(1,1) NREM.reward(1,1)+1800]);
                    TS.run.run1 = movement.reward;
                    TS.run.run2 = movement.aversive;
                else
%                     TS.pre = REM.baseline;
%                     TS.post = REM.reward;   
                    TS.pre = Restrict(NREM.baseline,[NREM.baseline(end,2)-1800 NREM.baseline(end,2)]);
                    TS.post = Restrict(NREM.reward,[NREM.reward(1,1) NREM.reward(1,1)+1800]);   
                    TS.run.run1 = movement.reward;
                    TS.run.run2 = movement.aversive;                  
                end
                
                if numberV>0
                    pks.vHPC.reward = assemblies_peaks([Bins' SpksTrains.vHPC] , patterns.reward.vHPC , th);
                    [pre1 , post1] = Reactivation_Rate(pks.vHPC.reward,TS);
                    [pre2 , post2] = Reactivation_Mean([Bins' SpksTrains.vHPC] , patterns.reward.vHPC , TS);

                    [run1 run2] = Activation_Run_Rate(pks.vHPC.reward,TS.run);
                    [run3 run4] = Activation_Run_Mean([Bins' SpksTrains.vHPC] , patterns.reward.vHPC , TS.run);  
                    
                    Pre.vHPC.reward = [Pre.vHPC.reward ; pre1 pre2]; clear pre1 pre2
                    Post.vHPC.reward = [Post.vHPC.reward ; post1 post2]; clear post1 post2
                    Run.vHPC.reward = [Run.vHPC.reward ; run1 run2 run3 run4]; clear run1 run2 run3 run4
                end
                
                if numberD>0
                    pks.dHPC.reward = assemblies_peaks([Bins' SpksTrains.dHPC] , patterns.reward.dHPC , th);
                    [pre1 , post1] = Reactivation_Rate(pks.dHPC.reward,TS);
                    [pre2 , post2] = Reactivation_Mean([Bins' SpksTrains.dHPC] , patterns.reward.dHPC , TS);
                    
                    [run1 run2] = Activation_Run_Rate(pks.dHPC.reward,TS.run);
                    [run3 run4] = Activation_Run_Mean([Bins' SpksTrains.dHPC] , patterns.reward.dHPC , TS.run);  
                    
                    Pre.dHPC.reward = [Pre.dHPC.reward ; pre1 pre2]; clear pre1 pre2
                    Post.dHPC.reward = [Post.dHPC.reward ; post1 post2]; clear post1 post2
                    Run.dHPC.reward = [Run.dHPC.reward ; run1 run2 run3 run4]; clear run1 run2 run3 run4                
                end
            end            
            
        end
        disp(' ')
        clear A aversiveTS aversiveTS_run baselineTS rewardTS rewardTS_run
        clear behavior bins Cell_type_classification cellulartype cond
        clear is K Kinfo group_dHPC group_vHPC matrixC matrixD matrixV
        clear NREM REM WAKE segmentation tmp cond config
        clear spiketrains_dHPC spiketrains_vHPC opts MUA
        clear patterns Thresholded i  ii numberD numberV movement cross crossN
        clear Spikes bins Clusters Shocks_filt Rewards_filt config n_SU_D n_SU_V
        clear clusters coordinated coordinated_ripple_bursts coordinatedV
        clear cooridnated_event coordinatedV_refined coordinatedV_refined
        clear ripple_bursts ripple_event ripplesD ripplesV
        clear spks spks_dHPC spks_vHPC ripples cooridnated_event
        clear cooridnated_eventDV cooridnated_eventVD segments movement
    end
    clear num_assembliesA num_assembliesR
end


figure,
subplot(121)
x = [ones(size(Pre.dHPC.aversive(:,1))) ; ones(size(Post.dHPC.reward(:,1)))*2];
DIA = (Post.dHPC.aversive(:,1) ./ Pre.dHPC.aversive(:,1));% ./ (Post.dHPC.aversive(:,2) + Pre.dHPC.aversive(:,2));
DIR = (Post.dHPC.reward(:,1) ./ Pre.dHPC.reward(:,1));% ./ (Post.dHPC.reward(:,2) + Pre.dHPC.reward(:,2));
y = [DIA ; DIR];
scatter(x,y,"filled",'jitter','on', 'jitterAmount',0.1'),hold on,xlim([0 3])
scatter([1 2] , [nanmedian(DIA) nanmedian(DIR)],'filled'),set(gca, 'YScale', 'log')
% boxplot(y,x)
yline(1,'--')
signrank(DIA,0)
signrank(DIR,0)
ranksum(DIA,DIR)
%  ylim([-1 1])

subplot(122)
x = [ones(size(Pre.vHPC.aversive(:,1))) ; ones(size(Post.vHPC.reward(:,1)))*2];
DIA1 = (Post.vHPC.aversive(:,1) ./ Pre.vHPC.aversive(:,1));% ./ (Post.vHPC.aversive(:,2) + Pre.vHPC.aversive(:,2));
DIR1 = (Post.vHPC.reward(:,1) ./ Pre.vHPC.reward(:,1));% ./ (Post.vHPC.reward(:,2) + Pre.vHPC.reward(:,2));
y = [DIA1 ; DIR1];
scatter(x,y,"filled",'jitter','on', 'jitterAmount',0.1'),hold on,xlim([0 3])
scatter([1 2] , [nanmedian(DIA) nanmedian(DIR)],'filled'),set(gca, 'YScale', 'log')
% boxplot(y,x)
yline(1,'--')
signrank(DIA1,0)
signrank(DIR1,0)
ranksum(DIA1,DIR1)
% ylim([-1 1])



figure
subplot(221)
fitlm(Run.dHPC.aversive(:,1),Post.dHPC.aversive(:,1))
plot(ans),ylim([0 2]),xlim([0 3])

subplot(222)
fitlm(Run.dHPC.reward(:,1),Post.dHPC.reward(:,1))
plot(ans),ylim([0 2]),xlim([0 3])

[p,t,df,slope1,slope2,sem1,sem2] = CompareSlopes(Post.dHPC.aversive(:,1),Run.dHPC.aversive(:,1),Post.dHPC.reward(:,1),Run.dHPC.reward(:,1))

subplot(223)
fitlm(Run.vHPC.aversive(:,1),Post.vHPC.aversive(:,1))
plot(ans),ylim([0 2]),xlim([0 3])

subplot(224)
fitlm(Run.vHPC.reward(:,1),Post.vHPC.reward(:,1))
plot(ans),ylim([0 2]),xlim([0 3])

[p,t,df,slope1,slope2,sem1,sem2] = CompareSlopes(Post.vHPC.aversive(:,1),Run.vHPC.aversive(:,1),Post.vHPC.reward(:,1),Run.vHPC.reward(:,1))


end


