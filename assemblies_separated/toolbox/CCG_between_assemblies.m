function [Pre , Post , Run , Iterator , T] = CCG_between_assemblies(path)
% This function calculates a CCG between dorsal and ventral assemblies
% in Pre, Post sleep, and Run. Also, it defines which assemblies are
% correlated (mean activation > 99percentile from a surrogate)
%
% INPUTS
% path: cell, inside each cell you should put the path of each session you
%       want to analyze
%
% OUTPUT
% Pre and Post: Structure, it store the CCG Assemblies between assemblies.
%
% Run: Structure, it store the CCG between assemblies during run.
%
% Iterator: Structure, it contains a logical array (1: correlated / 0: not)
%
% T: time vector.
%
%               Architecture of each output:
%                   Pre.aversive / reward
%                   Post.aversive / reward
%                   Run.aversive / reward
%                   Iterator.aversive / reward
%
% Morci Juan Facundo 04/2024

%variables to use in the script
criteria_fr = 0;

%variables for CCG construction
th = 5;           % threshold for peak detection
sm = 1;           % smooth
dur = 0.5;        % window duration
b = 0.005;        % binsize
iterations = 500; % number of iterations for Run surrogate

% storage variables
Pre.aversive = [];       Post.aversive = [];
Pre.reward = [];         Post.reward = [];
Run.aversive = [] ;      Run.reward = [];
Iterator.aversive = [];  Iterator.reward = [];


T = [];

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
        REM.all = ToIntervals(states==5);    NREM.all = ToIntervals(states==3);    WAKE.all = ToIntervals(states==1);
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
        if and(numberD > 3 , numberV > 3)
            % --- Aversive ---
            disp('Lets go for the assemblies')
            if isfile('separated_assemblies.mat')
                load('separated_assemblies.mat')
            end
            
            % Definition of limits for assemblie activity strength
            limits = [0 segments.Var1(end)/1000];
            events = [];
            
            % Spiketrains construction
            [SpksTrains.dHPC , Bins , Cluster] = spike_train_construction(spks_dHPC, clusters.dHPC, cellulartype, 0.005, limits, events, false,true);
            [SpksTrains.vHPC , Bins , Cluster] = spike_train_construction(spks_vHPC, clusters.vHPC, cellulartype, 0.005, limits, events, false,true);
            
            if isfield(patterns,'aversive') %checking if there are aversive assemblies
                if aversiveTS_run(1)<rewardTS_run(1)
                    TS.pre = Restrict(NREM.baseline,[NREM.baseline(end,2)-1800 NREM.baseline(end,2)]);
                    TS.post = Restrict(NREM.aversive,[NREM.aversive(1,1) NREM.aversive(1,1)+1800]);
%                     TS.pre = NREM.baseline;
%                     TS.post = NREM.aversive;
                else
                    TS.pre = Restrict(NREM.reward,[NREM.reward(end,2)-1800 NREM.reward(end,2)]);
                    TS.post = Restrict(NREM.aversive,[NREM.aversive(1,1) NREM.aversive(1,1)+1800]);
%                     TS.pre = NREM.reward;
%                     TS.post = NREM.aversive;
                end
                
                % Peaks detection
                pks.vHPC.aversive = assemblies_peaks([Bins' SpksTrains.vHPC] , patterns.aversive.vHPC , th);
                pks.dHPC.aversive = assemblies_peaks([Bins' SpksTrains.dHPC] , patterns.aversive.dHPC , th);
                
                % Run cooridnation
                for i = 1 : size(patterns.aversive.dHPC,2)
                    x = Restrict(pks.dHPC.aversive{i}(:,1) , movement.aversive);
                    for ii = 1 : size(patterns.aversive.vHPC,2)
                        y = Restrict(pks.vHPC.aversive{ii}(:,1) ,  movement.aversive);
                        if and(size(x,1)>5,size(y,1)>5)
                            [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                            [ccg,T] = CCG(s,ids,'binSize',b,'duration',dur,'smooth',sm,'mode','ccg');
                            Run.aversive = [Run.aversive , (ccg(:,1,2)./b)];
                            
                            % timepoints to define the mean value
                            [~ , up] = min(abs(T-0.1));
                            [~ , down] = min(abs(T-(-0.1)));
                            
                            % calculation of mean value
                            M = nanmean(ccg(down:up,1,2)./b);
                           
                            % surrogate construction
                            surrogate = [];
                            for i = 1 : iterations
                                [s,ids,groups] = CCGParameters(x,ones(length(x),1),ShuffleSpks(y),ones(length(y),1)*2);
                                [ccg,T] = CCG(s,ids,'binSize',b,'duration',dur,'smooth',sm,'mode','ccg');
                                surrogate = [surrogate ; nanmean(ccg(down:up,1,2)./b)];
                            end
                            
                            % Storage of result
                            if M > prctile(surrogate , 75)
                                Iterator.aversive = [Iterator.aversive ; 1];
                            else
                                Iterator.aversive = [Iterator.aversive ; 0];
                            end
                            clear y s ids groups ccg surrogate M up down
                        else
                            Iterator.aversive = [Iterator.aversive ; nan];
                             Run.aversive = [Run.aversive , NaN(size(Run.aversive , 1),1)];
                        end
                    end
                    clear x
                end
                
                    
                % Sleep coordination
                for i = 1 : size(patterns.aversive.dHPC,2)
                    x = Restrict(pks.dHPC.aversive{i}(:,1) , TS.pre);
                    for ii = 1 : size(patterns.aversive.vHPC,2)
                        y = Restrict(pks.vHPC.aversive{ii}(:,1) , TS.pre);
                        if and(size(x,1)>5,size(y,1)>5)
                            [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                            [ccg,T] = CCG(s,ids,'binSize',b,'duration',dur,'smooth',sm,'mode','ccg');
                            Pre.aversive = [Pre.aversive , (ccg(:,1,2)./b)]; clear y s ids groups ccg
                        else
                            Pre.aversive = [Pre.aversive , NaN(size(Pre.aversive,1),1)];
                        end
                    end
                    clear x
                end
                
                for i = 1 : size(patterns.aversive.dHPC,2)
                    x = Restrict(pks.dHPC.aversive{i}(:,1) , TS.post);
                    for ii = 1 : size(patterns.aversive.vHPC,2)
                        y = Restrict(pks.vHPC.aversive{ii}(:,1) , TS.post);
                        if and(size(x,1)>5,size(y,1)>5)
                            [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                            [ccg,T] = CCG(s,ids,'binSize',b,'duration',dur,'smooth',sm,'mode','ccg');
                            Post.aversive = [Post.aversive , (ccg(:,1,2)./b)]; clear y s ids groups ccg
                        else
                            Post.aversive = [Post.aversive , NaN(size(Post.aversive,1),1)];
                        end
                    end
                    clear x
                end
            end
            
            if isfield(patterns,'reward') %checking if there are reward assemblies
                if aversiveTS_run(1)<rewardTS_run(1)
                    TS.pre = Restrict(NREM.aversive,[NREM.aversive(end,2)-1800 NREM.aversive(end,2)]);
                    TS.post = Restrict(NREM.reward,[NREM.reward(1,1) NREM.reward(1,1)+1800]);
%                     TS.pre = NREM.aversive;
%                     TS.post = NREM.reward;
                else
                    TS.pre = Restrict(NREM.baseline,[NREM.baseline(end,2)-1800 NREM.baseline(end,2)]);
                    TS.post = Restrict(NREM.reward,[NREM.reward(1,1) NREM.reward(1,1)+1800]);
%                     TS.pre = NREM.baseline;
%                     TS.post = NREM.reward;
                    
                end
                
                % Peaks detection
                pks.vHPC.reward = assemblies_peaks([Bins' SpksTrains.vHPC] , patterns.reward.vHPC , th);
                pks.dHPC.reward = assemblies_peaks([Bins' SpksTrains.dHPC] , patterns.reward.dHPC , th);
                
                % Run cooridnation
                for i = 1 : size(patterns.reward.dHPC,2)
                    x = Restrict(pks.dHPC.reward{i}(:,1) , movement.reward);
                    for ii = 1 : size(patterns.reward.vHPC,2)
                        y = Restrict(pks.vHPC.reward{ii}(:,1) ,  movement.reward);
                        if and(size(x,1)>5,size(y,1)>5)
                            [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                            [ccg,T] = CCG(s,ids,'binSize',b,'duration',dur,'smooth',sm,'mode','ccg');
                            Run.reward = [Run.reward , (ccg(:,1,2)./b)];
                            
                            % timepoints to define the mean value
                            [~ , up] = min(abs(T-0.1));
                            [~ , down] = min(abs(T-(-0.1)));
                            
                            % calculation of mean value
                            M = nanmean(ccg(down:up,1,2)./b);
                           
                            % surrogate construction
                            surrogate = [];
                            for i = 1 : iterations
                                [s,ids,groups] = CCGParameters(x,ones(length(x),1),ShuffleSpks(y),ones(length(y),1)*2);
                                [ccg,T] = CCG(s,ids,'binSize',b,'duration',dur,'smooth',sm,'mode','ccg');
                                surrogate = [surrogate ; nanmean(ccg(down:up,1,2)./b)];
                            end
                            
                            % Storage of result
                            if M > prctile(surrogate , 75)
                                Iterator.reward = [Iterator.reward ; 1];
                            else
                                Iterator.reward = [Iterator.reward ; 0];
                            end
                            clear y s ids groups ccg surrogate M up down
                        else
                            Iterator.reward = [Iterator.aversive ; nan];
                            Run.reward = [Run.reward , NaN(size(Run.reward , 1),1)];
                        end
                    end
                    clear x
                end      
                
                % Sleep Coordination
                for i = 1 : size(patterns.reward.dHPC,2)
                    x = Restrict(pks.dHPC.reward{i}(:,1) , TS.pre);
                    for ii = 1 : size(patterns.reward.vHPC,2)
                        y = Restrict(pks.vHPC.reward{ii}(:,1) , TS.pre);
                        if and(size(x,1)>5,size(y,1)>5)
                            [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                            [ccg,T] = CCG(s,ids,'binSize',b,'duration',dur,'smooth',sm,'mode','ccg');
                            Pre.reward = [Pre.reward , (ccg(:,1,2)./b)]; clear y s ids groups ccg
                        else
                            Pre.reward = [Pre.reward ,NaN(size(Pre.reward,1),1)];
                        end
                    end
                    clear x
                end
                
                for i = 1 : size(patterns.reward.dHPC,2)
                    x = Restrict(pks.dHPC.reward{i}(:,1) , TS.post);
                    for ii = 1 : size(patterns.reward.vHPC,2)
                        y = Restrict(pks.vHPC.reward{ii}(:,1) , TS.post);
                        if and(size(x,1)>5,size(y,1)>5)
                            [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                            [ccg,T] = CCG(s,ids,'binSize',b,'duration',dur,'smooth',sm,'mode','ccg');
                            Post.reward = [Post.reward , (ccg(:,1,2)./b)]; clear y s ids groups ccg
                        else
                            Post.reward = [Post.reward , NaN(size(Post.reward,1),1)];
                        end
                    end
                    clear x
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
        clear ripple_bursts ripple_event ripplesD ripplesV Bins
        clear spks spks_dHPC spks_vHPC ripples cooridnated_event SpksTrains
        clear cooridnated_eventDV cooridnated_eventVD segments pks
    end
    clear num_assembliesA num_assembliesR
end


%  IA = Iterator.aversive;
%  IA(isnan(IA)) = 0;
%  IA = logical(IA);
%  
%  IR = Iterator.reward;
%  IR(isnan(IR)) = 0;
%  IR = logical(IR);
%  
% figure,
% subplot(221)
% x1 = zscore(Pre.aversive(:,IA))';
% [i ii] = max(x1');
% [i ii] = sort(ii);
% imagesc(T , [1:1:size(x1,2)] , x1(ii,:)),caxis([-3 3]),colormap 'jet'
% 
% subplot(222)
% x2 = zscore(Post.aversive(:,IA))';
% [i ii] = max(x2');
% [i ii] = sort(ii);
% imagesc(T , [1:1:size(x2,2)] , x2(ii,:)),caxis([-3 3]),colormap 'jet'
% 
% subplot(223)
% x3 = zscore(Pre.reward(:,IR))';
% [i ii] = max(x3');
% [i ii] = sort(ii);
% imagesc(T , [1:1:size(x3,2)] , x3(ii,:)),caxis([-3 3]),colormap 'jet'
% 
% 
% subplot(224)
% x4 = zscore(Post.reward(:,IR))';
% [i ii] = max(x4');
% [i ii] = sort(ii);
% imagesc(T , [1:1:size(x4,2)] , x4(ii,:)),caxis([-3 3]),colormap 'jet'
% 
% 
% figure
% plot(T,nanmean(x1),'k'),hold on
% ciplot(nanmean(x1)-nansem(x1) , nanmean(x1)+nansem(x1) , T , 'k'),alpha 0.5
% plot(T,nanmean(x2),'r'),hold on
% ciplot(nanmean(x2)-nansem(x2) , nanmean(x2)+nansem(x2) , T , 'r'),alpha 0.5
% 
% 
% figure
% plot(T,nanmean(x3),'k'),hold on
% ciplot(nanmean(x3)-nansem(x3) , nanmean(x3)+nansem(x3) , T , 'k'),alpha 0.5
% plot(T,nanmean(x4),'b'),hold on
% ciplot(nanmean(x4)-nansem(x4) , nanmean(x4)+nansem(x4) , T , 'b'),alpha 0.5
% 
% 
% figure,
% subplot(121)
% x = zscore(Run.aversive)';
% [i ii] = max(x');
% [i ii] = sort(ii);
% imagesc(T , [1:1:size(Pre.aversive,2)] , x(ii,:)),caxis([-3 3]),colormap 'jet'
% 
% subplot(122)
% x = zscore(Run.reward)';
% [i ii] = max(x');
% [i ii] = sort(ii);
% imagesc(T , [1:1:size(Post.aversive,2)] , x(ii,:)),caxis([-3 3]),colormap 'jet'
% 
% figure
% plot(T,nanmean(zscore(Run.reward(:,IR))')),hold on
% plot(T,nanmean(zscore(Run.aversive(:,IA))')),hold on

end
