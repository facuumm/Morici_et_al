function [C] = assemblies_crossCorrelation(path)
% This function calculates a CCG between dorsal and ventral assemblies
% in Pre, Post sleep, and Run. Also, it defines which assemblies are
% correlated (mean activation > 99percentile from a surrogate)
%
% INPUTS
% path: cell, inside each cell you should put the path of each session you
%       want to analyze
%
% OUTPUT
% C: Structure, it store the CCG Assemblies between assemblies.
%
% Morci Juan Facundo 04/2024

%variables to use in the script
criteria_fr = 0;

%variables for CCG construction
th = 3;           % threshold for peak detection
b = 0.1;        % binsize

% storage variables
C.aversive.pre = [];  C.reward.pre = [];
C.aversive.post = [];  C.reward.post = [];

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
                else
                    TS.pre = Restrict(NREM.reward,[NREM.reward(end,2)-1800 NREM.reward(end,2)]);
                    TS.post = Restrict(NREM.aversive,[NREM.aversive(1,1) NREM.aversive(1,1)+1800]);
                end
                
                Activity.vHPC = assembly_activity(patterns.aversive.vHPC , SpksTrains.vHPC');
                Activity.dHPC = assembly_activity(patterns.aversive.dHPC , SpksTrains.dHPC');
                
                TS.post = InIntervals(Bins,TS.post);
                TS.pre = InIntervals(Bins,TS.pre);
                for i = 1 : size( Activity.dHPC,1)
                    tmp =  Activity.dHPC(i,TS.post);
                    for ii = 1 : size(Activity.vHPC,1)
                        tmp1 =  Activity.vHPC(ii,TS.post);
                        c = xcorr(tmp,tmp1,500/5,'normalized'); clear tmp1
                        C.aversive.post = [C.aversive.post ; c]; clear c
                    end
                    
                    tmp =  Activity.dHPC(i,TS.pre);
                    for ii = 1 : size(Activity.vHPC,1)
                        tmp1 =  Activity.vHPC(ii,TS.pre);
                        c = xcorr(tmp,tmp1,500/5,'normalized'); clear tmp1
                        C.aversive.pre = [C.aversive.pre ; c]; clear c
                    end
                end
                clear TS Activity
            end

            
            if isfield(patterns,'reward') %checking if there are reward assemblies
                if aversiveTS_run(1)<rewardTS_run(1)
                    TS.pre = Restrict(NREM.aversive,[NREM.aversive(end,2)-1800 NREM.aversive(end,2)]);
                    TS.post = Restrict(NREM.reward,[NREM.reward(1,1) NREM.reward(1,1)+1800]);
                else
                    TS.pre = Restrict(NREM.baseline,[NREM.baseline(end,2)-1800 NREM.baseline(end,2)]);
                    TS.post = Restrict(NREM.reward,[NREM.reward(1,1) NREM.reward(1,1)+1800]);
                    
                end
                
                Activity.vHPC = assembly_activity(patterns.reward.vHPC , SpksTrains.vHPC');
                Activity.dHPC = assembly_activity(patterns.reward.dHPC , SpksTrains.dHPC');
                
                TS.post = InIntervals(Bins,TS.post);
                TS.pre = InIntervals(Bins,TS.pre);
                for i = 1 : size( Activity.dHPC,1)
                    tmp =  Activity.dHPC(i,TS.post);
                    for ii = 1 : size(Activity.vHPC,1)
                        tmp1 =  Activity.vHPC(ii,TS.post);
                        c = xcorr(tmp,tmp1,500/5,'normalized'); clear tmp1
                        C.reward.post = [C.reward.post ; c]; clear c
                    end
                    
                    tmp =  Activity.dHPC(i,TS.pre);
                    for ii = 1 : size(Activity.vHPC,1)
                        tmp1 =  Activity.vHPC(ii,TS.pre);
                        c = xcorr(tmp,tmp1,500/5,'normalized'); clear tmp1
                        C.reward.pre = [C.reward.pre ; c]; clear c
                    end
                end
                clear TS Activity
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

% T = [-200 : 5 : 200];
% 
% Pre = [];
% for i = 1:size(C.aversive.pre,1)
%     Pre = [Pre , Smooth(C.aversive.pre(i,:),2)];
% end
% Post = [];
% for i = 1:size(C.aversive.post,1)
%     Post = [Post , Smooth(C.aversive.post(i,:),2)];
% end
% 
% [~ , up] = min(abs(T-100));
% [~ , down] = min(abs(T-(-100)));
% [~ , iterator] = sort(nanmean(Post(down:up,:)))
% 
% % [~ , up] = max(Pre);
% % [~ , iterator] = sort(up)
% 
% figure,
% subplot(121),imagesc(T,[1:1:size(Pre,2)],Pre(:,iterator)'),caxis([-0.1 0.1]),colormap 'jet'
% subplot(122),imagesc(T,[1:1:size(Post,2)],Post(:,iterator)'),caxis([-0.1 0.1]),colormap 'jet'
% 
% 
% 
% Pre = [];
% for i = 1:size(C.reward.pre,1)
%     Pre = [Pre , Smooth(C.reward.pre(i,:),1)];
% end
% Post = [];
% for i = 1:size(C.reward.post,1)
%     Post = [Post , Smooth(C.reward.post(i,:),1)];
% end
% 
% [~ , up] = min(abs(T-50));
% [~ , down] = min(abs(T-(-50)));
% [~ , iterator] = sort(nanmean(Post(down:up,:)))
% 
% figure,
% 
% subplot(121),imagesc(T,[1:1:size(Pre,2)],Pre(:,iterator)'),caxis([-0.01 0.01])
% subplot(122),imagesc(T,[1:1:size(Post,2)],Post(:,iterator)'),caxis([-0.01 0.01])
% 
% 
% 
% 
% figure
% plot(T,nanmean(Pre'),'k'), hold on
% plot(T,nanmean(Post'),'r')
% 
% figure
% plot(nanmean(C.reward.pre),'k'), hold on
% plot(nanmean(C.reward.post),'b')


end
