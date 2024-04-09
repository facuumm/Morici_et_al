function [cross time] = cross_assemblies(path)
% This function a CCG between timestamps from assemblies peaks.
%
% INPUTS
% path: cell, inside each cell you should put the path of each session you
%       want to analyze
%
% OUTPUT
% cross: structure, it contains the CCG for aversive and reward
%                .dHPC_joint.Pre
%                           .Post
%                .vHPC_joint.Pre
%                           .Post
%                .dHPC_vHPC.Pre
%                          .Post
%
% time: column vector, it contains the time axis for the CCG plot.
%
% other functions: CCG from FMA toolbox
% Morci Juan Facundo 01/2024

% variables for CCG
sm = 1;
dur = 1;

% Initialization of structures that wil contain the outputs
cross.dHPC_joint.aversive.Pre = [];   cross.dHPC_joint.aversive.Post = [];
cross.vHPC_joint.aversive.Pre = [];   cross.vHPC_joint.aversive.Post = [];
cross.dHPC_vHPC.aversive.Pre = [];    cross.dHPC_vHPC.aversive.Post = [];

cross.dHPC_joint.reward.Pre = [];   cross.dHPC_joint.reward.Post = [];
cross.vHPC_joint.reward.Pre = [];   cross.vHPC_joint.reward.Post = [];
cross.dHPC_vHPC.reward.Pre = [];    cross.dHPC_vHPC.reward.Post = [];

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
        clear y A R
        
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
        
        if aversiveTS_run(1) < rewardTS_run(1)
            Is.aversive.pre = NREM.baseline;
            Is.aversive.post = NREM.aversive;
            Is.reward.pre = NREM.aversive;
            Is.reward.post = NREM.reward;
        else
            Is.aversive.pre = NREM.reward;
            Is.aversive.post = NREM.aversive;
            Is.reward.pre = NREM.baseline;
            Is.reward.post = NREM.reward;
        end
        
        
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
        
        if or(numberD >2 , numberV > 2)
            %% --- Aversive ---
            disp('Lets go for the assemblies')
            if isfile('dorsalventral_assemblies_aversive.mat')
                disp('Loading Aversive template')
                load('dorsalventral_assemblies_aversive.mat')
            else
                Th = [];
                pat = [];
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
                cond1 =  logical(0); %checking of dHPC SU
                cond2 =  logical(0); %checking of vHPC SU
                cond.dHPC.aversive = and(cond1 , not(cond2));
                cond.vHPC.aversive = and(cond2 , not(cond1));
                cond.both.aversive = and(cond1 , cond2); clear cond1 cond2
            end
            
            %% --- Reward ---
            disp('Loading Reward template')
            if isfile('dorsalventral_assemblies_reward.mat')
                load('dorsalventral_assemblies_reward.mat')
            else
                Th = [];
                pat = [];
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
            
            
            %% SpikeTrain
            limits = [0 segments.Var1(end)/1000];
            [Spikes , bins , Clusters] = spike_train_construction([spks_dHPC;spks_vHPC], clusters.all, cellulartype, 0.025, limits, [] , false, true);
            Activities.aversive = assembly_activity(patterns.all.aversive,Spikes');
            Activities.reward = assembly_activity(patterns.all.reward,Spikes');
            clear limits Spikes Clusters
            
            %% Aversive
            if and(sum(cond.both.aversive)>0 , sum(cond.dHPC.aversive)>0)
                templates1 = patterns.all.aversive(:,cond.both.aversive);
                templates2 = patterns.all.aversive(:,cond.dHPC.aversive);
                
                Act1 = (Activities.aversive(cond.both.aversive,:));
                Act2 = (Activities.aversive(cond.dHPC.aversive,:));
                
                for i = 1:size(templates1,2)
                    [pks1,loc1] = findpeaks(zscore(Act1(i,:)),bins,'MinPeakHeight',5);
                    for ii = 1:size(templates2,2)
                        [pks2,loc2] = findpeaks(zscore(Act2(ii,:)),bins,'MinPeakHeight',5);
                        % for Pre
                        x = Restrict(loc1,Is.aversive.pre);
                        y = Restrict(loc2,Is.aversive.pre);
                        if and(length(x)>5 , length(y)>5)
                            [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                            [ccg,T] = CCG(s,ids,'binSize',0.025,'duration',dur,'smooth',sm,'mode','ccg');
                            ccg = ccg(:,1,2)./sum(ccg(:,1,2));
                            cross.dHPC_joint.aversive.Pre = [cross.dHPC_joint.aversive.Pre , ccg];
                            time = T; clear ccg T x y s ids groups
                        end
                        
                        % for Post
                        x = Restrict(loc1,Is.aversive.post);
                        y = Restrict(loc2,Is.aversive.post);
                        if and(length(x)>5 , length(y)>5)
                            [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                            [ccg,T] = CCG(s,ids,'binSize',0.025,'duration',dur,'smooth',sm,'mode','ccg');
                            ccg = ccg(:,1,2)./sum(ccg(:,1,2));
                            cross.dHPC_joint.aversive.Post = [cross.dHPC_joint.aversive.Post , ccg];
                            time = T; clear ccg T x y s ids groups
                        end
                    end
                end
                clear  templates1 templates2 Act1 Act2
            end
            
            
            
            if and(sum(cond.both.aversive)>0 , sum(cond.vHPC.aversive)>0)
                templates1 = patterns.all.aversive(:,cond.both.aversive);
                templates2 = patterns.all.aversive(:,cond.vHPC.aversive);
                
                Act1 = (Activities.aversive(cond.both.aversive,:));
                Act2 = (Activities.aversive(cond.vHPC.aversive,:));
                
                for i = 1:size(templates1,2)
                    [pks1,loc1] = findpeaks(zscore(Act1(i,:)),bins,'MinPeakHeight',5);
                    for ii = 1:size(templates2,2)
                        [pks2,loc2] = findpeaks(zscore(Act2(ii,:)),bins,'MinPeakHeight',5);
                        % for Pre
                        x = Restrict(loc1,Is.aversive.pre);
                        y = Restrict(loc2,Is.aversive.pre);
                        if and(length(x)>5 , length(y)>5)
                            [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                            [ccg,T] = CCG(s,ids,'binSize',0.025,'duration',dur,'smooth',sm,'mode','ccg');
                            ccg = ccg(:,1,2)./sum(ccg(:,1,2));
                            cross.vHPC_joint.aversive.Pre = [cross.vHPC_joint.aversive.Pre , ccg];
                            time = T; clear ccg T x y s ids groups
                        end
                        
                        % for Post
                        x = Restrict(loc1,Is.aversive.post);
                        y = Restrict(loc2,Is.aversive.post);
                        if and(length(x)>5 , length(y)>5)
                            [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                            [ccg,T] = CCG(s,ids,'binSize',0.025,'duration',dur,'smooth',sm,'mode','ccg');
                            ccg = ccg(:,1,2)./sum(ccg(:,1,2));
                            cross.vHPC_joint.aversive.Post = [cross.vHPC_joint.aversive.Post , ccg];
                            time = T; clear ccg T x y s ids groups
                        end
                    end
                end
                clear  templates1 templates2 Act1 Act2
            end
            
            
            if and(sum(cond.dHPC.aversive)>0 , sum(cond.vHPC.aversive)>0)
                
                templates1 = patterns.all.aversive(:,cond.dHPC.aversive);
                templates2 = patterns.all.aversive(:,cond.vHPC.aversive);
                
                Act1 = (Activities.aversive(cond.dHPC.aversive,:));
                Act2 = (Activities.aversive(cond.vHPC.aversive,:));
                
                for i = 1:size(templates1,2)
                    [pks1,loc1] = findpeaks(zscore(Act1(i,:)),bins,'MinPeakHeight',5);
                    for ii = 1:size(templates2,2)
                        [pks2,loc2] = findpeaks(zscore(Act2(ii,:)),bins,'MinPeakHeight',5);
                        % for Pre
                        x = Restrict(loc1,Is.aversive.pre);
                        y = Restrict(loc2,Is.aversive.pre);
                        if and(length(x)>5 , length(y)>5)
                            [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                            [ccg,T] = CCG(s,ids,'binSize',0.025,'duration',dur,'smooth',sm,'mode','ccg');
                            ccg = ccg(:,1,2)./sum(ccg(:,1,2));
                            cross.dHPC_vHPC.aversive.Pre = [cross.dHPC_vHPC.aversive.Pre , ccg];
                            time = T; clear ccg T x y s ids groups
                        end
                        
                        % for Post
                        x = Restrict(loc1,Is.aversive.post);
                        y = Restrict(loc2,Is.aversive.post);
                        if and(length(x)>5 , length(y)>5)
                            [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                            [ccg,T] = CCG(s,ids,'binSize',0.025,'duration',dur,'smooth',sm,'mode','ccg');
                            ccg = ccg(:,1,2)./sum(ccg(:,1,2));
                            cross.dHPC_vHPC.aversive.Post = [cross.dHPC_vHPC.aversive.Post , ccg];
                            time = T; clear ccg T x y s ids groups
                        end
                    end
                end
                clear  templates1 templates2 Act1 Act2
            end
            
            
            %% Reward
                if and(sum(cond.both.reward)>0 , sum(cond.dHPC.reward)>0)
                templates1 = patterns.all.reward(:,cond.both.reward);
                templates2 = patterns.all.reward(:,cond.dHPC.reward);
                
                Act1 = (Activities.reward(cond.both.reward,:));
                Act2 = (Activities.reward(cond.dHPC.reward,:));
                
                for i = 1:size(templates1,2)
                    [pks1,loc1] = findpeaks(zscore(Act1(i,:)),bins,'MinPeakHeight',5);
                    for ii = 1:size(templates2,2)
                        [pks2,loc2] = findpeaks(zscore(Act2(ii,:)),bins,'MinPeakHeight',5);
                        % for Pre
                        x = Restrict(loc1,Is.reward.pre);
                        y = Restrict(loc2,Is.reward.pre);
                        if and(length(x)>5 , length(y)>5)
                            [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                            [ccg,T] = CCG(s,ids,'binSize',0.025,'duration',dur,'smooth',sm,'mode','ccg');
                            ccg = ccg(:,1,2)./sum(ccg(:,1,2));
                            cross.dHPC_joint.reward.Pre = [cross.dHPC_joint.reward.Pre , ccg];
                            time = T; clear ccg T x y s ids groups
                        end
                        
                        % for Post
                        x = Restrict(loc1,Is.reward.post);
                        y = Restrict(loc2,Is.reward.post);
                        if and(length(x)>5 , length(y)>5)
                            [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                            [ccg,T] = CCG(s,ids,'binSize',0.025,'duration',dur,'smooth',sm,'mode','ccg');
                            ccg = ccg(:,1,2)./sum(ccg(:,1,2));
                            cross.dHPC_joint.reward.Post = [cross.dHPC_joint.reward.Post , ccg];
                            time = T; clear ccg T x y s ids groups
                        end
                    end
                end
                clear  templates1 templates2 Act1 Act2
            end
            
            
            
            if and(sum(cond.both.reward)>0 , sum(cond.vHPC.reward)>0)
                templates1 = patterns.all.reward(:,cond.both.reward);
                templates2 = patterns.all.reward(:,cond.vHPC.reward);
                
                Act1 = (Activities.reward(cond.both.reward,:));
                Act2 = (Activities.reward(cond.vHPC.reward,:));
                
                for i = 1:size(templates1,2)
                    [pks1,loc1] = findpeaks(zscore(Act1(i,:)),bins,'MinPeakHeight',5);
                    for ii = 1:size(templates2,2)
                        [pks2,loc2] = findpeaks(zscore(Act2(ii,:)),bins,'MinPeakHeight',5);
                        % for Pre
                        x = Restrict(loc1,Is.reward.pre);
                        y = Restrict(loc2,Is.reward.pre);
                        if and(length(x)>5 , length(y)>5)
                            [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                            [ccg,T] = CCG(s,ids,'binSize',0.025,'duration',dur,'smooth',sm,'mode','ccg');
                            ccg = ccg(:,1,2)./sum(ccg(:,1,2));
                            cross.vHPC_joint.reward.Pre = [cross.vHPC_joint.reward.Pre , ccg];
                            time = T; clear ccg T x y s ids groups
                        end
                        
                        % for Post
                        x = Restrict(loc1,Is.reward.post);
                        y = Restrict(loc2,Is.reward.post);
                        if and(length(x)>5 , length(y)>5)
                            [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                            [ccg,T] = CCG(s,ids,'binSize',0.025,'duration',dur,'smooth',sm,'mode','ccg');
                            ccg = ccg(:,1,2)./sum(ccg(:,1,2));
                            cross.vHPC_joint.reward.Post = [cross.vHPC_joint.reward.Post , ccg];
                            time = T; clear ccg T x y s ids groups
                        end
                    end
                end
                clear  templates1 templates2 Act1 Act2
            end
            
            
            if and(sum(cond.dHPC.reward)>0 , sum(cond.vHPC.reward)>0)
                
                templates1 = patterns.all.reward(:,cond.dHPC.reward);
                templates2 = patterns.all.reward(:,cond.vHPC.reward);
                
                Act1 = (Activities.reward(cond.dHPC.reward,:));
                Act2 = (Activities.reward(cond.vHPC.reward,:));
                
                for i = 1:size(templates1,2)
                    [pks1,loc1] = findpeaks(zscore(Act1(i,:)),bins,'MinPeakHeight',5);
                    for ii = 1:size(templates2,2)
                        [pks2,loc2] = findpeaks(zscore(Act2(ii,:)),bins,'MinPeakHeight',5);
                        % for Pre
                        x = Restrict(loc1,Is.reward.pre);
                        y = Restrict(loc2,Is.reward.pre);
                        if and(length(x)>5 , length(y)>5)
                            [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                            [ccg,T] = CCG(s,ids,'binSize',0.025,'duration',dur,'smooth',sm,'mode','ccg');
                            ccg = ccg(:,1,2)./sum(ccg(:,1,2));
                            cross.dHPC_vHPC.reward.Pre = [cross.dHPC_vHPC.reward.Pre , ccg];
                            time = T; clear ccg T x y s ids groups
                        end
                        
                        % for Post
                        x = Restrict(loc1,Is.reward.post);
                        y = Restrict(loc2,Is.reward.post);
                        if and(length(x)>5 , length(y)>5)
                            [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                            [ccg,T] = CCG(s,ids,'binSize',0.025,'duration',dur,'smooth',sm,'mode','ccg');
                            ccg = ccg(:,1,2)./sum(ccg(:,1,2));
                            cross.dHPC_vHPC.reward.Post = [cross.dHPC_vHPC.reward.Post , ccg];
                            time = T; clear ccg T x y s ids groups
                        end
                    end
                end
                clear  templates1 templates2 Act1 Act2
            end
            
            clear aversiveTS aversiveTS_run baselineTS bins Cell_type_classification
            clear cellulartype clusters Events group_dHPC group_vHPC i ii iii K
            clear Kinfo NREM REM numberD numberV p parameter patterns post pre
            clear SpikesD SpikesV spks spks_dHPC spks_vHPC Thresholded WAKE
            clear Kinfo K ii iii ripples ripplesD ripplesV ripple_event Spikes
            clear x y loc1 loc2 pks1 pks2 cond Is rewardTS rewardTS_run aversiveTS aversiveTS_run
            clear segments c
        end
    end
end

% %Detection of lag where maximal value was detected
% %Aversive
% [h p] = max(cross_members.aversive.pre);
% lags.aversive.pre = [];
% for i = 1 : size(p,2)
%     lags.aversive.pre = [lags.aversive.pre ; time(p(i))];
% end
%
% [h p] = max(cross_members.aversive.post);
% lags.aversive.post = [];
% for i = 1 : size(p,2)
%     lags.aversive.post = [lags.aversive.post ; time(p(i))];
% end
%
% % Reward
% [h p] = max(cross_members.reward.pre);
% lags.reward.pre = [];
% for i = 1 : size(p,2)
%     lags.reward.pre = [lags.reward.pre ; time(p(i))];
% end
%
% [h p] = max(cross_members.reward.post);
% lags.reward.post = [];
% for i = 1 : size(p)
%     lags.reward.post = [lags.reward.post ; time(p(i))];
% end
%

figure,
plot(time,nanmean(cross.dHPC_joint.aversive.Post,2),'r')
hold on
plot(time,nanmean(cross.dHPC_joint.aversive.Pre,2),'k')

figure
plot(time,nanmean(cross.vHPC_joint.aversive.Post,2),'r')
hold on
plot(time,nanmean(cross.vHPC_joint.aversive.Pre,2),'k')

figure
plot(time,nanmean(cross.dHPC_vHPC.aversive.Post,2),'r')
hold on
plot(time,nanmean(cross.dHPC_vHPC.aversive.Pre,2),'k')

figure,
plot(time,nanmean(cross.dHPC_joint.reward.Post,2),'b')
hold on
plot(time,nanmean(cross.dHPC_joint.reward.Pre,2),'k')

figure,
plot(time,nanmean(cross.vHPC_joint.reward.Post,2),'b')
hold on
plot(time,nanmean(cross.vHPC_joint.reward.Pre,2),'k')

figure,
plot(time,nanmean(cross.dHPC_vHPC.reward.Post,2),'b')
hold on
plot(time,nanmean(cross.dHPC_vHPC.reward.Pre,2),'k')

end