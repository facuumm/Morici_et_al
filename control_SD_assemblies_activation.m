clc
clear

%% Parameters
path = {'E:\Rat126\Ephys\in_Pyr';'E:\Rat103\usable';'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat132\recordings\in_pyr';'E:\Rat165\in_pyr\'};%List of folders from the path

percentages.aversive = [];
percentages.reward = [];
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
            %% --- Aversive ---
            disp('Lets go for the assemblies')
            if isfile('dorsalventral_assemblies_aversiveVF.mat')
                disp('Loading Aversive template')
                load('dorsalventral_assemblies_aversiveVF.mat')
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
            if isfile('dorsalventral_assemblies_rewardVF.mat')
                load('dorsalventral_assemblies_rewardVF.mat')
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
            cellulartype = [K(:,1) , K(:,4)];
            [Spikes , bins , Clusters] = spike_train_construction([spks_dHPC;spks_vHPC], clusters.all, cellulartype, 0.025, limits, [] , false, true);
            clear limits
            
            % Restrict data to NREM sleep
            IN = InIntervals(bins,NREM.all);
            Spikes = Spikes(IN,:);
            bins = bins(IN); clear IN
            
            if sum(cond.both.aversive)>0
                [P] = bin_proportions_in_function_of_Th(patterns.all.aversive , cond.both.aversive , [bins' Spikes], Thresholded.aversive.all, [numberD numberV]);
                percentages.aversive = [percentages.aversive ; P]; clear P
            end
            
            if sum(cond.both.reward)>0
                [P] = bin_proportions_in_function_of_Th(patterns.all.reward , cond.both.reward , [bins' Spikes], Thresholded.reward.all, [numberD numberV]);
                percentages.reward = [percentages.reward ; P]; clear P
            end
            
            clear aversiveTS aversiveTS_run rewardTS rewardTS_run baselineTS Cell_type_classification cellulartype
            clear clusters Clusters cond group_dHPC group_vHPC K Kinfo numberD numberV
            clear Spikes spks spks_dHPC spks_vHPC Thresholded segments patterns bins ans NREM
        end
    end
end

figure
subplot(4,4,1)
histogram(percentages.aversive(:,1),'BinEdges',[0:5:100],'Normalization','probability'),ylim([0 0.4]), hold on, xline(50)
subplot(4,4,2)
histogram(percentages.aversive(:,2),'BinEdges',[0:5:100],'Normalization','probability'),ylim([0 0.4]), hold on, xline(50)
subplot(4,4,3)
histogram(percentages.aversive(:,3),'BinEdges',[0:5:100],'Normalization','probability'),ylim([0 0.4]), hold on, xline(50)
subplot(4,4,4)
histogram(percentages.aversive(:,4),'BinEdges',[0:5:100],'Normalization','probability'),ylim([0 0.4]), hold on, xline(50)
subplot(4,4,5)
histogram(percentages.aversive(:,5),'BinEdges',[0:5:100],'Normalization','probability'),ylim([0 0.4]), hold on, xline(50)
subplot(4,4,6)
histogram(percentages.aversive(:,6),'BinEdges',[0:5:100],'Normalization','probability'),ylim([0 0.4]), hold on, xline(50)
subplot(4,4,7)
histogram(percentages.aversive(:,7),'BinEdges',[0:5:100],'Normalization','probability'),ylim([0 0.4]), hold on, xline(50)
subplot(4,4,8)
histogram(percentages.aversive(:,8),'BinEdges',[0:5:100],'Normalization','probability'),ylim([0 0.4]), hold on, xline(50)
subplot(4,4,9)
histogram(percentages.aversive(:,9),'BinEdges',[0:5:100],'Normalization','probability'),ylim([0 0.4]), hold on, xline(50)
subplot(4,4,10)
histogram(percentages.aversive(:,10),'BinEdges',[0:5:100],'Normalization','probability'),ylim([0 0.4]), hold on, xline(50)
subplot(4,4,11)
histogram(percentages.aversive(:,11),'BinEdges',[0:5:100],'Normalization','probability'),ylim([0 0.4]), hold on, xline(50)
subplot(4,4,12)
histogram(percentages.aversive(:,12),'BinEdges',[0:5:100],'Normalization','probability'),ylim([0 0.4]), hold on, xline(50)
subplot(4,4,13)
histogram(percentages.aversive(:,13),'BinEdges',[0:5:100],'Normalization','probability'),ylim([0 0.4]), hold on, xline(50)
subplot(4,4,14)
histogram(percentages.aversive(:,14),'BinEdges',[0:5:100],'Normalization','probability'),ylim([0 0.4]), hold on, xline(50)
subplot(4,4,15)
histogram(percentages.aversive(:,15),'BinEdges',[0:5:100],'Normalization','probability'),ylim([0 0.4]), hold on, xline(50)
subplot(4,4,16)
histogram(percentages.aversive(:,16),'BinEdges',[0:5:100],'Normalization','probability'),ylim([0 0.4]), hold on, xline(50)

figure
subplot(4,4,1)
histogram(percentages.reward(:,1),'BinEdges',[0:5:100],'Normalization','probability'),ylim([0 0.4]), hold on, xline(50)
subplot(4,4,2)
histogram(percentages.reward(:,2),'BinEdges',[0:5:100],'Normalization','probability'),ylim([0 0.4]), hold on, xline(50)
subplot(4,4,3)
histogram(percentages.reward(:,3),'BinEdges',[0:5:100],'Normalization','probability'),ylim([0 0.4]), hold on, xline(50)
subplot(4,4,4)
histogram(percentages.reward(:,4),'BinEdges',[0:5:100],'Normalization','probability'),ylim([0 0.4]), hold on, xline(50)
subplot(4,4,5)
histogram(percentages.reward(:,5),'BinEdges',[0:5:100],'Normalization','probability'),ylim([0 0.4]), hold on, xline(50)
subplot(4,4,6)
histogram(percentages.reward(:,6),'BinEdges',[0:5:100],'Normalization','probability'),ylim([0 0.4]), hold on, xline(50)
subplot(4,4,7)
histogram(percentages.reward(:,7),'BinEdges',[0:5:100],'Normalization','probability'),ylim([0 0.4]), hold on, xline(50)
subplot(4,4,8)
histogram(percentages.reward(:,8),'BinEdges',[0:5:100],'Normalization','probability'),ylim([0 0.4]), hold on, xline(50)
subplot(4,4,9)
histogram(percentages.reward(:,9),'BinEdges',[0:5:100],'Normalization','probability'),ylim([0 0.4]), hold on, xline(50)
subplot(4,4,10)
histogram(percentages.reward(:,10),'BinEdges',[0:5:100],'Normalization','probability'),ylim([0 0.4]), hold on, xline(50)
subplot(4,4,11)
histogram(percentages.reward(:,11),'BinEdges',[0:5:100],'Normalization','probability'),ylim([0 0.4]), hold on, xline(50)
subplot(4,4,12)
histogram(percentages.reward(:,12),'BinEdges',[0:5:100],'Normalization','probability'),ylim([0 0.4]), hold on, xline(50)
subplot(4,4,13)
histogram(percentages.reward(:,13),'BinEdges',[0:5:100],'Normalization','probability'),ylim([0 0.4]), hold on, xline(50)
subplot(4,4,14)
histogram(percentages.reward(:,14),'BinEdges',[0:5:100],'Normalization','probability'),ylim([0 0.4]), hold on, xline(50)
subplot(4,4,15)
histogram(percentages.reward(:,15),'BinEdges',[0:5:100],'Normalization','probability'),ylim([0 0.4]), hold on, xline(50)
subplot(4,4,16)
histogram(percentages.reward(:,16),'BinEdges',[0:5:100],'Normalization','probability'),ylim([0 0.4]), hold on, xline(50)