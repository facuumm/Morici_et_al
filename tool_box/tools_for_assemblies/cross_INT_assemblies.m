function [cross time] = cross_INT_assemblies(path)
% This function a CCG between putative-interneuron spiking activity and the
% timestamps from assemblies peaks.
%
% INPUTS
% path: cell, inside each cell you should put the path of each session you
%       want to analyze
%
% OUTPUT
% cross_members: structure, it contains the CCG for aversive and reward
%                joint assemblies, during the Pre and Post sleep.
%                cross_members.aversive.pre = [];    cross_members.aversive.post = [];
%                cross_members.reward.pre = [];    cross_members.reward.post = [];
%
% cross_nonmembers: structure, same as the output described above.
%
% time: column vector, it contains the time axis for the CCG plot.
%
% lags: structure, it contains the max value of each CCG for each
%       condition. SImilar structure as thoe outputs described above.
%
%
% other functions: CCG from FMA toolbox
% Morci Juan Facundo 01/2024


% Initialization of structures that wil contain the outputs
cross.dInt.dHPCA = [];   cross.vInt.dHPCA = [];
cross.dInt.vHPCA = [];   cross.vInt.vHPCA = [];
cross.dInt.jointA = [];  cross.vInt.jointA = [];
cross.dInt.dHPCR = [];   cross.vInt.dHPCR = [];
cross.dInt.vHPCR = [];   cross.vInt.vHPCR = [];
cross.dInt.jointR = [];  cross.vInt.jointR = [];

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
            cellulartype = [K(:,1) , K(:,4)];
            [Spikes , bins , Clusters] = spike_train_construction([spks_dHPC;spks_vHPC], clusters.all, cellulartype, 0.025, limits, [] , true, true);
            clear limits
            
            % Joint assemblies
            if sum(cond.both.aversive)>0
                if not(isempty(clusters.int.vHPC))
                    for i = 1:size(clusters.int.vHPC,1)
                        tmp = Restrict(spks(spks(:,1)==clusters.int.vHPC(i),2),aversiveTS_run./1000);
                        [R , ttttt] = triggered_CCG(patterns.all.aversive , cond.both.aversive , [bins' , Spikes], 1 , 0.025 , tmp , 5 , [],'Probability'); clear b
                        cross.vInt.jointA = [cross.vInt.jointA , R]; clear R
                        time = ttttt; clear ttttt
                    end
                end
                
                if not(isempty(clusters.int.dHPC))
                    for i = 1:size(clusters.int.dHPC,1)
                        tmp = Restrict(spks(spks(:,1)==clusters.int.dHPC(i),2),aversiveTS_run./1000);
                        [R , ttttt] = triggered_CCG(patterns.all.aversive , cond.both.aversive , [bins' , Spikes], 1 , 0.025 , tmp , 5 , [],'Probability'); clear b
                        cross.dInt.jointA = [cross.dInt.jointA , R]; clear R
                        time = ttttt; clear ttttt
                    end
                end
            end
            
            
            if sum(cond.both.reward)>0
                if not(isempty(clusters.int.vHPC))
                    for i = 1:size(clusters.int.vHPC,1)
                        tmp = Restrict(spks(spks(:,1)==clusters.int.vHPC(i),2),rewardTS_run./1000);
                        [R , ttttt] = triggered_CCG(patterns.all.reward , cond.both.reward , [bins' , Spikes], 1 , 0.025 , tmp , 5 , [],'Probability'); clear b
                        cross.vInt.jointR = [cross.vInt.jointR , R]; clear R
                        time = ttttt; clear ttttt
                    end
                end
                
                if not(isempty(clusters.int.dHPC))
                    for i = 1:size(clusters.int.dHPC,1)
                        tmp = Restrict(spks(spks(:,1)==clusters.int.dHPC(i),2),rewardTS_run./1000);
                        [R , ttttt] = triggered_CCG(patterns.all.reward , cond.both.reward , [bins' , Spikes], 1 , 0.025 , tmp , 5 , [],'Probability'); clear b
                        cross.dInt.jointR = [cross.dInt.jointR , R]; clear R
                        time = ttttt; clear ttttt
                    end
                end
            end
            
            % dHPC assemblies
            if sum(cond.dHPC.aversive)>0
                if not(isempty(clusters.int.vHPC))
                    for i = 1:size(clusters.int.vHPC,1)
                        tmp = Restrict(spks(spks(:,1)==clusters.int.vHPC(i),2),aversiveTS_run./1000);
                        [R , ttttt] = triggered_CCG(patterns.all.aversive , cond.dHPC.aversive , [bins' , Spikes], 1 , 0.025 , tmp , 5 , [],'Probability'); clear b
                        cross.vInt.dHPCA = [cross.vInt.dHPCA , R]; clear R
                        time = ttttt; clear ttttt
                    end
                end
                
                if not(isempty(clusters.int.dHPC))
                    for i = 1:size(clusters.int.dHPC,1)
                        tmp = Restrict(spks(spks(:,1)==clusters.int.dHPC(i),2),aversiveTS_run./1000);
                        [R , ttttt] = triggered_CCG(patterns.all.aversive , cond.dHPC.aversive , [bins' , Spikes], 1 , 0.025 , tmp , 5 , [],'Probability'); clear b
                        cross.dInt.dHPCA = [cross.dInt.dHPCA , R]; clear R
                        time = ttttt; clear ttttt
                    end
                end
            end
            
            
            if sum(cond.dHPC.reward)>0
                if not(isempty(clusters.int.vHPC))
                    for i = 1:size(clusters.int.vHPC,1)
                        tmp = Restrict(spks(spks(:,1)==clusters.int.vHPC(i),2),rewardTS_run./1000);
                        [R , ttttt] = triggered_CCG(patterns.all.reward , cond.dHPC.reward , [bins' , Spikes], 1 , 0.025 , tmp , 5 , [],'Probability'); clear b
                        cross.vInt.dHPCR = [cross.vInt.dHPCR , R]; clear R
                        time = ttttt; clear ttttt
                    end
                end
                
                if not(isempty(clusters.int.dHPC))
                    for i = 1:size(clusters.int.dHPC,1)
                        tmp = Restrict(spks(spks(:,1)==clusters.int.dHPC(i),2),rewardTS_run./1000);
                        [R , ttttt] = triggered_CCG(patterns.all.reward , cond.dHPC.reward , [bins' , Spikes], 1 , 0.025 , tmp , 5 , [],'Probability'); clear b
                        cross.dInt.dHPCR = [cross.dInt.dHPCR , R]; clear R
                        time = ttttt; clear ttttt
                    end
                end
            end
            
            
            % vHPC assemblies
            if sum(cond.vHPC.aversive)>0
                if not(isempty(clusters.int.vHPC))
                    for i = 1:size(clusters.int.vHPC,1)
                        tmp = Restrict(spks(spks(:,1)==clusters.int.vHPC(i),2),aversiveTS_run./1000);
                        [R , ttttt] = triggered_CCG(patterns.all.aversive , cond.vHPC.aversive , [bins' , Spikes], 1 , 0.025 , tmp , 5 , [],'Probability'); clear b
                        cross.vInt.vHPCA = [cross.vInt.vHPCA , R]; clear R
                        time = ttttt; clear ttttt
                    end
                end
                
                if not(isempty(clusters.int.dHPC))
                    for i = 1:size(clusters.int.dHPC,1)
                        tmp = Restrict(spks(spks(:,1)==clusters.int.dHPC(i),2),aversiveTS_run./1000);
                        [R , ttttt] = triggered_CCG(patterns.all.aversive , cond.vHPC.aversive , [bins' , Spikes], 1 , 0.025 , tmp , 5 , [],'Probability'); clear b
                        cross.dInt.dHPCA = [cross.dInt.dHPCA , R]; clear R
                        time = ttttt; clear ttttt
                    end
                end
            end
            
            
            if sum(cond.vHPC.reward)>0
                if not(isempty(clusters.int.vHPC))
                    for i = 1:size(clusters.int.vHPC,1)
                        tmp = Restrict(spks(spks(:,1)==clusters.int.vHPC(i),2),rewardTS_run./1000);
                        [R , ttttt] = triggered_CCG(patterns.all.reward , cond.vHPC.reward , [bins' , Spikes], 1 , 0.025 , tmp , 5 , [],'Probability'); clear b
                        cross.vInt.vHPCR = [cross.vInt.vHPCR , R]; clear R
                        time = ttttt; clear ttttt
                    end
                end
                
                if not(isempty(clusters.int.dHPC))
                    for i = 1:size(clusters.int.dHPC,1)
                        tmp = Restrict(spks(spks(:,1)==clusters.int.dHPC(i),2),rewardTS_run./1000);
                        [R , ttttt] = triggered_CCG(patterns.all.reward , cond.vHPC.reward , [bins' , Spikes], 1 , 0.025 , tmp , 5 , [],'Probability'); clear b
                        cross.dInt.vHPCR = [cross.dInt.vHPCR , R]; clear R
                        time = ttttt; clear ttttt
                    end
                end
            end
            
            clear aversiveTS aversiveTS_run baselineTS bins Cell_type_classification
            clear cellulartype clusters Events group_dHPC group_vHPC i ii iii K
            clear Kinfo NREM REM numberD numberV p parameter patterns post pre
            clear SpikesD SpikesV spks spks_dHPC spks_vHPC Thresholded WAKE
            clear Kinfo K ii iii ripples ripplesD ripplesV ripple_event Spikes
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


end