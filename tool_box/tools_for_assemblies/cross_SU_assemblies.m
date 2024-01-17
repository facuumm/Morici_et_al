function [cross_members cross_nonmembers time lags] = cross_SU_assemblies(path,events,type_assemblies)
% HACER EL HELP DE ESTA FUNCION

% Initialization of structures that wil contain the outputs
cross_members.aversive.pre = [];    cross_members.aversive.post = [];
cross_nonmembers.aversive.pre = [];    cross_nonmembers.aversive.post = [];
cross_members.reward.pre = [];    cross_members.reward.post = [];
cross_nonmembers.reward.pre = [];    cross_nonmembers.reward.post = [];

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
        
        % Defining what condition was first
        if aversiveTS_run(1) < rewardTS_run(1)
            config = 1;
        else
            config = 2;
        end
                
        %% load sleep states
        disp('Uploading sleep scoring')
        x = dir([cd,'\*-states.mat']);    states = load([cd,'\',x.name]);    states = states.states;
        REM.all = ToIntervals(states==5);    NREM.all = ToIntervals(states==3);    WAKE.all = ToIntervals(states==1);
        %         REM.all(REM.all(:,2)-REM.all(:,1)<60,:) = [];
        clear x states
        
        NREM.all(NREM.all(:,2)-NREM.all(:,1)<60*1,:)=[];
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
%         if criteria_type == 0 %pyr
            cellulartype = [K(:,1) , K(:,3)];
%         elseif criteria_type == 1 % int
%             cellulartype = [K(:,1) , K(:,4)];
%         elseif criteria_type == 2 % all
%             cellulartype = [K(:,1) , ones(length(K),1)];
%         end
        
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
                if or(a > 0 ,  r > 0)
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
                if or(a > 0 ,  r > 0)
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
        
        
            %% --- Aversive ---
            disp('Lets go for the assemblies')
            if isfile('dorsalventral_assemblies_aversive3.mat')
                disp('Loading Aversive template')
                load('dorsalventral_assemblies_aversive3.mat')
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
            if isfile('dorsalventral_assemblies_reward3.mat')
                load('dorsalventral_assemblies_reward3.mat')
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
                    
        if strcmp(type_assemblies,'Both')
            if strcmp(events,'NREM')
                Events = NREM;
            elseif strcmp(events,'REM')
                Events = REM;
            end
            
            if sum(cond.both.aversive)>0
                %Defining the assemblies to work with
                p = Thresholded.aversive.all(:,cond.both.aversive);
                
                %Definint Events to study for this assembly
                if config==1
                    pre = Events.baseline;
                    post = Events.aversive;
                else
                    pre = Events.reward;
                    post = Events.aversive;                 
                end
                
                % Here I will iterate across dHPC members and perform a CCG
                % between members from the vHPC.
                for i = 1 : sum(cond.both.aversive)
                    parameter = p(:,i);
                    for ii = 1 :size(clusters.dHPC,1)
                        if parameter(ii)
                            SpikesD = spks_dHPC(spks_dHPC(:,1)==clusters.dHPC(ii),2);
                            for iii = 1 :size(clusters.vHPC,1)
                                if parameter(size(clusters.dHPC,1)+iii)
                                    SpikesV = spks_vHPC(spks_vHPC(:,1)==clusters.vHPC(iii),2);
                                    
                                    % For Pre
                                    x = Restrict(SpikesD,pre);
                                    y = Restrict(SpikesV,pre);
                                    if and(not(isempty(x)),not(isempty(y)))
                                        [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                                        [ccg,tttt] = CCG(s,ids,'binSize',0.003,'duration',0.5,'smooth',0,'mode','ccg'); %ccg calculation
                                        time = tttt;
                                        cross_members.aversive.pre = [cross_members.aversive.pre , ccg(:,1,2)];
                                        clear x y s ids groups ccg tttt
                                    end
                                    
                                    % for post
                                    x = Restrict(SpikesD,post);
                                    y = Restrict(SpikesV,post);
                                    if and(not(isempty(x)),not(isempty(y)))
                                        [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                                        [ccg,tttt] = CCG(s,ids,'binSize',0.003,'duration',0.5,'smooth',0,'mode','ccg'); %ccg calculation
                                        time = tttt;
                                        cross_members.aversive.post = [cross_members.aversive.post , ccg(:,1,2)];
                                        clear x y s ids groups ccg tttt
                                    end
                                else
                                    SpikesV = spks_vHPC(spks_vHPC(:,1)==clusters.vHPC(iii),2);
                                    
                                    % For Pre
                                    x = Restrict(SpikesD,pre);
                                    y = Restrict(SpikesV,pre);
                                    if and(not(isempty(x)),not(isempty(y)))
                                        [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                                        [ccg,tttt] = CCG(s,ids,'binSize',0.003,'duration',0.5,'smooth',0,'mode','ccg'); %ccg calculation
                                        time = tttt;
                                        cross_nonmembers.aversive.pre = [cross_nonmembers.aversive.pre , ccg(:,1,2)];
                                        clear x y s ids groups ccg tttt
                                    end
                                    
                                    % for post
                                    x = Restrict(SpikesD,post);
                                    y = Restrict(SpikesV,post);
                                    if and(not(isempty(x)),not(isempty(y)))
                                        [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                                        [ccg,tttt] = CCG(s,ids,'binSize',0.003,'duration',0.5,'smooth',0,'mode','ccg'); %ccg calculation
                                        time = tttt;
                                        cross_nonmembers.aversive.post = [cross_nonmembers.aversive.post , ccg(:,1,2)];
                                        clear x y s ids groups ccg tttt
                                    end
                                end
                            end
                        end
                    end
                end
            end
            
            if sum(cond.both.reward)>0
               %Defining the assemblies to work with
                p = Thresholded.reward.all(:,cond.both.reward);
                
                %Definint Events to study for this assembly
                if config==1
                    pre = Events.aversive;
                    post = Events.reward;
                else
                    pre = Events.baseline;
                    post = Events.reward;                 
                end
                
                % Here I will iterate across dHPC members and perform a CCG
                % between members from the vHPC.
                for i = 1 : sum(cond.both.reward)
                    parameter = p(:,i);
                    for ii = 1 :size(clusters.dHPC,1)
                        if parameter(ii)
                            SpikesD = spks_dHPC(spks_dHPC(:,1)==clusters.dHPC(ii),2);
                            for iii = 1 :size(clusters.vHPC,1)
                                if parameter(size(clusters.dHPC,1)+iii)
                                    SpikesV = spks_vHPC(spks_vHPC(:,1)==clusters.vHPC(iii),2);
                                    
                                    % For Pre
                                    x = Restrict(SpikesD,pre);
                                    y = Restrict(SpikesV,pre);
                                    if and(not(isempty(x)),not(isempty(y)))
                                    [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                                    [ccg,tttt] = CCG(s,ids,'binSize',0.003,'duration',0.5,'smooth',0,'mode','ccg'); %ccg calculation
                                    time = tttt;
                                    cross_members.reward.pre = [cross_members.reward.pre , ccg(:,1,2)];
                                    clear x y s ids groups ccg tttt
                                    end
                                    
                                    % for post
                                    x = Restrict(SpikesD,post);
                                    y = Restrict(SpikesV,post);
                                    if and(not(isempty(x)),not(isempty(y)))
                                    [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                                    [ccg,tttt] = CCG(s,ids,'binSize',0.003,'duration',0.5,'smooth',0,'mode','ccg'); %ccg calculation
                                    time = tttt;
                                    cross_members.reward.post = [cross_members.reward.post , ccg(:,1,2)];
                                    clear x y s ids groups ccg tttt SpikesV  
                                    end
                                else
                                    SpikesV = spks_vHPC(spks_vHPC(:,1)==clusters.vHPC(iii),2);
                                    
                                    % For Pre
                                    x = Restrict(SpikesD,pre);
                                    y = Restrict(SpikesV,pre);
                                    if and(not(isempty(x)),not(isempty(y)))
                                    [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                                    [ccg,tttt] = CCG(s,ids,'binSize',0.003,'duration',0.5,'smooth',0,'mode','ccg'); %ccg calculation
                                    time = tttt;
                                    cross_nonmembers.reward.pre = [cross_nonmembers.reward.pre , ccg(:,1,2)];
                                    clear x y s ids groups ccg tttt
                                    end
                                    
                                    % for post
                                    x = Restrict(SpikesD,post);
                                    y = Restrict(SpikesV,post);
                                    if and(not(isempty(x)),not(isempty(y)))
                                    [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                                    [ccg,tttt] = CCG(s,ids,'binSize',0.003,'duration',0.5,'smooth',0,'mode','ccg'); %ccg calculation
                                    time = tttt;
                                    cross_nonmembers.reward.post = [cross_nonmembers.reward.post , ccg(:,1,2)];
                                    clear x y s ids groups ccg tttt SpikesV  
                                    end
                                end
                            end
                            clear SpikesD
                        end
                    end
                end                
            end
        end
        clear aversiveTS aversiveTS_run baselineTS bins Cell_type_classification
        clear cellulartype clusters Events group_dHPC group_vHPC i ii iii K
        clear Kinfo NREM REM numberD numberV p parameter patterns post pre
        clear SpikesD SpikesV spks spks_dHPC spks_vHPC Thresholded WAKE
    end
end

%Detection of lag where maximal value was detected
%Aversive
[h p] = max(cross_members.aversive.pre);
lags.aversive.pre = [];
for i = 1 : size(p,2)
    lags.aversive.pre = [lags.aversive.pre ; time(p(i))];
end

[h p] = max(cross_members.aversive.post);
lags.aversive.post = [];
for i = 1 : size(p,2)
    lags.aversive.post = [lags.aversive.post ; time(p(i))];
end

% Reward
[h p] = max(cross_members.reward.pre);
lags.reward.pre = [];
for i = 1 : size(p,2)
    lags.reward.pre = [lags.reward.pre ; time(p(i))];
end

[h p] = max(cross_members.reward.post);
lags.reward.post = [];
for i = 1 : size(p)
    lags.reward.post = [lags.reward.post ; time(p(i))];
end



end