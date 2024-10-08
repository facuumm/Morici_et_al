function [Pre Post T] = AssembliesPeaks_SU_CCG(path)
% This function calculates the Pre and Post sleep assemblies rate.
%
% INPUTS
% path: cell, inside each cell you should put the path of each session you
%       want to analyze
%
% OUTPUT
% Pre, Post and Run: Structure, it store the Assemblies Rate for each type
%                    of assemblies.
%
%               Architecture of each output:
%                   Pre.dHPC
%                      .vHPC.aversive
%                           .reward
%                   Post.dHPC
%                       .vHPC.aversive
%                            .reward
%
%
% Morci Juan Facundo 04/2024

% variables to use in the script
criteria_fr = 0;
criteria_type = 0; % criteria for celltype (0:pyr, 1:int, 2:all)
th = 5; % threshold for detecting peak assemblies
sm = 2;
dur = 0.6;

% storage variables
Pre.aversive.ccg = [] ;         Post.aversive.ccg = [] ;
Pre.reward.ccg = [] ;           Post.reward.ccg = [] ;

Pre.aversive.iterator = [] ;    Post.aversive.iterator = [] ;
Pre.reward.iterator = [] ;      Post.reward.iterator = [] ;

Pre.aversive.ccgM = [] ;        Post.aversive.ccgM = [] ;
Pre.reward.ccgM = [] ;          Post.reward.ccgM = [] ;

Pre.aversive.ccgM1 = [] ;        Post.aversive.ccgM1 = [] ;
Pre.reward.ccgM1 = [] ;          Post.reward.ccgM1 = [] ;

I.aversive.dHPC = [];                I.reward.dHPC = [];
I.aversive.vHPC = [];                I.reward.vHPC = [];

Post.aversive.members.dHPC = [];         Post.aversive.members.vHPC = [];
Post.reward.members.dHPC = [];           Post.reward.members.vHPC = [];
Pre.aversive.members.dHPC = [];          Pre.aversive.members.vHPC = [];
Pre.reward.members.dHPC = [];            Pre.reward.members.vHPC = [];

Post.aversive.nomembers.dHPC = [];       Post.aversive.nomembers.vHPC = [];
Post.reward.nomembers.dHPC = [];         Post.reward.nomembers.vHPC = [];
Pre.aversive.nomembers.dHPC = [];        Pre.aversive.nomembers.vHPC = [];
Pre.reward.nomembers.dHPC = [];          Pre.reward.nomembers.vHPC = [];

durations.dHPC = [];             durations.vHPC = [];
durations.event = [];            durations.bursts = [];

Pre.aversive.ccgD = [] ;        Post.aversive.ccgD = [] ;
Pre.reward.ccgD = [] ;          Post.reward.ccgD = [] ;

Pre.aversive.ccgV = [] ;        Post.aversive.ccgV = [] ;
Pre.reward.ccgV = [] ;          Post.reward.ccgV = [] ;

Pre.aversive.ccgB = [] ;        Post.aversive.ccgB = [] ;
Pre.reward.ccgB = [] ;          Post.reward.ccgB = [] ;
countA = 1;                     countR = 1;

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
        
        baselineTS = baselineTS./1000;
        aversiveTS = aversiveTS./1000;
        rewardTS = rewardTS./1000;
        aversiveTS_run = aversiveTS_run./1000;
        rewardTS_run = rewardTS_run./1000;
        
        NREM.baseline = Restrict(NREM.all,baselineTS);
        NREM.aversive = Restrict(NREM.all,aversiveTS);
        NREM.reward = Restrict(NREM.all,rewardTS);
        
        REM.baseline = Restrict(REM.all,baselineTS);
        REM.aversive = Restrict(REM.all,aversiveTS);
        REM.reward = Restrict(REM.all,rewardTS);
        
        %% Load ripples
        if exist('ripplesD_customized2.csv')
            ripplesD = table2array(readtable('ripplesD_customized2.csv'));
            RD = true;
        else
            RD = false;
        end
        
        if exist('ripplesV_customized2.csv')
            ripplesV = table2array(readtable('ripplesV_customized2.csv'));
            RV = true;
        else
            RV = false;
        end
        
        if and(RD,RV)
            % coordination
            coordinated = [];
            coordinatedV = [];
            cooridnated_event = [];
            for i = 1:length(ripplesD)
                r = ripplesD(i,:);
                tmp = sum(and(ripplesV(:,2)>= r(1,2)-0.1, ripplesV(:,2)<= r(1,2)+0.1));
                if tmp>0
                    z = ripplesV(and(ripplesV(:,2)>= r(1,2)-0.1, ripplesV(:,2)<= r(1,2)+0.1),:);
                    coordinatedV = [coordinatedV ; z];
                    [p,indice] = min(abs(r(2)-z(:,2)));
                    coordinated = [coordinated ; r];
                    
                    %                     temporal1 = sum(and(ripplesD(:,2) < r(2) + 0.4 , ripplesD(:,2) > r(2) - 0.4))>1;
                    %                     temporal2 = sum(and(ripplesV(:,2) < z(indice,2) + 0.4 , ripplesV(:,2) > z(indice,2) - 0.4))>1;
                    
                    %                     if and(not(temporal1),not(temporal1))
                    peak = min(min(z(indice,2),r(2)))+abs(z(indice,2)-r(2))/2;
                    low = min([r(1) , z(indice,1)]);
                    up = max([r(3) , z(indice,3)]);
                    cooridnated_event = [cooridnated_event ; low , peak , up];
                    durations.dHPC = [durations.dHPC ; r(3)-r(1)];
                    durations.vHPC = [durations.vHPC ; z(indice,3)-z(indice,1)];
                    durations.event = [durations.event ; up-low];
                    %                     end
                    clear tmp2 tmp1 p indice z peak low up
                end
                clear r
            end
            clear x tmp i
            
            [C,IA,IC] = unique(coordinatedV(:,1));
            coordinatedV  = coordinatedV(IA,:); clear C IA IC
            
            % Store events time stamps
            % dRipples
            ripples.dHPC.baseline = Restrict(ripplesD , NREM.baseline);
            ripples.dHPC.reward = Restrict(ripplesD , NREM.reward);
            ripples.dHPC.aversive = Restrict(ripplesD , NREM.aversive);
            % vRipples
            ripples.vHPC.baseline = Restrict(ripplesV , NREM.baseline);
            ripples.vHPC.reward = Restrict(ripplesV , NREM.reward);
            ripples.vHPC.aversive = Restrict(ripplesV , NREM.aversive);
            % coordinated dRipples
            ripples.dHPC.coordinated.all = coordinated;
            ripples.dHPC.uncoordinated.all = ripplesD(not(ismember(ripplesD(:,2) , coordinated(:,2))),:);
            ripples.dHPC.coordinated.baseline = Restrict(coordinated , NREM.baseline);
            ripples.dHPC.coordinated.reward = Restrict(coordinated , NREM.reward);
            ripples.dHPC.coordinated.aversive = Restrict(coordinated , NREM.aversive);
            % coordinated vRipples
            ripples.vHPC.coordinated.all = coordinatedV;
            ripples.vHPC.uncoordinated.all = ripplesV(not(ismember(ripplesV(:,2) , coordinatedV(:,2))),:);
            ripples.vHPC.coordinated.baseline = Restrict(coordinatedV , NREM.baseline);
            ripples.vHPC.coordinated.reward = Restrict(coordinatedV , NREM.reward);
            ripples.vHPC.coordinated.aversive = Restrict(coordinatedV , NREM.aversive);
            %coordinated event
            cooridnated_event((cooridnated_event(:,3)-cooridnated_event(:,1)<0.04),:) = [];
            ripple_event.baseline = Restrict(cooridnated_event,baselineTS);
            ripple_event.reward = Restrict(cooridnated_event,rewardTS);
            ripple_event.aversive = Restrict(cooridnated_event,aversiveTS);
            ripple_event.all = cooridnated_event;
            
            load([cd , '\coordinated_ripple_bursts.mat'])
            
            % Separation of bursts across emotional condition
            bursts.coordinated.DV(bursts.coordinated.DV(:,2)-bursts.coordinated.DV(:,1)<0.06, : ) = [];
            bursts.coordinated.DV(bursts.coordinated.DV(:,2)-bursts.coordinated.DV(:,1)>0.2, : ) = [];
            
            bursts.baseline = Restrict(bursts.coordinated.DV,NREM.baseline);
            bursts.aversive = Restrict(bursts.coordinated.DV,NREM.aversive);
            bursts.reward = Restrict(bursts.coordinated.DV,NREM.reward);
            durations.bursts = [durations.bursts ; bursts.coordinated.DV(:,2)-bursts.coordinated.DV(:,1)];
            
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
                a = length(Restrict(spks(spks(:,1)==cluster,2),aversiveTS_run)) / ((aversiveTS_run(2)-aversiveTS_run(1)));
                r = length(Restrict(spks(spks(:,1)==cluster,2),rewardTS_run)) / ((rewardTS_run(2)-rewardTS_run(1)));
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
                a = length(Restrict(spks(spks(:,1)==cluster,2),aversiveTS_run)) / ((aversiveTS_run(2)-aversiveTS_run(1)));
                r = length(Restrict(spks(spks(:,1)==cluster,2),rewardTS_run)) / ((rewardTS_run(2)-rewardTS_run(1)));
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
            
            limits = [0 segments.Var1(end)/1000];
            events = [];
            [Spikes , bins , Clusters] = spike_train_construction([spks_dHPC;spks_vHPC], clusters.all, cellulartype, 0.025, limits, events, false, true);
            clear limits events
            
            % --- Aversive ---
            disp('Lets go for the assemblies')
            if isfile('dorsalventral_assemblies_aversiveVF.mat')
                disp('Loading Aversive template')
                load('dorsalventral_assemblies_aversiveVF.mat')
                
                Thresholded.aversive = Th;
                patterns.aversive = pat;
                clear cond Th pat
                
                patterns.aversive = patterns.aversive .* Thresholded.aversive;
                
                % Detection of members
                if not(isempty(patterns.aversive))
                    if numberD>0
                        cond1 =  sum(Thresholded.aversive(1:size(clusters.dHPC,1),:),1)>0; %checking of dHPC SU
                        cond2 =  sum(Thresholded.aversive(size(clusters.dHPC,1)+1:end,:),1)>0; %checking of vHPC SU
                        cond.dHPC = and(cond1 , not(cond2));
                        cond.vHPC = and(cond2 , not(cond1));
                        cond.both = and(cond1 , cond2); clear cond1 cond2
                    else
                        cond1 =  logical(zeros(1,size(Thresholded.aversive,2))); %checking of dHPC SU
                        cond2 =  logical(ones(1,size(Thresholded.aversive,2))); %checking of vHPC SU
                        cond.dHPC = and(cond1 , not(cond2));
                        cond.vHPC = and(cond2 , not(cond1));
                        cond.both = and(cond1 , cond2); clear cond1 cond2
                    end
                else
                    cond1 =  false; %checking of dHPC SU
                    cond2 =  logical(0); %checking of vHPC SU
                    cond.dHPC = and(cond1 , not(cond2));
                    cond.vHPC = and(cond2 , not(cond1));
                    cond.both = and(cond1 , cond2); clear cond1 cond2
                end
                
                if sum(cond.both)>0
                    
                    patterns.aversive = patterns.aversive(:,cond.both);
                    
                    if aversiveTS_run(1)<rewardTS_run(1)
                        TS.pre = bursts.baseline;
                        TS.post = bursts.aversive;
                        TS.run.run1 = movement.aversive;
                        TS.run.run2 = movement.reward;
                    else
                        TS.pre = bursts.reward;
                        TS.post = bursts.aversive;
                        TS.run.run1 = movement.aversive;
                        TS.run.run2 = movement.reward;
                    end
%                     pks.aversive = assemblies_peaks_joint(patterns.aversive , cond.both , [bins' Spikes], Thresholded.aversive, [numberD numberV]);
                    pks.aversive = assemblies_peaks([bins' Spikes] , patterns.aversive , th);
                    
                    members = Thresholded.aversive(:,cond.both);
                    
                    % Iteration across assemblies
                    for i = 1 : size(pks.aversive,1)
                        
                        %% Pre aversive members
                        for ii = 1 : numberD
                            if (members(ii,i))
                                member1 = clusters.dHPC(ii);
                                y = spks_dHPC(spks_dHPC(:,1) == member1 , 2);
                                y = Restrict(y,TS.pre);
                                x = Restrict(pks.aversive{i}(:,1),TS.pre);
                                if and(size(y,1)>10 , size(x,1)>5)
                                    [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                                    [ccg,T] = CCG(s,ids,'binSize',0.005,'duration',1,'smooth',0,'mode','ccg');
                                    ccg = ccg(:,1,2)./length(x); ccg = ccg./0.005;
                                    
                                    Pre.aversive.members.dHPC = [Pre.aversive.members.dHPC , ccg]; clear ccg x y m s ids groups
                                    I.aversive.dHPC = [I.aversive.dHPC ; countA];
                                end
                            end
                        end
                        
                        for iii = 1 : numberV
                            if (members(iii+numberD,i))
                                member2 = clusters.vHPC(iii);
                                y = spks_vHPC(spks_vHPC(:,1) == member2 , 2);
                                y = Restrict(y,TS.pre);
                                x = Restrict(pks.aversive{i}(:,1),TS.pre);
                                if and(size(y,1)>10 , size(x,1)>5)
                                    [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                                    [ccg,T] = CCG(s,ids,'binSize',0.005,'duration',1,'smooth',0,'mode','ccg');
                                    ccg = ccg(:,1,2)./length(x); ccg = ccg./0.005;
                                    
                                    Pre.aversive.members.vHPC = [Pre.aversive.members.vHPC , ccg]; clear ccg x y m s ids groups
                                    I.aversive.vHPC = [I.aversive.vHPC ; countA];
                                end
                            end
                        end
                        
                        %% Post  aversive members
                        for ii = 1 : numberD
                            if (members(ii,i))
                                member1 = clusters.dHPC(ii);
                                y = spks_dHPC(spks_dHPC(:,1) == member1 , 2);
                                y = Restrict(y,TS.post);
                                x = Restrict(pks.aversive{i}(:,1),TS.post);
                                if and(size(y,1)>10 , size(x,1)>5)
                                    [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                                    [ccg,T] = CCG(s,ids,'binSize',0.005,'duration',1,'smooth',0,'mode','ccg');
                                    ccg = ccg(:,1,2)./length(x); ccg = ccg./0.005;
                                    ;
                                    Post.aversive.members.dHPC = [Post.aversive.members.dHPC , ccg]; clear ccg x y m s ids groups
                                    I.aversive.dHPC = [I.aversive.dHPC ; countA];
                                end
                            end
                        end
                        
                        for iii = 1 : numberV
                            if (members(iii+numberD,i))
                                member2 = clusters.vHPC(iii);
                                y = spks_vHPC(spks_vHPC(:,1) == member2 , 2);
                                y = Restrict(y,TS.post);
                                x = Restrict(pks.aversive{i}(:,1),TS.post);
                                if and(size(y,1)>10 , size(x,1)>5)
                                    [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                                    [ccg,T] = CCG(s,ids,'binSize',0.005,'duration',1,'smooth',0,'mode','ccg');
                                    ccg = ccg(:,1,2)./length(x); ccg = ccg./0.005;
                                    
                                    Post.aversive.members.vHPC = [Post.aversive.members.vHPC , ccg]; clear ccg x y m s ids groups
                                    I.aversive.vHPC = [I.aversive.vHPC ; countA];
                                end
                            end
                        end
                        
                    end
                    
                    clear pks TS
                    
                    
                    % --- Reward ---
                    if isfile('dorsalventral_assemblies_rewardVF.mat')
                        disp('Loading Reward template')
                        load('dorsalventral_assemblies_rewardVF.mat')
                        
                        Thresholded.reward = Th;
                        patterns.reward = pat;
                        clear cond Th pat
                        
                        patterns.reward = patterns.reward .* Thresholded.reward;
                        
                        % Detection of members
                        if not(isempty(patterns.reward))
                            if numberD>0
                                cond1 =  sum(Thresholded.reward(1:size(clusters.dHPC,1),:),1)>0; %checking of dHPC SU
                                cond2 =  sum(Thresholded.reward(size(clusters.dHPC,1)+1:end,:),1)>0; %checking of vHPC SU
                                cond.dHPC = and(cond1 , not(cond2));
                                cond.vHPC = and(cond2 , not(cond1));
                                cond.both = and(cond1 , cond2); clear cond1 cond2
                            else
                                cond1 =  logical(zeros(1,size(Thresholded.reward,2))); %checking of dHPC SU
                                cond2 =  logical(ones(1,size(Thresholded.reward,2))); %checking of vHPC SU
                                cond.dHPC = and(cond1 , not(cond2));
                                cond.vHPC = and(cond2 , not(cond1));
                                cond.both = and(cond1 , cond2); clear cond1 cond2
                            end
                        else
                            cond1 =  false; %checking of dHPC SU
                            cond2 =  logical(0); %checking of vHPC SU
                            cond.dHPC = and(cond1 , not(cond2));
                            cond.vHPC = and(cond2 , not(cond1));
                            cond.both = and(cond1 , cond2); clear cond1 cond2
                        end
                        
                        
                        if sum(cond.both)>0
                            
                            patterns.reward = patterns.reward(:,cond.both);
                            
                            if aversiveTS_run(1)<rewardTS_run(1)
                                TS.pre = bursts.aversive;
                                TS.post = bursts.reward;
                                TS.run.run1 = movement.reward;
                                TS.run.run2 = movement.aversive;
                            else
                                TS.pre = bursts.baseline;
                                TS.post = bursts.reward;
                                TS.run.run2 = movement.aversive;
                            end
%                             pks.reward = assemblies_peaks_joint(patterns.reward , cond.both , [bins' Spikes], Thresholded.reward, [numberD numberV]);
                            pks.reward = assemblies_peaks([bins' Spikes] , patterns.reward , th);
                            members = Thresholded.reward(:,cond.both);
                            
                            % Iteration across assemblies
                            for i = 1 : size(pks.reward,1)
                                
                                %% Pre reward members
                                for ii = 1 : numberD
                                    if (members(ii,i))
                                        member1 = clusters.dHPC(ii);
                                        y = spks_dHPC(spks_dHPC(:,1) == member1 , 2);
                                        y = Restrict(y,TS.pre);
                                        x = Restrict(pks.reward{i}(:,1),TS.pre);
                                        if and(size(y,1)>10 , size(x,1)>5)
                                            [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                                            [ccg,T] = CCG(s,ids,'binSize',0.005,'duration',1,'smooth',0,'mode','ccg');
                                            ccg = ccg(:,1,2)./length(x); ccg = ccg./0.005;
                                            
                                            Pre.reward.members.dHPC = [Pre.reward.members.dHPC , ccg]; clear ccg x y m s ids groups
                                            I.reward.dHPC = [I.reward.dHPC ; countR];
                                        end
                                    end
                                end
                                
                                for iii = 1 : numberV
                                    if (members(iii+numberD,i))
                                        member2 = clusters.vHPC(iii);
                                        y = spks_vHPC(spks_vHPC(:,1) == member2 , 2);
                                        y = Restrict(y,TS.pre);
                                        x = Restrict(pks.reward{i}(:,1),TS.pre);
                                        if and(size(y,1)>10 , size(x,1)>5)
                                            [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                                            [ccg,T] = CCG(s,ids,'binSize',0.005,'duration',1,'smooth',0,'mode','ccg');
                                            ccg = ccg(:,1,2)./length(x); ccg = ccg./0.005;
                                            
                                            Pre.reward.members.vHPC = [Pre.reward.members.vHPC , ccg]; clear ccg x y m s ids groups
                                            I.reward.vHPC = [I.reward.vHPC ; countR];
                                        end
                                    end
                                end
                                
                                %% Post  reward members
                                for ii = 1 : numberD
                                    if (members(ii,i))
                                        member1 = clusters.dHPC(ii);
                                        y = spks_dHPC(spks_dHPC(:,1) == member1 , 2);
                                        y = Restrict(y,TS.post);
                                        x = Restrict(pks.reward{i}(:,1),TS.post);
                                        if and(size(y,1)>10 , size(x,1)>5)
                                            [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                                            [ccg,T] = CCG(s,ids,'binSize',0.005,'duration',1,'smooth',0,'mode','ccg');
                                            ccg = ccg(:,1,2)./length(x); ccg = ccg./0.005;
                                            ;
                                            Post.reward.members.dHPC = [Post.reward.members.dHPC , ccg]; clear ccg x y m s ids groups
                                            I.reward.dHPC = [I.reward.dHPC ; countR];
                                        end
                                    end
                                end
                                
                                for iii = 1 : numberV
                                    if (members(iii+numberD,i))
                                        member2 = clusters.vHPC(iii);
                                        y = spks_vHPC(spks_vHPC(:,1) == member2 , 2);
                                        y = Restrict(y,TS.post);
                                        x = Restrict(pks.reward{i}(:,1),TS.post);
                                        if and(size(y,1)>10 , size(x,1)>5)
                                            [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                                            [ccg,T] = CCG(s,ids,'binSize',0.005,'duration',1,'smooth',0,'mode','ccg');
                                            ccg = ccg(:,1,2)./length(x); ccg = ccg./0.005;
                                            
                                            Post.reward.members.vHPC = [Post.reward.members.vHPC , ccg]; clear ccg x y m s ids groups
                                            I.reward.vHPC = [I.reward.vHPC ; countR];
                                        end
                                    end
                                end
                                
                                clear RBR
                            end
                            clear pks TS
                        end
                    end
                    
                    
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
end

%% cleaning data
% Aversive
Pre.aversive.ccgM(:,sum(isnan(Post.aversive.ccgM))>0) = [];
% I.aversive(sum(isnan(Post.aversive.ccgM))>0,:) = [];
Post.aversive.ccgM(:,sum(isnan(Post.aversive.ccgM))>0) = [];

% Reward
Pre.reward.ccgM(:,sum(isnan(Post.reward.ccgM))>0) = [];
% I.reward(sum(isnan(Post.reward.ccgM))>0,:) = [];
Post.reward.ccgM(:,sum(isnan(Post.reward.ccgM))>0) = [];


end
