function [Pre Post T I durations] = Ripple_Assemblies_CCG(path)
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

I.aversive.dHPC = [];                I.reward = [];
I.aversive.vHPC = [];

Post.aversive.dHPC = [];         Post.aversive.vHPC = [];
durations.dHPC = [];             durations.vHPC = [];
durations.event = [];            durations.bursts = [];

Pre.aversive.ccgD = [] ;        Post.aversive.ccgD = [] ;
Pre.reward.ccgD = [] ;          Post.reward.ccgD = [] ;

Pre.aversive.ccgV = [] ;        Post.aversive.ccgV = [] ;
Pre.reward.ccgV = [] ;          Post.reward.ccgV = [] ;

Pre.aversive.ccgB = [] ;        Post.aversive.ccgB = [] ;
Pre.reward.ccgB = [] ;          Post.reward.ccgB = [] ;
countA = 1;

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
            bursts.coordinated.DV(bursts.coordinated.DV(:,2)-bursts.coordinated.DV(:,1)<0.05, : ) = [];
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
%                     patterns.aversive = patterns.aversive(:,cond.both);

                    if aversiveTS_run(1)<rewardTS_run(1)
                        TS.pre = NREM.baseline;
                        TS.post = NREM.aversive;
                        TS.run.run1 = movement.aversive;
                        TS.run.run2 = movement.reward;
                    else
                        TS.pre = NREM.reward;
                        TS.post = NREM.aversive;
                        TS.run.run1 = movement.aversive;
                        TS.run.run2 = movement.reward;
                    end
                    
                    pks.aversive = assemblies_peaks([bins' Spikes] , patterns.aversive , th);
                    
                    % Events for baseline calculation
                    x =[ripplesD(:,1)-0.05 ripplesD(:,3)+0.05 ; ripplesV(:,1)-0.05 ripplesV(:,3)+0.05];
%                     x =[ripple_event.all-0.05 ripple_event.all+0.05];
                    [c cc] = sort(x(:,1)); x = x(cc,:); x = ConsolidateIntervals(x);
                    
                    iterador = [];
                    % Iteration across assemblies
                    for i = 1 : size(pks.aversive,1)
                        % --- Pre ---
                        y = Restrict(pks.aversive{i}(:,1),TS.pre);
                        %                             x1 = Restrict(ripple_event.all(:,2),TS.pre);
                        x1 = Restrict(bursts.coordinated.DV(:,1),TS.pre);
                        if and(length(y)>5 , length(x1)>2)
                            [pInc1, pDec1 surp1] = RippleModulation_assemblies(ripple_event.all,y,TS.pre);
                            


                            [s,ids,groups] = CCGParameters(x1,ones(length(x1),1),y,ones(length(y),1)*2);
                            [ccg1,T] = CCG(s,ids,'binSize',0.02,'duration',4,'smooth',0,'mode','ccg');
                            cond1 = true;
%                             x1 = Restrict(x,TS.pre);
%                             m1 = InvertIntervals(x1,TS.pre(1,1) , TS.pre(end,2));
                            m1 = InvertIntervals(x,NREM.all(1,1),NREM.all(end,2));
                            m1 = length(Restrict(y,m1))/(sum(m1(:,2)-m1(:,1))); clear c cc
                        else
                            cond1 = false;
                        end
                        
                        % --- Post ---
                        y = Restrict(pks.aversive{i}(:,1),TS.post);
                        %                             x2 = Restrict(ripple_event.all(:,2),TS.post);
                        x2 = Restrict(bursts.coordinated.DV(:,1),TS.post);
                        
                        if and(length(y)>5 , length(x2)>2)
                            [pInc2 pDec2 surp2] = RippleModulation_assemblies(ripple_event.all,y,TS.post);
                            
                            [s,ids,groups] = CCGParameters(x2,ones(length(x2),1),y,ones(length(y),1)*2);
                            [ccg2,T] = CCG(s,ids,'binSize',0.02,'duration',4,'smooth',0,'mode','ccg');
                            cond2 = true;
%                             x2 = Restrict(x,TS.post);
%                             m2 = InvertIntervals(x2,TS.post(1,1) , TS.post(end,2));
                            m2 = InvertIntervals(x,NREM.all(1,1),NREM.all(end,2));
                            m2 = length(Restrict(y,m2))/(sum(m2(:,2)-m2(:,1)));                            
                        else
                            cond2 = false;
                        end
                        
                        if and(cond1,cond2)
%                             x1 = Restrict(ripple_event.all(:,2),TS.pre);
%                             x2 = Restrict(ripple_event.all(:,2),TS.post);
                            x1 = Restrict(bursts.coordinated.DV(:,1),TS.pre);
                            x2 = Restrict(bursts.coordinated.DV(:,1),TS.post);

                            ccg1 = (ccg1(:,1,2)./length(x1)); ccg1 = ccg1./0.02; 
                            ccg2 = (ccg2(:,1,2)./length(x2)); ccg2 = ccg2./0.02; 
                            Pre.aversive.ccg = [Pre.aversive.ccg , (ccg1)];%./m1)];%
                            Pre.aversive.iterator = [Pre.aversive.iterator ; pInc1 pDec1] ;
                            Post.aversive.ccg = [Post.aversive.ccg , (ccg2)];%./m2)];
                            Post.aversive.iterator = [Post.aversive.iterator ; pInc2 pDec2] ;
                            
                            % all
                            [average1 , time] = triggered_average_Ripples(x1,patterns.aversive(:,i),logical(1),[bins' Spikes],[-0.2 0.6],'all',[numberD numberV]);
                            [average2 , time] = triggered_average_Ripples(x2,patterns.aversive(:,i),logical(1),[bins' Spikes],[-0.2 0.6],'all',[numberD numberV]);
                            Pre.aversive.ccgM = [Pre.aversive.ccgM  , average1];
                            Post.aversive.ccgM = [Post.aversive.ccgM  , average2]; clear average1 average2
                            
                            % Just dHPC
                            [average1 , time] = triggered_average_Ripples(x1,patterns.aversive(:,i),logical(1),[bins' Spikes],[-0.2 0.6],'dHPC',[numberD numberV]);
                            [average2 , time] = triggered_average_Ripples(x2,patterns.aversive(:,i),logical(1),[bins' Spikes],[-0.2 0.6],'dHPC',[numberD numberV]);
                            Pre.aversive.ccgD = [Pre.aversive.ccgD  , average1];
                            Post.aversive.ccgD = [Post.aversive.ccgD  , average2]; clear average1 average2
                            
                            % Just vHPC
                            [average1 , time] = triggered_average_Ripples(x1,patterns.aversive(:,i),logical(1),[bins' Spikes],[-0.2 0.6],'vHPC',[numberD numberV]);
                            [average2 , time] = triggered_average_Ripples(x2,patterns.aversive(:,i),logical(1),[bins' Spikes],[-0.2 0.6],'vHPC',[numberD numberV]);
                            Pre.aversive.ccgV = [Pre.aversive.ccgV  , average1];
                            Post.aversive.ccgV = [Post.aversive.ccgV  , average2]; clear average1 average2        
                            
                            % Just Both
                            [average1 , time] = triggered_average_Ripples(x1,patterns.aversive(:,i),logical(1),[bins' Spikes],[-0.2 0.6],'both',[numberD numberV]);
                            [average2 , time] = triggered_average_Ripples(x2,patterns.aversive(:,i),logical(1),[bins' Spikes],[-0.2 0.6],'both',[numberD numberV]);
                            Pre.aversive.ccgB = [Pre.aversive.ccgB  , average1];
                            Post.aversive.ccgB = [Post.aversive.ccgB  , average2]; clear average1 average2                                
                            
                            [c cc] = min(abs(T-(-0.2)));
                            [c ccc] = min(abs(T-(0.2)));
                            c = nanmean((ccg2(cc:ccc,:)));
                            C = nanmean((ccg1(cc:ccc,:)));
                            
                            p = (c>C);
                            
                            if c>C
                                iterador = [iterador , true];
                            else
                                iterador = [iterador , false];
                            end
                        else
                            iterador = [iterador , nan];
                        end
                        clear y pInc1 pInc2 pDec1 pDec2 surp1 surp2
                        clear s ids groups ccg m1 m2 x1 x2 ccg1 ccg2
                        clear average1 average2
                    end
                    clear pks
                    
                    
                    x1 =[ripple_event.all-0.05 ripple_event.all+0.05];
                    [c cc] = sort(x(:,1));% x = x(cc,:); x = ConsolidateIntervals(x);
                    
                    members = Thresholded.aversive(:,cond.both);
                    
                    for i = 1 : sum(cond.both)
                        if iterador(i)
                            for ii = 1 : numberD
                                if not(members(ii,i))
                                    member1 = clusters.dHPC(ii);
                                    y = spks_dHPC(spks_dHPC(:,1) == member1 , 2);
                                    y = Restrict(y,TS.post);
                                    x = Restrict(bursts.coordinated.DV(:,1),TS.post);
                                    if and(size(y,1)>5 , size(x,1)>5)
                                        [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                                        [ccg,T] = CCG(s,ids,'binSize',0.005,'duration',1,'smooth',0,'mode','ccg');
                                        ccg = ccg(:,1,2)./length(x); ccg = ccg./0.005;
                                        
                                        m = InvertIntervals(x1,NREM.all(1,1),NREM.all(end,2));
                                        m = length(Restrict(y,m))/(sum(m(:,2)-m(:,1)));
                                        ccg = ccg./m;
                                        Post.aversive.dHPC = [Post.aversive.dHPC , ccg]; clear ccg x y m s ids groups
                                        I.aversive.dHPC = [I.aversive.dHPC ; countA];
                                    end
                                end
                            end
                            
                            for iii = 1 : numberV
                                if not(members(iii+numberD,i))
                                    member2 = clusters.vHPC(iii);
                                    y = spks_vHPC(spks_vHPC(:,1) == member2 , 2);
                                    y = Restrict(y,TS.post);
                                    x = Restrict(bursts.coordinated.DV(:,1),TS.post);
                                    if and(size(y,1)>5 , size(x,1)>5)
                                        [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                                        [ccg,T] = CCG(s,ids,'binSize',0.005,'duration',1,'smooth',0,'mode','ccg');
                                        ccg = ccg(:,1,2)./length(x); ccg = ccg./0.005;
                                        
                                        m = InvertIntervals(x1,NREM.all(1,1),NREM.all(end,2));
                                        m = length(Restrict(y,m))/(sum(m(:,2)-m(:,1)));
                                        ccg = ccg./m;
                                        Post.aversive.vHPC = [Post.aversive.vHPC , ccg]; clear ccg x y m s ids groups
                                        I.aversive.vHPC = [I.aversive.vHPC ; countA];
                                    end                                    
                                end
                            end
                        end
                        countA = countA+1;
                    end
                    
%                                         % ccg between members
%                     if aversiveTS_run(1)<rewardTS_run(1)
%                         TS.pre = [ripple_event.baseline(:,1) ripple_event.baseline(:,3)];
%                         TS.post = [ripple_event.aversive(:,1) ripple_event.aversive(:,3)];
%                     else
%                         TS.pre = [ripple_event.reward(:,1) ripple_event.reward(:,3)];
%                         TS.post = [ripple_event.aversive(:,1) ripple_event.aversive(:,3)];                      
%                     end
%                     
%                     load('react_coord_ripples_aversive.mat')
%                     % between members
%                     members = Thresholded.aversive(:,cond.both);
%                     for i = 1 : sum(cond.both)
%                         if iterador(i)
%                         for ii = 1 : numberD
%                                 if members(ii,i)
%                                     member1 = clusters.dHPC(ii);
%                                     x = spks_dHPC(spks_dHPC(:,1) == member1 , 2);
%                                     for iii = 1 : numberV
%                                         if members(iii+numberD,i)
%                                             member2 = clusters.vHPC(iii);
%                                             y = spks_vHPC(spks_vHPC(:,1) == member2 , 2);
%                                             
%                                             % --- Pre ---
%                                             xx = Restrict(x,TS.pre);
%                                             yy = Restrict(y,TS.pre);
% %                                             yy = y;
%                                             if and(length(xx)>5 , length(yy)>5)
%                                                 [s,ids,groups] = CCGParameters(xx,ones(length(xx),1),yy,ones(length(yy),1)*2);
%                                                 [ccg1,T] = CCG(s,ids,'binSize',0.003,'duration',dur,'smooth',sm,'mode','ccg');
%                                                 cond1 = true;
%                                                 
%                                                 surrogate = surrogate_ccg(x,y,0.1,dur,sm,TS.pre);
%                                                 surrogate = quantile(surrogate,0.90);
%                                                 
%                                                 ccg1 = zscore(ccg1(:,1,2)./sum(ccg1(:,1,2)));
%                                                 [~,up] = min(abs(T-0.1));    [~,down] = min(abs(T-(-0.1)));
%                                                 
%                                                 increase1 = max(ccg1(down:up))>surrogate;
%                                                 clear surrogate up down s ids groups
%                                             else
%                                                 cond1 = false;
%                                             end
%                                             
%                                             % --- Post ---
%                                             xx = Restrict(x,TS.post);
%                                             yy = Restrict(y,TS.post);
% %                                             yy = y;
%                                             if and(length(xx)>5 , length(yy)>5)
%                                                 [s,ids,groups] = CCGParameters(xx,ones(length(xx),1),yy,ones(length(yy),1)*2);
%                                                 [ccg2,T] = CCG(s,ids,'binSize',0.003,'duration',dur,'smooth',sm,'mode','ccg');
%                                                 cond2 = true;
%                                                 
%                                                 surrogate = surrogate_ccg(x,y,0.1,dur,sm,TS.post);
%                                                 surrogate = quantile(surrogate,0.90);
%                                                 
%                                                 ccg2 = (ccg2(:,1,2)./sum(ccg2(:,1,2)));
%                                                 [~,up] = min(abs(T-0.1));    [~,down] = min(abs(T-(-0.1)));
%                                                 
%                                                 increase2 = max(ccg2(down:up))>surrogate;
%                                                 clear surrogate up down s ids groups
%                                             else
%                                                 cond2 = false;
%                                             end
%                                             
%                                             if and(cond1,cond2)
%                                                 Pre.aversive.ccgM = [Pre.aversive.ccgM , (ccg1)];
%                                                 Post.aversive.ccgM = [Post.aversive.ccgM , (ccg2)];
%                                                 if and(increase2,not(increase1))
%                                                     I.aversive = [I.aversive ; true , RBA(i,1)];
%                                                 else
%                                                     I.aversive = [I.aversive ; false , RBA(i,1)];
%                                                 end
%                                             end
%                                             clear ccg1 ccg2 cond1 cond2
%                                         end
%                                     end
%                                 end
%                             end
%                         end
%                     end
                    
%                     % between non-members
%                     for i = 1 : sum(cond.both)
%                         for ii = 1 : numberD
%                             if (members(ii,i))
%                                 member1 = clusters.dHPC(ii);
%                                 x = spks_dHPC(spks_dHPC(:,1) == member1 , 2);
%                                 for iii = 1 : numberV
%                                     if not(members(iii+numberD,i))
%                                         member2 = clusters.vHPC(iii);
%                                         y = spks_vHPC(spks_vHPC(:,1) == member2 , 2);
%                                         
%                                         % --- Pre ---
%                                         xx = Restrict(x,TS.pre);
%                                         yy = Restrict(y,TS.pre);
%                                         if and(length(xx)>5 , length(yy)>5)
%                                             [s,ids,groups] = CCGParameters(xx,ones(length(xx),1),yy,ones(length(yy),1)*2);
%                                             [ccg1,T] = CCG(s,ids,'binSize',0.003,'duration',dur,'smooth',sm,'mode','ccg');
%                                             cond1 = true;
%                                         else
%                                             cond1 = false;
%                                         end
%                                         
%                                         % --- Post ---
%                                         xx = Restrict(x,TS.post);
%                                         yy = Restrict(y,TS.post);
%                                         if and(length(xx)>5 , length(yy)>5)
%                                             [s,ids,groups] = CCGParameters(xx,ones(length(xx),1),yy,ones(length(yy),1)*2);
%                                             [ccg2,T] = CCG(s,ids,'binSize',0.003,'duration',dur,'smooth',sm,'mode','ccg');
%                                             cond2 = true;
%                                         else
%                                             cond2 = false;
%                                         end
%                                         
%                                         if and(cond1,cond2)
%                                             Pre.aversive.ccgM1 = [Pre.aversive.ccgM1 , zscore(ccg1(:,1,2)./sum(ccg1(:,1,2)))];
%                                             Post.aversive.ccgM1 = [Post.aversive.ccgM1 , zscore(ccg2(:,1,2)./sum(ccg2(:,1,2)))];
% %                                             if conditional(i)
% %                                                 I.aversive = [I.aversive ; true , RBA(i,1)];
% %                                             else
% %                                                 I.aversive = [I.aversive ; false , RBA(i,1)];
% %                                             end
%                                         end
%                                         clear ccg1 ccg2 cond1 cond2
%                                         
%                                     end
%                                 end
%                             end
%                         end
%                     end
                                        
                    clear RBA
                    clear pks TS
                end
            end
            
            
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
                        TS.pre = NREM.aversive;
                        TS.post = NREM.reward;
                        TS.run.run1 = movement.reward;
                        TS.run.run2 = movement.aversive;
                    else
                        TS.pre = NREM.baseline;
                        TS.post = NREM.reward;
                        TS.run.run2 = movement.aversive;
                    end
%                     pks.reward = assemblies_peaks_joint(patterns.reward , cond.both , [bins' Spikes], Thresholded.reward, [numberD numberV]);
                    pks.reward = assemblies_peaks([bins' Spikes] , patterns.reward , th);
                    
                    % Events for baseline calculation
                    x =[ripplesD(:,1)-0.05 ripplesD(:,3)+0.05 ; ripplesV(:,1)-0.05 ripplesV(:,3)+0.05];
%                     x =[ripple_event.all-0.05 ripple_event.all+0.05];
                    [c cc] = sort(x(:,1)); x = x(cc,:); x = ConsolidateIntervals(x);
                    
                    % Iteration across assemblies
                    conditional = [];
                    for i = 1 : size(pks.reward,1)
                        % --- Pre ---
                        y = Restrict(pks.reward{i}(:,1),TS.pre);
                        %                             x1 = Restrict(ripple_event.all(:,2),TS.pre);
                        x1 = Restrict(bursts.coordinated.DV(:,1),TS.pre);
                        if and(length(y)>5 , length(x1)>2)
                            [pInc1 pDec1 surp1] = RippleModulation_assemblies(ripple_event.all,y,TS.pre);
                            
                            
                            
                            [s,ids,groups] = CCGParameters(x1,ones(length(x1),1),y,ones(length(y),1)*2);
                            [ccg1,T] = CCG(s,ids,'binSize',0.02,'duration',4,'smooth',0,'mode','ccg');
                            cond1 = true;
                            %                             x1 = Restrict(x,TS.pre);
                            %                             m1 = InvertIntervals(x1,TS.pre(1,1) , TS.pre(end,2));
                            m1 = InvertIntervals(x,NREM.all(1,1),NREM.all(end,2));
                            m1 = length(Restrict(y,m1))/(sum(m1(:,2)-m1(:,1))); clear c cc
                            cond1 = true;
                        else
                            cond1 = false;
                        end
                        
                        % --- Post ---
                        y = Restrict(pks.reward{i}(:,1),TS.post);
                        x2 = Restrict(bursts.coordinated.DV(:,1),TS.post);
%                             x2 = Restrict(ripple_event.all(:,2),TS.post);
                        if and(length(y)>5 , length(x2)>2)
                            [pInc2 pDec2 surp2] = RippleModulation_assemblies(ripple_event.all,y,TS.post);
  
                            [s,ids,groups] = CCGParameters(x2,ones(length(x2),1),y,ones(length(y),1)*2);
                            [ccg2,T] = CCG(s,ids,'binSize',0.02,'duration',4,'smooth',0,'mode','ccg');
                            cond2 = true;
%                             x2 = Restrict(x,TS.post);
%                             m2 = InvertIntervals(x2,TS.post(1,1) , TS.post(end,2));
                            m2 = InvertIntervals(x,NREM.all(1,1),NREM.all(end,2));
                            m2 = length(Restrict(y,m2))/(sum(m2(:,2)-m2(:,1))); 
                        else
                            cond2 = false;
                        end
                        
                        if and(cond1,cond2)                           
%                             x1 = Restrict(ripple_event.all(:,2),TS.pre);
%                             x2 = Restrict(ripple_event.all(:,2),TS.post);
                            x1 = Restrict(bursts.coordinated.DV(:,1),TS.pre);
                            x2 = Restrict(bursts.coordinated.DV(:,1),TS.post);

                            ccg1 = (ccg1(:,1,2)./length(x1)); ccg1 = ccg1./0.02; 
                            ccg2 = (ccg2(:,1,2)./length(x2)); ccg2 = ccg2./0.02; 
                            Pre.reward.ccg = [Pre.reward.ccg , (ccg1)];%./m1)];%
                            Pre.reward.iterator = [Pre.reward.iterator ; pInc1 pDec1] ;
                            Post.reward.ccg = [Post.reward.ccg , (ccg2)];%./m2)];
                            Post.reward.iterator = [Post.reward.iterator ; pInc2 pDec2] ;                            
                        
                            % all
                            [average1 , time] = triggered_average_Ripples(x1,patterns.reward(:,i),logical(1),[bins' Spikes],[-0.2 0.6],'all',[numberD numberV]);
                            [average2 , time] = triggered_average_Ripples(x2,patterns.reward(:,i),logical(1),[bins' Spikes],[-0.2 0.6],'all',[numberD numberV]);
                            Pre.reward.ccgM = [Pre.reward.ccgM  , average1];
                            Post.reward.ccgM = [Post.reward.ccgM  , average2]; clear average1 average2
                            
                            % Just dHPC
                            [average1 , time] = triggered_average_Ripples(x1,patterns.reward(:,i),logical(1),[bins' Spikes],[-0.2 0.6],'dHPC',[numberD numberV]);
                            [average2 , time] = triggered_average_Ripples(x2,patterns.reward(:,i),logical(1),[bins' Spikes],[-0.2 0.6],'dHPC',[numberD numberV]);
                            Pre.reward.ccgD = [Pre.reward.ccgD  , average1];
                            Post.reward.ccgD = [Post.reward.ccgD  , average2]; clear average1 average2
                            
                            % Just vHPC
                            [average1 , time] = triggered_average_Ripples(x1,patterns.reward(:,i),logical(1),[bins' Spikes],[-0.2 0.6],'vHPC',[numberD numberV]);
                            [average2 , time] = triggered_average_Ripples(x2,patterns.reward(:,i),logical(1),[bins' Spikes],[-0.2 0.6],'vHPC',[numberD numberV]);
                            Pre.reward.ccgV = [Pre.reward.ccgV  , average1];
                            Post.reward.ccgV = [Post.reward.ccgV  , average2]; clear average1 average2        
                            
                            % Just Both
                            [average1 , time] = triggered_average_Ripples(x1,patterns.reward(:,i),logical(1),[bins' Spikes],[-0.2 0.6],'both',[numberD numberV]);
                            [average2 , time] = triggered_average_Ripples(x2,patterns.reward(:,i),logical(1),[bins' Spikes],[-0.2 0.6],'both',[numberD numberV]);
                            Pre.reward.ccgB = [Pre.reward.ccgB  , average1];
                            Post.reward.ccgB = [Post.reward.ccgB  , average2]; clear average1 average2  
                            
                        end
                        clear y pInc1 pInc2 pDec1 pDec2 surp1 surp2
                        clear s ids groups ccg x1 x2 m1 m2 ccg ccg1 ccg2
                        clear average1 average2
                    end
                    clear pks TS
                    
                    
                    % ccg between members
                    if aversiveTS_run(1)<rewardTS_run(1)
                        TS.pre = [ripple_event.aversive(:,1) ripple_event.aversive(:,3)];
                        TS.post = [ripple_event.reward(:,1) ripple_event.reward(:,3)];
                    else
                        TS.pre = [ripple_event.baseline(:,1) ripple_event.baseline(:,3)];
                        TS.post = [ripple_event.reward(:,1) ripple_event.reward(:,3)];
                    end
                    
                    
                    load('react_coord_ripples_reward.mat')
                    members = Thresholded.reward(:,cond.both);
                    % between members
%                     for i = 1 : sum(cond.both)
%                         for ii = 1 : numberD
%                             if members(ii,i)
%                                 member1 = clusters.dHPC(ii);
%                                 x = spks_dHPC(spks_dHPC(:,1) == member1 , 2);
%                                 for iii = 1 : numberV
%                                     if members(iii+numberD,i)
%                                         member2 = clusters.vHPC(iii);
%                                         y = spks_vHPC(spks_vHPC(:,1) == member2 , 2);
%                                         
%                                         % --- Pre ---
%                                         xx = Restrict(x,TS.pre);
%                                         yy = Restrict(y,TS.pre);
%                                         if and(length(xx)>5 , length(yy)>5)
%                                             [s,ids,groups] = CCGParameters(xx,ones(length(xx),1),yy,ones(length(yy),1)*2);
%                                             [ccg1,T] = CCG(s,ids,'binSize',0.003,'duration',dur,'smooth',sm,'mode','ccg');
%                                             cond1 = true;
%                                         else
%                                             cond1 = false;
%                                         end
%                                         
%                                         % --- Post ---
%                                         xx = Restrict(x,TS.post);
%                                         yy = Restrict(y,TS.post);
%                                         if and(length(xx)>5 , length(yy)>5)
%                                             [s,ids,groups] = CCGParameters(xx,ones(length(xx),1),yy,ones(length(yy),1)*2);
%                                             [ccg2,T] = CCG(s,ids,'binSize',0.003,'duration',dur,'smooth',sm,'mode','ccg');
%                                             cond2 = true;
%                                         else
%                                             cond2 = false;
%                                         end
%                                         
%                                         if and(cond1,cond2)
%                                             Pre.reward.ccgM = [Pre.reward.ccgM , zscore(ccg1(:,1,2)./sum(ccg1(:,1,2)))];
%                                             Post.reward.ccgM = [Post.reward.ccgM , zscore(ccg2(:,1,2)./sum(ccg2(:,1,2)))];
%                                             if conditional(i)
%                                                 I.reward = [I.reward ; true , RBR(i,1)];
%                                             else
%                                                 I.reward = [I.reward ; false , RBR(i,1)];
%                                             end
%                                         end
%                                         clear ccg1 ccg2 cond1 cond2
%                                         
%                                     end
%                                 end
%                             end
%                         end
%                     end
%                     
%                     % between non-members
%                     for i = 1 : sum(cond.both)
%                         for ii = 1 : numberD
%                             if (members(ii,i))
%                                 member1 = clusters.dHPC(ii);
%                                 x = spks_dHPC(spks_dHPC(:,1) == member1 , 2);
%                                 for iii = 1 : numberV
%                                     if not(members(iii+numberD,i))
%                                         member2 = clusters.vHPC(iii);
%                                         y = spks_vHPC(spks_vHPC(:,1) == member2 , 2);
%                                         
%                                         % --- Pre ---
%                                         xx = Restrict(x,TS.pre);
%                                         yy = Restrict(y,TS.pre);
%                                         if and(length(xx)>5 , length(yy)>5)
%                                             [s,ids,groups] = CCGParameters(xx,ones(length(xx),1),yy,ones(length(yy),1)*2);
%                                             [ccg1,T] = CCG(s,ids,'binSize',0.003,'duration',dur,'smooth',sm,'mode','ccg');
%                                             cond1 = true;
%                                         else
%                                             cond1 = false;
%                                         end
%                                         
%                                         % --- Post ---
%                                         xx = Restrict(x,TS.post);
%                                         yy = Restrict(y,TS.post);
%                                         if and(length(xx)>5 , length(yy)>5)
%                                             [s,ids,groups] = CCGParameters(xx,ones(length(xx),1),yy,ones(length(yy),1)*2);
%                                             [ccg2,T] = CCG(s,ids,'binSize',0.003,'duration',dur,'smooth',sm,'mode','ccg');
%                                             cond2 = true;
%                                         else
%                                             cond2 = false;
%                                         end
%                                         
%                                         if and(cond1,cond2)
%                                             Pre.reward.ccgM1 = [Pre.reward.ccgM1 , zscore(ccg1(:,1,2)./sum(ccg1(:,1,2)))];
%                                             Post.reward.ccgM1 = [Post.reward.ccgM1 , zscore(ccg2(:,1,2)./sum(ccg2(:,1,2)))];
% %                                             if conditional(i)
% %                                                 I.reward = [I.reward ; true , RBR(i,1)];
% %                                             else
% %                                                 I.reward = [I.reward ; false , RBR(i,1)];
% %                                             end
%                                         end
%                                         clear ccg1 ccg2 cond1 cond2
%                                         
%                                     end
%                                 end
%                             end
%                         end
%                     end
                    clear RBR
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
% % % % 
% % %% Plot aversive assemblies activity
% % tmp1 = [];
% % tmp2 = [];
% % for j = 1:size(Post.aversive.ccg,2)
% %     
% %     [i ii] = min(abs(T-(-1)));
% %     [i iii] = min(abs(T-(1)));
% %     I = [Pre.aversive.ccg(1:ii,j) ; Pre.aversive.ccg(iii:end,j)];
% %     II = [Post.aversive.ccg(1:ii,j) ; Post.aversive.ccg(iii:end,j)];
% %     I = nanmean(I);
% %     II = nanmean(II);
% %     
% %     tmp1 = [tmp1 , Smooth(Pre.aversive.ccg(:,j)./I,0,'kernel','gaussian','type','l')];
% %     tmp2 = [tmp2 , Smooth(Post.aversive.ccg(:,j)./II,0,'kernel','gaussian','type','l')];
% %     
% %     clear I II i ii iii
% % end
% % 
% % [i ii] = min(abs(T-(-0.2)));
% % [i iii] = min(abs(T-(0.2)));
% % iA = nanmean((tmp2(ii:iii,:)));
% % IA = nanmean((tmp1(ii:iii,:)));
% % % % i=quantile(i,[0.33 0.66]);
% % % 
% % % p = not(i>=I);
% % % C1 = nanmean(tmp1(:,p)');
% % % S1 = nansem(tmp1(:,p)');
% % % 
% % % figure
% % % C2 = nanmean(tmp2(:,p)');
% % % S2 = nansem(tmp2(:,p)');
% % % subplot(121)
% % % plot(T,(C1),'k')
% % % hold on
% % % ciplot(C1-S1 , C1+S1 , T , 'k'),alpha 0.5
% % % plot(T,(C2),'r')
% % % ciplot(C2-S2 , C2+S2 , T , 'r'),alpha 0.5
% % % xlim([-0.6 0.6]),ylim([0 10])
% % % 
% % % subplot(122)
% % p = iA>=IA;
% % % C1 = nanmean(tmp1(:,p)');
% % % S1 = nansem(tmp1(:,p)');
% % % 
% % % % x = [I(p)',i(p)']
% % % % figure,subplot(121),boxplot(x),ylim([0 8])
% % % 
% % % C2 = nanmean(tmp2(:,p)');
% % % S2 = nansem(tmp2(:,p)');
% % % plot(T,(C1),'k')
% % % hold on
% % % ciplot(C1-S1 , C1+S1 , T , 'k'),alpha 0.5
% % % plot(T,(C2),'r')
% % % ciplot(C2-S2 , C2+S2 , T , 'r'),alpha 0.5
% % % xlim([-0.5 0.5]),ylim([0 6])
% % 
% % % Plot Colormaps
% % [i ii] = min(abs(T-(0)));
% % [i iii] = min(abs(T-(0.25)));
% % i = nanmean((tmp2(ii:iii,:)));
% % [ii iii] = sort(i,'descend');
% % 
% % figure,
% % subplot(211),imagesc(T,[1:size(tmp1,2)],tmp1(:,iii)'),colormap 'jet', caxis([0 10]),xlim([-0.6 0.6])
% % subplot(212),imagesc(T,[1:size(tmp2,2)],tmp2(:,iii)'),colormap 'jet', caxis([0 10]),xlim([-0.6 0.6])
% % 
% % figure
% % C1 = nanmean(tmp1');
% % S1 = nansem(tmp1');
% % 
% % C2 = nanmean(tmp2');
% % S2 = nansem(tmp2');
% % 
% % plot(T,(C1),'k')
% % hold on
% % ciplot(C1-S1 , C1+S1 , T , 'k'),alpha 0.5
% % plot(T,(C2),'r')
% % ciplot(C2-S2 , C2+S2 , T , 'r'),alpha 0.5
% % xlim([-0.6 0.6]),ylim([0 8])
% % 
% % %% Plot reward assemblies activity
% % tmp1 = [];
% % tmp2 = [];
% % for j = 1:size(Post.reward.ccg,2)
% %     
% %     [i ii] = min(abs(T-(-1)));
% %     [i iii] = min(abs(T-(1)));
% %     I = [Pre.reward.ccg(1:ii,j) ; Pre.reward.ccg(iii:end,j)];
% %     II = [Post.reward.ccg(1:ii,j) ; Post.reward.ccg(iii:end,j)];
% %     I = nanmean(I);
% %     II = nanmean(II);
% %     
% %     tmp1 = [tmp1 , Smooth(Pre.reward.ccg(:,j)./I,2,'kernel','gaussian','type','l')];
% %     tmp2 = [tmp2 , Smooth(Post.reward.ccg(:,j)./II,2,'kernel','gaussian','type','l')];
% %     
% %     clear I II i ii iii
% % end
% % 
% % [i ii] = min(abs(T-(-0.2)));
% % [i iii] = min(abs(T-(0.2)));
% % iR = nanmean((tmp2(ii:iii,:)));
% % IR = nanmean((tmp1(ii:iii,:)));
% % % i=quantile(i,[0.33 0.66]);
% % % 
% % % p = not(i>I);
% % % C1 = nanmean(tmp1(:,p)');
% % % S1 = nansem(tmp1(:,p)');
% % % 
% % % figure
% % % C2 = nanmean(tmp2(:,p)');
% % % S2 = nansem(tmp2(:,p)');
% % % subplot(211)
% % % plot(T,(C1),'k')
% % % hold on
% % % ciplot(C1-S1 , C1+S1 , T , 'k'),alpha 0.5
% % % plot(T,(C2),'r')
% % % ciplot(C2-S2 , C2+S2 , T , 'r'),alpha 0.5
% % % xlim([-0.6 0.6]),ylim([0 10])
% % % 
% % % subplot(212)
% % pp = iR>=IR;
% % % C1 = nanmean(tmp1(:,pp)');
% % % S1 = nansem(tmp1(:,pp)');
% % % 
% % % 
% % % x = [I(pp)',i(pp)']
% % % subplot(122),boxplot(x),ylim([0 8])
% % % 
% % % 
% % % C2 = nanmean(tmp2(:,pp)');
% % % S2 = nansem(tmp2(:,pp)');
% % % plot(T,(C1),'k')
% % % hold on
% % % ciplot(C1-S1 , C1+S1 , T , 'k'),alpha 0.5
% % % plot(T,(C2),'b')
% % % ciplot(C2-S2 , C2+S2 , T , 'b'),alpha 0.5
% % % xlim([-0.5 0.5]),ylim([0 6])
% % 
% % 
% % % Plot colormap
% % [i ii] = min(abs(T-(0)));
% % [i iii] = min(abs(T-(0.2)));
% % i = nanmean((tmp2(ii:iii,:)));
% % [ii iii] = sort(i,'descend');
% % 
% % figure,
% % subplot(211),imagesc(T,[1:size(tmp1,2)],tmp1(:,iii)'),colormap 'jet', caxis([0 10]),xlim([-0.6 0.6])
% % subplot(212),imagesc(T,[1:size(tmp2,2)],tmp2(:,iii)'),colormap 'jet', caxis([0 10]),xlim([-0.6 0.6])
% % 
% % figure
% % C1 = nanmean(tmp1');
% % S1 = nansem(tmp1');
% % 
% % C2 = nanmean(tmp2');
% % S2 = nansem(tmp2');
% % 
% % plot(T,(C1),'k')
% % hold on
% % ciplot(C1-S1 , C1+S1 , T , 'k'),alpha 0.5
% % plot(T,(C2),'r')
% % ciplot(C2-S2 , C2+S2 , T , 'r'),alpha 0.5
% % xlim([-0.6 0.6]),ylim([0 8])
% % 
% % 
% % %% Plot Percentages
% % p = (sum(p)/size(Pre.aversive.ccg,2))*100;
% % pp = (sum(pp)/size(Pre.reward.ccg,2))*100;
% % figure
% % subplot(121),pie([p , 100-p])
% % subplot(122),pie([pp , 100-pp])
% % 
% % % Plot Gain across conditions
% % x = [IA , iA , IR , iR];
% % x = [x;[ones(size(IA)) , ones(size(iA)) , ones(size(IR))*2 , ones(size(iR))*2]];
% % x = [x;[ones(size(IA)) , ones(size(iA))*2 , ones(size(IR)) , ones(size(iR))*2]];
% % 
% % grps = [ones(size(IA)) , ones(size(iA))*2 , ones(size(IR))*3 , ones(size(iR))*4];
% % figure,scatter(grps,x(1,:),"filled",'jitter','on', 'jitterAmount',0.1),ylim([0 8]) , xlim([0 5]),hold on
% % scatter([1 2 3 4],[nanmean(IA) , nanmean(iA) , nanmean(IR) , nanmean(iR)] ,'filled')
% % 
% % [~,~,stats] = anovan(x(1,:)' , {x(2,:)' , x(3,:)'} , 'model','interaction','varnames',{'Structure','phase'})
% % [results,~,~,gnames] = multcompare(stats,"Dimension",[1 2]);
% % 
% % 
% % % 
% % grps = [ones(size(iA)) , ones(size(iR))*2];
% % x = [iA , iR];
% % figure,scatter(grps,x(1,:),"filled",'jitter','on', 'jitterAmount',0.1),ylim([0 8]) , xlim([0 3]),hold on
% % scatter([1 2],[nanmean(iA) , nanmean(iR)] ,'filled')
% % [h p] = ranksum(iA , iR)
% % 
% %% Plot SU activity ACA FACU ---------------------------------------
% Post.aversive.dHPC(:,sum(isnan(Post.aversive.dHPC))>0) = [];
% Post.aversive.vHPC(:,sum(isnan(Post.aversive.vHPC))>0) = [];
% 
% tmp = [];
% for i = 1 : size(Post.aversive.dHPC,2)
%     tmp = [tmp , Smooth((Post.aversive.dHPC(:,i)),2,'kernel','gaussian','type','l')];
% end
% 
% [i ii] = max(tmp);
% [i ii] = sort(ii);
% figure,
% imagesc([-0.5 : 0.01 : 0.5] , [i:size(tmp,2)],tmp(:,ii)'),clim([0 100])
% xline(0)
% 
% tmp1 = [];
% for i = 1 : size(Post.aversive.vHPC,2)
%     tmp1 = [tmp1 , Smooth((Post.aversive.vHPC(:,i)),2,'kernel','gaussian','type','l')];
% end
% 
% [i ii] = max(Post.aversive.vHPC);
% [i ii] = sort(ii);
% figure,
% imagesc([-0.5 : 0.01 : 0.5] , [i:size(tmp1,2)],tmp1(:,ii)'),clim([0 100])
% xline(0)
% 
% 
% figure,plot([-0.5 : 0.005 : 0.5] ,nanmean(tmp'),'g'),hold on
% ciplot(nanmean(tmp')-nansem(tmp') , nanmean(tmp')+nansem(tmp'),[-0.5 : 0.005 : 0.5],'g'),alpha 0.2
% plot([-0.5 : 0.005 : 0.5] ,nanmean(tmp1'),'b')
% ciplot(nanmean(tmp1')-nansem(tmp1') , nanmean(tmp1')+nansem(tmp1'),[-0.5 : 0.005 : 0.5],'b'),alpha 0.2
% xline(0,'--'),xlim([-0.05 0.15])
% 
% 
% groups = unique(I.aversive.dHPC);
% Corr = [];
% for i = 1 : length(groups)
%     g = groups(i);
%     template1 = I.aversive.dHPC==g;
%     template1 = tmp(:,template1');
%     
%     template2 = I.aversive.vHPC==g;
%     template2 = tmp1(:,template2');
%     
%     for ii = 1 : size(template1,2)
%         for iii = 1 : size(template2,2)
%             A = template1(:,ii)-nanmean(template1(:,ii));
%             B = template2(:,iii)-nanmean(template2(:,iii));
%             P = nanmean((A*B'/(norm(A)*norm(B))));
% %             P = xcorr(B,A,100,'normalized');
% %             P = A.*B;
%             Corr = [Corr , Smooth(P,2,'kernel','gaussian','type','l')];
%         end
%     end
% end
% 
% figure,plot([-0.5 : 0.005 : 0.5],nanmean(Corr')),xline(0,'--'),xlim([-0.05 0.15])
% 
% % % 
% % % 
% % % x = [];
% % % for i = 1 : size(Post.aversive.dHPC,2)
% % %     tmp = Post.aversive.dHPC(:,i);
% % %     tmp = tmp-min(tmp);
% % %     tmp = tmp./max(tmp);
% % %     x = [x , Smooth(tmp,1)]; clear tmp
% % % end
% % % 
% % % y = [];
% % % for i = 1 : size(Post.aversive.vHPC,2)
% % %     tmp = Post.aversive.vHPC(:,i);
% % %     tmp = tmp-min(tmp);
% % %     tmp = tmp./max(tmp);
% % %     y = [y , Smooth(tmp,1)]; clear tmp
% % % end
% % % 
% % % [i ii] = max(x);
% % % [i ii] = sort(ii);
% % % figure,imagesc([-0.5 : 0.01 : 0.5],[1:size(x,2)],x(:,ii)'),xline(0,'--')
% % % 
% % % [i ii] = max(y);
% % % [i ii] = sort(ii);
% % % figure,imagesc([-0.5 : 0.01 : 0.5],[1:size(y,2)],y(:,ii)'),xline(0,'--')
% % % 
% % % figure,plot([-0.5 : 0.01 : 0.5] ,nanmean(x'),'g'),hold on
% % % plot([-0.5 : 0.01 : 0.5] ,nanmean(y'),'b')
% % % xline(0,'--')
% % % 
% % % figure,plot([-0.5 : 0.01 : 0.5] ,nanmean(x'),'g'),hold on
% % % plot([-0.5 : 0.01 : 0.5] ,nanmean(y'),'b')
% % % xline(0,'--')
% % % xlim([-0.1 0.1])
% % % 
% % % [~ , down] = min(abs([-0.1 : 0.01 : 0.1]-(-0.05)));
% % % [~ , up] = min(abs([-0.1 : 0.01 : 0.1]-(0.05)));
% % % 
% % % [~ , ii] = max(x(down:up,:));
% % % [~ , iii] = max(y(down:up,:));
% % % 
% % % time = [-0.05 : 0.01 : 0.05];
% % % lagsD = [];
% % % for i = 1 : size(ii,2)
% % %     lagsD = [lagsD , time(ii(i))];
% % % end
% % % 
% % % lagsV = [];
% % % for i = 1 : size(iii,2)
% % %     lagsV = [lagsV , time(iii(i))];
% % % end
% % % 
% % % 
% % % 
% % % x = Post.aversive.ccgM(:,logical(I.aversive(:,1)));
% % % [i ii] = max(x);
% % % [i ii] = sort(ii);
% % % 
% % % figure,imagesc([-0.3:0.003:0.3],[1:size(Post.aversive.ccgM,2)],x(:,ii)'),xline(0,'--')
% % % figure, plot([-0.3:0.003:0.3],nanmean(x')),xline(0,'--')
% % % 
% % % 
% % % x = Pre.aversive.ccgM(:,logical(I.aversive(:,1)));
% % % [i ii] = max(x);
% % % [i ii] = sort(ii);
% % % 
% % % figure,imagesc([-0.3:0.003:0.3],[1:size(Post.aversive.ccgM,2)],x(:,ii)'),xline(0,'--')
% % % figure, plot([-0.3:0.003:0.3],nanmean(x')),xline(0,'--')
% % % 
% % % 
% % % 
% % % Plot aversive assemblies activity
% % tmp1 = [];
% % tmp2 = [];
% % tmp3 = [];
% % tmp4 = [];
% % for j = 1:size(Post.aversive.ccg,2)
% %     t = [-4 : 0.025 : 4];
% %     t1 = Smooth(zscore(Pre.aversive.ccgM(:,j)),1,'kernel','gaussian','type','l');
% %     t2 = Smooth(zscore(Post.aversive.ccgM(:,j)),1,'kernel','gaussian','type','l');
% %     tmp1 = [tmp1 , t1];
% %     tmp2 = [tmp2 , t2];
% %     
% %     
% %     [i ii] = min(abs(t-(0)));
% %     [i iii] = min(abs(t-(0.2)));
% % 
% %    tmp3 = [tmp3 ; nanmean(t1(ii:iii))];
% %    tmp4 = [tmp4 ; nanmean(t2(ii:iii))];
% % 
% %    
% %     clear I II i ii iii t1 t2
% % end
% % 
% % figure
% % time = [-4 : 0.025 : 4];
% % m = nanmean(tmp1(:,p)');
% % s = nansem(tmp1(:,p)');
% % plot(time,m,'k'),hold on
% % ciplot(m-s,m+s,time,'k'),alpha 0.5
% % 
% % mm = nanmean(tmp2(:,p)');
% % ss = nansem(tmp2(:,p)');
% % plot(time,mm,'r'),hold on
% % ciplot(mm-ss,mm+ss,time,'r'),alpha 0.5
% % xlim([-0.6 0.6])
% % 
% % figure
% % time = [-4 : 0.025 : 4];
% % m = nanmean(tmp1');
% % s = nansem(tmp1');
% % plot(time,m,'k'),hold on
% % ciplot(m-s,m+s,time,'k'),alpha 0.5
% % 
% % 
% % 
% % mm = nanmean(tmp2');
% % ss = nansem(tmp2');
% % plot(time,mm,'r'),hold on
% % ciplot(mm-ss,mm+ss,time,'r'),alpha 0.5
% % xlim([-0.6 0.6])
% % 
% % figure,
% % grps = [ones(sum(p),1) ; ones(sum(p),1)*2];
% % x = [tmp3(p) ; tmp4(p)];
% % scatter(grps,x,"filled",'jitter','on', 'jitterAmount',0.1), hold on
% % scatter([1 2],[nanmean(tmp3(p)) , nanmean(tmp4(p))],"filled"), hold on
% % xlim([0 3])
% % 
% % [h hh] = ttest2(tmp3(p) , tmp4(p))
% % 
% % 
% % % Plot reward assemblies activity
% % tmp1 = [];
% % tmp2 = [];
% % tmp3 = [];
% % tmp4 = [];
% % for j = 1:size(Post.reward.ccg,2)
% %     t = [-4 : 0.025 : 4];
% %     t1 = Smooth((Pre.reward.ccgM(:,j)),1,'kernel','gaussian','type','l');
% %     t2 = Smooth((Post.reward.ccgM(:,j)),1,'kernel','gaussian','type','l');
% %     tmp1 = [tmp1 , t1];
% %     tmp2 = [tmp2 , t2];
% %     
% %     
% %     [i ii] = min(abs(t-(0.1)));
% %     [i iii] = min(abs(t-(0.5)));
% % 
% %    tmp3 = [tmp3 ; nanmean(t1(ii:iii))];
% %    tmp4 = [tmp4 ; nanmean(t2(ii:iii))];
% % 
% %    
% %     clear I II i ii iii t1 t2
% % end
% % 
% % figure
% % time = [-4 : 0.025 : 4];
% % m = nanmean(tmp1(:,pp)');
% % s = nansem(tmp1(:,pp)');
% % plot(time,m,'k'),hold on
% % ciplot(m-s,m+s,time,'k'),alpha 0.5
% % 
% % mm = nanmean(tmp2(:,pp)');
% % ss = nansem(tmp2(:,pp)');
% % plot(time,mm,'r'),hold on
% % ciplot(mm-ss,mm+ss,time,'r'),alpha 0.5
% % xlim([-0.6 0.6])
% % 
% % figure
% % time = [-4 : 0.025 : 4];
% % m = nanmean(tmp1');
% % s = nansem(tmp1');
% % plot(time,m,'k'),hold on
% % ciplot(m-s,m+s,time,'k'),alpha 0.5
% % 
% % mm = nanmean(tmp2');
% % ss = nansem(tmp2');
% % plot(time,mm,'r'),hold on
% % ciplot(mm-ss,mm+ss,time,'r'),alpha 0.5
% % xlim([-0.6 0.6])
% % ylim([0 0.7])
% % 
% % 
% % 
% % figure,
% % grps = [ones(sum(pp),1) ; ones(sum(pp),1)*2];
% % x = [tmp3(pp) ; tmp4(pp)];
% % scatter(grps,x,"filled",'jitter','on', 'jitterAmount',0.1), hold on
% % scatter([1 2],[nanmean(tmp3(pp)) , nanmean(tmp4(pp))],"filled"), hold on
% % xlim([0 3])
% % [h hh] = ttest2(tmp3(pp) , tmp4(pp))
% % 
% % 
% % 
% % 
% % %% --- Plot components separated ---
% % % Plot dHPC assemblies activity
% % tmp1 = [];
% % tmp2 = [];
% % for j = 1:size(Post.aversive.ccgD,2)
% %     t = [-0.2 : 0.025 : 0.6];
% %     t1 = Smooth((Pre.aversive.ccgD(:,j)),1,'kernel','gaussian','type','l');
% %     t2 = Smooth((Post.aversive.ccgD(:,j)),1,'kernel','gaussian','type','l');
% %     
% % %     t1 = t1 - min(t1); t1 = t1./max(t1);
% % %     t2 = t2 - min(t2); t2 = t2./max(t2);
% %     
% %     tmp1 = [tmp1 , t1];
% %     tmp2 = [tmp2 , t2];
% %     clear I II i ii iii t1 t2
% % end
% % 
% % figure
% % time = [-0.2 : 0.025 : 0.6];
% % m = nanmean(tmp1');
% % s = nansem(tmp1');
% % plot(time,m,'k'),hold on
% % ciplot(m-s,m+s,time,'k'),alpha 0.5
% % 
% % mm = nanmean(tmp2');
% % ss = nansem(tmp2');
% % plot(time,mm,'r'),hold on
% % ciplot(mm-ss,mm+ss,time,'r'),alpha 0.5
% % xlim([-0.2 0.6])
% % 
% % 
% % % Plot vHPC assemblies activity
% % tmp1 = [];
% % tmp2 = [];
% % for j = 1:size(Post.aversive.ccgV,2)
% %     t = [-0.2 : 0.025 : 0.6];
% %     t1 = Smooth((Pre.aversive.ccgV(:,j)),1,'kernel','gaussian','type','l');
% %     t2 = Smooth((Post.aversive.ccgV(:,j)),1,'kernel','gaussian','type','l');
% %     
% % %     t1 = t1 - min(t1); t1 = t1./max(t1);
% % %     t2 = t2 - min(t2); t2 = t2./max(t2);
% %     
% %     tmp1 = [tmp1 , t1];
% %     tmp2 = [tmp2 , t2];
% %     clear I II i ii iii t1 t2
% % end
% % 
% % figure
% % time = [-0.2 : 0.025 : 0.6];
% % m = nanmean(tmp1');
% % s = nansem(tmp1');
% % plot(time,m,'k'),hold on
% % ciplot(m-s,m+s,time,'k'),alpha 0.5
% % 
% % mm = nanmean(tmp2');
% % ss = nansem(tmp2');
% % plot(time,mm,'r'),hold on
% % ciplot(mm-ss,mm+ss,time,'r'),alpha 0.5
% % xlim([-0.2 0.6])
% % 
% % 
% % % Plot both assemblies activity
% % tmp1 = [];
% % tmp2 = [];
% % for j = 1:size(Post.aversive.ccgB,2)
% %     t = [-4 : 0.025 : 4];
% %     t1 = Smooth((Pre.aversive.ccgB(:,j)),1,'kernel','gaussian','type','l');
% %     t2 = Smooth((Post.aversive.ccgB(:,j)),1,'kernel','gaussian','type','l');
% %     
% % %     t1 = t1 - min(t1); t1 = t1./max(t1);
% % %     t2 = t2 - min(t2); t2 = t2./max(t2);
% %     
% %     tmp1 = [tmp1 , t1];
% %     tmp2 = [tmp2 , t2];
% %     clear I II i ii iii t1 t2
% % end
% % 
% % figure
% % time = [-0.2 : 0.025 : 0.6];
% % m = nanmean(tmp1');
% % s = nansem(tmp1');
% % plot(time,m,'k'),hold on
% % ciplot(m-s,m+s,time,'k'),alpha 0.5
% % 
% % mm = nanmean(tmp2');
% % ss = nansem(tmp2');
% % plot(time,mm,'r'),hold on
% % ciplot(mm-ss,mm+ss,time,'r'),alpha 0.5
% % xlim([-0.2 0.6])
% % 
% % % % 
% % % % 