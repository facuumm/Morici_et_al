function [Pre Post T I] = Ripple_Assemblies_CCG(path)
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
% 3 SD funciona
% 5 SD funciona
sm = 1;
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

I.aversive = [];                I.reward = [];

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
                    
                    cooridnated_event = [cooridnated_event ; min([r(1) , z(indice,1)]) , min(min(z(indice,2),r(2)))+abs(z(indice,2)-r(2))/2 , max([r(3) , z(indice,3)])];
                    
                    clear tmp2 tmp1 p indice z
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
            bursts.baseline = Restrict(bursts.coordinated,NREM.baseline);
            bursts.aversive = Restrict(bursts.coordinated,NREM.aversive);
            bursts.reward = Restrict(bursts.coordinated,NREM.reward);
            
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
                    pks.aversive = assemblies_peaks([bins' Spikes] , patterns.aversive(:,cond.both) , th);
                    
                    % Iteration across assemblies
                    conditional = [];
                    for i = 1 : size(pks.aversive,1)
                        % --- Pre ---
                        y = Restrict(pks.aversive{i}(:,1),TS.pre);
                        if length(y)>5
                            [pInc1 pDec1 surp1] = RippleModulation_assemblies(ripple_event.all,y,TS.pre);
                            
                            x1 = Restrict(ripple_event.all(:,2),TS.pre);
                            [s,ids,groups] = CCGParameters(x1,ones(length(x1),1),y,ones(length(y),1)*2);
                            [ccg1,T] = CCG(s,ids,'binSize',0.025,'duration',0.4,'smooth',1,'mode','ccg');
                            cond1 = true;
                        else
                            cond1 = false;
                        end
                        
                        % --- Post ---
                        y = Restrict(pks.aversive{i}(:,1),TS.post);
                        if length(y)>5
                            [pInc2 pDec2 surp2] = RippleModulation_assemblies(ripple_event.all,y,TS.post);
                            
                            x2 = Restrict(ripple_event.all(:,2),TS.post);
                            [s,ids,groups] = CCGParameters(x2,ones(length(x2),1),y,ones(length(y),1)*2);
                            [ccg2,T] = CCG(s,ids,'binSize',0.025,'duration',0.4,'smooth',1,'mode','ccg');
                            cond2 = true;
                        else
                            cond2 = false;
                        end
                        
                        if and(cond1,cond2)
                            Pre.aversive.ccg = [Pre.aversive.ccg , zscore(ccg1(:,1,2)./0.025)./length(x1)] ;
                            Pre.aversive.iterator = [Pre.aversive.iterator ; pInc1 pDec1] ;
                            Post.aversive.ccg = [Post.aversive.ccg , zscore(ccg2(:,1,2)./0.025)./length(x2)] ;
                            Post.aversive.iterator = [Post.aversive.iterator ; pInc2 pDec2] ;
                            conditional = [conditional or(and(pInc1<0.05 , pInc2<0.05) , or(pInc1<0.05 , pInc2<0.05))];
                        end
                        clear y x pInc1 pInc2 pDec1 pDec2 surp1 surp2
                        clear s ids groups ccg
                    end
                    clear pks TS
                    
                    % ccg between members
                    if aversiveTS_run(1)<rewardTS_run(1)
                        TS.pre = [ripple_event.baseline(:,1) ripple_event.baseline(:,3)];
                        TS.post = [ripple_event.aversive(:,1) ripple_event.aversive(:,3)];
                        
%                         TS.pre = bursts.baseline;
%                         TS.post = bursts.aversive;                        
                    else
                        TS.pre = [ripple_event.reward(:,1) ripple_event.reward(:,3)];
                        TS.post = [ripple_event.aversive(:,1) ripple_event.aversive(:,3)];
                        
%                         TS.pre = bursts.reward;
%                         TS.post = bursts.aversive;                        
                    end
                    
                    load('react_coord_ripples_aversive.mat')
                    % between members
                    members = Thresholded.aversive(:,cond.both);
                    for i = 1 : sum(cond.both)
                        for ii = 1 : numberD
                            if members(ii,i)
%                                 break
                                member1 = clusters.dHPC(ii);
                                x = spks_dHPC(spks_dHPC(:,1) == member1 , 2);
                                for iii = 1 : numberV
                                    if members(iii+numberD,i)
%                                         break
                                        member2 = clusters.vHPC(iii);
                                        y = spks_vHPC(spks_vHPC(:,1) == member2 , 2);
                                        
                                        % --- Pre ---
                                        xx = Restrict(x,TS.pre);
                                        yy = Restrict(y,TS.pre);
                                        if and(length(xx)>5 , length(yy)>5)
                                            [s,ids,groups] = CCGParameters(xx,ones(length(xx),1),yy,ones(length(yy),1)*2);
                                            [ccg1,T] = CCG(s,ids,'binSize',0.003,'duration',dur,'smooth',sm,'mode','ccg');
                                            cond1 = true;
                                            
                                            surrogate = surrogate_ccg(xx,yy);
                                            
                                        else
                                            cond1 = false;
                                        end
                                        
                                        % --- Post ---
                                        xx = Restrict(x,TS.post);
                                        yy = Restrict(y,TS.post);
                                        if and(length(xx)>5 , length(yy)>5)
                                            [s,ids,groups] = CCGParameters(xx,ones(length(xx),1),yy,ones(length(yy),1)*2);
                                            [ccg2,T] = CCG(s,ids,'binSize',0.003,'duration',dur,'smooth',sm,'mode','ccg');
                                            cond2 = true;
                                        else
                                            cond2 = false;
                                        end
                                        
                                        if and(cond1,cond2)
                                            Pre.aversive.ccgM = [Pre.aversive.ccgM , zscore(ccg1(:,1,2)./sum(ccg1(:,1,2)))];
                                            Post.aversive.ccgM = [Post.aversive.ccgM , zscore(ccg2(:,1,2)./sum(ccg2(:,1,2)))];
                                            if conditional(i)
                                                I.aversive = [I.aversive ; true , RBA(i,1)];
                                            else
                                                I.aversive = [I.aversive ; false , RBA(i,1)];
                                            end
                                        end
                                        clear ccg1 ccg2 cond1 cond2
                                        
                                    end
                                end
                            end
                        end
                    end
                    
                    % between non-members
                    for i = 1 : sum(cond.both)
                        for ii = 1 : numberD
                            if (members(ii,i))
                                member1 = clusters.dHPC(ii);
                                x = spks_dHPC(spks_dHPC(:,1) == member1 , 2);
                                for iii = 1 : numberV
                                    if not(members(iii+numberD,i))
                                        member2 = clusters.vHPC(iii);
                                        y = spks_vHPC(spks_vHPC(:,1) == member2 , 2);
                                        
                                        % --- Pre ---
                                        xx = Restrict(x,TS.pre);
                                        yy = Restrict(y,TS.pre);
                                        if and(length(xx)>5 , length(yy)>5)
                                            [s,ids,groups] = CCGParameters(xx,ones(length(xx),1),yy,ones(length(yy),1)*2);
                                            [ccg1,T] = CCG(s,ids,'binSize',0.003,'duration',dur,'smooth',sm,'mode','ccg');
                                            cond1 = true;
                                        else
                                            cond1 = false;
                                        end
                                        
                                        % --- Post ---
                                        xx = Restrict(x,TS.post);
                                        yy = Restrict(y,TS.post);
                                        if and(length(xx)>5 , length(yy)>5)
                                            [s,ids,groups] = CCGParameters(xx,ones(length(xx),1),yy,ones(length(yy),1)*2);
                                            [ccg2,T] = CCG(s,ids,'binSize',0.003,'duration',dur,'smooth',sm,'mode','ccg');
                                            cond2 = true;
                                        else
                                            cond2 = false;
                                        end
                                        
                                        if and(cond1,cond2)
                                            Pre.aversive.ccgM1 = [Pre.aversive.ccgM1 , zscore(ccg1(:,1,2)./sum(ccg1(:,1,2)))];
                                            Post.aversive.ccgM1 = [Post.aversive.ccgM1 , zscore(ccg2(:,1,2)./sum(ccg2(:,1,2)))];
%                                             if conditional(i)
%                                                 I.aversive = [I.aversive ; true , RBA(i,1)];
%                                             else
%                                                 I.aversive = [I.aversive ; false , RBA(i,1)];
%                                             end
                                        end
                                        clear ccg1 ccg2 cond1 cond2
                                        
                                    end
                                end
                            end
                        end
                    end
                                        
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
                    %                     patterns.reward = patterns.reward(:,cond.both);
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
                    pks.reward = assemblies_peaks([bins' Spikes] , patterns.reward(:,cond.both) , th);
                    
                    
                    % Iteration across assemblies
                    conditional = [];
                    for i = 1 : size(pks.reward,1)
                        % --- Pre ---
                        y = Restrict(pks.reward{i}(:,1),TS.pre);
                        if length(y)>5
                            [pInc1 pDec1 surp1] = RippleModulation_assemblies(ripple_event.all,y,TS.pre);
                            
                            x1 = Restrict(ripple_event.all(:,2),TS.pre);
                            [s,ids,groups] = CCGParameters(x1,ones(length(x1),1),y,ones(length(y),1)*2);
                            [ccg1,T] = CCG(s,ids,'binSize',0.025,'duration',0.4,'smooth',1,'mode','ccg');
                            cond1 = true;
                        else
                            cond1 = false;
                        end
                        
                        % --- Post ---
                        y = Restrict(pks.reward{i}(:,1),TS.post);
                        if length(y)>5
                            [pInc2 pDec2 surp2] = RippleModulation_assemblies(ripple_event.all,y,TS.post);
                            
                            x2 = Restrict(ripple_event.all(:,2),TS.post);
                            [s,ids,groups] = CCGParameters(x2,ones(length(x2),1),y,ones(length(y),1)*2);
                            [ccg2,T] = CCG(s,ids,'binSize',0.025,'duration',0.4,'smooth',1,'mode','ccg');
                            cond2 = true;
                        else
                            cond2 = false;
                        end
                        
                        if and(cond1,cond2)
                            Pre.reward.ccg = [Pre.reward.ccg , zscore(ccg1(:,1,2)./0.025)./length(x1)] ;
                            Pre.reward.iterator = [Pre.reward.iterator ; pInc1 pDec1] ;
                            Post.reward.ccg = [Post.reward.ccg , zscore(ccg2(:,1,2)./0.025)./length(x2)] ;
                            Post.reward.iterator = [Post.reward.iterator ; pInc2 pDec2] ;
                            conditional = [conditional or(and(pInc1<0.05 , pInc2<0.05) , or(pInc1<0.05 , pInc2<0.05))];
                        end
                        clear y x pInc1 pInc2 pDec1 pDec2 surp1 surp2
                        clear s ids groups ccg
                    end
                    clear pks TS
                    
                    
                    % ccg between members
                    if aversiveTS_run(1)<rewardTS_run(1)
                        TS.pre = [ripple_event.aversive(:,1) ripple_event.aversive(:,3)];
                        TS.post = [ripple_event.reward(:,1) ripple_event.reward(:,3)];

%                         TS.pre = bursts.aversive;
%                         TS.post = bursts.reward;
                    else
                        TS.pre = [ripple_event.baseline(:,1) ripple_event.baseline(:,3)];
                        TS.post = [ripple_event.reward(:,1) ripple_event.reward(:,3)];

%                         TS.pre = bursts.baseline;
%                         TS.post = bursts.reward;
                    end
                    
                    
                    load('react_coord_ripples_reward.mat')
                    members = Thresholded.reward(:,cond.both);
                    % between members
                    for i = 1 : sum(cond.both)
                        for ii = 1 : numberD
                            if members(ii,i)
                                member1 = clusters.dHPC(ii);
                                x = spks_dHPC(spks_dHPC(:,1) == member1 , 2);
                                for iii = 1 : numberV
                                    if members(iii+numberD,i)
                                        member2 = clusters.vHPC(iii);
                                        y = spks_vHPC(spks_vHPC(:,1) == member2 , 2);
                                        
                                        % --- Pre ---
                                        xx = Restrict(x,TS.pre);
                                        yy = Restrict(y,TS.pre);
                                        if and(length(xx)>5 , length(yy)>5)
                                            [s,ids,groups] = CCGParameters(xx,ones(length(xx),1),yy,ones(length(yy),1)*2);
                                            [ccg1,T] = CCG(s,ids,'binSize',0.003,'duration',dur,'smooth',sm,'mode','ccg');
                                            cond1 = true;
                                        else
                                            cond1 = false;
                                        end
                                        
                                        % --- Post ---
                                        xx = Restrict(x,TS.post);
                                        yy = Restrict(y,TS.post);
                                        if and(length(xx)>5 , length(yy)>5)
                                            [s,ids,groups] = CCGParameters(xx,ones(length(xx),1),yy,ones(length(yy),1)*2);
                                            [ccg2,T] = CCG(s,ids,'binSize',0.003,'duration',dur,'smooth',sm,'mode','ccg');
                                            cond2 = true;
                                        else
                                            cond2 = false;
                                        end
                                        
                                        if and(cond1,cond2)
                                            Pre.reward.ccgM = [Pre.reward.ccgM , zscore(ccg1(:,1,2)./sum(ccg1(:,1,2)))];
                                            Post.reward.ccgM = [Post.reward.ccgM , zscore(ccg2(:,1,2)./sum(ccg2(:,1,2)))];
                                            if conditional(i)
                                                I.reward = [I.reward ; true , RBR(i,1)];
                                            else
                                                I.reward = [I.reward ; false , RBR(i,1)];
                                            end
                                        end
                                        clear ccg1 ccg2 cond1 cond2
                                        
                                    end
                                end
                            end
                        end
                    end
                    
                    % between non-members
                    for i = 1 : sum(cond.both)
                        for ii = 1 : numberD
                            if (members(ii,i))
                                member1 = clusters.dHPC(ii);
                                x = spks_dHPC(spks_dHPC(:,1) == member1 , 2);
                                for iii = 1 : numberV
                                    if not(members(iii+numberD,i))
                                        member2 = clusters.vHPC(iii);
                                        y = spks_vHPC(spks_vHPC(:,1) == member2 , 2);
                                        
                                        % --- Pre ---
                                        xx = Restrict(x,TS.pre);
                                        yy = Restrict(y,TS.pre);
                                        if and(length(xx)>5 , length(yy)>5)
                                            [s,ids,groups] = CCGParameters(xx,ones(length(xx),1),yy,ones(length(yy),1)*2);
                                            [ccg1,T] = CCG(s,ids,'binSize',0.003,'duration',dur,'smooth',sm,'mode','ccg');
                                            cond1 = true;
                                        else
                                            cond1 = false;
                                        end
                                        
                                        % --- Post ---
                                        xx = Restrict(x,TS.post);
                                        yy = Restrict(y,TS.post);
                                        if and(length(xx)>5 , length(yy)>5)
                                            [s,ids,groups] = CCGParameters(xx,ones(length(xx),1),yy,ones(length(yy),1)*2);
                                            [ccg2,T] = CCG(s,ids,'binSize',0.003,'duration',dur,'smooth',sm,'mode','ccg');
                                            cond2 = true;
                                        else
                                            cond2 = false;
                                        end
                                        
                                        if and(cond1,cond2)
                                            Pre.reward.ccgM1 = [Pre.reward.ccgM1 , zscore(ccg1(:,1,2)./sum(ccg1(:,1,2)))];
                                            Post.reward.ccgM1 = [Post.reward.ccgM1 , zscore(ccg2(:,1,2)./sum(ccg2(:,1,2)))];
%                                             if conditional(i)
%                                                 I.reward = [I.reward ; true , RBR(i,1)];
%                                             else
%                                                 I.reward = [I.reward ; false , RBR(i,1)];
%                                             end
                                        end
                                        clear ccg1 ccg2 cond1 cond2
                                        
                                    end
                                end
                            end
                        end
                    end
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
I.aversive(sum(isnan(Post.aversive.ccgM))>0,:) = [];
Post.aversive.ccgM(:,sum(isnan(Post.aversive.ccgM))>0) = [];
% 
% tmp = [];
% for i = 1 : size(Post.aversive.ccgM,2)
%     mm = Post.aversive.ccgM1(:,randperm(size(Post.aversive.ccgM1,2)));
%     mm(:,sum(isnan(mm))>0) = [];
%     mm = mm(:,1:size(Post.aversive.ccgM,2));
%     
%     tmp = [tmp , mm(:,23)];
% end
% Post.aversive.ccgM1 = tmp; clear tmp
% 
% tmp = [];
% for i = 1 : size(Pre.aversive.ccgM,2)
%     mm = Pre.aversive.ccgM1(:,randperm(size(Pre.aversive.ccgM1,2)));
%     mm(:,sum(isnan(mm))>0) = [];
%     mm = mm(:,1:size(Pre.aversive.ccgM,2));
%     
%     tmp = [tmp , mm(:,23)];
% end
% Pre.aversive.ccgM1 = tmp; clear tmp
% 
% 
% Reward
Pre.reward.ccgM(:,sum(isnan(Post.reward.ccgM))>0) = [];
I.reward(sum(isnan(Post.reward.ccgM))>0,:) = [];
Post.reward.ccgM(:,sum(isnan(Post.reward.ccgM))>0) = [];
% 
% tmp = [];
% for i = 1 : size(Post.reward.ccgM,2)
%     mm = Post.reward.ccgM1(:,randperm(size(Post.reward.ccgM1,2)));
%     mm(:,sum(isnan(mm))>0) = [];
%     mm = mm(:,1:size(Post.reward.ccgM,2));
%     
%     tmp = [tmp , mm(:,23)];
% end
% Post.reward.ccgM1 = tmp; clear tmp
% 
% tmp = [];
% for i = 1 : size(Pre.reward.ccgM,2)
%     mm = Pre.reward.ccgM1(:,randperm(size(Pre.reward.ccgM1,2)));
%     mm(:,sum(isnan(mm))>0) = [];
%     mm = mm(:,1:size(Pre.reward.ccgM,2));
%     
%     tmp = [tmp , mm(:,23)];
% end
% Pre.reward.ccgM1 = tmp; clear tmp


end
% 
% 
% % Cleaning the data
% 
% m = Post.aversive.ccgM;%(:,logical(I.aversive));
% mm = Pre.aversive.ccgM;%(:,logical(I.aversive));
% 
% % mm = Post.aversive.ccgM1(:,randperm(size(Post.aversive.ccgM1,2)));
% % mm(:,sum(isnan(mm))>0) = [];
% % mm = mm(:,1:size(m,2));
% 
% 
% 
% % Separation deppending on their Reactivation Strength
% React = unique(I.aversive(:,2));
% q = quantile(React,[0.3333 0.6666]);
% 
% % R1 = React<=q(1);
% % R2 = and(React>q(1),React<=q(2));
% % R3 = React>q(2);
% 
% R1 = not(isnan(React));%<=nanmean(React);
% R2 = React>=nanmean(React);
% 
% %First Q
% R1 = React(R1);
% R1 = ismember(I.aversive(:,2),R1);
% m1 = m(:,R1);
% mm1 = mm(:,R1);
% 
% figure,
% subplot(131)
% me1 = nanmean(mm1');
% se1 = nansem(mm1');
% plot(T,me1,'k'),hold on
% ciplot(me1-se1 , me1+se1, T,'k'),alpha 0.5,hold on
% 
% 
% me2 = nanmean(m1');
% se2 = nansem(m1');
% plot(T,me2,'r')
% ciplot(me2-se2 , me2+se2, T,'r'),alpha 0.5
% ylim([-0.55 2.2])
% % xlim([-0.2 0.2])
% xline(0,'--')
% 
% 
% %Second Q
% R2 = React(R2);
% R2 = ismember(I.aversive(:,2),R2);
% m2 = m(:,R2);
% mm2 = mm(:,R2);
% 
% subplot(132)
% me1 = nanmean(mm2');
% se1 = nansem(mm2');
% plot(T,me1,'k'),hold on
% ciplot(me1-se1 , me1+se1, T,'k'),alpha 0.5,hold on
% 
% me2 = nanmean(m2');
% se2 = nansem(m2');
% plot(T,me2,'r')
% ciplot(me2-se2 , me2+se2, T,'r'),alpha 0.5
% ylim([-0.55 2.2])
% % xlim([-0.2 0.2])
% xline(0,'--')
% 
% 
% %Thrid Q
% R3 = React(R3);
% R3 = ismember(I.aversive(:,2),R3);
% m2 = m(:,R3);
% mm2 = mm(:,R3);
% 
% subplot(133)
% me1 = nanmean(mm2');
% se1 = nansem(mm2');
% plot(T,me1,'k'),hold on
% ciplot(me1-se1 , me1+se1, T,'k'),alpha 0.5,hold on
% 
% me2 = nanmean(m2');
% se2 = nansem(m2');
% plot(T,me2,'r')
% ciplot(me2-se2 , me2+se2, T,'r'),alpha 0.5
% % ylim([0.2 0.6])
% xlim([-0.2 0.2])
% xline(0,'--')
% 
% 
% 
% 
% 
% % 
% % All
% 
% m = Post.aversive.ccgM;%(:,logical(I.aversive));
% 
% mm = Post.aversive.ccgM1(:,randperm(size(Post.aversive.ccgM1,2)));
% mm(:,sum(isnan(mm))>0) = [];
% mm = mm(:,1:size(m,2));
% 
% 
% m1 = m;
% mm1 = mm;
% 
% figure,
% subplot(121)
% me1 = nanmean(mm1');
% se1 = nansem(mm1');
% plot(T,me1,'k'),hold on
% ciplot(me1-se1 , me1+se1, T,'k'),alpha 0.5,hold on
% 
% 
% me2 = nanmean(m1');
% se2 = nansem(m1');
% plot(T,me2,'r')
% ciplot(me2-se2 , me2+se2, T,'r'),alpha 0.5
% xlim([-0.3 0.3])
% % ylim([-0.7 2.2])
% xline(0,'--')
% 
% 
% 
% % Cleaning the data
% m = Pre.aversive.ccgM;%(:,logical(I.aversive));
% 
% mm = Pre.aversive.ccgM1(:,randperm(size(Pre.aversive.ccgM1,2)));
% mm(:,sum(isnan(mm))>0) = [];
% mm = mm(:,1:size(m,2));
% 
% m1 = m;
% mm1 = mm;
% 
% subplot(122)
% me1 = nanmean(mm1');
% se1 = nansem(mm1');
% plot(T,me1,'k'),hold on
% ciplot(me1-se1 , me1+se1, T,'k'),alpha 0.5,hold on
% 
% 
% me2 = nanmean(m1');
% se2 = nansem(m1');
% plot(T,me2,'r')
% ciplot(me2-se2 , me2+se2, T,'r'),alpha 0.5
% xlim([-0.3 0.3])
% % ylim([-0.7 2.2])
% xline(0,'--')
% 
% 
% 
% 
% % All just members
% 
% m = Post.aversive.ccgM;%(:,logical(I.aversive));
% 
% m1 = m;
% 
% figure,
% me2 = nanmean(m1');
% se2 = nansem(m1');
% plot(T,me2,'r')
% ciplot(me2-se2 , me2+se2, T,'k'),alpha 0.5
% xlim([-0.3 0.3])
% % ylim([-0.7 2.2])
% xline(0,'--')
% hold on
% 
% % Cleaning the data
% m = Pre.aversive.ccgM;%(:,logical(I.aversive));
% m1 = m;
% 
% me2 = nanmean(m1');
% se2 = nansem(m1');
% plot(T,me2,'r')
% ciplot(me2-se2 , me2+se2, T,'r'),alpha 0.5
% xlim([-0.3 0.3])
% % ylim([-0.7 2.2])
% xline(0,'--')