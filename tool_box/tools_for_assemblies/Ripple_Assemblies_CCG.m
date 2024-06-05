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

% storage variables
Pre.aversive.ccg = [] ;         Post.aversive.ccg = [] ;
Pre.reward.ccg = [] ;           Post.reward.ccg = [] ;

Pre.aversive.iterator = [] ;    Post.aversive.iterator = [] ;
Pre.reward.iterator = [] ;      Post.reward.iterator = [] ;

Pre.aversive.ccgM = [] ;        Post.aversive.ccgM = [] ;
Pre.reward.ccgM = [] ;          Post.reward.ccgM = [] ;

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
                            [ccg1,T] = CCG(s,ids,'binSize',0.025,'duration',2,'smooth',2,'mode','ccg');
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
                            [ccg2,T] = CCG(s,ids,'binSize',0.025,'duration',2,'smooth',2,'mode','ccg');
                            cond2 = true;
                        else
                            cond2 = false;
                        end
                        
                        if and(cond1,cond2)
                            Pre.aversive.ccg = [Pre.aversive.ccg , (ccg1(:,1,2)./0.025)./length(x1)] ;
                            Pre.aversive.iterator = [Pre.aversive.iterator ; pInc1 pDec1] ;
                            Post.aversive.ccg = [Post.aversive.ccg , (ccg2(:,1,2)./0.025)./length(x2)] ;
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
                    else
                        TS.pre = [ripple_event.reward(:,1) ripple_event.reward(:,3)];
                        TS.post = [ripple_event.aversive(:,1) ripple_event.aversive(:,3)];
                    end
                    
                    
                    members = Thresholded.aversive(:,cond.both);
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
                                            [ccg1,T] = CCG(s,ids,'binSize',0.005,'duration',2,'smooth',1,'mode','ccg');
                                            cond1 = true;
                                        else
                                            cond1 = false;
                                        end
                                        
                                        % --- Post ---
                                        xx = Restrict(x,TS.post);
                                        yy = Restrict(y,TS.post);
                                        if and(length(xx)>5 , length(yy)>5)
                                            [s,ids,groups] = CCGParameters(xx,ones(length(xx),1),yy,ones(length(yy),1)*2);
                                            [ccg2,T] = CCG(s,ids,'binSize',0.005,'duration',2,'smooth',1,'mode','ccg');
                                            cond2 = true;
                                        else
                                            cond2 = false;
                                        end
                                        
                                        if and(cond1,cond2)
                                            Pre.aversive.ccgM = [Pre.aversive.ccgM , zscore(ccg1(:,1,2)./sum(ccg1(:,1,2)))];
                                            Post.aversive.ccgM = [Post.aversive.ccgM , zscore(ccg2(:,1,2)./sum(ccg2(:,1,2)))];
                                            if conditional(i)
                                                I.aversive = [I.aversive , true];
                                            else
                                                I.aversive = [I.aversive , false];
                                            end
                                        end
                                        clear ccg1 ccg2 cond1 cond2
                                        
                                    end
                                end
                            end
                        end
                    end
                    
                    
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
                            [ccg1,T] = CCG(s,ids,'binSize',0.025,'duration',2,'smooth',2,'mode','ccg');
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
                            [ccg2,T] = CCG(s,ids,'binSize',0.025,'duration',2,'smooth',2,'mode','ccg');
                            cond2 = true;
                        else
                            cond2 = false;
                        end
                        
                        if and(cond1,cond2)
                            Pre.reward.ccg = [Pre.reward.ccg , (ccg1(:,1,2)./0.025)./length(x1)] ;
                            Pre.reward.iterator = [Pre.reward.iterator ; pInc1 pDec1] ;
                            Post.reward.ccg = [Post.reward.ccg , (ccg2(:,1,2)./0.025)./length(x2)] ;
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
                    else
                        TS.pre = [ripple_event.baseline(:,1) ripple_event.baseline(:,3)];
                        TS.post = [ripple_event.reward(:,1) ripple_event.reward(:,3)];
                    end
                    
                    

                    members = Thresholded.reward(:,cond.both);
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
                                            [ccg1,T] = CCG(s,ids,'binSize',2,'duration',0.2,'smooth',1,'mode','ccg');
                                            cond1 = true;
                                        else
                                            cond1 = false;
                                        end
                                        
                                        % --- Post ---
                                        xx = Restrict(x,TS.post);
                                        yy = Restrict(y,TS.post);
                                        if and(length(xx)>5 , length(yy)>5)
                                            [s,ids,groups] = CCGParameters(xx,ones(length(xx),1),yy,ones(length(yy),1)*2);
                                            [ccg2,T] = CCG(s,ids,'binSize',2,'duration',0.2,'smooth',1,'mode','ccg');
                                            cond2 = true;
                                        else
                                            cond2 = false;
                                        end
                                        
                                        if and(cond1,cond2)
                                            Pre.reward.ccgM = [Pre.reward.ccgM , zscore(ccg1(:,1,2)./sum(ccg1(:,1,2)))];
                                            Post.reward.ccgM = [Post.reward.ccgM , zscore(ccg2(:,1,2)./sum(ccg2(:,1,2)))];
                                            if conditional(i)
                                                I.reward = [I.reward , true];
                                            else
                                                I.reward = [I.reward , false];
                                            end
                                        end
                                        clear ccg1 ccg2 cond1 cond2
                                        
                                    end
                                end
                            end
                        end
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
%
% x = and((Pre.aversive.iterator(:,1)<0.001) , Post.aversive.iterator(:,1)<0.001);
% y.pre = Pre.aversive.ccg(:,x);
% y.post = Post.aversive.ccg(:,x);
%
% tmp = [];
% for i = 1 : sum(x)
%     tmp = [tmp ; max(y.post(35:47,i)) - max(y.pre(35:47,i))];
% end
%
%
%
% x = and((Pre.reward.iterator(:,1)<0.001) , Post.reward.iterator(:,1)<0.001);
% y.pre = Pre.reward.ccg(:,x);
% y.post = Post.reward.ccg(:,x);
%
% tmp1 = [];
% for i = 1 : sum(x)
%     tmp1 = [tmp1 ; max(y.post(35:47,i)) - max(y.pre(35:47,i))];
% end
%
%
%
%
%
%
%
% figure
% tmp = and((Pre.aversive.iterator(:,1)<0.05) , Post.aversive.iterator(:,1)<0.05);
% ccg1 = Pre.aversive.ccg;
% ccg2 = Post.aversive.ccg;
% 
% plot(T,nanmean(ccg1(:,tmp)')),hold on
% plot(T,nanmean(ccg2(:,tmp)')),hold on
% 
% 
% figure
% tmp = and(not(Pre.reward.iterator(:,1)<0.05) , Post.reward.iterator(:,1)<0.05)
% ccg1 = Pre.reward.ccg;
% ccg2 = Post.reward.ccg;
% 
% plot(nanmean(ccg1(:,tmp)')),hold on
% plot(nanmean(ccg2(:,tmp)')),hold on
% % 
% % Pre.aversive.ccgM(:,sum(isnan(Post.aversive.ccgM))>0) = [];
% % I.aversive(:,sum(isnan(Post.aversive.ccgM))>0) = [];
% % Post.aversive.ccgM(:,sum(isnan(Post.aversive.ccgM))>0) = [];
% % 
% % m = Post.aversive.ccgM(:,logical(I.aversive));
% % mm = Pre.aversive.ccgM(:,logical(I.aversive));
% % 
% % %% Plot of Aversive assemblies
% % figure
% % [i ii] = max(m);
% % [i ii] = sort(ii);
% % subplot(122),imagesc(T,[1:1:size(m,2)],m(:,ii)')
% % subplot(121),imagesc(T,[1:1:size(mm,2)],mm(:,ii)')
% % 
% % 
% % figure,
% % m1 = nanmean(mm');
% % s1 = nansem(mm');
% % plot(T,m1,'k'),hold on
% % ciplot(m1-s1 , m1+s1, T,'k'),alpha 0.5
% % 
% % 
% % m2 = nanmean(m');
% % s2 = nansem(m');
% % plot(T,m2,'r')
% % ciplot(m2-s2 , m2+s2, T,'r'),alpha 0.5
% % 
% % xline(0,'--')
% % xlim([-0.05 0.05])
% % %% Plot of Reward assemblies
% % figure
% % [i ii] = max(Post.reward.ccgM);
% % [i ii] = sort(ii);
% % subplot(122),imagesc(T,[1:1:size(Post.reward.ccgM,2)],Post.reward.ccgM(:,ii)')
% % 
% % % [i ii] = max(Pre.reward.ccgM);
% % % [i ii] = sort(ii);
% % subplot(121),imagesc(T,[1:1:size(Post.reward.ccgM,2)],Pre.reward.ccgM(:,ii)')
% % 
% % 
% % figure,
% % m = nanmean(Pre.reward.ccgM');
% % s = nansem(Pre.reward.ccgM');
% % plot(T,m,'k'),hold on
% % ciplot(m-s , m+s, T,'k'),alpha 0.5
% % 
% % 
% % m = nanmean(Post.reward.ccgM');
% % s = nansem(Post.reward.ccgM');
% % plot(T,m,'r')
% % ciplot(m-s , m+s, T,'r'),alpha 0.5
% % 
% % xline(0,'--')