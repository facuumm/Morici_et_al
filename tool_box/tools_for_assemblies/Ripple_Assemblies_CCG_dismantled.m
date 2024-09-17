function [Pre Post T I durations Mod shock] = Ripple_Assemblies_CCG_dismantled(path)
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
Pre.aversive.iterator = [] ;            Post.aversive.iterator = [] ;
Pre.reward.iterator = [] ;              Post.reward.iterator = [] ;

I.aversive.dHPC = [];                    I.reward.dHPC = [];
I.aversive.vHPC = [];                    I.reward.vHPC = [];

Post.aversive.members.dHPC = [];         Post.aversive.members.vHPC = [];
Post.reward.members.dHPC = [];           Post.reward.members.vHPC = [];
Pre.aversive.members.dHPC = [];          Pre.aversive.members.vHPC = [];
Pre.reward.members.dHPC = [];            Pre.reward.members.vHPC = [];

Post.aversive.nomembers.dHPC = [];       Post.aversive.nomembers.vHPC = [];
Post.reward.nomembers.dHPC = [];         Post.reward.nomembers.vHPC = [];
Pre.aversive.nomembers.dHPC = [];        Pre.aversive.nomembers.vHPC = [];
Pre.reward.nomembers.dHPC = [];          Pre.reward.nomembers.vHPC = [];

durations.dHPC = [];                     durations.vHPC = [];
durations.event = [];                    durations.bursts = [];

Mod.pre.aversive.members.dHPC = [];      Mod.pre.aversive.members.vHPC = [];
Mod.pre.reward.members.dHPC = [];        Mod.pre.reward.members.vHPC = [];
Mod.post.aversive.members.dHPC = [];     Mod.post.aversive.members.vHPC = [];
Mod.post.reward.members.dHPC = [];       Mod.post.reward.members.vHPC = [];

shock.aversive.members.dHPC  = [];       shock.aversive.members.vHPC  = [];

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
%             bursts.coordinated.DV(bursts.coordinated.DV(:,2)-bursts.coordinated.DV(:,1)>0.3, : ) = [];
            bursts.baseline = Restrict(bursts.coordinated.DV,NREM.baseline);
            bursts.aversive = Restrict(bursts.coordinated.DV,NREM.aversive);
            bursts.reward = Restrict(bursts.coordinated.DV,NREM.reward);
            durations.bursts = [durations.bursts ; bursts.coordinated.DV(:,2)-bursts.coordinated.DV(:,1)];
            
            %% Detection of first event of the bursts
            tmp1 = []; % dorsal
            tmp2 = []; % ventral
            tmp3 = []; % all dorsal
            tmp4 = []; % all ventral
            for i = 1 : size(bursts.coordinated.DV,1)
                index = bursts.coordinated.DV(i,:);
                D = Restrict(ripplesD(:,1:3),index);
                V = Restrict(ripplesV(:,1:3),index);
                
%                 if and(size(D,1)>1 , size(V,1)>1)
                tmp1 = [tmp1 ; D(1,:)];
                tmp2 = [tmp2 ; V(1,:)];
%                 end
                
                tmp3 = [tmp3 ; D];
                tmp4 = [tmp4 ; V];
                
                clear D V index
            end
            
            bursts.FirstRipple.dorsal  = tmp1;
            bursts.FirstRipple.ventral = tmp2; clear tmp1 tmp2
            bursts.allRipples.dorsal   = tmp3;
            bursts.allRipples.ventral  = tmp4; clear tmp3 tmp4
            
        end
        
        %% Spikes
        % Load Units
        disp('Uploading Spiking activity')
        cd 'Spikesorting'
        [clusters , numberD , numberV , spks , spks_dHPC , spks_vHPC , cellulartype] = load_SU_FM(cd,criteria_type,criteria_fr,aversiveTS_run,rewardTS_run);
        
        %% Assemblies detection
        if and(numberD > 3 , numberV > 3)
            
            % Load data regarding the shock
            load('vHPC_responsivness_all.mat')
            load('dHPC_responsivness_all.mat')
            
            % --- Aversive ---
            disp('Lets go for the assemblies')
            if isfile('dorsalventral_assemblies_aversiveVF.mat')
                disp('Loading Aversive template')
                load('dorsalventral_assemblies_aversiveVF.mat')
                
                Thresholded.aversive = Th;    patterns.aversive = pat;     clear cond Th pat
                patterns.aversive = patterns.aversive .* Thresholded.aversive;
                
                cond  = classification_of_asselblies(Thresholded.aversive,clusters.dHPC);% Detection of members
                
                % Check if I have joint assemblies
                if sum(cond.both)>0
                    patterns.aversive = patterns.aversive(:,cond.both);
                    
                    if aversiveTS_run(1)<rewardTS_run(1)
                        TS.pre = NREM.baseline;
                        TS.post = NREM.aversive;
                        TS.run.run1 = movement.aversive;
                        TS.run.run2 = movement.reward;
                        TS.sleep.pre = baselineTS;
                        TS.sleep.post = aversiveTS;
                    else
                        TS.pre = NREM.reward;
                        TS.post = NREM.aversive;
                        TS.run.run1 = movement.aversive;
                        TS.run.run2 = movement.reward;
                        TS.sleep.pre = rewardTS;
                        TS.sleep.post = aversiveTS;
                    end
                    
                    % Members definition
                    members = Thresholded.aversive(:,cond.both);
                    members = sum(members,2);
                    members = members>=1;
                    
                    if not(isempty(bursts.coordinated.DV))
                        i = 1;
                        %% Pre aversive members
                        iteratorD = [];
                        for ii = 1 : numberD
                            if (members(ii,i))
                                member1 = clusters.dHPC(ii);
                                y = spks_dHPC(spks_dHPC(:,1) == member1 , 2);
                                y = Restrict(y,TS.pre);
%                                 x = Restrict(bursts.coordinated.DV(:,1:3),TS.pre);
                                x = Restrict(bursts.coordinated.DV,TS.pre);
                                if and(size(x,1)>5 , size(y,1)>5)
%                                     [pInc pDec surp] = RippleModulation(ripplesV(:,1:3),spks,member1,TS.pre);
                                    [s,ids,groups] = CCGParameters(x(:,1),ones(length(x(:,2)),1),y,ones(length(y),1)*2);
                                    [ccg,T] = CCG(s,ids,'binSize',0.005,'duration',2,'smooth',0,'mode','ccg');
                                    ccg = ccg(:,1,2)./length(x); ccg = ccg./0.005;
                                    
%                                     y = spks_dHPC(spks_dHPC(:,1) == member1 , 2);
%                                     [m] = meanFR_outside_ripples(ripplesD(:,1:3) , NREM.all , y);
%                                     ccg = ccg./m; clear m
                                    
                                    Pre.aversive.members.dHPC = [Pre.aversive.members.dHPC , ccg]; clear ccg x y m s ids groups
                                    I.aversive.dHPC = [I.aversive.dHPC ; countA];
                                    iteratorD = [iteratorD ; true];
%                                     Mod.pre.aversive.members.dHPC = [Mod.pre.aversive.members.dHPC ; pInc , pDec]; clear pDec pInc surp
                                    
                                    % Shock information storage
                                     shock.aversive.members.dHPC  = [shock.aversive.members.dHPC ; dHPC_resp.resp_ave(dHPC_resp.id == member1)];
                                else
                                    iteratorD = [iteratorD ; false];
                                end
                            else
                                iteratorD = [iteratorD ; false];
                            end
                        end
                        
                        iteratorV = [];
                        for iii = 1 : numberV
                            if (members(iii+numberD,i))
                                member2 = clusters.vHPC(iii);
                                y = spks_vHPC(spks_vHPC(:,1) == member2 , 2);
                                y = Restrict(y,TS.pre);
%                                 x = Restrict(bursts.coordinated.DV(:,1:3),TS.pre);
                                x = Restrict(bursts.coordinated.DV,TS.pre);
                                if and(size(x,1)>5 , size(y,1)>5)
%                                     [pInc pDec surp] = RippleModulation(ripplesD(:,1:3),spks,member1,TS.pre);
                                    [s,ids,groups] = CCGParameters(x(:,1),ones(length(x(:,2)),1),y,ones(length(y),1)*2);
                                    [ccg,T] = CCG(s,ids,'binSize',0.005,'duration',2,'smooth',0,'mode','ccg');
                                    ccg = ccg(:,1,2)./length(x); ccg = ccg./0.005;
                                    
%                                     y = spks_vHPC(spks_vHPC(:,1) == member2 , 2);
%                                     [m] = meanFR_outside_ripples(ripplesV(:,1:3) , NREM.all , y);
%                                     ccg = ccg./m; clear m
                                    
                                    Pre.aversive.members.vHPC = [Pre.aversive.members.vHPC , ccg]; clear ccg x y m s ids groups
                                    I.aversive.vHPC = [I.aversive.vHPC ; countA];
                                    iteratorV = [iteratorV ; true];
%                                     Mod.pre.aversive.members.vHPC = [Mod.pre.aversive.members.vHPC ; pInc , pDec]; clear pDec pInc surp
                                    
                                    % Shock information storage
                                    shock.aversive.members.vHPC  = [shock.aversive.members.vHPC ; vHPC_resp.resp_ave(vHPC_resp.id == member2)];
                                else
                                    iteratorV = [iteratorV ; false];
                                end
                            else
                                iteratorV = [iteratorV ; false];
                            end
                        end
                        
                        %% Post  aversive members
                        for ii = 1 : numberD
                            if (members(ii,i))
                                member1 = clusters.dHPC(ii);
                                y = spks_dHPC(spks_dHPC(:,1) == member1 , 2);
                                y = Restrict(y,TS.post);
%                                 x = Restrict(bursts.coordinated.DV(:,1:3),TS.post);
                                x = Restrict(bursts.coordinated.DV,TS.post);
                                if and(size(x,1)>5 , size(y,1)>5)
                                    if iteratorD(ii)
%                                         [pInc pDec surp] = RippleModulation(ripplesV(:,1:3),spks,member1,TS.post);
                                        [s,ids,groups] = CCGParameters(x(:,1),ones(length(x(:,2)),1),y,ones(length(y),1)*2);
                                        [ccg,T] = CCG(s,ids,'binSize',0.005,'duration',2,'smooth',0,'mode','ccg');
                                        ccg = ccg(:,1,2)./length(x); ccg = ccg./0.005;
                                        
%                                         y = spks_dHPC(spks_dHPC(:,1) == member1 , 2);
%                                         [m] = meanFR_outside_ripples(ripplesD(:,1:3) , NREM.all , y);
%                                         ccg = ccg./m; clear m
                                        
                                        Post.aversive.members.dHPC = [Post.aversive.members.dHPC , ccg]; clear ccg x y m s ids groups
                                        I.aversive.dHPC = [I.aversive.dHPC ; countA];
%                                         Mod.post.aversive.members.dHPC = [Mod.post.aversive.members.dHPC ; pInc , pDec]; clear pDec pInc surp
                                    end
                                end
                            end
                        end
                        
                        for iii = 1 : numberV
                            if (members(iii+numberD,i))
                                member2 = clusters.vHPC(iii);
                                y = spks_vHPC(spks_vHPC(:,1) == member2 , 2);
                                y = Restrict(y,TS.post);
%                                 x = Restrict(bursts.coordinated.DV(:,1:3),TS.post);
                                x = Restrict(bursts.coordinated.DV,TS.post);
                                if and(size(x,1)>5 , size(y,1)>5)
                                    if iteratorV(iii)
%                                         [pInc pDec surp] = RippleModulation(ripplesD(:,1:3),spks,member1,TS.post);
                                        [s,ids,groups] = CCGParameters(x(:,1),ones(length(x(:,2)),1),y,ones(length(y),1)*2);
                                        [ccg,T] = CCG(s,ids,'binSize',0.005,'duration',2,'smooth',0,'mode','ccg');
                                        ccg = ccg(:,1,2)./length(x); ccg = ccg./0.005;
                                        
%                                         y = spks_vHPC(spks_vHPC(:,1) == member2 , 2);
%                                         [m] = meanFR_outside_ripples(ripplesV(:,1:3) , NREM.all , y);
%                                         ccg = ccg./m; clear m                                        
                                        
                                        Post.aversive.members.vHPC = [Post.aversive.members.vHPC , ccg]; clear ccg x y m s ids groups
                                        I.aversive.vHPC = [I.aversive.vHPC ; countA];
%                                         Mod.post.aversive.members.vHPC = [Mod.post.aversive.members.vHPC ; pInc , pDec]; clear pDec pInc surp
                                    end
                                end
                            end
                        end
                    end
                    countA = countA+1;
                    clear RBA
                    clear pks TS
                end
            end
            
            
            % --- Reward ---
            if isfile('dorsalventral_assemblies_rewardVF.mat')
                disp('Loading Reward template')
                load('dorsalventral_assemblies_rewardVF.mat')
                
                Thresholded.reward = Th;   patterns.reward = pat;    clear cond Th pat
                
                patterns.reward = patterns.reward .* Thresholded.reward;
                cond  = classification_of_asselblies(Thresholded.reward,clusters.dHPC);% Detection of members
                
                % Check if I have joint assemblies
                if sum(cond.both)>0
                    
                    patterns.reward = patterns.reward(:,cond.both);
                    
                    if aversiveTS_run(1)<rewardTS_run(1)
                        TS.pre = NREM.aversive;
                        TS.post = NREM.reward;
                        TS.run.run1 = movement.reward;
                        TS.run.run2 = movement.aversive;
                        TS.sleep.pre = aversiveTS;
                        TS.sleep.post = rewardTS;
                    else
                        TS.pre = NREM.baseline;
                        TS.post = NREM.reward;
                        TS.run.run2 = movement.aversive;
                        TS.sleep.pre = baselineTS;
                        TS.sleep.post = rewardTS;
                    end
                    
                    % Members definition
                    members = Thresholded.reward(:,cond.both);
                    members = sum(members,2);
                    members = members>=1;
                    
                    if not(isempty(bursts.coordinated.DV))
                        i = 1;
                        %% Pre reward members
                        iteratorD = [];
                        for ii = 1 : numberD
                            if (members(ii,i))
                                member1 = clusters.dHPC(ii);
                                y = spks_dHPC(spks_dHPC(:,1) == member1 , 2);
                                y = Restrict(y,TS.pre);
%                                 x = Restrict(bursts.coordinated.DV(:,1:3),TS.pre);
                                x = Restrict(bursts.coordinated.DV,TS.pre);
                                if and(size(x,1)>5 , size(y,1)>5)
%                                     [pInc pDec surp] = RippleModulation(ripplesV(:,1:3),spks,member1,TS.pre);
                                    [s,ids,groups] = CCGParameters(x(:,1),ones(length(x(:,2)),1),y,ones(length(y),1)*2);
                                    [ccg,T] = CCG(s,ids,'binSize',0.005,'duration',2,'smooth',0,'mode','ccg');
                                    ccg = ccg(:,1,2)./length(x); ccg = ccg./0.005;
                                    
%                                     y = spks_dHPC(spks_dHPC(:,1) == member1 , 2);
%                                     [m] = meanFR_outside_ripples(ripplesD(:,1:3) , NREM.all , y);
%                                     ccg = ccg./m; clear m
                                        
                                    Pre.reward.members.dHPC = [Pre.reward.members.dHPC , ccg]; clear ccg x y m s ids groups
                                    I.reward.dHPC = [I.reward.dHPC ; countR];
                                    iteratorD = [iteratorD , true];
%                                     Mod.pre.reward.members.dHPC = [Mod.pre.reward.members.dHPC ; pInc , pDec]; clear pDec pInc surp
                                else
                                    iteratorD = [iteratorD , false];
                                end
                            else
                                iteratorD = [iteratorD , false];
                            end
                        end
                        
                        iteratorV = [];
                        for iii = 1 : numberV
                            if (members(iii+numberD,i))
                                member2 = clusters.vHPC(iii);
                                y = spks_vHPC(spks_vHPC(:,1) == member2 , 2);
                                y = Restrict(y,TS.pre);
%                                 x = Restrict(bursts.coordinated.DV(:,1:3),TS.pre);
                                x = Restrict(bursts.coordinated.DV,TS.pre);
                                if and(size(x,1)>5 , size(y,1)>5)
%                                     [pInc pDec surp] = RippleModulation(ripplesD(:,1:3),spks,member1,TS.pre);
                                    [s,ids,groups] = CCGParameters(x(:,1),ones(length(x(:,2)),1),y,ones(length(y),1)*2);
                                    [ccg,T] = CCG(s,ids,'binSize',0.005,'duration',2,'smooth',0,'mode','ccg');
                                    ccg = ccg(:,1,2)./length(x); ccg = ccg./0.005;
                                    
%                                     y = spks_vHPC(spks_vHPC(:,1) == member2 , 2);
%                                     [m] = meanFR_outside_ripples(ripplesV(:,1:3) , NREM.all , y);
%                                     ccg = ccg./m; clear m
                                    
                                    Pre.reward.members.vHPC = [Pre.reward.members.vHPC , ccg]; clear ccg x y m s ids groups
                                    I.reward.vHPC = [I.reward.vHPC ; countR];
                                    iteratorV = [iteratorV , true];
%                                     Mod.pre.reward.members.vHPC = [Mod.pre.reward.members.vHPC ; pInc , pDec]; clear pDec pInc surp
                                else
                                    iteratorV = [iteratorV , false];
                                end
                            else
                                iteratorV = [iteratorV , false];
                            end
                        end
                        
                        %% Post reward members
                        for ii = 1 : numberD
                            if (members(ii,i))
                                member1 = clusters.dHPC(ii);
                                y = spks_dHPC(spks_dHPC(:,1) == member1 , 2);
                                y = Restrict(y,TS.post);
%                                 x = Restrict(bursts.coordinated.DV(:,1:3),TS.post);
                                x = Restrict(bursts.coordinated.DV,TS.post);
                                if and(size(x,1)>5 , size(y,1)>5)
                                    if iteratorD(ii)
%                                         [pInc pDec surp] = RippleModulation(ripplesV(:,1:3),spks,member1,TS.post);
                                        [s,ids,groups] = CCGParameters(x(:,1),ones(length(x(:,2)),1),y,ones(length(y),1)*2);
                                        [ccg,T] = CCG(s,ids,'binSize',0.005,'duration',2,'smooth',0,'mode','ccg');
                                        ccg = ccg(:,1,2)./length(x); ccg = ccg./0.005;
                                        
%                                         y = spks_dHPC(spks_dHPC(:,1) == member1 , 2);
%                                         [m] = meanFR_outside_ripples(ripplesD(:,1:3) , NREM.all , y);
%                                         ccg = ccg./m; clear m
                                        
                                        Post.reward.members.dHPC = [Post.reward.members.dHPC , ccg]; clear ccg x y m s ids groups
                                        I.reward.dHPC = [I.reward.dHPC ; countR];
%                                         Mod.post.reward.members.dHPC = [Mod.post.reward.members.dHPC ; pInc , pDec]; clear pDec pInc surp
                                    end
                                end
                            end
                        end
                        
                        for iii = 1 : numberV
                            if (members(iii+numberD,i))
                                member2 = clusters.vHPC(iii);
                                y = spks_vHPC(spks_vHPC(:,1) == member2 , 2);
                                y = Restrict(y,TS.post);
%                                 x = Restrict(bursts.coordinated.DV(:,1:3),TS.post);
                                x = Restrict(bursts.coordinated.DV,TS.post);
                                if and(size(x,1)>5 , size(y,1)>5)
                                    if iteratorV(iii)
%                                         [pInc pDec surp] = RippleModulation(ripplesD(:,1:3),spks,member1,TS.post);
                                        [s,ids,groups] = CCGParameters(x(:,1),ones(length(x(:,2)),1),y,ones(length(y),1)*2);
                                        [ccg,T] = CCG(s,ids,'binSize',0.005,'duration',2,'smooth',0,'mode','ccg');
                                        ccg = ccg(:,1,2)./length(x); ccg = ccg./0.005;
                                        
%                                         y = spks_vHPC(spks_vHPC(:,1) == member2 , 2);
%                                         [m] = meanFR_outside_ripples(ripplesV(:,1:3) , NREM.all , y);
%                                         ccg = ccg./m; clear m
                                        
                                        Post.reward.members.vHPC = [Post.reward.members.vHPC , ccg]; clear ccg x y m s ids groups
                                        I.reward.vHPC = [I.reward.vHPC ; countR];
%                                         Mod.post.reward.members.vHPC = [Mod.post.reward.members.vHPC ; pInc , pDec]; clear pDec pInc surp
                                    end
                                end
                            end
                        end
                    end
                    countR = countR+1;
                    
                    clear RBR
                    clear pks TS
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
    clear vHPC_resp dHPC_resp
end

end
