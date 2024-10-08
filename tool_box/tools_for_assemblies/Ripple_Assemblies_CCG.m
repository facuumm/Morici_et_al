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

I.aversive.dHPC = [];            I.reward.dHPC = [];
I.aversive.vHPC = [];            I.reward.vHPC = [];

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
                    
                    if r(2)<z(2) % keep only when dorsal happen first
                        cooridnated_event = [cooridnated_event ; r];
                    end
                    
%                     peak = min(min(z(indice,2),r(2)))+abs(z(indice,2)-r(2))/2;
%                     low = min([r(1) , z(indice,1)]);
%                     up = max([r(3) , z(indice,3)]);
%                     cooridnated_event = [cooridnated_event ; low , peak , up];    
%                     
%                     durations.dHPC = [durations.dHPC ; r(3)-r(1)];
%                     durations.vHPC = [durations.vHPC ; z(indice,3)-z(indice,1)];
%                     durations.event = [durations.event ; up-low];
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
            
        end
        
        %% Spikes
        % Load Units
        disp('Uploading Spiking activity')
        cd 'Spikesorting'
        [clusters , numberD , numberV , spks , spks_dHPC , spks_vHPC , cellulartype] = load_SU_FM(cd,criteria_type,criteria_fr,aversiveTS_run,rewardTS_run);
        
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
                    
                    % Detection of peaks
                    pks.aversive = assemblies_peaks([bins' Spikes] , patterns.aversive , th);
                                        
                    iterador = [];
                    % Iteration across assemblies
                    for i = 1 : size(pks.aversive,1)
                        % --- Pre ---
                        y = Restrict(pks.aversive{i}(:,1),TS.pre);
                          x1 = Restrict(ripple_event.all(:,2),TS.pre);
%                         x1 = Restrict(bursts.coordinated.DV(:,1),TS.pre);
                        if and(length(y)>5 , length(x1)>2)
                            [pInc1, pDec1 surp1] = RippleModulation_assemblies(ripple_event.all,y,TS.pre);

                            [s,ids,groups] = CCGParameters(x1,ones(length(x1),1),y,ones(length(y),1)*2);
                            [ccg1,T] = CCG(s,ids,'binSize',0.01,'duration',4,'smooth',0,'mode','ccg');
                            cond1 = true;
                            
                            x1 = Restrict([ripple_event.all(:,1)+0.05 ripple_event.all(:,3)+0.05],TS.pre);
                            m1 = InvertIntervals(x1,TS.pre(1,1) , TS.pre(end,2));
                            m1 = length(Restrict(y,m1))/(sum(m1(:,2)-m1(:,1))); clear c cc
                        else
                            cond1 = false;
                        end
                        
                        % --- Post ---
                        y = Restrict(pks.aversive{i}(:,1),TS.post);
                        x2 = Restrict(ripple_event.all(:,2),TS.post);
%                         x2 = Restrict(bursts.coordinated.DV(:,1),TS.post);
                        
                        if and(length(y)>5 , length(x2)>2)
                            [pInc2 pDec2 surp2] = RippleModulation_assemblies(ripple_event.all,y,TS.post);
                            
                            [s,ids,groups] = CCGParameters(x2,ones(length(x2),1),y,ones(length(y),1)*2);
                            [ccg2,T] = CCG(s,ids,'binSize',0.01,'duration',4,'smooth',0,'mode','ccg');
                            cond2 = true;
                            
                            x2 = Restrict([ripple_event.all(:,1)+0.05 ripple_event.all(:,3)+0.05],TS.post);
                            m2 = InvertIntervals(x2,TS.post(1,1) , TS.post(end,2));
                            m2 = length(Restrict(y,m2))/(sum(m2(:,2)-m2(:,1)));
                        else
                            cond2 = false;
                        end
                        
                        if and(cond1,cond2)
                           x1 = Restrict(ripple_event.all(:,2),TS.pre);
                           x2 = Restrict(ripple_event.all(:,2),TS.post);
%                             x1 = Restrict(bursts.coordinated.DV(:,1),TS.pre);
%                             x2 = Restrict(bursts.coordinated.DV(:,1),TS.post);
                            
                            ccg1 = (ccg1(:,1,2)./length(x1)); ccg1 = ccg1./0.01;
                            ccg2 = (ccg2(:,1,2)./length(x2)); ccg2 = ccg2./0.01;
                            Pre.aversive.ccg = [Pre.aversive.ccg , (ccg1)];%./m1)];%
                            Pre.aversive.iterator = [Pre.aversive.iterator ; pInc1 pDec1] ;
                            Post.aversive.ccg = [Post.aversive.ccg , (ccg2)];%./m2)];
                            Post.aversive.iterator = [Post.aversive.iterator ; pInc2 pDec2] ;
                            
%                             % all
%                             [average1 , time] = triggered_average_Ripples(x1,patterns.aversive(:,i),logical(1),[bins' Spikes],[-0.2 0.6],'all',[numberD numberV]);
%                             [average2 , time] = triggered_average_Ripples(x2,patterns.aversive(:,i),logical(1),[bins' Spikes],[-0.2 0.6],'all',[numberD numberV]);
%                             Pre.aversive.ccgM = [Pre.aversive.ccgM  , average1];
%                             Post.aversive.ccgM = [Post.aversive.ccgM  , average2]; clear average1 average2
%                             
%                             % Just dHPC
%                             [average1 , time] = triggered_average_Ripples(x1,patterns.aversive(:,i),logical(1),[bins' Spikes],[-0.2 0.6],'dHPC',[numberD numberV]);
%                             [average2 , time] = triggered_average_Ripples(x2,patterns.aversive(:,i),logical(1),[bins' Spikes],[-0.2 0.6],'dHPC',[numberD numberV]);
%                             Pre.aversive.ccgD = [Pre.aversive.ccgD  , average1];
%                             Post.aversive.ccgD = [Post.aversive.ccgD  , average2]; clear average1 average2
%                             
%                             % Just vHPC
%                             [average1 , time] = triggered_average_Ripples(x1,patterns.aversive(:,i),logical(1),[bins' Spikes],[-0.2 0.6],'vHPC',[numberD numberV]);
%                             [average2 , time] = triggered_average_Ripples(x2,patterns.aversive(:,i),logical(1),[bins' Spikes],[-0.2 0.6],'vHPC',[numberD numberV]);
%                             Pre.aversive.ccgV = [Pre.aversive.ccgV  , average1];
%                             Post.aversive.ccgV = [Post.aversive.ccgV  , average2]; clear average1 average2
%                             
%                             % Just Both
%                             [average1 , time] = triggered_average_Ripples(x1,patterns.aversive(:,i),logical(1),[bins' Spikes],[-0.2 0.6],'both',[numberD numberV]);
%                             [average2 , time] = triggered_average_Ripples(x2,patterns.aversive(:,i),logical(1),[bins' Spikes],[-0.2 0.6],'both',[numberD numberV]);
%                             Pre.aversive.ccgB = [Pre.aversive.ccgB  , average1];
%                             Post.aversive.ccgB = [Post.aversive.ccgB  , average2]; clear average1 average2
                            
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
                        x1 = Restrict(ripple_event.all(:,2),TS.pre);
%                         x1 = Restrict(bursts.coordinated.DV(:,1),TS.pre);
                        if and(length(y)>5 , length(x1)>2)
                            [pInc1 pDec1 surp1] = RippleModulation_assemblies(ripple_event.all,y,TS.pre);
                            
                            [s,ids,groups] = CCGParameters(x1,ones(length(x1),1),y,ones(length(y),1)*2);
                            [ccg1,T] = CCG(s,ids,'binSize',0.01,'duration',4,'smooth',0,'mode','ccg');
                            cond1 = true;
                            
                            x1 = Restrict([ripple_event.all(:,1)+0.05 ripple_event.all(:,3)+0.05],TS.pre);
                            m1 = InvertIntervals(x1,TS.pre(1,1) , TS.pre(end,2));
                            m1 = length(Restrict(y,m1))/(sum(m1(:,2)-m1(:,1))); clear c cc
                            cond1 = true;
                        else
                            cond1 = false;
                        end
                        
                        % --- Post ---
                        y = Restrict(pks.reward{i}(:,1),TS.post);
%                         x2 = Restrict(bursts.coordinated.DV(:,1),TS.post);
                        x2 = Restrict(ripple_event.all(:,2),TS.post);
                        if and(length(y)>5 , length(x2)>2)
                            [pInc2 pDec2 surp2] = RippleModulation_assemblies(ripple_event.all,y,TS.post);
                            
                            [s,ids,groups] = CCGParameters(x2,ones(length(x2),1),y,ones(length(y),1)*2);
                            [ccg2,T] = CCG(s,ids,'binSize',0.01,'duration',4,'smooth',0,'mode','ccg');
                            cond2 = true;
                            
                            x2 = Restrict([ripple_event.all(:,1)+0.05 ripple_event.all(:,3)+0.05],TS.post);
                            m2 = InvertIntervals(x2,TS.post(1,1) , TS.post(end,2));
                            m2 = length(Restrict(y,m2))/(sum(m2(:,2)-m2(:,1)));
                        else
                            cond2 = false;
                        end
                        
                        if and(cond1,cond2)
                            x1 = Restrict(ripple_event.all(:,2),TS.pre);
                            x2 = Restrict(ripple_event.all(:,2),TS.post);
%                             x1 = Restrict(bursts.coordinated.DV(:,1),TS.pre);
%                             x2 = Restrict(bursts.coordinated.DV(:,1),TS.post);
                            
                            ccg1 = (ccg1(:,1,2)./length(x1)); ccg1 = ccg1./0.01;
                            ccg2 = (ccg2(:,1,2)./length(x2)); ccg2 = ccg2./0.01;
                            Pre.reward.ccg = [Pre.reward.ccg , (ccg1)];%./m1)];%
                            Pre.reward.iterator = [Pre.reward.iterator ; pInc1 pDec1] ;
                            Post.reward.ccg = [Post.reward.ccg , (ccg2)];%./m2)];
                            Post.reward.iterator = [Post.reward.iterator ; pInc2 pDec2] ;
                            
%                             % all
%                             [average1 , time] = triggered_average_Ripples(x1,patterns.reward(:,i),logical(1),[bins' Spikes],[-0.2 0.6],'all',[numberD numberV]);
%                             [average2 , time] = triggered_average_Ripples(x2,patterns.reward(:,i),logical(1),[bins' Spikes],[-0.2 0.6],'all',[numberD numberV]);
%                             Pre.reward.ccgM = [Pre.reward.ccgM  , average1];
%                             Post.reward.ccgM = [Post.reward.ccgM  , average2]; clear average1 average2
%                             
%                             % Just dHPC
%                             [average1 , time] = triggered_average_Ripples(x1,patterns.reward(:,i),logical(1),[bins' Spikes],[-0.2 0.6],'dHPC',[numberD numberV]);
%                             [average2 , time] = triggered_average_Ripples(x2,patterns.reward(:,i),logical(1),[bins' Spikes],[-0.2 0.6],'dHPC',[numberD numberV]);
%                             Pre.reward.ccgD = [Pre.reward.ccgD  , average1];
%                             Post.reward.ccgD = [Post.reward.ccgD  , average2]; clear average1 average2
%                             
%                             % Just vHPC
%                             [average1 , time] = triggered_average_Ripples(x1,patterns.reward(:,i),logical(1),[bins' Spikes],[-0.2 0.6],'vHPC',[numberD numberV]);
%                             [average2 , time] = triggered_average_Ripples(x2,patterns.reward(:,i),logical(1),[bins' Spikes],[-0.2 0.6],'vHPC',[numberD numberV]);
%                             Pre.reward.ccgV = [Pre.reward.ccgV  , average1];
%                             Post.reward.ccgV = [Post.reward.ccgV  , average2]; clear average1 average2
%                             
%                             % Just Both
%                             [average1 , time] = triggered_average_Ripples(x1,patterns.reward(:,i),logical(1),[bins' Spikes],[-0.2 0.6],'both',[numberD numberV]);
%                             [average2 , time] = triggered_average_Ripples(x2,patterns.reward(:,i),logical(1),[bins' Spikes],[-0.2 0.6],'both',[numberD numberV]);
%                             Pre.reward.ccgB = [Pre.reward.ccgB  , average1];
%                             Post.reward.ccgB = [Post.reward.ccgB  , average2]; clear average1 average2
                            
                        end
                        clear y pInc1 pInc2 pDec1 pDec2 surp1 surp2
                        clear s ids groups ccg x1 x2 m1 m2 ccg ccg1 ccg2
                        clear average1 average2
                    end
                    clear pks
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
% 

tmp = [];
for i = 1 : size(Post.aversive.ccg,2)
    t = (Pre.aversive.ccg(:,i));
    t = Smooth((t),2,'kernel','gaussian');
    
    [ii iii] = min(abs([-2 : 0.01 : 2]-(-0.3)));
    [ii iiii] = min(abs([-0.5 : 0.005 : 0.5]-0.3));
    C = (nanmean([t(1:iii) ; t(iiii:end)]));
    t = t./C;  
%     t = zscore(t);

    tmp = [tmp , t]; clear t
end

tmp1 = [];
for i = 1 : size(Post.aversive.ccg,2)
    t = (Post.aversive.ccg(:,i));
    t = Smooth((t),2,'kernel','gaussian');
    
    [ii iii] = min(abs([-2 : 0.01 : 2]-(-0.3)));
    [ii iiii] = min(abs([-0.5 : 0.005 : 0.5]-0.2));
    C = (nanmean([t(1:iii) ; t(iiii:end)]));
    t = t./C;
%     t = zscore(t);
    
    tmp1 = [tmp1 , t]; clear t
end

plot(T,nanmean(tmp'),'k'),hold on
ciplot(nanmean(tmp')-nansem(tmp') , nanmean(tmp')+nansem(tmp') , T , 'k'),alpha 0.5
plot(T,nanmean(tmp1'),'r'),hold on
ciplot(nanmean(tmp1')-nansem(tmp1') , nanmean(tmp1')+nansem(tmp1') , T , 'r'),alpha 0.5
xlim([-0.3 0.3])