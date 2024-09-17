function [Post T I shock] = CCG_members_function(path)
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
sm = 2;
dur = 0.6;

% storage variables
Post.aversive.members = [];             Post.reward.members = [];

durations.dHPC = [];                     durations.vHPC = [];
durations.event = [];                    durations.bursts = [];

shock.aversive.members.dHPC  = [];       shock.aversive.members.vHPC  = [];
shock.reward.members.dHPC  = [];         shock.reward.members.vHPC  = [];


I.reward.post = [];                      I.aversive.post = [];
    
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
                        TS.pre = Restrict(bursts.coordinated.DV,NREM.baseline);
                        TS.post = Restrict(bursts.coordinated.DV,NREM.aversive);
                        TS.run.run1 = movement.aversive;
                        TS.run.run2 = movement.reward;
                        TS.sleep.pre = baselineTS;
                        TS.sleep.post = aversiveTS;
                    else
                        TS.pre = Restrict(bursts.coordinated.DV,NREM.reward);
                        TS.post = Restrict(bursts.coordinated.DV,NREM.aversive);
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
                        %% Post aversive members
                        for ii = 1 : numberD
                            if (members(ii,i))
                                member1 = clusters.dHPC(ii);
                                x = spks_dHPC(spks_dHPC(:,1) == member1 , 2);
                                x = Restrict(x,TS.post);
                                if size(x,1)>5
                                    for iii = 1 : numberV
                                        if (members(iii+numberD,i))
                                            member2 = clusters.vHPC(iii);
                                            y = spks_vHPC(spks_vHPC(:,1) == member2 , 2);
                                            y = Restrict(y,TS.post);
                                            if size(y,1)>5
                                                [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                                                [ccg,T] = CCG(s,ids,'binSize',0.003,'duration',0.5,'smooth',sm,'mode','ccg');
                                                ccg = ccg(:,1,2)./sum(ccg(:,1,2));
                                                
                                                surrogate = surrogate_ccg(x,y,0.025,0.5,sm,[]);
                                                surrogate = quantile(surrogate,0.9);
                                                
                                                [TT TTT] = min(abs(T-(-0.025)));
                                                [TT TTTT] = min(abs(T-(0.025)));
                                                
                                                if max(ccg(TTT:TTTT)) > surrogate
                                                    I.aversive.post = [I.aversive.post ; true];
                                                else
                                                    I.aversive.post = [I.aversive.post ; false];
                                                end
                                                
                                                Post.aversive.members = [Post.aversive.members , ccg]; clear ccg y m s ids groups TT TTT TTTT
                                                
                                                % Shock information storage
                                                shock.aversive.members.dHPC  = [shock.aversive.members.dHPC ; dHPC_resp.resp_ave(dHPC_resp.id == member1)];
                                                shock.aversive.members.vHPC  = [shock.aversive.members.vHPC ; vHPC_resp.resp_ave(vHPC_resp.id == member2)];
                                            end
                                        end
                                    end
                                    clear x
                                end
                            end
                            
                        end
                        clear TS
                    end
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
                        TS.pre = Restrict(bursts.coordinated.DV,NREM.aversive);
                        TS.post = Restrict(bursts.coordinated.DV,NREM.reward);
                        TS.run.run1 = movement.reward;
                        TS.run.run2 = movement.aversive;
                        TS.sleep.pre = aversiveTS;
                        TS.sleep.post = rewardTS;
                    else
                        TS.pre = Restrict(bursts.coordinated.DV,NREM.baseline);
                        TS.post = Restrict(bursts.coordinated.DV,NREM.reward);
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
                        %% Post reward members
                        for ii = 1 : numberD
                            if (members(ii,i))
                                member1 = clusters.dHPC(ii);
                                x = spks_dHPC(spks_dHPC(:,1) == member1 , 2);
                                x = Restrict(x,TS.post);
                                if size(x,1)>5
                                    for iii = 1 : numberV
                                        if (members(iii+numberD,i))
                                            member2 = clusters.vHPC(iii);
                                            y = spks_vHPC(spks_vHPC(:,1) == member2 , 2);
                                            y = Restrict(y,TS.post);
                                            if size(y,1)>5
                                                [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                                                [ccg,T] = CCG(s,ids,'binSize',0.003,'duration',0.5,'smooth',sm,'mode','ccg');
                                                ccg = ccg(:,1,2)./sum(ccg(:,1,2));
                                                
                                                surrogate = surrogate_ccg(x,y,0.025,0.5,sm,[]);
                                                surrogate = quantile(surrogate,0.9);
                                                
                                                [TT TTT] = min(abs(T-(-0.025)));
                                                [TT TTTT] = min(abs(T-(0.025)));
                                                
                                                if max(ccg(TTT:TTTT)) > surrogate
                                                    I.reward.post = [I.reward.post ; true];
                                                else
                                                    I.reward.post = [I.reward.post ; false];
                                                end
                                                
                                                Post.reward.members = [Post.reward.members , ccg]; clear ccg y m s ids groups TT TTT TTTT
                                                
                                                % Shock information storage
                                                shock.reward.members.dHPC  = [shock.reward.members.dHPC ; dHPC_resp.resp_ave(dHPC_resp.id == member1)];
                                                shock.reward.members.vHPC  = [shock.reward.members.vHPC ; vHPC_resp.resp_ave(vHPC_resp.id == member2)];
                                            end
                                        end
                                    end
                                end
                            end
                            clear x
                        end
                        clear TS
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
end
