function [Pre , Post , T] = PHIST_SUxRipples_Coor_Uncoor(path)
% This function calculates the Pre and Post sleep SU activity from one
% region and take coordinated and uncoordinated ripples from the other
% region.
%
% INPUTS
% path: cell, inside each cell you should put the path of each session you
%       want to analyze
%
% OUTPUT
% Pre, Post: Structure, it store the SU response sourrounding the ripple.
%
% T: vector, it contains time vector for visualization
%
% I: Structure, it contains the tag to determine if the SU is or not
%    shock responssive
%
% curve: Structure, it contains the response of the neuron sourrounding the
%        shock. It is zscored.
%
% Morci Juan Facundo 07/2025

% variables to use in the script
criteria_fr = 0;
criteria_type = 0; % criteria for celltype (0:pyr, 1:int, 2:all)
sm = 0;
dur = 4;
bin = 0.01;

% storage variables
Pre.aversive.dHPC = [];      Post.aversive.dHPC = [];
Pre.reward.dHPC = [];        Post.reward.dHPC = [];


Pre.aversive.vHPC = [];      Post.aversive.vHPC = [];
Pre.reward.vHPC = [];        Post.reward.vHPC = [];

T = [];
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
        
        %% load sleep states
        disp('Uploading sleep scoring')
        x = dir([cd,'\*-states.mat']);    states = load([cd,'\',x.name]);    states = states.states;
        REM.all = ToIntervals(states==5);    NREM.all = ToIntervals(states==3);
        clear x states
        
        baselineTS = baselineTS./1000;                  aversiveTS = aversiveTS./1000;                    rewardTS = rewardTS./1000;
        aversiveTS_run = aversiveTS_run./1000;          rewardTS_run = rewardTS_run./1000;
        
        NREM.baseline = Restrict(NREM.all,baselineTS);   NREM.aversive = Restrict(NREM.all,aversiveTS);   NREM.reward = Restrict(NREM.all,rewardTS);
        REM.baseline = Restrict(REM.all,baselineTS);     REM.aversive = Restrict(REM.all,aversiveTS);     REM.reward = Restrict(REM.all,rewardTS);
        
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
            RB = true;
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
                    
                    peak = min(min(z(indice,2),r(2)))+abs(z(indice,2)-r(2))/2;
                    low = min([r(1) , z(indice,1)]);
                    up = max([r(3) , z(indice,3)]);
                    cooridnated_event = [cooridnated_event ; low , peak , up];
                    
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
        else
            RB = false;
        end
        
        %% Spikes
        % Load Units
        disp('Uploading Spiking activity')
        cd 'Spikesorting'
        [clusters , numberD , numberV , spks , spks_dHPC , spks_vHPC , cellulartype] = load_SU_FM(cd,criteria_type,criteria_fr,aversiveTS_run,rewardTS_run);
        
        if RB
            for i = 1 : 2
                if i == 1
                    if aversiveTS_run(1)<rewardTS_run(1)
                        TS.pre = NREM.baseline;       TS.post = NREM.aversive;
                    else
                        TS.pre = NREM.reward;         TS.post = NREM.aversive;
                    end
                else
                    if aversiveTS_run(1)<rewardTS_run(1)
                        TS.pre = NREM.aversive;        TS.post = NREM.reward;
                    else
                        TS.pre = NREM.baseline;        TS.post = NREM.reward;
                    end
                end
                
                for ii = 1 : numberD
                    % Pre
                    y = Restrict(spks(spks(:,1)==clusters.dHPC(ii),2),TS.pre);
                    x = Restrict(ripples.vHPC.coordinated.all(:,1:3),TS.pre);
%                     [m] = meanFR_outside_ripples(ripplesV(:,1:3) , NREM.all , y , []);
                    if and(size(x,1)>5 , size(y,1)>5)
                        [s,ids,groups] = CCGParameters(x(:,2),ones(length(x(:,2)),1),y,ones(length(y),1)*2);
                        [ccg1,T] = CCG(s,ids,'binSize',bin,'duration',dur,'smooth',sm,'mode','ccg');
                        ccg1 = ccg1(:,1,2)./size(x,1);    ccg1 = ccg1./bin;% ccg1 = ccg1./m;
                    else
                        ccg1 = [];
                    end
                    clear x y s ids groups
                    
                    % Post
                    y = Restrict(spks(spks(:,1)==clusters.dHPC(ii),2),TS.post);
                    x = Restrict(ripples.vHPC.coordinated.all(:,1:3),TS.post);
                    if and(size(x,1)>5 , size(y,1)>5)
                        [s,ids,groups] = CCGParameters(x(:,2),ones(length(x(:,2)),1),y,ones(length(y),1)*2);
                        [ccg2,T] = CCG(s,ids,'binSize',bin,'duration',dur,'smooth',sm,'mode','ccg');
                        ccg2 = ccg2(:,1,2)./size(x,1);    ccg2 = ccg2./bin;% ccg2 = ccg2./m;
                    else
                        ccg2 = [];
                    end
                    clear x y s ids groups
                    
                    if and(not(isempty(ccg1)) ,  not(isempty(ccg2)))
                        if i == 1
                            Pre.aversive.dHPC = [Pre.aversive.dHPC , ccg1];
                            Post.aversive.dHPC = [Post.aversive.dHPC , ccg2];
                        else
                            Pre.reward.dHPC = [Pre.reward.dHPC , ccg1];
                            Post.reward.dHPC = [Post.reward.dHPC , ccg2];
                        end
                    end
                    clear ccg1 ccg2
                end
                clear TS
            end
            clear i ii
        end
        
        if RB
            for i = 1 : 2
                if i == 1
                    if aversiveTS_run(1)<rewardTS_run(1)
                        TS.pre = NREM.baseline;       TS.post = NREM.aversive;
                    else
                        TS.pre = NREM.reward;         TS.post = NREM.aversive;
                    end
                else
                    if aversiveTS_run(1)<rewardTS_run(1)
                        TS.pre = NREM.aversive;        TS.post = NREM.reward;
                    else
                        TS.pre = NREM.baseline;        TS.post = NREM.reward;
                    end
                end
                
                for ii = 1 : numberV
                    % Pre
                    y = Restrict(spks(spks(:,1)==clusters.vHPC(ii),2),TS.pre);
                    x = Restrict(ripples.dHPC.coordinated.all(:,1:3),TS.pre);
%                     [m] = meanFR_outside_ripples(ripplesD(:,1:3) , NREM.all , y , []);
                    if and(size(x,1)>5 , size(y,1)>5)
                        [s,ids,groups] = CCGParameters(x(:,2),ones(length(x(:,2)),1),y,ones(length(y),1)*2);
                        [ccg1,T] = CCG(s,ids,'binSize',bin,'duration',dur,'smooth',sm,'mode','ccg');
                        ccg1 = ccg1(:,1,2)./size(x,1);    ccg1 = ccg1./bin;% ccg1 = ccg1./m;
                    else
                        ccg1 = [];
                    end
                    clear x y s ids groups
                    
                    % Post
                    y = Restrict(spks(spks(:,1)==clusters.vHPC(ii),2),TS.post);
                    x = Restrict(ripples.dHPC.coordinated.all(:,1:3),TS.post);
                    
                    if and(size(x,1)>5 , size(y,1)>5)
                        [s,ids,groups] = CCGParameters(x(:,2),ones(length(x(:,2)),1),y,ones(length(y),1)*2);
                        [ccg2,T] = CCG(s,ids,'binSize',bin,'duration',dur,'smooth',sm,'mode','ccg');
                        ccg2 = ccg2(:,1,2)./size(x,1);    ccg2 = ccg2./bin;% ccg2 = ccg2./m;
                    else
                        ccg2 = [];
                    end
                    clear x y s ids groups
                    
                    if and(not(isempty(ccg1)) ,  not(isempty(ccg2)))
                        if i == 1
                            Pre.aversive.vHPC = [Pre.aversive.vHPC , ccg1];
                            Post.aversive.vHPC = [Post.aversive.vHPC , ccg2];
                        else
                            Pre.reward.vHPC = [Pre.reward.vHPC , ccg1];
                            Post.reward.vHPC = [Post.reward.vHPC , ccg2];
                        end
                    end
                    clear ccg1 ccg2
                end
                clear TS
            end
            clear i ii
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
    clear cooridnated_eventDV cooridnated_eventVD segments movement RD RV RB
end
end
% figure
% % --- dHPC ---
% subplot(221)
% s = nansem(zscore(Pre.aversive.dHPC)');
% m = nanmean(zscore(Pre.aversive.dHPC)');
% plot(T,m,'k'),hold on
% ciplot(m-s,m+s,T,'k'), alpha 0.5
% s = nansem(zscore(Post.aversive.dHPC)');
% m = nanmean(zscore(Post.aversive.dHPC)');
% plot(T,m,'r'),hold on
% ciplot(m-s,m+s,T,'r'), alpha 0.5
% xlim([-0.6 0.6])
% subplot(222)
% s = nansem(zscore(Pre.reward.dHPC)');
% m = nanmean(zscore(Pre.reward.dHPC)');
% plot(T,m,'k'),hold on
% ciplot(m-s,m+s,T,'k'), alpha 0.5
% s = nansem(zscore(Post.aversive.dHPC)');
% m = nanmean(zscore(Post.aversive.dHPC)');
% plot(T,m,'b'),hold on
% ciplot(m-s,m+s,T,'b'), alpha 0.5
% xlim([-0.6 0.6])
% 
% % --- vHPC ---
% subplot(223)
% s = nansem(zscore(Pre.aversive.vHPC)');
% m = nanmean(zscore(Pre.aversive.vHPC)');
% plot(T,m,'k'),hold on
% ciplot(m-s,m+s,T,'k'), alpha 0.5
% s = nansem(zscore(Post.aversive.vHPC)');
% m = nanmean(zscore(Post.aversive.vHPC)');
% plot(T,m,'r'),hold on
% ciplot(m-s,m+s,T,'r'), alpha 0.5
% xlim([-0.6 0.6])
% 
% subplot(224)
% s = nansem(zscore(Pre.reward.vHPC)');
% m = nanmean(zscore(Pre.reward.vHPC)');
% plot(T,m,'k'),hold on
% ciplot(m-s,m+s,T,'k'), alpha 0.5
% s = nansem(zscore(Post.aversive.vHPC)');
% m = nanmean(zscore(Post.aversive.vHPC)');
% plot(T,m,'b'),hold on
% ciplot(m-s,m+s,T,'b'), alpha 0.5
% xlim([-0.6 0.6])
% 