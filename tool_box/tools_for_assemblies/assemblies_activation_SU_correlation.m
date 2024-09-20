function [Pre Post I] = assemblies_activation_SU_correlation(path)
% This function calculates the Pre and Post sleep assemblies rate.
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
% Morci Juan Facundo 08/2024

% variables to use in the script
criteria_fr = 0;
criteria_type = 0; % criteria for celltype (0:pyr, 1:int, 2:all)
sm = 0;
dur = 0.1;
bin = 0.003;

% storage variables
Pre.aversive.members = [];      Post.aversive.members = [];
Pre.reward.members = [];        Post.reward.members = [];

I.aversive = [];
I.reward = [];


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
            if or(numberD > 3 , numberV > 3)
                
                limits = [0 segments.Var1(end)/1000];
                events = [];
                [Spikes , bins , Clusters] = spike_train_construction([spks_dHPC;spks_vHPC], clusters.all, cellulartype, 0.025, limits, events, false, true);
                clear limits events
                
                % timestamps for reactivation calculation
                is.sws.baseline = InIntervals(bins,[ripple_event.baseline(:,1) ripple_event.baseline(:,3)]);
                is.sws.reward = InIntervals(bins,[ripple_event.reward(:,1) ripple_event.reward(:,3)]);
                is.sws.aversive = InIntervals(bins,[ripple_event.aversive(:,1) ripple_event.aversive(:,3)]);
                
                ripple_event.filtered.baseline = ripple_event.baseline;
                ripple_event.filtered.aversive = ripple_event.aversive;
                ripple_event.filtered.reward = ripple_event.reward;
                
                is.sws.timestamps.sleep.aversive = NREM.aversive;
                is.sws.timestamps.sleep.reward = NREM.reward;
                is.sws.timestamps.sleep.baseline = NREM.baseline;
                
                
                is.sws.runaversive = InIntervals(bins,movement.aversive);
                is.sws.runreward = InIntervals(bins,movement.reward);
                
                is.sws.timestamps.run.aversive = InIntervals(bins,aversiveTS_run./1000);
                is.sws.timestamps.run.reward = InIntervals(bins,rewardTS_run./1000);
                
                is.sws.timestamps.aversiveSleep = aversiveTS./1000;
                is.sws.timestamps.aversiveRun = aversiveTS_run./1000;
                is.sws.timestamps.rewardSleep = rewardTS./1000;
                is.sws.timestamps.rewardRun = rewardTS_run./1000;
                is.sws.timestamps.baselineSleep = baselineTS./1000;
                
                % --- Aversive ---
                disp('Lets go for the assemblies')
                if isfile('dorsalventral_assemblies_aversiveVF.mat')
                    disp('Loading Aversive template')
                    [patterns , cond , Thresholded] = load_assemblies(cd , 'dorsalventral_assemblies_aversiveVF.mat', clusters, numberD , 'aversive');
                    th = logical(Thresholded.aversive(:,cond.both));
                    
                    % TS
                    if aversiveTS_run(1)<rewardTS_run(1)
                        TS.pre = [ripple_event.baseline(:,1) ripple_event.baseline(:,3)];
                        TS.post = [ripple_event.aversive(:,1) ripple_event.aversive(:,3)];
                    else
                        TS.pre = [ripple_event.reward(:,1) ripple_event.reward(:,3)];
                        TS.post = [ripple_event.aversive(:,1) ripple_event.aversive(:,3)];
                    end
                    
                    % loop for ccg
                    for j = 1 : size(th,2)
                        [R] = reactivation_strength(patterns.aversive , cond.both , [bins' , Spikes] , is.sws , 5 , 'A' , config , true , []); clear templates
                        
                        members = clusters.all(th(:,j));
                        limits = [0 segments.Var1(end)/1000];
                        events = [];
                        [x , y , c] = spike_train_construction(spks, members, cellulartype, 0.005, limits, events, false, true);
                        clear limits events
                        
                        % Pre
                        x1 = InIntervals(y,TS.pre);
                        xp = x(:,not(sum(x(x1,:))==0));
                        x1 = corr(xp(x1,:),'rows','complete');     x1(logical(diag(diag(x1)))) = nan;
                        x1 = mean(mean(x1,'omitnan'));
                        
                        % Post
                        x2 = InIntervals(y,TS.post);
                        xp = x(:,not(sum(x(x2,:))==0));
                        x2 = corr(xp(x2,:),'rows','complete');     x2(logical(diag(diag(x2)))) = nan;
                        x2 = mean(mean(x2,'omitnan'));
                        
                        Pre.aversive.members = [Pre.aversive.members , x1];
                        Post.aversive.members = [Post.aversive.members , x2];
                        
                        I.aversive = [I.aversive ; R(j,1)]; clear p x1 x2 x y c members xp
                    end
                end
                
                % --- Reward ---
                disp('Lets go for the assemblies')
                if isfile('dorsalventral_assemblies_rewardVF.mat')
                    disp('Loading Reward template')
                    [patterns , cond , Thresholded] = load_assemblies(cd , 'dorsalventral_assemblies_rewardVF.mat', clusters, numberD , 'reward');
                    th = logical(Thresholded.reward(:,cond.both));
                    
                    % TS
                    if aversiveTS_run(1)<rewardTS_run(1)
                        TS.pre = [ripple_event.aversive(:,1) ripple_event.aversive(:,3)];
                        TS.post = [ripple_event.reward(:,1) ripple_event.reward(:,3)];
                    else
                        TS.pre = [ripple_event.baseline(:,1) ripple_event.baseline(:,3)];
                        TS.post = [ripple_event.reward(:,1) ripple_event.reward(:,3)];
                    end
                    
                    % loop for ccg
                    for j = 1 : size(th,2)
                        [R] = reactivation_strength(patterns.reward , cond.both , [bins' , Spikes] , is.sws , 5 , 'R' , config , true , []); clear templates
                        
                        members = clusters.all(th(:,j));
                        limits = [0 segments.Var1(end)/1000];
                        events = [];
                        [x , y , c] = spike_train_construction(spks, members, cellulartype, 0.005, limits, events, false, true);
                        clear limits events
                        
                        % Pre
                        x1 = InIntervals(y,TS.pre);
                        xp = x(:,not(sum(x(x1,:))==0));
                        x1 = corr(xp(x1,:),'rows','complete');     x1(logical(diag(diag(x1)))) = nan;
                        x1 = mean(mean(x1,'omitnan'));
                        
                        % Post
                        x2 = InIntervals(y,TS.post);
                        xp = x(:,not(sum(x(x2,:))==0));
                        x2 = corr(xp(x2,:),'rows','complete');     x2(logical(diag(diag(x2)))) = nan;
                        x2 = mean(mean(x2,'omitnan'));
                        
                        Pre.reward.members = [Pre.reward.members , x1];
                        Post.reward.members = [Post.reward.members , x2];
                        
                        I.reward = [I.reward ; R(j,1)]; clear p x1 x2 x y c members xp
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
        clear cooridnated_eventDV cooridnated_eventVD segments movement RD RV RB
    end
end



end