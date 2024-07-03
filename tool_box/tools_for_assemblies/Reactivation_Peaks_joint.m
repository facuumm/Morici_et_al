function [Pre Post Run Peaks] = Reactivation_Peaks_joint(path)
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
Pre.aversive.thirty = [];   Post.aversive.thirty = [];
Pre.reward.thirty = [];     Post.reward.thirty = [];
Run.aversive = [];          Run.reward = [];

Peaks.aversive.Post = {};   Peaks.reward.Post = {};
Peaks.aversive.Pre = {};    Peaks.reward.Pre = {};
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
        if or(numberD > 3 , numberV > 3)
            
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
                
                if sum(cond.dHPC)>0
%                     if aversiveTS_run(1)<rewardTS_run(1)
%                     patterns.aversive = patterns.aversive(:,cond.both);
                    if aversiveTS_run(1)<rewardTS_run(1)
                        TS.pre = NREM.baseline;
                        TS.post = NREM.aversive;
%                         TS.pre = Restrict(NREM.baseline,[aversiveTS_run(1)-3600 aversiveTS_run(1)]);
%                         TS.post = Restrict(NREM.aversive,[aversiveTS_run(2) aversiveTS_run(2)+3600]);
%                         TS.pre = REM.baseline;
%                         TS.post = REM.aversive;
%                         TS.pre = Restrict([ripple_event.baseline(:,1) , ripple_event.baseline(:,3)],[aversiveTS_run(1)-3600 aversiveTS_run(1)]);
%                         TS.post = Restrict([ripple_event.aversive(:,1) , ripple_event.aversive(:,3)],[aversiveTS_run(2) aversiveTS_run(2)+3600]);
%                         TS.pre = Restrict(bursts.coordinated,NREM.baseline);
%                         TS.post = Restrict(bursts.coordinated,NREM.aversive);
                        TS.run.run1 = movement.aversive;
                        TS.run.run2 = movement.reward;
                    else
                        TS.pre = NREM.reward;
                        TS.post = NREM.aversive;
%                         TS.pre = Restrict(NREM.reward,[aversiveTS_run(1)-3600 aversiveTS_run(1)]);
%                         TS.post = Restrict(NREM.aversive,[aversiveTS_run(2) aversiveTS_run(2)+3600]);
%                         TS.pre = REM.reward;
%                         TS.post = REM.aversive;
%                         TS.pre = Restrict([ripple_event.reward(:,1) , ripple_event.reward(:,3)],[aversiveTS_run(1)-3600 aversiveTS_run(1)]);
%                         TS.post = Restrict([ripple_event.aversive(:,1) , ripple_event.aversive(:,3)],[aversiveTS_run(2) aversiveTS_run(2)+3600]);
%                         TS.pre = Restrict(bursts.coordinated,NREM.reward);
%                         TS.post = Restrict(bursts.coordinated,NREM.aversive);
                        TS.run.run1 = movement.aversive;
                        TS.run.run2 = movement.reward;
                    end
%                     pks.aversive = assemblies_peaks_joint(patterns.aversive , cond.both , [bins' Spikes], Thresholded.aversive, [numberD numberV]);
                    pks.aversive = assemblies_peaks([bins' Spikes] , patterns.aversive(:,cond.dHPC) , th);
                    [pre1 , post1] = Reactivation_Rate(pks.aversive,TS);
%                     [pre2 , post2] = Reactivation_Mean([bins' Spikes] , patterns.aversive(:,cond.both) , TS);
                    [pre2 , post2] = Reactivation_Mean_peaks(pks.aversive(:,1),TS);
                    
                    [run1 run2] = Activation_Run_Rate(pks.aversive,TS.run);
                    [run3 run4] = Activation_Run_Mean([bins' Spikes] , patterns.aversive(:,cond.dHPC) , TS.run);
                    
                    Pre.aversive.thirty = [Pre.aversive.thirty ; pre1  pre2]; clear pre1 pre2
                    Post.aversive.thirty = [Post.aversive.thirty ; post1 post2]; clear post1 post2
                    Run.aversive = [Run.aversive ; run1 run2 run3 run4]; clear run1 run2 run3 run4
                    
                    for i = 1 : size(pks.aversive,1)
                        tmp = [];
                        for ii = 1 : size(TS.post,1)
                            if ii == 1
                                tmp = Restrict(pks.aversive{i}(:,1),TS.post(ii,:)) - TS.post(ii,1);
                                dt = TS.post(ii,2) - TS.post(ii,1);
                            else
                                tmp = [tmp ; Restrict(pks.aversive{i}(:,1),TS.post(ii,:)) - TS.post(ii,1) + dt;];
                                dt = (TS.post(ii,2) - TS.post(ii,1))+dt;
                            end
                        end
                        if isempty(Peaks.aversive.Post)
                            Peaks.aversive.Post{1} =  tmp;
                        else
                            Peaks.aversive.Post{end+1} =  tmp;
                        end
                        
                        clear dt
                        
                        tmp = [];
                        for ii = 1 : size(TS.pre,1)
                            if ii == 1
                                tmp = Restrict(pks.aversive{i}(:,1),TS.pre(ii,:)) - TS.pre(ii,1);
                                dt = TS.pre(ii,2) - TS.pre(ii,1);
                            else
                                tmp = [tmp ; Restrict(pks.aversive{i}(:,1),TS.pre(ii,:)) - TS.pre(ii,1) + dt;];
                                dt = (TS.pre(ii,2) - TS.pre(ii,1))+dt;
                            end
                        end
                        if isempty(Peaks.aversive.Pre)
                            Peaks.aversive.Pre{1} =  tmp;
                        else
                            Peaks.aversive.Pre{end+1} =  tmp;
                        end
                        
                        clear tmp dt ii m n
                    end
%                 end
                end
            end
            
            
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
                
                
                if sum(cond.dHPC)>0
%                     if aversiveTS_run(1)>rewardTS_run(1)
%                     patterns.reward = patterns.reward(:,cond.both);
                    if aversiveTS_run(1)<rewardTS_run(1)
                        TS.pre = NREM.aversive;
                        TS.post = NREM.reward;                        
%                         TS.pre = Restrict(NREM.aversive,[rewardTS_run(1)-3600 rewardTS_run(1)]);
%                         TS.post = Restrict(NREM.reward,[rewardTS_run(2) rewardTS_run(2)+3600]);
%                         TS.pre = REM.aversive;
%                         TS.post = REM.reward;
%                         TS.pre = Restrict([ripple_event.aversive(:,1) , ripple_event.aversive(:,3)],[rewardTS_run(1)-3600 rewardTS_run(1)]);
%                         TS.post = Restrict([ripple_event.reward(:,1) , ripple_event.reward(:,3)],[rewardTS_run(2) rewardTS_run(2)+3600]);
%                         TS.pre = Restrict(bursts.coordinated,NREM.aversive);
%                         TS.post = Restrict(bursts.coordinated,NREM.reward);
                        TS.run.run1 = movement.reward;
                        TS.run.run2 = movement.aversive;
                    else
%                         TS.pre = NREM.baseline;
%                         TS.post = NREM.reward;                        
%                         TS.pre = Restrict(NREM.baseline,[rewardTS_run(1)-3600 rewardTS_run(1)]);
%                         TS.post = Restrict(NREM.reward,[rewardTS_run(2) rewardTS_run(2)+3600]);
                        TS.pre = REM.baseline;
                        TS.post = REM.reward;
%                         TS.pre = Restrict([ripple_event.baseline(:,1) , ripple_event.baseline(:,3)],[rewardTS_run(1)-3600 rewardTS_run(1)]);
%                         TS.post = Restrict([ripple_event.reward(:,1) , ripple_event.reward(:,3)],[rewardTS_run(2) rewardTS_run(2)+3600]);
%                         TS.pre = Restrict(bursts.coordinated,NREM.baseline);
%                         TS.post = Restrict(bursts.coordinated,NREM.reward);
                        TS.run.run1 = movement.reward;
                        TS.run.run2 = movement.aversive;
                    end
                    
%                     pks.reward = assemblies_peaks_joint(patterns.reward , cond.both , [bins' Spikes], Thresholded.reward, [numberD numberV]);
                    pks.reward = assemblies_peaks([bins' Spikes] , patterns.reward(:,cond.dHPC) , th);
                    [pre1 , post1] = Reactivation_Rate(pks.reward,TS);
%                     [pre2 , post2] = Reactivation_Mean([bins' Spikes] , patterns.reward(:,cond.both) , TS);
                    [pre2 , post2] = Reactivation_Mean_peaks(pks.reward,TS);
                    [run1 run2] = Activation_Run_Rate(pks.reward,TS.run);
                    [run3 run4] = Activation_Run_Mean([bins' Spikes] , patterns.reward(:,cond.dHPC) , TS.run);                    
                    
                    Pre.reward.thirty = [Pre.reward.thirty ; pre1  pre2]; clear pre1 pre2
                    Post.reward.thirty = [Post.reward.thirty ; post1 post2]; clear post1 post2
                    Run.reward = [Run.reward ; run1 run2 run3 run4]; clear run1 run2 run3 run4
                    

                    for i = 1 : size(pks.reward,1)
                        tmp = [];
                        for ii = 1 : size(TS.post,1)
                            if ii == 1
                                tmp = Restrict(pks.reward{i}(:,1),TS.post(ii,:)) - TS.post(ii,1);
                                dt = (TS.post(ii,2) - TS.post(ii,1));
                            else
                                tmp = [tmp ; Restrict(pks.reward{i}(:,1),TS.post(ii,:)) - TS.post(ii,1) + dt;];
                                dt = (TS.post(ii,2) - TS.post(ii,1))+dt;
                            end
                        end
                        if isempty(Peaks.reward.Post)
                            Peaks.reward.Post{1} =  tmp;
                        else
                            Peaks.reward.Post{end+1} =  tmp;
                        end
                        clear dt
                        
                        tmp = [];
                        for ii = 1 : size(TS.pre,1)
                            if ii == 1
                                tmp = Restrict(pks.reward{i}(:,1),TS.pre(ii,:)) - TS.pre(ii,1);
                                dt = (TS.pre(ii,2) - TS.pre(ii,1));
                            else
                                tmp = [tmp ; Restrict(pks.reward{i}(:,1),TS.pre(ii,:)) - TS.pre(ii,1) + dt;];
                                dt = (TS.pre(ii,2) - TS.pre(ii,1))+dt;
                            end
                        end
                        if isempty(Peaks.reward.Pre)
                            Peaks.reward.Pre{1} =  tmp;
                        else
                            Peaks.reward.Pre{end+1} =  tmp;
                        end
                        clear tmp dt ii m n
                    end
%                 end
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
    clear num_assembliesA num_assembliesR
end



%% NREM

cumA = [];    cumAP = [];
cumA2 = [];  cumA2P = [];
DurationA = [];
ReactA.post = nan(40,50);
ReactA.pre = nan(40,50);
for i = 1 : size(Peaks.aversive.Post,2)
    %% Initialization of variables
    tmp = Peaks.aversive.Post{i};
    tmp1 = Peaks.aversive.Pre{i};
    
    % Calculation of duration of each assemblie activation
    AA = tmp1(end)-tmp1(1);    A = tmp(end)-tmp(1);
    
    %% Reactivation Trains
    [tmp,bins]=binspikes(tmp,1/60,[0 round(tmp(end))]);    % Post
    [tmp1,bins]=binspikes(tmp1,1/60,[0 round(tmp1(end))]); % Pre
    
    %% Rate
    tmp1 = flip(tmp1./60); % Pre
    tmp = tmp./60;         % Post
    
    % Definition of maximal extension
    if length(tmp) > length(tmp1)
        m = length(tmp1);
    else
        m = length(tmp);
    end
    
    % Keeping same quantity of data
    tmp = tmp(1:m);    tmp1 = tmp1(1:m);
    
    
    %% Saving data
    % Post
    p = [];
    if length(tmp)>=20
        cumA = [cumA ; tmp(1:20)];
        %          p = [p , nanmedian(tmp(1:20))];
    end
    
    if length(tmp)<20
        cumA = [cumA ; tmp(1:end)];
%         p = [p , nanmedian( tmp(1:20))];
    end
    
    if length(tmp)>=40
        cumA2 = [cumA2 ; tmp(20:40)];
    end
    
    if and(length(tmp)>20 , length(tmp)<40)
        cumA2 = [cumA2 ; tmp(20:end)];
    end
    
    % Pre
    pp = [];
    if length(tmp1)>=20
        cumAP = [cumAP ; tmp1(1:20)];
%         pp = [pp ; nanmedian(tmp1(1:20))];
    end
    
    if length(tmp1)<20
        cumAP = [cumAP ; tmp1(1:end)];
%         pp = [pp ; nanmedian(tmp1(1:20))];
    end
    
    if length(tmp1)>=40
        cumA2P = [cumA2P ; tmp1(20:40)];
    end
    
    if and(length(tmp1)>20 , length(tmp1)<40)
        cumA2P = [cumA2P ; tmp1(20:end)];
    end
    DurationA = [DurationA ; AA , A]; clear AA A
    
    
    if length(tmp)>=40
        ReactA.pre(:,i) = tmp1(1:40);
        ReactA.post(:,i) = tmp(1:40);
    else
        ReactA.pre(1:m,i) = tmp1(1:end);
        ReactA.post(1:m,i) = tmp(1:end);
    end
    
    
    %      ReactA = [ReactA ; p , pp]; clear p pp
end


figure,

ii = nanmean(ReactA.post(1:20,:));
[i ii] = sort(ii,'descend');

subplot(121)
imagesc([1:40],[1:size(ReactA.pre,2)],ReactA.pre(:,ii)'),caxis([0 0.5]),xlim([1 40])
subplot(122)
imagesc([1:40],[1:size(ReactA.post,2)],ReactA.post(:,ii)'),caxis([0 0.5]),xlim([1 40])


cumR = [];    cumRP = [];
cumR2 = [];  cumR2P = [];
DurationR = [];
ReactA.post = nan(40,51);
ReactA.pre = nan(40,51);
for i = 1 : size(Peaks.reward.Post,2)
    %% Initialization of variables
    tmp = Peaks.reward.Post{i};      
    tmp1 = Peaks.reward.Pre{i};      
    if and(not(isempty(tmp)),not(isempty(tmp1)))
    % Calculation of duration of each assemblie activation
    AA = tmp1(end)-tmp1(1);    A = tmp(end)-tmp(1);    
    
    %% Reactivation Trains 
    [tmp,bins]=binspikes(tmp,1/60,[0 round(tmp(end))]);    % Post
    [tmp1,bins]=binspikes(tmp1,1/60,[0 round(tmp1(end))]); % Pre    
    
    %% Rate
    tmp1 = flip(tmp1./60); % Pre
    tmp = tmp./60;         % Post

    % Definition of maximal extension
    if length(tmp) > length(tmp1)
        m = length(tmp1);
    else
        m = length(tmp);
    end
    
    tmp = tmp(1:m);
    tmp1 = tmp1(1:m);

    
    %% Saving data
    % Post
     if length(tmp)>=20
         cumR = [cumR ; tmp(1:20)];
     end
     
     if length(tmp)<20
         cumR = [cumR ; tmp(1:end)];
     end
     
     if length(tmp)>=40
         cumR2 = [cumR2 ; tmp(20:40)];
     end
     
     if and(length(tmp)>20 , length(tmp)<40)
         cumR2 = [cumR2 ; tmp(20:end)];
     end          
    
    % Pre    
     if length(tmp1)>=20
         cumRP = [cumRP ; tmp1(1:20)];
     end
     
     if length(tmp1)<20
         cumRP = [cumRP ; tmp1(1:end)];
     end     

     if length(tmp1)>=40
         cumR2P = [cumR2P ; tmp1(20:40)];
     end
     
     if and(length(tmp1)>20 , length(tmp1)<40)
         cumR2P = [cumR2P ; tmp1(20:end)];
     end     
     DurationR = [DurationR ; AA , A]; clear AA A
     
    if length(tmp)>=40
        ReactR.pre(:,i) = tmp1(1:40);
        ReactR.post(:,i) = tmp(1:40);
    else
        ReactR.pre(1:m,i) = tmp1(1:end);
        ReactR.post(1:m,i) = tmp(1:end);
    end
    end
end


% 
% figure,
% ii = nanmean(ReactR.post(1:20,:));
% [i ii] = sort(ii,'descend');
% 
% subplot(121)
% imagesc([1:40],[1:size(ReactR.pre,2)],ReactR.pre(:,ii)'),caxis([0 0.5]),xlim([1 40])
% subplot(122)
% imagesc([1:40],[1:size(ReactR.post,2)],ReactR.post(:,ii)'),caxis([0 0.5]),xlim([1 40])
% 
% 
% 
% %% Plor cumulatives
% % Aversive
% figure,
% subplot(121)
% cdfplot(cumA),hold on
% cdfplot(cumAP),
% xlim([0 1.2]),
% [h p] = kstest2(cumA,cumAP,'Tail','smaller')
% 
% subplot(122)
% cdfplot(cumA2),hold on
% cdfplot(cumA2P),
% xlim([0 1.2]),
% [h p] = kstest2(cumA2,cumA2P,'Tail','smaller')
% 
% 
% 
% % Reward
% figure,
% subplot(121)
% cdfplot(cumR),hold on
% cdfplot(cumRP),
% xlim([0 1.2]),
% [h, p] = kstest2(cumR,cumRP,'Tail','smaller')
% 
% subplot(122)
% cdfplot(cumR2),hold on
% cdfplot(cumR2P),
% xlim([0 1.2]),
% [h p] = kstest2(cumR2,cumR2P,'Tail','smaller')
% 
% 
% 
% [h p] = kstest2(cumA,cumR,'Tail','smaller')
% figure
% subplot(121)
% cdfplot(cumA),hold on
% cdfplot(cumR),
% xlim([0 0.6]),
% 
% [h p] = kstest2(cumA2,cumR2,'Tail','smaller')
% subplot(122)
% cdfplot(cumA2),hold on
% cdfplot(cumR2),
% xlim([0 0.6]),


% 
% figure,
% subplot(121)
% pd = fitdist(cumA,'Kernel','Kernel','epanechnikov')
% x = min(cumA):0.01:max(cumA);
% p = cdf(pd,x);
% plot(x,p)
% hold on
% pd = fitdist(cumAP,'Kernel','Kernel','epanechnikov')
% x = min(cumAP):0.01:max(cumAP);
% p = cdf(pd,x);
% plot(x,p)
% xlim([0 1]),
% 
% subplot(122)
% pd = fitdist(cumR,'Kernel','Kernel','epanechnikov')
% x = min(cumR):0.01:max(cumR);
% p = cdf(pd,x);
% plot(x,p)
% hold on
% pd = fitdist(cumRP,'Kernel','Kernel','epanechnikov')
% x = min(cumRP):0.01:max(cumRP);
% p = cdf(pd,x);
% plot(x,p)
% xlim([0 1]),
% 
% 
% figure,
% subplot(121)
% pd = fitdist(cumA2,'Kernel','Kernel','epanechnikov')
% x = min(cumA2):0.01:max(cumA2);
% p = cdf(pd,x);
% plot(x,p)
% hold on
% pd = fitdist(cumA2P,'Kernel','Kernel','epanechnikov')
% x = min(cumA2P):0.01:max(cumA2P);
% p = cdf(pd,x);
% plot(x,p)
% xlim([0 1]),
% 
% subplot(122)
% pd = fitdist(cumR2,'Kernel','Kernel','epanechnikov')
% x = min(cumR2):0.01:max(cumR2);
% p = cdf(pd,x);
% plot(x,p)
% hold on
% pd = fitdist(cumR2P,'Kernel','Kernel','epanechnikov')
% x = min(cumR2P):0.01:max(cumR2P);
% p = cdf(pd,x);
% plot(x,p)
% xlim([0 1]),


% %% REM
% 
% cumA = [];    cumAP = [];
% for i = 1 : size(Peaks.aversive.Post,2)
%     %% Initialization of variables
%     tmp = Peaks.aversive.Post{i};
%     tmp1 = Peaks.aversive.Pre{i};
%     
%     if and(not(isempty(tmp)),not(isempty(tmp1)))    % Calculation of duration of each assemblie activation
%         AA = tmp1(end)-tmp1(1);    A = tmp(end)-tmp(1);
%         
%         %% Reactivation Trains
%         [tmp,bins]=binspikes(tmp,1/60,[0 round(tmp(end))]);    % Post
%         [tmp1,bins]=binspikes(tmp1,1/60,[0 round(tmp1(end))]); % Pre
%         
%         %% Rate
%         tmp1 = flip(tmp1./60); % Pre
%         tmp = tmp./60;         % Post
%         
%         
%         %% Saving data
%         % Post
%         cumA = [cumA ; tmp];
%         % Pre
%         cumAP = [cumAP ; tmp1];
%     end
% end
% 
% 
% cumR = [];    cumRP = [];
% for i = 1 : size(Peaks.reward.Post,2)
%     %% Initialization of variables
%     tmp = Peaks.reward.Post{i};
%     tmp1 = Peaks.reward.Pre{i};
%     
%     if and(not(isempty(tmp)),not(isempty(tmp1)))    % Calculation of duration of each assemblie activation
%         AA = tmp1(end)-tmp1(1);    A = tmp(end)-tmp(1);
%         
%         %% Reactivation Trains
%         [tmp,bins]=binspikes(tmp,1/60,[0 round(tmp(end))]);    % Post
%         [tmp1,bins]=binspikes(tmp1,1/60,[0 round(tmp1(end))]); % Pre
%         
%         %% Rate
%         tmp1 = flip(tmp1./60); % Pre
%         tmp = tmp./60;         % Post
%         
%         
%         %% Saving data
%         % Post
%         cumR = [cumR ; tmp];
%         % Pre
%         cumRP = [cumRP ; tmp1];
%     end
% end
% 
% 
% %% Plor cumulatives
% % Aversive
% figure,
% subplot(121)
% cdfplot(cumA),hold on
% cdfplot(cumAP),
% xlim([0 2.6]),
% [h p] = kstest2(cumA,cumAP,'Tail','smaller')
% 
% % Reward
% subplot(122)
% cdfplot(cumR),hold on
% cdfplot(cumRP),
% xlim([0 2.6]),
% [h p] = kstest2(cumR,cumRP,'Tail','smaller')

end


