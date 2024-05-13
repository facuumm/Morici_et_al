function [Pre Post Run Peaks] = Reactivation_Peaks_joint(path)
% This function calculates the Pre and Post sleep assemblies rate.
%
% INPUTS
% path: cell, inside each cell you should put the path of each session you
%       want to analyze
%
% OUTPUT
% Pre and Post: Structure, it store the Assemblies Rate for each type of
%               assemblies.
%
%               Architecture of each output:
%                   Pre.dHPC
%                      .vHPC.aversive
%                           .reward
%                   Post.dHPC
%                       .vHPC.aversive
%                            .reward
%
% Morci Juan Facundo 04/2024

% variables to use in the script
criteria_fr = 0;
criteria_type = 0; % criteria for celltype (0:pyr, 1:int, 2:all)
th = 3; % threshold for detecting peak assemblies
% 3 SD funciona
% 5 SD funciona

% storage variables
Pre.aversive.thirty = [];   Post.aversive.thirty = [];
Pre.reward.thirty = [];     Post.reward.thirty = [];
Run.aversive = [];          Run.reward = [];

Peaks.aversive.Post = {};   Peaks.reward.Post = {};

Pre.aversive.sixty = [];   Post.aversive.sixty = [];

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
                
%                 patterns.aversive = patterns.aversive .* Thresholded.aversive;
                
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
                    pks.aversive = assemblies_peaks_joint(patterns.aversive , cond.both , [bins' Spikes], Thresholded.aversive, [numberD numberV]);
%                     pks.aversive = assemblies_peaks([bins' Spikes] , patterns.aversive(:,cond.both) , th);
                    [pre1 , post1] = Reactivation_Rate(pks.aversive,TS);
%                     [pre2 , post2] = Reactivation_Mean([bins' Spikes] , patterns.aversive(:,cond.both) , TS);
                    [pre2 , post2] = Reactivation_Mean_peaks(pks.aversive(:,1),TS);
                    
                    [run1 run2] = Activation_Run_Rate(pks.aversive,TS.run);
                    [run3 run4] = Activation_Run_Mean([bins' Spikes] , patterns.aversive(:,cond.both) , TS.run);
                    
                    Pre.aversive.thirty = [Pre.aversive.thirty ; pre1  pre2]; clear pre1 pre2
                    Post.aversive.thirty = [Post.aversive.thirty ; post1 post2]; clear post1 post2
                    Run.aversive = [Run.aversive ; run1 run2 run3 run4]; clear run1 run2 run3 run4
                    
                    for i = 1 : size(pks.aversive,1)
                        if isempty(Peaks.aversive.Post)
                            Peaks.aversive.Post{1} =  Restrict(pks.aversive{i}(:,1),TS.post) - TS.post(1,1);
                        else
                            Peaks.aversive.Post{end+1} =  Restrict(pks.aversive{i}(:,1),TS.post) - TS.post(1,1);
                        end
                    end
                end
            end
            
            
            if isfile('dorsalventral_assemblies_rewardVF.mat')
                disp('Loading Reward template')
                load('dorsalventral_assemblies_rewardVF.mat')
                
                Thresholded.reward = Th;
                patterns.reward = pat;
                clear cond Th pat
                
%                 patterns.reward = patterns.reward .* Thresholded.reward;
                
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
                        TS.pre = NREM.baseline;
                        TS.post = NREM.reward;                        
%                         TS.pre = Restrict(NREM.baseline,[rewardTS_run(1)-3600 rewardTS_run(1)]);
%                         TS.post = Restrict(NREM.reward,[rewardTS_run(2) rewardTS_run(2)+3600]);
%                         TS.pre = REM.baseline;
%                         TS.post = REM.reward;
%                         TS.pre = Restrict([ripple_event.baseline(:,1) , ripple_event.baseline(:,3)],[rewardTS_run(1)-3600 rewardTS_run(1)]);
%                         TS.post = Restrict([ripple_event.reward(:,1) , ripple_event.reward(:,3)],[rewardTS_run(2) rewardTS_run(2)+3600]);
%                         TS.pre = Restrict(bursts.coordinated,NREM.baseline);
%                         TS.post = Restrict(bursts.coordinated,NREM.reward);
                        TS.run.run1 = movement.reward;
                        TS.run.run2 = movement.aversive;
                    end
                    
                    pks.reward = assemblies_peaks_joint(patterns.reward , cond.both , [bins' Spikes], Thresholded.reward, [numberD numberV]);
%                     pks.reward = assemblies_peaks([bins' Spikes] , patterns.reward(:,cond.both) , th);
                    [pre1 , post1] = Reactivation_Rate(pks.reward,TS);
%                     [pre2 , post2] = Reactivation_Mean([bins' Spikes] , patterns.reward(:,cond.both) , TS);
                    [pre2 , post2] = Reactivation_Mean_peaks(pks.reward,TS);
                    [run1 run2] = Activation_Run_Rate(pks.reward,TS.run);
                    [run3 run4] = Activation_Run_Mean([bins' Spikes] , patterns.reward(:,cond.both) , TS.run);                    
                    
                    Pre.reward.thirty = [Pre.reward.thirty ; pre1  pre2]; clear pre1 pre2
                    Post.reward.thirty = [Post.reward.thirty ; post1 post2]; clear post1 post2
                    Run.reward = [Run.reward ; run1 run2 run3 run4]; clear run1 run2 run3 run4
                    
                    for i = 1 : size(pks.reward,1)
                        if isempty(Peaks.reward.Post)
                            Peaks.reward.Post{1} =  Restrict(pks.reward{i}(:,1),TS.post) - TS.post(1,1);
                        else
                            Peaks.reward.Post{end+1} =  Restrict(pks.reward{i}(:,1),TS.post) - TS.post(1,1);
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
    clear num_assembliesA num_assembliesR
end
% 
% figure,
% subplot(121)
% scatter(Run.aversive(:,1),Run.aversive(:,2),'filled','r'),hold on
% scatter(Run.reward(:,1),Run.reward(:,2),'filled','b'),hold on
% subplot(122)
% scatter(Run.aversive(:,3),Run.aversive(:,4),'filled','r'),hold on
% scatter(Run.reward(:,3),Run.reward(:,4),'filled','b'),hold on
% 
% figure,
% subplot(121)
% x = [ones(size(Post.reward.thirty(:,1))) ; ones(size(Pre.aversive.thirty(:,1)))*2];
% DIA = (Post.aversive.thirty(:,1) - Pre.aversive.thirty(:,1)) ./ (Post.aversive.thirty(:,1) + Pre.aversive.thirty(:,1));
% DIR = (Post.reward.thirty(:,1) - Pre.reward.thirty(:,1)) ./ (Post.reward.thirty(:,1) + Pre.reward.thirty(:,1));
% y = [DIR ; DIA];
% scatter(x,y,"filled",'jitter','on', 'jitterAmount',0.1'),hold on,xlim([0 3])
% scatter([1 2] , [nanmedian(DIR) nanmedian(DIA)],'filled'),%set(gca, 'YScale', 'log')
% yline(0,'--')
% signrank(DIA,0)
% signrank(DIR,0)
% ranksum(DIA,DIR)
% ylim([-1 1])
% 
% subplot(222)
% fitlm(Post.reward.thirty(:,1),Run.reward(:,1))
% plot(ans), xlim([0 1]) , ylim([0 2.5])
% 
% subplot(224)
% fitlm(Post.aversive.thirty(:,1),Run.aversive(:,1))
% plot(ans), xlim([0 1]) , ylim([0 2.5])
% 
% % 
% % subplot(122)
% % x = [ones(size(Post.aversive.sixty(:,1))) ; ones(size(Pre.reward.sixty(:,1)))*2];
% % DIA = (Post.aversive.sixty(:,1) ./ Pre.aversive.sixty(:,1));% ./ (Post.aversive.sixty(:,1) + Pre.aversive.sixty(:,1));
% % DIR = (Post.reward.sixty(:,1) ./ Pre.reward.sixty(:,1));% ./ (Post.reward.sixty(:,1) + Pre.reward.sixty(:,1));
% % y = [DIA ; DIR];
% % scatter(x,y,"filled",'jitter','on', 'jitterAmount',0.1'),hold on,xlim([0 3])
% % scatter([1 2] , [nanmedian(DIA) nanmedian(DIR)],'filled'),set(gca, 'YScale', 'log')
% % yline(1,'--')
% % signrank(DIA,1)
% % signrank(DIR,1)
% % ranksum(DIA,DIR)
% % ylim([0 10])
% % boxplot(y,x),set(gca, 'YScale', 'log')
% % 
% % 
% % 
% % figure,
% % subplot(121)
% % x = [ones(size(Pre.dHPC.aversive(:,2))) ; ones(size(Post.dHPC.reward(:,2)))*2];
% % DIA = (Post.dHPC.aversive(:,2) ./ Pre.dHPC.aversive(:,2));% ./ (Post.dHPC.aversive(:,2) + Pre.dHPC.aversive(:,2));
% % DIR = (Post.dHPC.reward(:,2) ./ Pre.dHPC.reward(:,2));% ./ (Post.dHPC.reward(:,2) + Pre.dHPC.reward(:,2));
% % y = [DIA ; DIR];
% % scatter(x,y,"filled",'jitter','on', 'jitterAmount',0.1'),hold on,xlim([0 3])
% % scatter([1 2] , [nanmedian(DIA) nanmedian(DIR)],'filled'),set(gca, 'YScale', 'log')
% % yline(1,'--')
% % signrank(DIA,1)
% % signrank(DIR,1)
% % % ylim([0 10])
% % 
% % subplot(122)
% % x = [ones(size(Pre.vHPC.aversive(:,2))) ; ones(size(Post.vHPC.reward(:,2)))*2];
% % DIA = (Post.vHPC.aversive(:,2) ./ Pre.vHPC.aversive(:,2));% ./ (Post.vHPC.aversive(:,2) + Pre.vHPC.aversive(:,2));
% % DIR = (Post.vHPC.reward(:,2) ./ Pre.vHPC.reward(:,2));% ./ (Post.vHPC.reward(:,2) + Pre.vHPC.reward(:,2));
% % y = [DIA ; DIR];
% % scatter(x,y,"filled",'jitter','on', 'jitterAmount',0.1'),hold on,xlim([0 3])
% % scatter([1 2] , [nanmedian(DIA) nanmedian(DIR)],'filled'),set(gca, 'YScale', 'log')
% % yline(1,'--')
% % signrank(DIA,1)
% % signrank(DIR,1)
% % % ylim([0 10])
% % 
% % 
% % 
% % 
% % 
% % 
% % figure,
% % subplot(121)
% % x = [ones(size(Pre.aversive(:,2))) ; ones(size(Post.aversive(:,2)))*2];
% % DIA = (Pre.aversive(:,1));% ./ (Post.dHPC.aversive(:,2) + Pre.dHPC.aversive(:,2));
% % DIR = (Post.aversive(:,1));% ./ (Post.dHPC.reward(:,2) + Pre.dHPC.reward(:,2));
% % y = [DIA ; DIR];
% % scatter(x,y,"filled",'jitter','on', 'jitterAmount',0.1'),hold on,xlim([0 3])
% % scatter([1 2] , [nanmedian(DIA) nanmedian(DIR)],'filled')%,set(gca, 'YScale', 'log')
% % ranksum(DIA,DIR)
% % 
% % 
% % subplot(122)
% % x = [ones(size(Pre.reward(:,2))) ; ones(size(Post.reward(:,2)))*2];
% % DIA = (Pre.reward(:,1));% ./ (Post.dHPC.aversive(:,2) + Pre.dHPC.aversive(:,2));
% % DIR = (Post.reward(:,1));% ./ (Post.dHPC.reward(:,2) + Pre.dHPC.reward(:,2));
% % y = [DIA ; DIR];
% % scatter(x,y,"filled",'jitter','on', 'jitterAmount',0.1'),hold on,xlim([0 3])
% % scatter([1 2] , [nanmedian(DIA) nanmedian(DIR)],'filled')%,set(gca, 'YScale', 'log')
% % ranksum(DIA,DIR)
% 
% figure
result = {};
result1 = nan(200,50);
for i = 1 : size(Peaks.aversive.Post,2)
    tmp = Peaks.aversive.Post{i};
    [tmp,bins]=binspikes(tmp,1/60,[0 round(tmp(end))]);
    result1(1:size(tmp,1),i) = Smooth(tmp./60,2);
    if isempty(result)
%         plot(tmp./60,'k'),hold on
        result{1} = [bins' tmp./60]; clear tmp bins
    else
%         plot(tmp./60,'k'),hold on
        result{end+1} = [bins' tmp./60]; clear tmp bins
    end
end

figure
plot(nanmean(result1'))




% figure
result = {};
result2 = nan(200,50);
for i = 1 : size(Peaks.reward.Post,2)
    tmp = Peaks.reward.Post{i};
    [tmp,bins]=binspikes(tmp,1/60,[0 round(tmp(end))]);
    result2(1:size(tmp,1),i) = Smooth(tmp./60,2);
    if isempty(result)
%         plot(tmp./60,'k'),hold on
        result{1} = [bins' tmp./60]; clear tmp bins
    else
%         plot(tmp./60,'k'),hold on
        result{end+1} = [bins' tmp./60]; clear tmp bins
    end
end

figure
plot([1:200],nanmean(result2'),'b'),hold on
ciplot(nanmean(result2')-nansem(result2') , nanmean(result2')+nansem(result2') , [1:200] , 'b') , alpha 0.5
plot(nanmean(result1'),'r'),hold on
ciplot(nanmean(result1')-nansem(result1') , nanmean(result1')+nansem(result1') , [1:200] , 'r') , alpha 0.5


end


