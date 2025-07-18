clear
clc
% close all

%% Parameters
path = {'E:\Rat126\Ephys\in_Pyr';'E:\Rat103\usable';'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat132\recordings\in_pyr';'E:\Rat165\in_pyr\'};%List of folders from the path

% for SU
criteria_fr = 0; %criteria to include or not a SU into the analysis
criteria_n = [3 3]; % minimal number of neurons from each structure [vHPC dHPC]
criteria_type = 0; %criteria for celltype (0:pyr, 1:int, 2:all)
binSize = 0.025; %for qssemblie detection qnd qxctivity strength
normalization = true; % to define if normalization over Reactivation Strength is applied or not
th = 5; % threshold for peak detection
iterations = 5; % number of iterations for downsampling

% Storage variables
% For Reactivation measure
% Coordinated
reactivation.aversive.dvHPC.coordinated = []; reactivation.reward.dvHPC.coordinated = [];
reactivation.aversive.dHPC.coordinated = []; reactivation.reward.dHPC.coordinated = [];
reactivation.aversive.vHPC.coordinated = []; reactivation.reward.vHPC.coordinated = [];
SAVED.aversive = {};
count = 0;

% Uncooridnated
reactivation.aversive.dvHPC.uncoordinated = []; reactivation.reward.dvHPC.uncoordinated = [];
reactivation.aversive.dHPC.uncoordinated = []; reactivation.reward.dHPC.uncoordinated = [];
reactivation.aversive.vHPC.uncoordinated = []; reactivation.reward.vHPC.uncoordinated = [];

% General metrics
percentages = [];
Number_of_assemblies.aversive = [];
Number_of_assemblies.reward = [];

%% Main loop, to iterate across sessions
for tt = 2:length(path)
    %List of folders from the path
    files = dir(path{tt});
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
    clear files dirFlags
    num_assembliesR = [];
    num_assembliesA = [];
    for t = 1 : length(subFolders)-2
        disp(['-- Initiating analysis of folder #' , num2str(t) , ' from rat #',num2str(tt) , ' --'])
        session = [subFolders(t+2).folder,'\',subFolders(t+2).name];
        cd(session)
        
        %Loading TS of the sessions
        disp('Uploading session data')
        load('session_organization.mat')
        
        % Awake
        disp('Uploading Behavioral Data')
        load('behavioral_data.mat','movement')
        
        %% load sleep states
        disp('Uploading sleep scoring')
        x = dir([cd,'\*-states.mat']);    states = load([cd,'\',x.name]);    states = states.states;
        REM.all = ToIntervals(states==5);    NREM.all = ToIntervals(states==3);    WAKE.all = ToIntervals(states==1);
        %         REM.all(REM.all(:,2)-REM.all(:,1)<60,:) = [];
        clear x states
        
        %         NREM.all(NREM.all(:,2)-NREM.all(:,1)<60*time_criteria,:)=[];
        NREM.baseline = Restrict(NREM.all,baselineTS./1000);
        NREM.aversive = Restrict(NREM.all,aversiveTS./1000);
        NREM.reward = Restrict(NREM.all,rewardTS./1000);
        
        REM.baseline = Restrict(REM.all,baselineTS./1000);
        REM.aversive = Restrict(REM.all,aversiveTS./1000);
        REM.reward = Restrict(REM.all,rewardTS./1000);
        
        WAKE.baseline = Restrict(WAKE.all,baselineTS./1000);
        WAKE.aversive = Restrict(WAKE.all,aversiveTS./1000);
        WAKE.reward = Restrict(WAKE.all,rewardTS./1000);
        
        
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
        
        if and(RV, RD)
            % coordination
            coordinated = [];
            coordinatedV = [];
            coordinatedV_refined = [];
            cooridnated_event = [];
            cooridnated_eventDV = [];
            cooridnated_eventVD = [];
            for i = 1:length(ripplesD)
                r = ripplesD(i,:);
                tmp = sum(and(ripplesV(:,2)>= r(1,2)-0.1, ripplesV(:,2)<= r(1,2)+0.1));
                if tmp>0
                    z = ripplesV(and(ripplesV(:,2)>= r(1,2)-0.1, ripplesV(:,2)<= r(1,2)+0.1),:);
                    coordinatedV = [coordinatedV ; z];
                    [p,indice] = min(abs(r(2)-z(:,2)));
                    coordinatedV_refined = [coordinatedV_refined ; z(indice,:)];
                    coordinated = [coordinated ; r];
                    
                    cooridnated_event = [cooridnated_event ; min([r(1) , z(indice,1)]) , min(min(z(indice,2),r(2)))+abs(z(indice,2)-r(2))/2 , max([r(3) , z(indice,3)])];
                    
                    if r(2)<z(indice,2)
                        cooridnated_eventDV = [cooridnated_eventDV ; min([r(1) , z(indice,1)]) , min(min(z(indice,2),r(2)))+abs(z(indice,2)-r(2))/2 , max([r(3) , z(indice,3)])];
                    else
                        cooridnated_eventVD = [cooridnated_eventVD ; min([r(1) , z(indice,1)]) , min(min(z(indice,2),r(2)))+abs(z(indice,2)-r(2))/2 , max([r(3) , z(indice,3)])];
                    end
                    
                    clear tmp2 tmp1 p indice z
                end
                clear r
            end
            clear x tmp i
            
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
            %eliminate ripples <50 ms
%             ri = (ripples.dHPC.uncoordinated.all(:,3) - ripples.dHPC.uncoordinated.all(:,1))<0.05;
%             ripples.dHPC.uncoordinated.all(ri,:) = []; clear ri
            ripples.dHPC.coordinated.baseline = Restrict(coordinated , NREM.baseline);
            ripples.dHPC.coordinated.reward = Restrict(coordinated , NREM.reward);
            ripples.dHPC.coordinated.aversive = Restrict(coordinated , NREM.aversive);
            % coordinated vRipples
            ripples.vHPC.coordinated.all = coordinatedV_refined;
            ripples.vHPC.uncoordinated.all = ripplesV(not(ismember(ripplesV(:,2) , coordinatedV_refined(:,2))),:);
            %eliminate ripples <50 ms
%             ri = (ripples.vHPC.uncoordinated.all(:,3) - ripples.vHPC.uncoordinated.all(:,1))<0.05;
%             ripples.vHPC.uncoordinated.all(ri,:) = []; clear ri
            ripples.vHPC.coordinated.baseline = Restrict(coordinatedV_refined , NREM.baseline);
            ripples.vHPC.coordinated.reward = Restrict(coordinatedV_refined , NREM.reward);
            ripples.vHPC.coordinated.aversive = Restrict(coordinatedV_refined , NREM.aversive);
            % ALL Uncoordinated
            % downsampling to have same amount of ripples
            if size(ripples.dHPC.uncoordinated.all,1)>size(ripples.dHPC.coordinated.all,1)
                ran = randperm(size(ripples.dHPC.uncoordinated.all,1));
                ran = ripples.dHPC.uncoordinated.all(ran,:);
                ran = ran(1 : size(ripples.dHPC.coordinated.all,1) , :);
            else
                ran = ripples.dHPC.uncoordinated.all;
            end
            
            if size(ripples.vHPC.uncoordinated.all,1)>size(ripples.vHPC.coordinated.all,1)
                ran1 = randperm(size(ripples.vHPC.uncoordinated.all,1));
                ran1 = ripples.vHPC.uncoordinated.all(ran1,:);
                ran1 = ran1(1 : size(ripples.vHPC.coordinated.all,1) , :);
            else
                ran1 = ripples.vHPC.uncoordinated.all;
            end
            
            u = [ran ; ran1]; clear ran ran1
            
%             u = [ripples.dHPC.uncoordinated.all ripples.dHPC.uncoordinated.all ; ripples.vHPC.uncoordinated.all ripples.vHPC.uncoordinated.all];
            
            [B I] = sort(u(:,1));
            u = ConsolidateIntervals([u(I,1) u(I,3)]);
            
            ripple_event.uncoordinated = u; clear u
            %coordinated event
            cooridnated_event((cooridnated_event(:,3)-cooridnated_event(:,1)<0.04),:) = [];
            ripple_event.baseline = Restrict(cooridnated_event,baselineTS./1000);
            ripple_event.reward = Restrict(cooridnated_event,rewardTS./1000);
            ripple_event.aversive = Restrict(cooridnated_event,aversiveTS./1000);
            ripple_event.all = cooridnated_event;
        elseif RD
            % Store events time stamps
            % dRipples
            ripples.dHPC.baseline = Restrict(ripplesD , NREM.baseline);
            ripples.dHPC.reward = Restrict(ripplesD , NREM.reward);
            ripples.dHPC.aversive = Restrict(ripplesD , NREM.aversive);
        elseif RV
            % vRipples
            ripples.vHPC.baseline = Restrict(ripplesV , NREM.baseline);
            ripples.vHPC.reward = Restrict(ripplesV , NREM.reward);
            ripples.vHPC.aversive = Restrict(ripplesV , NREM.aversive);
        end
        
        %% Spikes
        % Load Units
        disp('Uploading Spiking activity')
        cd 'Spikesorting'
        [clusters , numberD , numberV , spks , spks_dHPC , spks_vHPC , cellulartype] = load_SU_FM(cd,criteria_type,criteria_fr,aversiveTS_run./1000,rewardTS_run./1000);
        
        %% Assemblies detection
        if or(numberD > 3 , numberV > 3)
            % --- Aversive ---
            disp('Lets go for the assemblies')
            if isfile('dorsalventral_assemblies_aversive_quiet_epochs.mat')
                disp('Loading Aversive template')
                load('dorsalventral_assemblies_aversive_quiet_epochs.mat')
            end
            
            Thresholded.aversive.all = Th;
            patterns.all.aversive = pat;
            clear cond Th pat
            
            patterns.all.aversive = patterns.all.aversive .* Thresholded.aversive.all;
            
            % Detection of members
            if not(isempty(Thresholded.aversive.all))
                if numberD>0
                    cond1 =  sum(Thresholded.aversive.all(1:size(clusters.dHPC,1),:),1)>0; %checking of dHPC SU
                    cond2 =  sum(Thresholded.aversive.all(size(clusters.dHPC,1)+1:end,:),1)>0; %checking of vHPC SU
                    cond.dHPC.aversive = and(cond1 , not(cond2));
                    cond.vHPC.aversive = and(cond2 , not(cond1));
                    cond.both.aversive = and(cond1 , cond2); clear cond1 cond2
                else
                    cond1 =  logical(zeros(1,size(Thresholded.aversive.all,2))); %checking of dHPC SU
                    cond2 =  logical(ones(1,size(Thresholded.aversive.all,2))); %checking of vHPC SU
                    cond.dHPC.aversive = and(cond1 , not(cond2));
                    cond.vHPC.aversive = and(cond2 , not(cond1));
                    cond.both.aversive = and(cond1 , cond2); clear cond1 cond2
                end
            else
                cond1 =  false; %checking of dHPC SU
                cond2 =  logical(0); %checking of vHPC SU
                cond.dHPC.aversive = and(cond1 , not(cond2));
                cond.vHPC.aversive = and(cond2 , not(cond1));
                cond.both.aversive = and(cond1 , cond2); clear cond1 cond2
            end
            num_assembliesA = [num_assembliesA ; sum(cond.both.aversive) sum(cond.dHPC.aversive) sum(cond.vHPC.aversive)];
            
            % --- Reward ---
            disp('Loading Reward template')
            if isfile('dorsalventral_assemblies_reward_quiet_epochs.mat')
                load('dorsalventral_assemblies_reward_quiet_epochs.mat')
            end
            
            Thresholded.reward.all = Th;
            patterns.all.reward = pat;
            clear Th pat
            
            patterns.all.reward = patterns.all.reward .* Thresholded.reward.all;
            
            % Detection of members using
            if not(isempty(Thresholded.reward.all))
                if numberD>0
                    cond1 =  sum(Thresholded.reward.all(1:size(clusters.dHPC,1),:),1)>0; %checking of dHPC SU
                    cond2 =  sum(Thresholded.reward.all(size(clusters.dHPC,1)+1:end,:),1)>0; %checking of vHPC SU
                    cond.dHPC.reward = and(cond1 , not(cond2));
                    cond.vHPC.reward = and(cond2 , not(cond1));
                    cond.both.reward = and(cond1 , cond2); clear cond1 cond2
                else
                    cond1 =  logical(zeros(1,size(Thresholded.reward.all,2))); %checking of dHPC SU
                    cond2 =  logical(ones(1,size(Thresholded.reward.all,2))); %checking of vHPC SU
                    cond.dHPC.reward = and(cond1 , not(cond2));
                    cond.vHPC.reward = and(cond2 , not(cond1));
                    cond.both.reward = and(cond1 , cond2); clear cond1 cond2
                end
            else
                cond1 =  logical(0); %checking of dHPC SU
                cond2 =  logical(0); %checking of vHPC SU
                cond.dHPC.reward = and(cond1 , not(cond2));
                cond.vHPC.reward = and(cond2 , not(cond1));
                cond.both.reward = and(cond1 , cond2); clear cond1 cond2
            end
            num_assembliesR = [num_assembliesR ; sum(cond.both.reward) sum(cond.dHPC.reward) sum(cond.vHPC.reward)];
            
            %% SpikeTrains construction
            limits = [0 segments.Var1(end)/1000];
            events = [];
            [Spikes , bins , Clusters] = spike_train_construction([spks_dHPC;spks_vHPC], clusters.all, cellulartype, binSize, limits, events, false, true);
            clear limits events
        
            % NREM sleep
            if and(RD,RV)
                % coordinated event
                Coor.is.sws.baseline = InIntervals(bins,[ripple_event.baseline(:,1) ripple_event.baseline(:,3)]);
                Coor.is.sws.reward = InIntervals(bins,[ripple_event.reward(:,1) ripple_event.reward(:,3)]);
                Coor.is.sws.aversive = InIntervals(bins,[ripple_event.aversive(:,1) ripple_event.aversive(:,3)]);
                
                Coor.is.sws.timestamps.sleep.aversive = NREM.aversive;
                Coor.is.sws.timestamps.sleep.reward = NREM.reward;
                Coor.is.sws.timestamps.sleep.baseline = NREM.baseline;
                
                Coor.is.sws.runaversive = InIntervals(bins,movement.aversive);
                Coor.is.sws.runreward = InIntervals(bins,movement.reward);
                
                Coor.is.sws.timestamps.run.aversive = InIntervals(bins,aversiveTS_run./1000);
                Coor.is.sws.timestamps.run.reward = InIntervals(bins,rewardTS_run./1000);
                
                Coor.is.sws.timestamps.aversiveSleep = aversiveTS./1000;
                Coor.is.sws.timestamps.aversiveRun = aversiveTS_run./1000;
                Coor.is.sws.timestamps.rewardSleep = rewardTS./1000;
                Coor.is.sws.timestamps.rewardRun = rewardTS_run./1000;
                Coor.is.sws.timestamps.baselineSleep = baselineTS./1000;
                
                % Uncoordinated
                % Baseline
                c = ripple_event.uncoordinated;
                UnCoor.is.sws.baseline = InIntervals(bins,Restrict(c,NREM.baseline)); clear tmp
                
                % Reward
                UnCoor.is.sws.reward = InIntervals(bins,Restrict(c,NREM.reward)); clear tmp
                
                % Aversive
                UnCoor.is.sws.aversive = InIntervals(bins,Restrict(c,NREM.aversive)); clear tmp
                
                UnCoor.is.sws.timestamps.sleep.aversive = NREM.aversive;
                UnCoor.is.sws.timestamps.sleep.reward = NREM.reward;
                UnCoor.is.sws.timestamps.sleep.baseline = NREM.baseline;
                
                UnCoor.is.sws.runaversive = InIntervals(bins,movement.aversive);
                UnCoor.is.sws.runreward = InIntervals(bins,movement.reward);
                
                UnCoor.is.sws.timestamps.run.aversive = InIntervals(bins,aversiveTS_run./1000);
                UnCoor.is.sws.timestamps.run.reward = InIntervals(bins,rewardTS_run./1000);
                
                UnCoor.is.sws.timestamps.aversiveSleep = aversiveTS./1000;
                UnCoor.is.sws.timestamps.aversiveRun = aversiveTS_run./1000;
                UnCoor.is.sws.timestamps.rewardSleep = rewardTS./1000;
                UnCoor.is.sws.timestamps.rewardRun = rewardTS_run./1000;
                UnCoor.is.sws.timestamps.baselineSleep = baselineTS./1000;
            end
            
            %% Reactivation Strenght
            % Joint Aversive Assemblies
            if sum(cond.both.aversive)>=1
                [R] = reactivation_strength(patterns.all.aversive , cond.both.aversive , [bins' , Spikes] , Coor.is.sws , th , 'A' , config , normalization , []); clear templates
                reactivation.aversive.dvHPC.coordinated = [reactivation.aversive.dvHPC.coordinated ; R];
                
                for i = 1 : size(R,1)
                    count = count + 1;
                    x = R(i,1);
                    y = [tt t];
                    z = Thresholded.aversive.all(:,cond.both.aversive); z = z(:,i);
                    z = clusters.all(z);
                    z = struct('Reactivation',x,'Session',y,'ids',z);
                    SAVED.aversive{count} = z; clear x y z
                end
                
%                 R = Shuffle_and_Reactivation(ripple_event.uncoordinated, size(ripple_event.all,1), iterations, UnCoor.is.sws, patterns.all.aversive, cond.both.aversive, th , [bins' , Spikes], 'A' , config, normalization);
                [R] = reactivation_strength(patterns.all.aversive , cond.both.aversive , [bins' , Spikes] , UnCoor.is.sws , th , 'A' , config , normalization , []); clear templates
                reactivation.aversive.dvHPC.uncoordinated = [reactivation.aversive.dvHPC.uncoordinated ; R];
            end
            
            %             if sum(cond.dHPC.aversive)>=1
            %                 [R] = reactivation_strength(patterns.all.aversive , cond.dHPC.aversive , [bins' , Spikes] , Coor.is.sws , th , 'A' , config , normalization , []); clear templates
            %                 reactivation.aversive.dHPC.coordinated = [reactivation.aversive.dHPC.coordinated ; R];
            %
            %                 R = Shuffle_and_Reactivation(ripple_event.uncoordinated, size(ripple_event.all,1), iterations, UnCoor.is.sws, patterns.all.aversive, cond.dHPC.aversive, th , [bins' , Spikes], 'A' , config, normalization);
            % %                 [R] = reactivation_strength(patterns.all.aversive , cond.dHPC.aversive , [bins' , Spikes] , UnCoor.is.sws , th , 'A' , config , normalization , []); clear templates
            %                 reactivation.aversive.dHPC.uncoordinated = [reactivation.aversive.dHPC.uncoordinated ; R];
            %             end
            %
            %             if sum(cond.vHPC.aversive)>=1
            %                 [R] = reactivation_strength(patterns.all.aversive , cond.vHPC.aversive , [bins' , Spikes] , Coor.is.sws , th , 'A' , config , normalization , []); clear templates
            %                 reactivation.aversive.vHPC.coordinated = [reactivation.aversive.vHPC.coordinated ; R];
            %
            %                 R = Shuffle_and_Reactivation(ripple_event.uncoordinated, size(ripple_event.all,1), iterations, UnCoor.is.sws, patterns.all.aversive, cond.vHPC.aversive, th , [bins' , Spikes], 'A' , config, normalization);
            % %                 [R] = reactivation_strength(patterns.all.aversive , cond.vHPC.aversive , [bins' , Spikes] , UnCoor.is.sws , th , 'A' , config , normalization , []); clear templates
            %                 reactivation.aversive.vHPC.uncoordinated = [reactivation.aversive.vHPC.uncoordinated ; R];
            %             end
            
            %% Same for Reward assemblies
            % Joint Reward Assemblies
            if sum(cond.both.reward)>=1
                [R] = reactivation_strength(patterns.all.reward , cond.both.reward , [bins' , Spikes] , Coor.is.sws , th , 'R' , config , normalization , []); clear templates
                reactivation.reward.dvHPC.coordinated = [reactivation.reward.dvHPC.coordinated ; R];
                
                %                 R = Shuffle_and_Reactivation(ripple_event.uncoordinated, size(ripple_event.all,1), iterations, UnCoor.is.sws, patterns.all.reward, cond.both.reward, th , [bins' , Spikes], 'R' , config, normalization);
                [R] = reactivation_strength(patterns.all.reward , cond.both.reward , [bins' , Spikes] , UnCoor.is.sws , th , 'R' , config , normalization , []); clear templates
                reactivation.reward.dvHPC.uncoordinated = [reactivation.reward.dvHPC.uncoordinated ; R];
            end
            
            %             if sum(cond.dHPC.reward)>=1
            %                 [R] = reactivation_strength(patterns.all.reward , cond.dHPC.reward , [bins' , Spikes] , Coor.is.sws , th , 'R' , config , normalization , []); clear templates
            %                 reactivation.reward.dHPC.coordinated = [reactivation.reward.dHPC.coordinated ; R];
            %
            %                 R = Shuffle_and_Reactivation(ripple_event.uncoordinated, size(ripple_event.all,1), iterations, UnCoor.is.sws, patterns.all.reward, cond.dHPC.reward, th , [bins' , Spikes], 'R' , config, normalization);
            % %                 [R] = reactivation_strength(patterns.all.reward , cond.dHPC.reward , [bins' , Spikes] , UnCoor.is.sws , th , 'R' , config , normalization , []); clear templates
            %                 reactivation.reward.dHPC.uncoordinated = [reactivation.reward.dHPC.uncoordinated ; R];
            %             end
            %
            %             if sum(cond.vHPC.reward)>=1
            %                 [R] = reactivation_strength(patterns.all.reward , cond.vHPC.reward , [bins' , Spikes] , Coor.is.sws , th , 'R' , config , normalization , []); clear templates
            %                 reactivation.reward.vHPC.coordinated = [reactivation.reward.vHPC.coordinated ; R];
            %
            %                 R = Shuffle_and_Reactivation(ripple_event.uncoordinated, size(ripple_event.all,1), iterations, UnCoor.is.sws, patterns.all.reward, cond.vHPC.reward, th , [bins' , Spikes], 'R' , config, normalization);
            % %                 [R] = reactivation_strength(patterns.all.reward , cond.vHPC.reward , [bins' , Spikes] , UnCoor.is.sws , th , 'R' , config , normalization , []); clear templates
            %                 reactivation.reward.vHPC.uncoordinated = [reactivation.reward.vHPC.uncoordinated ; R];
            %             end
            
            
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
        clear cooridnated_eventDV cooridnated_eventVD segments
        clear RBA RBR RDA RDR RVA RVR Shocks_filt Rewards_filt
    end
    
    Number_of_assemblies.aversive = [Number_of_assemblies.aversive ; sum(num_assembliesA)];
    Number_of_assemblies.reward = [Number_of_assemblies.reward ; sum(num_assembliesR)];
    clear num_assembliesA num_assembliesR
end

% save([cd,'\Reactivation_Strength_Data_Normalized_eliminating_withingD_5SD.mat'],'reactivation')

%% Peaks mean for joint assemblies
figure
x = reactivation.reward.dvHPC.coordinated(:,1);
y = reactivation.aversive.dvHPC.coordinated(:,1);

kstest(x)
kstest(y)
[h, p] = ranksum(x,y,'tail','left')  
[h, p] = signrank(y,0,'tail','right')
[h, p] = signrank(x,0,'tail','right')

subplot(121),
grps = [ones(size(x,1),1) ; ones(size(y,1),1)*2];
Y = [x;y];
scatter(grps,Y,"filled",'jitter','on', 'jitterAmount',0.1),hold on
scatter([1 2] , [nanmean(x) nanmean(y)],'filled'),xlim([0 3]),ylim([-1 1])


x1 = reactivation.reward.dvHPC.uncoordinated(:,1);
y1 = reactivation.aversive.dvHPC.uncoordinated(:,1);

kstest(x)
kstest(y)
[h, p] = ranksum(x1,y1,'tail','left')  
[h, p] = signrank(y1,0,'tail','right')
[h, p] = signrank(x1,0,'tail','right')

subplot(122),
grps = [ones(size(x1,1),1) ; ones(size(y1,1),1)*2];
Y = [x1;y1];
scatter(grps,Y,"filled",'jitter','on', 'jitterAmount',0.1),hold on
scatter([1 2] , [nanmean(x1) nanmean(y1)],'filled'),xlim([0 3]),ylim([-1 1])

%% ANOVA
x1(isnan(x)) = [];
y1(isnan(y)) = [];
x(isnan(x)) = [];
y(isnan(y)) = [];

x(isnan(x1)) = [];
y(isnan(y1)) = [];
x1(isnan(x1)) = [];
y1(isnan(y1)) = [];




data = [x ; y ; x1 ; y1];
ref = [ones(length(x),1) ; ones(length(y),1) ; ones(length(x1),1)*2 ; ones(length(y1),1)*2];
ref = [ref , [ones(length(x),1) ; ones(length(y),1)*2 ; ones(length(x1),1) ; ones(length(y1),1)*2]];



perform2WayANOVA(data, ref(:,1), ref(:,2))
performFriedmanTest(data, ref(:,2), ref(:,1))
perform2WayPermutationTest(data, ref(:,2), ref(:,1))

[~,~,stats] = friedman(data,ref);

[~,~,stats] = anova2(data,2);
% [~,~,stats] = anovan(data,{ref(:,1) ref(:,2)},"Model","interaction","Varnames",["ripple","type"]);
[results,~,~,gnames] = multcompare(stats,"Dimension",[1 2]);

%  for dHPC assemblies
x = reactivation.reward.dHPC(:,1);
y = reactivation.aversive.dHPC(:,1);
kstest(x)
kstest(y)
[h, p] = ranksum(x,y,'tail','left')  
[h, p] = signrank(y,0,'tail','right')
[h, p] = signrank(x,0,'tail','right')

subplot(132),
grps = [ones(size(x,1),1) ; ones(size(y,1),1)*2];
Y = [x;y];
scatter(grps,Y,"filled",'jitter','on', 'jitterAmount',0.1),hold on
scatter([1 2] , [nanmean(x) nanmean(y)],'filled'),xlim([0 3]),ylim([-1 1])

%  for vHPC assemblies
x = reactivation.reward.vHPC(:,1);
y = reactivation.aversive.vHPC(:,1);
kstest(x)
kstest(y)
[h, p] = ranksum(x,y,'tail','left')  
[h, p] = signrank(y,0,'tail','right')
[h, p] = signrank(x,0,'tail','right')

subplot(133),
grps = [ones(size(x,1),1) ; ones(size(y,1),1)*2];
Y = [x;y];
scatter(grps,Y,"filled",'jitter','on', 'jitterAmount',0.1),hold on
scatter([1 2] , [nanmean(x) nanmean(y)],'filled'),xlim([0 3]),ylim([-1 1])