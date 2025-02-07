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

% Behavior
minimal_speed = 7; % minimal speed to detect quite periods
minimal_speed_time = 2; % minimal time to detect quite periods

% Storage variables
% For Reactivation measure
% Coordinated
reactivation.aversive.dvHPC.coordinated = []; reactivation.reward.dvHPC.coordinated = [];
reactivation.aversive.dHPC.coordinated = []; reactivation.reward.dHPC.coordinated = [];
reactivation.aversive.vHPC.coordinated = []; reactivation.reward.vHPC.coordinated = [];

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
        load('behavioral_data.mat')
        
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
            ripples.dHPC.coordinated.baseline = Restrict(coordinated , NREM.baseline);
            ripples.dHPC.coordinated.reward = Restrict(coordinated , NREM.reward);
            ripples.dHPC.coordinated.aversive = Restrict(coordinated , NREM.aversive);
            % coordinated vRipples
            ripples.vHPC.coordinated.all = coordinatedV_refined;
            ripples.vHPC.uncoordinated.all = ripplesV(not(ismember(ripplesV(:,2) , coordinatedV_refined(:,2))),:);
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
            if isfile('dorsalventral_assemblies_aversiveVF.mat')
                disp('Loading Aversive template')
                load('dorsalventral_assemblies_aversiveVF.mat')
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
            if isfile('dorsalventral_assemblies_rewardVF.mat')
                load('dorsalventral_assemblies_rewardVF.mat')
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
                tmp = not(ismember(ripples.dHPC.baseline(:,2) , ripples.dHPC.coordinated.baseline(:,2)));
                ripples.dHPC.uncoordinated.baseline = [ripples.dHPC.baseline(tmp,1) ripples.dHPC.baseline(tmp,3)];
                tmp = [ripples.dHPC.baseline(tmp,1) ripples.dHPC.baseline(tmp,3)];
                is.sws.baseline = InIntervals(bins,tmp); clear tmp
                
                tmp = not(ismember(ripples.vHPC.baseline(:,2) , ripples.vHPC.coordinated.baseline(:,2)));
                ripples.vHPC.uncoordinated.baseline = [ripples.vHPC.baseline(tmp,1) ripples.vHPC.baseline(tmp,3)];
                tmp = [ripples.vHPC.baseline(tmp,1) ripples.vHPC.baseline(tmp,3)];
                UnCoor.is.sws.baseline = (is.sws.baseline + InIntervals(bins,tmp)) >= 1; clear tmp
                
                % Reward
                tmp = not(ismember(ripples.dHPC.reward(:,2) , ripples.dHPC.coordinated.reward(:,2)));
                ripples.dHPC.uncoordinated.reward = [ripples.dHPC.reward(tmp,1) ripples.dHPC.reward(tmp,3)];
                tmp = [ripples.dHPC.reward(tmp,1) ripples.dHPC.reward(tmp,3)];
                is.sws.reward = InIntervals(bins,tmp); clear tmp
                
                tmp = not(ismember(ripples.vHPC.reward(:,2) , ripples.vHPC.coordinated.reward(:,2)));
                ripples.vHPC.uncoordinated.reward = [ripples.vHPC.reward(tmp,1) ripples.vHPC.reward(tmp,3)];
                tmp = [ripples.vHPC.reward(tmp,1) ripples.vHPC.reward(tmp,3)];
                UnCoor.is.sws.reward = (is.sws.reward + InIntervals(bins,tmp)) >= 1; clear tmp
                
                % Aversive
                tmp = not(ismember(ripples.dHPC.aversive(:,2) , ripples.dHPC.coordinated.aversive(:,2)));
                ripples.dHPC.uncoordinated.aversive = [ripples.dHPC.aversive(tmp,1) ripples.dHPC.aversive(tmp,3)];
                tmp = [ripples.dHPC.aversive(tmp,1) ripples.dHPC.aversive(tmp,3)];
                is.sws.aversive = InIntervals(bins,tmp); clear tmp
                
                tmp = not(ismember(ripples.vHPC.aversive(:,2) , ripples.vHPC.coordinated.aversive(:,2)));
                ripples.vHPC.uncoordinated.aversive = [ripples.vHPC.aversive(tmp,1) ripples.vHPC.aversive(tmp,3)];
                tmp = [ripples.vHPC.aversive(tmp,2)-0.2 ripples.vHPC.aversive(tmp,2)+0.2];
                UnCoor.is.sws.aversive = (is.sws.aversive + InIntervals(bins,tmp)) >= 1; clear tmp
                
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
                
                R = Shuffle_and_Reactivation(ripple_event.uncoordinated, size(ripple_event.all,1), iterations, UnCoor.is.sws, patterns.all.aversive, cond.both.aversive, th , [bins' , Spikes], 'A' , config, normalization);
%                 [R] = reactivation_strength(patterns.all.aversive , cond.both.aversive , [bins' , Spikes] , UnCoor.is.sws , th , 'A' , config , normalization , []); clear templates
                reactivation.aversive.dvHPC.uncoordinated = [reactivation.aversive.dvHPC.uncoordinated ; R];
            end
            
            if sum(cond.dHPC.aversive)>=1
                [R] = reactivation_strength(patterns.all.aversive , cond.dHPC.aversive , [bins' , Spikes] , Coor.is.sws , th , 'A' , config , normalization , []); clear templates
                reactivation.aversive.dHPC.coordinated = [reactivation.aversive.dHPC.coordinated ; R];
                
                R = Shuffle_and_Reactivation(ripple_event.uncoordinated, size(ripple_event.all,1), iterations, UnCoor.is.sws, patterns.all.aversive, cond.dHPC.aversive, th , [bins' , Spikes], 'A' , config, normalization);
%                 [R] = reactivation_strength(patterns.all.aversive , cond.dHPC.aversive , [bins' , Spikes] , UnCoor.is.sws , th , 'A' , config , normalization , []); clear templates
                reactivation.aversive.dHPC.uncoordinated = [reactivation.aversive.dHPC.uncoordinated ; R];
            end
            
            if sum(cond.vHPC.aversive)>=1
                [R] = reactivation_strength(patterns.all.aversive , cond.vHPC.aversive , [bins' , Spikes] , Coor.is.sws , th , 'A' , config , normalization , []); clear templates
                reactivation.aversive.vHPC.coordinated = [reactivation.aversive.vHPC.coordinated ; R];
                
                R = Shuffle_and_Reactivation(ripple_event.uncoordinated, size(ripple_event.all,1), iterations, UnCoor.is.sws, patterns.all.aversive, cond.vHPC.aversive, th , [bins' , Spikes], 'A' , config, normalization);
%                 [R] = reactivation_strength(patterns.all.aversive , cond.vHPC.aversive , [bins' , Spikes] , UnCoor.is.sws , th , 'A' , config , normalization , []); clear templates
                reactivation.aversive.vHPC.uncoordinated = [reactivation.aversive.vHPC.uncoordinated ; R];
            end
            
            %% Same for Reward assemblies
            % Joint Reward Assemblies
            if sum(cond.both.reward)>=1
                [R] = reactivation_strength(patterns.all.reward , cond.both.reward , [bins' , Spikes] , Coor.is.sws , th , 'R' , config , normalization , []); clear templates
                reactivation.reward.dvHPC.coordinated = [reactivation.reward.dvHPC.coordinated ; R];
                
                R = Shuffle_and_Reactivation(ripple_event.uncoordinated, size(ripple_event.all,1), iterations, UnCoor.is.sws, patterns.all.reward, cond.both.reward, th , [bins' , Spikes], 'R' , config, normalization);
%                 [R] = reactivation_strength(patterns.all.reward , cond.both.reward , [bins' , Spikes] , UnCoor.is.sws , th , 'R' , config , normalization , []); clear templates
                reactivation.reward.dvHPC.uncoordinated = [reactivation.reward.dvHPC.uncoordinated ; R];
            end
            
            if sum(cond.dHPC.reward)>=1
                [R] = reactivation_strength(patterns.all.reward , cond.dHPC.reward , [bins' , Spikes] , Coor.is.sws , th , 'R' , config , normalization , []); clear templates
                reactivation.reward.dHPC.coordinated = [reactivation.reward.dHPC.coordinated ; R];
                
                R = Shuffle_and_Reactivation(ripple_event.uncoordinated, size(ripple_event.all,1), iterations, UnCoor.is.sws, patterns.all.reward, cond.dHPC.reward, th , [bins' , Spikes], 'R' , config, normalization);
%                 [R] = reactivation_strength(patterns.all.reward , cond.dHPC.reward , [bins' , Spikes] , UnCoor.is.sws , th , 'R' , config , normalization , []); clear templates
                reactivation.reward.dHPC.uncoordinated = [reactivation.reward.dHPC.uncoordinated ; R];
            end
            
            if sum(cond.vHPC.reward)>=1
                [R] = reactivation_strength(patterns.all.reward , cond.vHPC.reward , [bins' , Spikes] , Coor.is.sws , th , 'R' , config , normalization , []); clear templates
                reactivation.reward.vHPC.coordinated = [reactivation.reward.vHPC.coordinated ; R];
                
                R = Shuffle_and_Reactivation(ripple_event.uncoordinated, size(ripple_event.all,1), iterations, UnCoor.is.sws, patterns.all.reward, cond.vHPC.reward, th , [bins' , Spikes], 'R' , config, normalization);
%                 [R] = reactivation_strength(patterns.all.reward , cond.vHPC.reward , [bins' , Spikes] , UnCoor.is.sws , th , 'R' , config , normalization , []); clear templates
                reactivation.reward.vHPC.uncoordinated = [reactivation.reward.vHPC.uncoordinated ; R];
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



data = [x ; y ; x1 ; y1];
ref = [ones(length(x),1) ; ones(length(y),1) ; ones(length(x1),1)*2 ; ones(length(y1),1)*2];
ref = [ref , [ones(length(x),1) ; ones(length(y),1)*2 ; ones(length(x1),1) ; ones(length(y1),1)*2]];



perform2WayANOVA(data, ref(:,2), ref(:,1))
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

%% Plot Peaks of activation
%  for joint assemblies
figure
x = (reactivation.reward.dvHPC(:,7));
y = (reactivation.aversive.dvHPC(:,7));
kstest(x)
kstest(y)
[h, p] = ranksum(x,y)  
[h, p] = signrank(y,0)
[h, p] = signrank(x,0)


subplot(131),
grps = [ones(size(x,1),1) ; ones(size(y,1),1)*2];
Y = [x;y];
scatter(grps,Y,"filled",'jitter','on', 'jitterAmount',0.1),hold on
scatter([1 2] , [nanmean(x) nanmean(y)],'filled'),xlim([0 3]),%ylim([-0.3 0.3])

%  for dHPC assemblies
x = reactivation.reward.dHPC(:,7);
y = reactivation.aversive.dHPC(:,7);
kstest(x)
kstest(y)
[h, p] = ranksum(x,y)  
[h, p] = signrank(y,0)
[h, p] = signrank(x,0)


subplot(132),
grps = [ones(size(x,1),1) ; ones(size(y,1),1)*2];
Y = [x;y];
scatter(grps,Y,"filled",'jitter','on', 'jitterAmount',0.1),hold on
scatter([1 2] , [nanmean(x) nanmean(y)],'filled'),xlim([0 3]),%ylim([-0.3 0.3])

%  for vHPC assemblies
x = reactivation.reward.vHPC(:,7);
y = reactivation.aversive.vHPC(:,7);
kstest(x)
kstest(y)
[h, p] = ranksum(x,y,'tail','left')  
[h, p] = signrank(y,0,'tail','right')
[h, p] = signrank(x,0,'tail','left')

subplot(133),
grps = [ones(size(x,1),1) ; ones(size(y,1),1)*2];
Y = [x;y];
scatter(grps,Y,"filled",'jitter','on', 'jitterAmount',0.1),hold on
scatter([1 2] , [nanmean(x) nanmean(y)],'filled'),xlim([0 3]),ylim([-0.3 0.3])

%% Plot Activity Strength surrounding the ripples
[h p] = min(abs([-2:binSize:2+binSize]-0.1))
[h pp] = min(abs([-2:binSize:2+binSize]-(-0.1)))

[h p] = sort(mean(gain.both.aversive.post(:,pp:p),2),'descend')

figure,
x = [];
y = gain.both.aversive.pre(p,:);
for i = 1:size(y,1)
    x = [x ; Smooth(y(i,:),0,'type','l')'];
end
subplot(121),imagesc([-2:binSize:2],[1:size(gain.both.aversive.pre,1)],x),xlim([-0.2 0.2]),caxis([-2 10]), colormap 'jet'

x = [];
y = gain.both.aversive.post(p,:);
for i = 1:size(y,1)
    x = [x ; Smooth(y(i,:),0,'type','l')'];
end
subplot(122),imagesc([-2:binSize:2],[1:size(gain.both.aversive.post,1)],x),xlim([-0.2 0.2]),caxis([-2 10]), colormap 'jet'

figure,
plot([-2:binSize:2],gain.both.aversive.pre,'k'),hold on
plot([-2:binSize:2],gain.both.aversive.post,'g')



[h p] = min(abs([-2:binSize:2+binSize]-0.1))
[h pp] = min(abs([-2:binSize:2+binSize]-(-0.1)))

[h p] = sort(mean(gain.both.reward.post(:,pp:p),2),'descend')

figure,
x = [];
y = gain.both.reward.pre(p,:);
for i = 1:size(y,1)
    x = [x ; Smooth(y(i,:),1,'type','l')'];
end
subplot(121),imagesc([-2:binSize:2],[1:size(gain.both.reward.pre,1)],x),xlim([-0.2 0.2]),caxis([-2 4]), colormap 'jet'

x = [];
y = gain.both.reward.post(p,:);
for i = 1:size(y,1)
    x = [x ; Smooth(y(i,:),1,'type','l')'];
end
subplot(122),imagesc([-2:binSize:2],[1:size(gain.both.reward.post,1)],x),xlim([-0.2 0.2]),caxis([-2 4]), colormap 'jet'


%% Plot cumulative distribution FR
%  for joint assemblies
figure
x = reactivation.aversive.dvHPC(:,2);
y = reactivation.aversive.dvHPC(:,3);
[h p] = kstest2(x,y)

subplot(321),
boxplot([x;y] , [ones(length(x),1) ; ones(length(x),1)*2]),ylim([0 1])

x = reactivation.reward.dvHPC(:,2);
y = reactivation.reward.dvHPC(:,3);
[h p] = kstest2(x,y)

subplot(322),
boxplot([x;y] , [ones(length(x),1) ; ones(length(x),1)*2]),ylim([0 1])

% for dHPC
x = reactivation.aversive.dHPC(:,2);
y = reactivation.aversive.dHPC(:,3);
[h p] = kstest2(x,y)

subplot(323),
boxplot([x;y] , [ones(length(x),1) ; ones(length(x),1)*2]),ylim([0 1])

x = reactivation.reward.dHPC(:,2);
y = reactivation.reward.dHPC(:,3);
[h p] = kstest2(x,y)

subplot(324),
boxplot([x;y] , [ones(length(x),1) ; ones(length(x),1)*2]),ylim([0 1])

% for vHPC
x = reactivation.aversive.vHPC(:,2);
y = reactivation.aversive.vHPC(:,3);
[h p] = kstest2(x,y)

subplot(325),
boxplot([x;y] , [ones(length(x),1) ; ones(length(x),1)*2]),ylim([0 1])

x = reactivation.reward.vHPC(:,2);
y = reactivation.reward.vHPC(:,3);
[h p] = kstest2(x,y)

subplot(326),
boxplot([x;y] , [ones(length(x),1) ; ones(length(x),1)*2]),ylim([0 1])

%% Plot Correlations
figure
subplot(131)
scatter(reactivation.reward.dvHPC(:,4) , reactivation.reward.dvHPC(:,1),'filled','b'),hold on,ylim([-0.5 0.5]),xlim([0 50]),ylim([-1 1])
% scatter(reactivation.aversive.dvHPC(:,4) , reactivation.aversive.dvHPC(:,1),'filled','r'),hold on,xlim([0 50]),ylim([-1 1])

subplot(132)
scatter(reactivation.reward.dHPC(:,4) , reactivation.reward.dHPC(:,1),'filled','b'),hold on,xlim([0 50]),ylim([-1 1])
% scatter(reactivation.aversive.dHPC(:,4) , reactivation.aversive.dHPC(:,1),'filled','r'),hold on,xlim([0 50]),ylim([-1 1])

subplot(133)
scatter(reactivation.reward.vHPC(:,4) , reactivation.reward.vHPC(:,1),'filled','b'),hold on,xlim([0 50]),ylim([-1 1])
% scatter(reactivation.aversive.vHPC(:,4) , reactivation.aversive.vHPC(:,1),'filled','r'),hold on,xlim([0 50]),ylim([-1 1])

figure
subplot(131)
% fitlm(reactivation.aversive.dvHPC(:,4) , reactivation.aversive.dvHPC(:,1))
fitlm(reactivation.reward.dvHPC(:,4) , reactivation.reward.dvHPC(:,1))
plot(ans),xlim([0 50]),ylim([-1 1])

subplot(132)
% fitlm(reactivation.aversive.dHPC(:,4) , reactivation.aversive.dHPC(:,1))
% plot(ans),xlim([0 50]),ylim([-1 1]),hold on
fitlm(reactivation.reward.dHPC(:,4) , reactivation.reward.dHPC(:,1))
plot(ans),xlim([0 50]),ylim([-1 1])

subplot(133)
% fitlm(reactivation.aversive.vHPC(:,4) , reactivation.aversive.vHPC(:,1))
% plot(ans),xlim([0 50]),ylim([-1 1]),hold on
fitlm(reactivation.reward.vHPC(:,4) , reactivation.reward.vHPC(:,1))
plot(ans),xlim([0 50]),ylim([-1 1])



%% Plot Activity Strength
% x = [reactivation.reward.dvHPC(:,4) ; reactivation.reward.dvHPC(:,5) ; reactivation.aversive.dvHPC(:,4) ; reactivation.aversive.dvHPC(:,5)];
% y = [ones(length(reactivation.reward.dvHPC(:,4)),1) ; ones(length(reactivation.reward.dvHPC(:,5)),1)*2 ; ones(length(reactivation.aversive.dvHPC(:,4)),1)*3 ; ones(length(reactivation.aversive.dvHPC(:,5)),1)*4];
% scatter(y,x,"filled",'jitter','on', 'jitterAmount',0.1),xlim([0 5]),hold on
% scatter([1 2 3 4],[nanmean(reactivation.reward.dvHPC(:,4)) nanmean(reactivation.reward.dvHPC(:,5)) nanmean(reactivation.aversive.dvHPC(:,4)) , nanmean(reactivation.aversive.dvHPC(:,5))],"filled"),xlim([0 5]),hold on)

% joint reward
x = reactivation.reward.dvHPC(:,4);
y = reactivation.reward.dvHPC(:,5);

kstest(x)
kstest(y)
[h, p] = ttest2(x,y,'Tail','right')  
[h, p] = ttest(y)
[h, p] = ttest(x)

xx = [1 2];
yy = [nanmean(x) nanmean(y)];
err = [nansem(x) nansem(y)];

subplot(321),
bar(xx,yy),hold on
er = errorbar(xx,yy,err)%;ylim([-3 3])
er.Color = [0 0 0];                            
er.LineStyle = 'none';
hold off

% joint aversive
x = reactivation.aversive.dvHPC(:,4);
y = reactivation.aversive.dvHPC(:,5);

kstest(x)
kstest(y)
[h, p] = ttest2(x,y,'Tail','right')  
[h, p] = ttest(y)
[h, p] = ttest(x)

xx = [1 2];
yy = [nanmean(x) nanmean(y)];
err = [nansem(x) nansem(y)];

subplot(322),
bar(xx,yy),hold on
er = errorbar(xx,yy,err)%;ylim([-3 3])
er.Color = [0 0 0];                            
er.LineStyle = 'none';
hold off


% dHPC reward
x = reactivation.reward.dHPC(:,4);
y = reactivation.reward.dHPC(:,5);

kstest(x)
kstest(y)
[h, p] = ttest2(x,y,'Tail','right')  
[h, p] = ttest(y)
[h, p] = ttest(x)

xx = [1 2];
yy = [nanmean(x) nanmean(y)];
err = [nansem(x) nansem(y)];

subplot(323),
bar(xx,yy),hold on
er = errorbar(xx,yy,err)%;ylim([-3 3])
er.Color = [0 0 0];                            
er.LineStyle = 'none';
hold off

% dHPC aversive
x = reactivation.aversive.dHPC(:,4);
y = reactivation.aversive.dHPC(:,5);

kstest(x)
kstest(y)
[h, p] = ttest2(x,y,'Tail','right')  
[h, p] = ttest(y)
[h, p] = ttest(x)

xx = [1 2];
yy = [nanmean(x) nanmean(y)];
err = [nansem(x) nansem(y)];

subplot(324),
bar(xx,yy),hold on
er = errorbar(xx,yy,err)%;ylim([-3 3])
er.Color = [0 0 0];                            
er.LineStyle = 'none';
hold off


% vHPC reward
x = reactivation.reward.vHPC(:,4);
y = reactivation.reward.vHPC(:,5);

kstest(x)
kstest(y)
[h, p] = ttest2(x,y,'Tail','right')  
[h, p] = ttest(y)
[h, p] = ttest(x)

xx = [1 2];
yy = [nanmean(x) nanmean(y)];
err = [nansem(x) nansem(y)];

subplot(325),
bar(xx,yy),hold on
er = errorbar(xx,yy,err)%;ylim([-3 3])
er.Color = [0 0 0];                            
er.LineStyle = 'none';
hold off

% vHPC aversive
x = reactivation.aversive.vHPC(:,4);
y = reactivation.aversive.vHPC(:,5);

kstest(x)
kstest(y)
[h, p] = ttest2(x,y,'Tail','right')  
[h, p] = ttest(y)
[h, p] = ttest(x)

xx = [1 2];
yy = [nanmean(x) nanmean(y)];
err = [nansem(x) nansem(y)];

subplot(326),
bar(xx,yy),hold on
er = errorbar(xx,yy,err)%;ylim([-3 3])
er.Color = [0 0 0];                            
er.LineStyle = 'none';
hold off

%% Plot Activity Strength
% x = [reactivation.reward.dvHPC(:,4) ; reactivation.reward.dvHPC(:,5) ; reactivation.aversive.dvHPC(:,4) ; reactivation.aversive.dvHPC(:,5)];
% y = [ones(length(reactivation.reward.dvHPC(:,4)),1) ; ones(length(reactivation.reward.dvHPC(:,5)),1)*2 ; ones(length(reactivation.aversive.dvHPC(:,4)),1)*3 ; ones(length(reactivation.aversive.dvHPC(:,5)),1)*4];
% scatter(y,x,"filled",'jitter','on', 'jitterAmount',0.1),xlim([0 5]),hold on
% scatter([1 2 3 4],[nanmean(reactivation.reward.dvHPC(:,4)) nanmean(reactivation.reward.dvHPC(:,5)) nanmean(reactivation.aversive.dvHPC(:,4)) , nanmean(reactivation.aversive.dvHPC(:,5))],"filled"),xlim([0 5]),hold on)

% joint reward
figure
x1 = (reactivation.reward.dvHPC(:,4)-reactivation.reward.dvHPC(:,5))./(reactivation.reward.dvHPC(:,4)+reactivation.reward.dvHPC(:,5));
y=ones(length(reactivation.reward.dvHPC(:,4)),1);
kstest(x1)
[h, p] = signrank(x1,1,'tail','left')

scatter(y,x1,"filled",'jitter','on', 'jitterAmount',0.1),xlim([0 5]),hold on
scatter(1,nanmean(x1),"filled",'jitter','on', 'jitterAmount',0.1),xlim([0 5]),hold on

% joint aversive

x2 = (reactivation.aversive.dvHPC(:,4)-reactivation.aversive.dvHPC(:,5))./(reactivation.aversive.dvHPC(:,4)+reactivation.aversive.dvHPC(:,5));
y=ones(length(reactivation.aversive.dvHPC(:,4)),1)*2;
kstest(x2)
[h, p] = signrank(x2,1,'tail','left')

scatter(y,x2,"filled",'jitter','on', 'jitterAmount',0.1),xlim([0 5]),hold on
scatter(2,nanmean(x2),"filled",'jitter','on', 'jitterAmount',0.1),xlim([0 5]),hold on

% dHPC reward
x3 = (reactivation.reward.dHPC(:,4)-reactivation.reward.dHPC(:,5))./(reactivation.reward.dHPC(:,4)+reactivation.reward.dHPC(:,5));
y=ones(length(reactivation.reward.dHPC(:,4)),1)*3;
kstest(x3)
[h, p] = signrank(x3,1,'tail','left')

scatter(y,x3,"filled",'jitter','on', 'jitterAmount',0.1),xlim([0 5]),hold on
scatter(3,nanmean(x3),"filled",'jitter','on', 'jitterAmount',0.1),xlim([0 5]),hold on

% dHPC aversive
x4 = (reactivation.aversive.dHPC(:,4)-reactivation.aversive.dHPC(:,5))./(reactivation.aversive.dHPC(:,4)+reactivation.aversive.dHPC(:,5));
y=ones(length(reactivation.aversive.dHPC(:,4)),1)*4;
kstest(x4)
[h, p] = signrank(x4,1,'tail','left')

scatter(y,x4,"filled",'jitter','on', 'jitterAmount',0.1),xlim([0 5]),hold on
scatter(4,nanmean(x4),"filled",'jitter','on', 'jitterAmount',0.1),xlim([0 5]),hold on


% vHPC reward
x5 = (reactivation.reward.vHPC(:,4)-reactivation.reward.vHPC(:,5))./(reactivation.reward.vHPC(:,4)+reactivation.reward.vHPC(:,5));
y=ones(length(reactivation.reward.vHPC(:,4)),1)*5;
kstest(x5)
[h, p] = signrank(x5,1,'tail','left')

scatter(y,x5,"filled",'jitter','on', 'jitterAmount',0.1),xlim([0 5]),hold on
scatter(5,nanmean(x5),"filled",'jitter','on', 'jitterAmount',0.1),xlim([0 5]),hold on

% vHPC aversive
x6 = (reactivation.aversive.vHPC(:,4)-reactivation.aversive.vHPC(:,5))./(reactivation.aversive.vHPC(:,4)+reactivation.aversive.vHPC(:,5));
y=ones(length(reactivation.aversive.vHPC(:,4)),1)*6;
kstest(x6)
[h, p] = signrank(x6,0,'tail','right')

scatter(y,x6,"filled",'jitter','on', 'jitterAmount',0.1),xlim([0 5]),hold on
scatter(6,nanmean(x6),"filled",'jitter','on', 'jitterAmount',0.1),xlim([0 7]),hold on
% ylim([-0.5 1])


x = [x1;x2;x3;x4;x5;x6];
grps = [ones(length(x1),1) ; ones(length(x2),1) ; ones(length(x3),1)*2 ; ones(length(x4),1)*2 ; ones(length(x5),1)*3 ; ones(length(x6),1)*3];
grps = [grps , [ones(length(x1),1) ; ones(length(x2),1)*2 ; ones(length(x3),1) ; ones(length(x4),1)*2 ; ones(length(x5),1) ; ones(length(x6),1)*2]];

figure,
anovan(x,grps,'model','interaction','varnames',{'type','condition'})


%% Shock repsonsive assemblies
figure
% both aversive
[h p] = min(abs([-10:binSize:10]-0));
[h pp] = min(abs([-10:binSize:10]-1));

tmp = [];
for i = 1 : size(Shock_response.both.aversive,1)
    tmp = [tmp ; Smooth(Shock_response.both.aversive(i,:),2,'type','l')'];
%     tmp = [tmp ; movmean(Shock_response.both.aversive(i,:),8)];

end
criteria = nanmean(Shock_response.both.aversive(:,p:pp),2);
% p = prctile(criteria,75);
% criteria = criteria >= p;
[h p] = sort(nanmean(tmp(:,p:pp),2),'descend');
subplot(321),imagesc([-10 : binSize : 10] , [1:size(Shock_response.both.aversive,1)] , tmp(p,:)),xlim([-2 3]),caxis([-1 1])%,colormap 'jet'
title('Both Averisve')
xline(0,'-w'), xline(1,'-w'), 

% plot([-10 : binSize : 10],tmp(criteria,:))
% plot(mean(gain.both.aversive.pre(:,criteria),2))

% both reward
[h p] = min(abs([-10:binSize:10]-0));
[h pp] = min(abs([-10:binSize:10]-0.5));

tmp = [];
for i = 1 : size(Shock_response.both.reward,1)
    tmp = [tmp ; Smooth(Shock_response.both.reward(i,:),2,'type','l')'];
%     tmp = [tmp ; movmean(Shock_response.both.reward(i,:),4)];

end
[h p] = sort(nanmean(tmp(:,p:pp),2),'descend');
subplot(322),imagesc([-10 : binSize : 10] , [1:size(Shock_response.both.reward,1)] , tmp(p,:)),xlim([-2 3]),caxis([-1 1])%,colormap 'jet'
title('Both Reward')
xline(0,'-w'), xline(1,'-w'), 


% dHPC aversive
[h p] = min(abs([-10:binSize:10]-0));
[h pp] = min(abs([-10:binSize:10]-0.5));

tmp = [];
for i = 1 : size(Shock_response.dHPC.aversive,1)
    tmp = [tmp ; Smooth(Shock_response.dHPC.aversive(i,:),2,'type','l')'];
%         tmp = [tmp ; movmean(Shock_response.dHPC.aversive(i,:),4)];
end
[h p] = sort(nanmean(tmp(:,p:pp),2),'descend');
subplot(323),imagesc([-10 : binSize : 10] , [1:size(Shock_response.dHPC.aversive,1)] , tmp(p,:)),xlim([-2 3]),caxis([-1 1])%,colormap 'jet'
title('dHPC Averisve')
xline(0,'-w'), xline(1,'-w'), 

% dHPC reward
[h p] = min(abs([-10:binSize:10]-0));
[h pp] = min(abs([-10:binSize:10]-0.5));

tmp = [];
for i = 1 : size(Shock_response.dHPC.reward,1)
    tmp = [tmp ; Smooth(Shock_response.dHPC.reward(i,:),2,'type','l')'];
%         tmp = [tmp ; movmean(Shock_response.dHPC.reward(i,:),4)];
end
[h p] = sort(nanmean(tmp(:,p:pp),2),'descend');
subplot(324),imagesc([-10 : binSize : 10] , [1:size(Shock_response.dHPC.reward,1)] , tmp(p,:)),xlim([-2 3]),caxis([-1 1])%,colormap 'jet'
title('dHPC Reward')
xline(0,'-w'), xline(1,'-w'), 

% vHPC aversive
[h p] = min(abs([-10:binSize:10]-0));
[h pp] = min(abs([-10:binSize:10]-0.5));

tmp = [];
for i = 1 : size(Shock_response.vHPC.aversive,1)
    tmp = [tmp ; Smooth(Shock_response.vHPC.aversive(i,:),2,'type','l')'];
%         tmp = [tmp ; movmean(Shock_response.vHPC.aversive(i,:),4)];

end
[h p] = sort(nanmean(tmp(:,p:pp),2),'descend');
subplot(325),imagesc([-10 : binSize : 10] , [1:size(Shock_response.vHPC.aversive,1)] , tmp(p,:)),xlim([-2 3]),caxis([-1 1])%,colormap 'jet'
title('vHPC Averisve')
xline(0,'-w'), xline(1,'-w'), 

% vHPC reward
[h p] = min(abs([-10:binSize:10]-0));
[h pp] = min(abs([-10:binSize:10]-0.5));

tmp = [];
for i = 1 : size(Shock_response.vHPC.reward,1)
    tmp = [tmp ; Smooth(Shock_response.vHPC.reward(i,:),2,'type','l')'];
%         tmp = [tmp ; movmean(Shock_response.vHPC.reward(i,:),4)];

end
[h p] = sort(nanmean(tmp(:,p:pp),2),'descend');
subplot(326),imagesc([-10 : binSize : 10] , [1:size(Shock_response.vHPC.reward,1)] , tmp(p,:)),xlim([-2 3]),caxis([-1 1])%,colormap 'jet'
title('vHPC Reward')
xline(0,'-w'), xline(1,'-w'), 

%% Shock repsonsive assemblies
figure
tmp = [];
for i = 1 : size(Shock_response.both.aversive,1)
    tmp = [tmp ; Smooth(Shock_response.both.aversive(i,:),1,'type','l')'];
end
subplot(131),plot([-10 : binSize : 10] , mean(tmp),'r'),xlim([-1 2]),hold on
ciplot(mean(tmp)-nansem(tmp) , mean(tmp)+nansem(tmp) , [-10 : binSize : 10] , 'r') , alpha 0.5
title('Both')

% both reward
tmp = [];
for i = 1 : size(Shock_response.both.reward,1)
    tmp = [tmp ; Smooth(Shock_response.both.reward(i,:),1,'type','l')'];
end
plot([-10 : binSize : 10] , mean(tmp),'b'),xlim([-1 2]),ylim([-0.1 2])
ciplot(mean(tmp)-nansem(tmp) , mean(tmp)+nansem(tmp) , [-10 : binSize : 10] , 'b') , alpha 0.5

% dHPC aversive
tmp = [];
for i = 1 : size(Shock_response.dHPC.aversive,1)
    tmp = [tmp ; Smooth(Shock_response.dHPC.aversive(i,:),1,'type','l')'];
end
subplot(132),plot([-10 : binSize : 10] , mean(tmp),'r'),xlim([-1 2]),hold on
ciplot(mean(tmp)-nansem(tmp) , mean(tmp)+nansem(tmp) , [-10 : binSize : 10] , 'r') , alpha 0.5
title('dHPC')

% dHPC reward
tmp = [];
for i = 1 : size(Shock_response.dHPC.reward,1)
    tmp = [tmp ; Smooth(Shock_response.dHPC.reward(i,:),1,'type','l')'];
end
plot([-10 : binSize : 10] , mean(tmp),'b'),xlim([-1 2]),ylim([-0.1 2])
ciplot(mean(tmp)-nansem(tmp) , mean(tmp)+nansem(tmp) , [-10 : binSize : 10] , 'b') , alpha 0.5

% vHPC aversive
tmp = [];
for i = 1 : size(Shock_response.vHPC.aversive,1)
    tmp = [tmp ; Smooth(Shock_response.vHPC.aversive(i,:),1,'type','l')'];
end
subplot(133),plot([-10 : binSize : 10] , mean(tmp),'r'),xlim([-1 2]),hold on
ciplot(mean(tmp)-nansem(tmp) , mean(tmp)+nansem(tmp) , [-10 : binSize : 10] , 'r') , alpha 0.5
title('vHPC')

% dHPC reward
tmp = [];
for i = 1 : size(Shock_response.vHPC.reward,1)
    tmp = [tmp ; Smooth(Shock_response.vHPC.reward(i,:),1,'type','l')'];
end
plot([-10 : binSize : 10] , mean(tmp),'b'),xlim([-1 2]),ylim([-0.1 2])
ciplot(mean(tmp)-nansem(tmp) , mean(tmp)+nansem(tmp) , [-10 : binSize : 10] , 'b') , alpha 0.5





%%

plot(gain.both.aversive.pre(:,criteria),'k'),hold on
plot(gain.both.aversive.post(:,criteria,:),'r'),
