% clear
clc
% close all

%% Parameters
path = {'E:\Rat126\Ephys\in_Pyr';'E:\Rat103\usable';'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat132\recordings\in_pyr';'E:\Rat165\in_pyr\'};%List of folders from the path

% What par of the code I want to run
S = logical(1);   % Reactivation Strength Calculation
W = 'N'; % to select what kind of ripples I want to check
% E= all coordinated ripples, DV dRipple-vRipple, VD vRipple-dRipple
% D= uncoordinated dorsal, V= uncoordinated ventral
% CB = cooridnated bursts
% N= NREM, R= REM
TA =  logical(0); % Trigger Reactivation Strength
Shock_Activity = logical(0); % Trigger activity during shock.

% for SU
criteria_fr = 0; %criteria to include or not a SU into the analysis
criteria_n = [3 3]; % minimal number of neurons from each structure [vHPC dHPC]
criteria_type = 0; %criteria for celltype (0:pyr, 1:int, 2:all)
binSize = 0.025; %for qssemblie detection qnd qxctivity strength
normalization = true; % to define if normalization over Reactivation Strength is applied or not
th = 5; % threshold for peak detection
iterations = 100; % number of iterations for downsampling

% Behavior
minimal_speed = 7; % minimal speed to detect quite periods
minimal_speed_time = 2; % minimal time to detect quite periods

% Storage variables
% For ripple repsonse
gain.both.reward.pre = [];     gain.both.reward.post = [];
gain.both.aversive.pre = [];   gain.both.aversive.post = [];

% For Reactivation measure
reactivation.aversive.dvHPC = []; reactivation.reward.dvHPC = [];
reactivation.aversive.dHPC = []; reactivation.reward.dHPC = [];
reactivation.aversive.vHPC = []; reactivation.reward.vHPC = [];

% Shock repsonse
Shock_response.dHPC.aversive = []; Shock_response.vHPC.aversive = []; Shock_response.both.aversive = [];
Shock_response.dHPC.reward = []; Shock_response.vHPC.reward = []; Shock_response.both.reward = [];

% General metrics
percentages = [];
Number_of_assemblies.aversive = [];
Number_of_assemblies.reward = [];

%% Main loop, to iterate across sessions
for tt = 1:length(path)
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
        load('behavioral_data.mat','movement','Shocks_filt','Rewards_filt')
        
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
            %coordinated event
            cooridnated_event((cooridnated_event(:,3)-cooridnated_event(:,1)<0.04),:) = [];
            ripple_event.baseline = Restrict(cooridnated_event,baselineTS./1000);
            ripple_event.reward = Restrict(cooridnated_event,rewardTS./1000);
            ripple_event.aversive = Restrict(cooridnated_event,aversiveTS./1000);
            ripple_event.all = cooridnated_event;
            % coordinated event when dRipple was first
            ripple_event.DV.baseline = Restrict(cooridnated_eventDV,baselineTS./1000);
            ripple_event.DV.reward = Restrict(cooridnated_eventDV,rewardTS./1000);
            ripple_event.DV.aversive = Restrict(cooridnated_eventDV,aversiveTS./1000);
            ripple_event.DV.all = cooridnated_eventDV;
            % coordinated event when vRipple was first
            ripple_event.VD.baseline = Restrict(cooridnated_eventVD,baselineTS./1000);
            ripple_event.VD.reward = Restrict(cooridnated_eventVD,rewardTS./1000);
            ripple_event.VD.aversive = Restrict(cooridnated_eventVD,aversiveTS./1000);
            ripple_event.VD.all = cooridnated_eventVD;
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
        
%         %% Ripple bursts upload
%         load('coordinated_ripple_bursts.mat')
%         bursts.coordinated.all(or((bursts.coordinated.all(:,2)-bursts.coordinated.all(:,1))<0.08 , (bursts.coordinated.all(:,2)-bursts.coordinated.all(:,1))>1),:) = [];
%         ripple_bursts.baseline = Restrict(bursts.coordinated.all,NREM.baseline);
%         ripple_bursts.reward = Restrict(bursts.coordinated.all,NREM.reward);
%         ripple_bursts.aversive = Restrict(bursts.coordinated.all,NREM.aversive);        
%         clear B bursts
        
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
                a = length(Restrict(spks(spks(:,1)==cluster,2),aversiveTS_run./1000)) / ((aversiveTS_run(2)-aversiveTS_run(1))/1000);
                r = length(Restrict(spks(spks(:,1)==cluster,2),rewardTS_run./1000)) / ((rewardTS_run(2)-rewardTS_run(1))/1000);
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
                a = length(Restrict(spks(spks(:,1)==cluster,2),aversiveTS_run./1000)) / ((aversiveTS_run(2)-aversiveTS_run(1))/1000);
                r = length(Restrict(spks(spks(:,1)==cluster,2),rewardTS_run./1000)) / ((rewardTS_run(2)-rewardTS_run(1))/1000);
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
        
%         % to save the clusters I used for further analysis
%         save([cd,'\SUclusters.mat'],'clusters')
        
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
            
            if S
                % NREM sleep
                if strcmp(W,'N')
                    is.sws.baseline = InIntervals(bins,NREM.baseline);
                    is.sws.reward = InIntervals(bins,NREM.reward);
                    is.sws.aversive = InIntervals(bins,NREM.aversive);
                    
                    is.sws.timestamps.sleep.aversive = NREM.aversive;
                    is.sws.timestamps.sleep.reward = NREM.reward;
                    is.sws.timestamps.sleep.baseline = NREM.baseline;
                    
                elseif strcmp(W,'R')
                    is.sws.baseline = InIntervals(bins,REM.baseline);
                    is.sws.reward = InIntervals(bins,REM.reward);
                    is.sws.aversive = InIntervals(bins,REM.aversive);
                    
                    is.sws.timestamps.sleep.aversive = REM.aversive;
                    is.sws.timestamps.sleep.reward = REM.reward;
                    is.sws.timestamps.sleep.baseline = REM.baseline;
                    
                elseif or(or(strcmp(W,'E'),strcmp(W,'DV')),strcmp(W,'VD'))
                    if and(RD,RV)
                        % coordinated event
                        is.sws.baseline = InIntervals(bins,[ripple_event.baseline(:,1) ripple_event.baseline(:,3)]);
                        is.sws.reward = InIntervals(bins,[ripple_event.reward(:,1) ripple_event.reward(:,3)]);
                        is.sws.aversive = InIntervals(bins,[ripple_event.aversive(:,1) ripple_event.aversive(:,3)]);
                        
                        ripple_event.filtered.baseline = ripple_event.baseline;
                        ripple_event.filtered.aversive = ripple_event.aversive;
                        ripple_event.filtered.reward = ripple_event.reward;
                        
                        is.sws.timestamps.sleep.aversive = NREM.aversive;
                        is.sws.timestamps.sleep.reward = NREM.reward;
                        is.sws.timestamps.sleep.baseline = NREM.baseline;
                    else
                        continue
                    end
                elseif strcmp(W,'CB')
                    % coordinated event
                    is.sws.baseline = InIntervals(bins,ripple_bursts.baseline);
                    is.sws.reward = InIntervals(bins,ripple_bursts.reward);
                    is.sws.aversive = InIntervals(bins,ripple_bursts.aversive);
                    
                    is.sws.timestamps.sleep.aversive = NREM.aversive;
                    is.sws.timestamps.sleep.reward = NREM.reward;
                    is.sws.timestamps.sleep.baseline = NREM.baseline;
                        
                elseif strcmp(W,'D')
                    tmp = not(ismember(ripples.dHPC.baseline(:,2) , ripples.dHPC.coordinated.baseline(:,2)));
                    ripples.dHPC.uncoordinated.baseline = [ripples.dHPC.baseline(tmp,1) ripples.dHPC.baseline(tmp,3)];
                    tmp = [ripples.dHPC.baseline(tmp,1) ripples.dHPC.baseline(tmp,3)];
                    is.sws.baseline = InIntervals(bins,tmp); clear tmp
                    
                    tmp = not(ismember(ripples.dHPC.reward(:,2) , ripples.dHPC.coordinated.reward(:,2)));
                    ripples.dHPC.uncoordinated.reward = [ripples.dHPC.reward(tmp,1) ripples.dHPC.reward(tmp,3)];
                    tmp = [ripples.dHPC.reward(tmp,1) ripples.dHPC.reward(tmp,3)];
                    is.sws.reward = InIntervals(bins,tmp); clear tmp
                    
                    tmp = not(ismember(ripples.dHPC.aversive(:,2) , ripples.dHPC.coordinated.aversive(:,2)));
                    ripples.dHPC.uncoordinated.aversive = [ripples.dHPC.aversive(tmp,1) ripples.dHPC.aversive(tmp,3)];
                    tmp = [ripples.dHPC.aversive(tmp,1) ripples.dHPC.aversive(tmp,3)];
                    is.sws.aversive = InIntervals(bins,tmp); clear tmp
                    
                    is.sws.timestamps.sleep.aversive = NREM.aversive;
                    is.sws.timestamps.sleep.reward = NREM.reward;
                    is.sws.timestamps.sleep.baseline = NREM.baseline;
                                        
                elseif strcmp(W,'V')
                    tmp = not(ismember(ripples.vHPC.baseline(:,2) , ripples.vHPC.coordinated.baseline(:,2)));
                    ripples.vHPC.uncoordinated.baseline = [ripples.vHPC.baseline(tmp,1) ripples.vHPC.baseline(tmp,3)];
                    tmp = [ripples.vHPC.baseline(tmp,1) ripples.vHPC.baseline(tmp,3)];
                    is.sws.baseline = InIntervals(bins,tmp); clear tmp
                    
                    tmp = not(ismember(ripples.vHPC.reward(:,2) , ripples.vHPC.coordinated.reward(:,2)));
                    ripples.vHPC.uncoordinated.reward = [ripples.vHPC.reward(tmp,1) ripples.vHPC.reward(tmp,3)];
                    tmp = [ripples.vHPC.reward(tmp,1) ripples.vHPC.reward(tmp,3)];
                    is.sws.reward = InIntervals(bins,tmp); clear tmp
                    
                    tmp = not(ismember(ripples.vHPC.aversive(:,2) , ripples.vHPC.coordinated.aversive(:,2)));
                    ripples.vHPC.uncoordinated.aversive = [ripples.vHPC.aversive(tmp,1) ripples.vHPC.aversive(tmp,3)];
                    tmp = [ripples.vHPC.aversive(tmp,2)-0.2 ripples.vHPC.aversive(tmp,2)+0.2];
                    is.sws.aversive = InIntervals(bins,tmp); clear tmp
                    
                    is.sws.timestamps.sleep.aversive = NREM.aversive;
                    is.sws.timestamps.sleep.reward = NREM.reward;
                    is.sws.timestamps.sleep.baseline = NREM.baseline;
                                    
                end
                
                is.sws.runaversive = InIntervals(bins,movement.aversive);
                is.sws.runreward = InIntervals(bins,movement.reward);
                
                is.sws.timestamps.run.aversive = InIntervals(bins,aversiveTS_run./1000);
                is.sws.timestamps.run.reward = InIntervals(bins,rewardTS_run./1000);
                
                %% Save timestamps
                is.sws.timestamps.aversiveSleep = aversiveTS./1000;
                is.sws.timestamps.aversiveRun = aversiveTS_run./1000;
                is.sws.timestamps.rewardSleep = rewardTS./1000;
                is.sws.timestamps.rewardRun = rewardTS_run./1000;
                is.sws.timestamps.baselineSleep = baselineTS./1000;
                
                %% Reactivation Strenght
%                 if and(numberD >= criteria_n(1),numberV >= criteria_n(2))
                    % Joint Aversive Assemblies
                    if sum(cond.both.aversive)>=1
                        if and(not(strcmp(W,'V')) , not(strcmp(W,'D'))) % check if is not uncoordinated events
                            
                            templates = ones(size(patterns.all.aversive,1),1);
                            templates(size(clusters.dHPC)+1:end) = 0;
                            templates = [templates ,  ones(size(patterns.all.aversive,1),1)];
                            templates(1:size(clusters.dHPC),2) = 0;
                            
                            [R] = reactivation_strength(patterns.all.aversive , cond.both.aversive , [bins' , Spikes] , is.sws , th , 'A' , config , normalization , []); clear templates
                            
                            tmp = [];
                            for i = 1:size(R,1)
                                
                                x = sum(movement.aversive(:,2)-movement.aversive(:,1));
                                
                                tmp = [tmp ; size(Shocks_filt,1) , x];
                            end
                            
                            reactivation.aversive.dvHPC = [reactivation.aversive.dvHPC ; R , tmp];
                            RBA = R; clear R tmp
                            
                        elseif strcmp(W,'D') % if is uncoordinated dorsal ripples
                            if length(ripples.dHPC.uncoordinated.all) > length(ripples.dHPC.coordinated.all)
                                tmp = [];
                                for i = 1:iterations % downsampling ripples
                                    r = randperm(length(ripples.dHPC.uncoordinated.all));
                                    r = ripples.dHPC.uncoordinated.all(r(1:length(ripples.dHPC.coordinated.all)),:);
                                    r = sort([r(:,1) r(:,3)]);
                                    is.sws.baseline = InIntervals(bins,Restrict(r,NREM.baseline));
                                    is.sws.reward = InIntervals(bins,Restrict(r,NREM.reward));
                                    is.sws.aversive = InIntervals(bins,Restrict(r,NREM.aversive));
                                    [R] = reactivation_strength(patterns.all.aversive , cond.both.aversive , [bins' , Spikes] , is.sws , th , 'A' , config , normalization , []); clear templates
                                    tmp(:,:,i) = R;
                                    clear R r
                                end
                                
                                R = [];
                                for i = 1 : size(tmp,1) % mean the output of the previous loop
                                    r = [];
                                    for ii = 1 : size(tmp,2)
                                        r = [r , nanmean(squeeze(tmp(i,ii,:)))];
                                    end
                                    R = [R ; r]; clear r
                                end
                                reactivation.aversive.dvHPC = [reactivation.aversive.dvHPC ; R];
                                RBA = R; clear R
                            else % if length(uncoordinated) <= length(coordinated)
                                [R] = reactivation_strength(patterns.all.aversive , cond.both.aversive , [bins' , Spikes] , is.sws , th , 'A' , config , normalization , []); clear templates
                                reactivation.aversive.dvHPC = [reactivation.aversive.dvHPC ; R];
                                RBA = R; clear R
                            end
                            
                        elseif strcmp(W,'V') % it does the same if I am taking unccordinated ventral ripples
                            if length(ripples.vHPC.uncoordinated.all) > length(ripples.vHPC.coordinated.all)
                                tmp = [];
                                for i = 1:iterations
                                    r = randperm(length(ripples.vHPC.uncoordinated.all));
                                    r = ripples.vHPC.uncoordinated.all(r(1:length(ripples.vHPC.coordinated.all)),:);
                                    r = sort([r(:,1) r(:,3)]);
                                    is.sws.baseline = InIntervals(bins,Restrict(r,NREM.baseline));
                                    is.sws.reward = InIntervals(bins,Restrict(r,NREM.reward));
                                    is.sws.aversive = InIntervals(bins,Restrict(r,NREM.aversive));
                                    [R] = reactivation_strength(patterns.all.aversive , cond.both.aversive , [bins' , Spikes] , is.sws , th , 'A' , config , normalization , []); clear templates
                                    tmp(:,:,i) = R;
                                    clear R r
                                end
                                
                                R = [];
                                for i = 1 : size(tmp,1)
                                    r = [];
                                    for ii = 1 : size(tmp,2)
                                        r = [r , nanmean(squeeze(tmp(i,ii,:)))];
                                    end
                                    R = [R ; r]; clear r
                                end
                                reactivation.aversive.dvHPC = [reactivation.aversive.dvHPC ; R];
                                RBA = R; clear R
                            else
                                [R] = reactivation_strength(patterns.all.aversive , cond.both.aversive , [bins' , Spikes] , is.sws , th , 'A' , config , normalization , []); clear templates
                                reactivation.aversive.dvHPC = [reactivation.aversive.dvHPC ; R];
                                RBA = R; clear R
                            end
                        end
                    end
%                 end
                
                % dHPC Aversive Assemblies
                if sum(cond.dHPC.aversive)>=1
                    if and(not(strcmp(W,'V')) , not(strcmp(W,'D')))
                        [R] = reactivation_strength(patterns.all.aversive , cond.dHPC.aversive , [bins' , Spikes] , is.sws , th , 'A' , config , normalization , []); clear templates
                        reactivation.aversive.dHPC = [reactivation.aversive.dHPC ; R];
                        RDA = R; clear R
                    elseif strcmp(W,'D')
                        if length(ripples.dHPC.uncoordinated.all) > length(ripples.dHPC.coordinated.all)
                            tmp = [];
                            for i = 1:iterations
                                r = randperm(length(ripples.dHPC.uncoordinated.all));
                                r = ripples.dHPC.uncoordinated.all(r(1:length(ripples.dHPC.coordinated.all)),:);
                                r = sort([r(:,1) r(:,3)]);
                                is.sws.baseline = InIntervals(bins,Restrict(r,NREM.baseline));
                                is.sws.reward = InIntervals(bins,Restrict(r,NREM.reward));
                                is.sws.aversive = InIntervals(bins,Restrict(r,NREM.aversive));
                                [R] = reactivation_strength(patterns.all.aversive , cond.dHPC.aversive , [bins' , Spikes] , is.sws , th , 'A' , config , normalization , []); clear templates
                                tmp(:,:,i) = R;
                                clear R r
                            end
                            
                            R = [];
                            for i = 1 : size(tmp,1)
                                r = [];
                                for ii = 1 : size(tmp,2)
                                    r = [r , nanmean(squeeze(tmp(i,ii,:)))];
                                end
                                R = [R ; r]; clear r
                            end
                            reactivation.aversive.dHPC = [reactivation.aversive.dHPC ; R];
                            RDA = R; clear R
                        else
                            [R] = reactivation_strength(patterns.all.aversive , cond.dHPC.aversive , [bins' , Spikes] , is.sws , th , 'A' , config , normalization , []); clear templates
                            reactivation.aversive.dHPC = [reactivation.aversive.dHPC ; R];
                            RDA = R; clear R
                        end
                    elseif strcmp(W,'V')
                        if length(ripples.vHPC.uncoordinated.all) > length(ripples.vHPC.coordinated.all)
                            tmp = [];
                            for i = 1:iterations
                                r = randperm(length(ripples.vHPC.uncoordinated.all));
                                r = ripples.vHPC.uncoordinated.all(r(1:length(ripples.vHPC.coordinated.all)),:);
                                r = sort([r(:,1) r(:,3)]);
                                is.sws.baseline = InIntervals(bins,Restrict(r,NREM.baseline));
                                is.sws.reward = InIntervals(bins,Restrict(r,NREM.reward));
                                is.sws.aversive = InIntervals(bins,Restrict(r,NREM.aversive));
                                [R] = reactivation_strength(patterns.all.aversive , cond.dHPC.aversive , [bins' , Spikes] , is.sws , th , 'A' , config , normalization , []); clear templates
                                tmp(:,:,i) = R;
                                clear R r
                            end
                            
                            R = [];
                            for i = 1 : size(tmp,1)
                                r = [];
                                for ii = 1 : size(tmp,2)
                                    r = [r , nanmean(squeeze(tmp(i,ii,:)))];
                                end
                                R = [R ; r]; clear r
                            end
                            reactivation.aversive.dHPC = [reactivation.aversive.dHPC ; R];
                            RDA = R; clear R
                        else
                            [R] = reactivation_strength(patterns.all.aversive , cond.dHPC.aversive , [bins' , Spikes] , is.sws , th , 'A' , config , normalization , []); clear templates
                            reactivation.aversive.dHPC = [reactivation.aversive.dHPC ; R];
                            RDA = R; clear R
                        end
                    end
                end
                
                % vHPC Aversive Assemblies
                if sum(cond.vHPC.aversive)>=1
                    if and(not(strcmp(W,'V')) , not(strcmp(W,'D')))
                        [R] = reactivation_strength(patterns.all.aversive , cond.vHPC.aversive , [bins' , Spikes] , is.sws , th , 'A' , config , normalization , []); clear templates
                        reactivation.aversive.vHPC = [reactivation.aversive.vHPC ; R];
                        RVA = R; clear R
                    elseif strcmp(W,'D')
                        if length(ripples.dHPC.uncoordinated.all) > length(ripples.dHPC.coordinated.all)
                            tmp = [];
                            for i = 1:iterations
                                r = randperm(length(ripples.dHPC.uncoordinated.all));
                                r = ripples.dHPC.uncoordinated.all(r(1:length(ripples.dHPC.coordinated.all)),:);
                                r = sort([r(:,1) r(:,3)]);
                                is.sws.baseline = InIntervals(bins,Restrict(r,NREM.baseline));
                                is.sws.reward = InIntervals(bins,Restrict(r,NREM.reward));
                                is.sws.aversive = InIntervals(bins,Restrict(r,NREM.aversive));
                                [R] = reactivation_strength(patterns.all.aversive , cond.vHPC.aversive , [bins' , Spikes] , is.sws , th , 'A' , config , normalization , []); clear templates
                                tmp(:,:,i) = R;
                                clear R r
                            end
                            
                            R = [];
                            for i = 1 : size(tmp,1)
                                r = [];
                                for ii = 1 : size(tmp,2)
                                    r = [r , nanmean(squeeze(tmp(i,ii,:)))];
                                end
                                R = [R ; r]; clear r
                            end
                            reactivation.aversive.vHPC = [reactivation.aversive.vHPC ; R];
                            RVA = R; clear R
                        else
                            [R] = reactivation_strength(patterns.all.aversive , cond.vHPC.aversive , [bins' , Spikes] , is.sws , th , 'A' , config , normalization , []); clear templates
                            reactivation.aversive.vHPC = [reactivation.aversive.vHPC ; R];
                            RVA = R; clear R
                        end
                    elseif strcmp(W,'V')
                        if length(ripples.vHPC.uncoordinated.all) > length(ripples.vHPC.coordinated.all)
                            tmp = [];
                            for i = 1:iterations
                                r = randperm(length(ripples.vHPC.uncoordinated.all));
                                r = ripples.vHPC.uncoordinated.all(r(1:length(ripples.vHPC.coordinated.all)),:);
                                r = sort([r(:,1) r(:,3)]);
                                is.sws.baseline = InIntervals(bins,Restrict(r,NREM.baseline));
                                is.sws.reward = InIntervals(bins,Restrict(r,NREM.reward));
                                is.sws.aversive = InIntervals(bins,Restrict(r,NREM.aversive));
                                [R] = reactivation_strength(patterns.all.aversive , cond.vHPC.aversive , [bins' , Spikes] , is.sws , th , 'A' , config , normalization , []); clear templates
                                tmp(:,:,i) = R;
                                clear R r
                            end
                            
                            R = [];
                            for i = 1 : size(tmp,1)
                                r = [];
                                for ii = 1 : size(tmp,2)
                                    r = [r , nanmean(squeeze(tmp(i,ii,:)))];
                                end
                                R = [R ; r]; clear r
                            end
                            reactivation.aversive.vHPC = [reactivation.aversive.vHPC ; R];
                            RVA = R; clear R
                        else
                            [R] = reactivation_strength(patterns.all.aversive , cond.vHPC.aversive , [bins' , Spikes] , is.sws , th , 'A' , config , normalization , []); clear templates
                            reactivation.aversive.vHPC = [reactivation.aversive.vHPC ; R];
                            RVA = R; clear R
                        end
                    end
                end
                
                %% Same for Reward assemblies
                % Joint Reward Assemblies
%                 if and(numberD >= criteria_n(1),numberV >= criteria_n(2))
                    if sum(cond.both.reward)>=1
                        if and(not(strcmp(W,'V')) , not(strcmp(W,'D')))
                            templates = ones(size(patterns.all.aversive,1),1);
                            templates(size(clusters.dHPC)+1:end) = 0;
                            templates = [templates ,  ones(size(patterns.all.aversive,1),1)];
                            templates(1:size(clusters.dHPC),2) = 0;
                            
                            [R] = reactivation_strength(patterns.all.reward , cond.both.reward , [bins' , Spikes] , is.sws , th , 'R' , config , normalization , []); clear templates
                            
                            tmp = [];
                            for i = 1:size(R,1)
                                
                                x = sum(movement.reward(:,2)-movement.reward(:,1));
                                
                                tmp = [tmp ; size(Rewards_filt,1) , x];
                            end                            
                            
                            reactivation.reward.dvHPC = [reactivation.reward.dvHPC ; R , tmp];
                            RBR = R; clear R tmp
                        elseif strcmp(W,'D')
                            if length(ripples.dHPC.uncoordinated.all) > length(ripples.dHPC.coordinated.all)
                                tmp = [];
                                for i = 1:iterations
                                    r = randperm(length(ripples.dHPC.uncoordinated.all));
                                    r = ripples.dHPC.uncoordinated.all(r(1:length(ripples.dHPC.coordinated.all)),:);
                                    r = sort([r(:,1) r(:,3)]);
                                    is.sws.baseline = InIntervals(bins,Restrict(r,NREM.baseline));
                                    is.sws.reward = InIntervals(bins,Restrict(r,NREM.reward));
                                    is.sws.aversive = InIntervals(bins,Restrict(r,NREM.aversive));
                                    [R] = reactivation_strength(patterns.all.reward , cond.both.reward , [bins' , Spikes] , is.sws , th , 'R' , config , normalization , []); clear templates
                                    tmp(:,:,i) = R;
                                    clear R r
                                end
                                
                                R = [];
                                for i = 1 : size(tmp,1)
                                    r = [];
                                    for ii = 1 : size(tmp,2)
                                        r = [r , nanmean(squeeze(tmp(i,ii,:)))];
                                    end
                                    R = [R ; r]; clear r
                                end
                                reactivation.reward.dvHPC = [reactivation.reward.dvHPC ; R];
                                RBR = R; clear R
                            else
                                [R] = reactivation_strength(patterns.all.reward , cond.both.reward , [bins' , Spikes] , is.sws , th , 'R' , config , normalization , []); clear templates
                                reactivation.reward.dvHPC = [reactivation.reward.dvHPC ; R];
                                RBR = R; clear R
                            end
                        elseif strcmp(W,'V')
                            if length(ripples.vHPC.uncoordinated.all) > length(ripples.vHPC.coordinated.all)
                                tmp = [];
                                for i = 1:iterations
                                    r = randperm(length(ripples.vHPC.uncoordinated.all));
                                    r = ripples.vHPC.uncoordinated.all(r(1:length(ripples.vHPC.coordinated.all)),:);
                                    r = sort([r(:,1) r(:,3)]);
                                    is.sws.baseline = InIntervals(bins,Restrict(r,NREM.baseline));
                                    is.sws.reward = InIntervals(bins,Restrict(r,NREM.reward));
                                    is.sws.aversive = InIntervals(bins,Restrict(r,NREM.aversive));
                                    [R] = reactivation_strength(patterns.all.reward , cond.both.reward , [bins' , Spikes] , is.sws , th , 'R' , config , normalization , []); clear templates
                                    tmp(:,:,i) = R;
                                    clear R r
                                end
                                
                                R = [];
                                for i = 1 : size(tmp,1)
                                    r = [];
                                    for ii = 1 : size(tmp,2)
                                        r = [r , nanmean(squeeze(tmp(i,ii,:)))];
                                    end
                                    R = [R ; r]; clear r
                                end
                                reactivation.reward.dvHPC = [reactivation.reward.dvHPC ; R];
                                RBR = R; clear R
                            else
                                [R] = reactivation_strength(patterns.all.reward , cond.both.reward , [bins' , Spikes] , is.sws , th , 'R' , config , normalization , []); clear templates
                                reactivation.reward.dvHPC = [reactivation.reward.dvHPC ; R];
                                RBR = R; clear R
                            end
                        end
                    end
%                 end
                
                % dHPC Reward Assemblies
                if sum(cond.dHPC.reward)>=1
                    if and(not(strcmp(W,'V')) , not(strcmp(W,'D')))
                        [R] = reactivation_strength(patterns.all.reward , cond.dHPC.reward , [bins' , Spikes] , is.sws , th , 'R' , config , normalization , []); clear templates
                        reactivation.reward.dHPC = [reactivation.reward.dHPC ; R];
                        RDR = R; clear R
                    elseif strcmp(W,'D')
                        if length(ripples.dHPC.uncoordinated.all) > length(ripples.dHPC.coordinated.all)
                            tmp = [];
                            for i = 1:iterations
                                r = randperm(length(ripples.dHPC.uncoordinated.all));
                                r = ripples.dHPC.uncoordinated.all(r(1:length(ripples.dHPC.coordinated.all)),:);
                                r = sort([r(:,1) r(:,3)]);
                                is.sws.baseline = InIntervals(bins,Restrict(r,NREM.baseline));
                                is.sws.reward = InIntervals(bins,Restrict(r,NREM.reward));
                                is.sws.aversive = InIntervals(bins,Restrict(r,NREM.aversive));
                                [R] = reactivation_strength(patterns.all.reward , cond.dHPC.reward , [bins' , Spikes] , is.sws , th , 'R' , config , normalization , []); clear templates
                                tmp(:,:,i) = R;
                                clear R r
                            end
                            
                            R = [];
                            for i = 1 : size(tmp,1)
                                r = [];
                                for ii = 1 : size(tmp,2)
                                    r = [r , nanmean(squeeze(tmp(i,ii,:)))];
                                end
                                R = [R ; r]; clear r
                            end
                            reactivation.reward.dHPC = [reactivation.reward.dHPC ; R];
                            RDR = R; clear R
                        else
                            [R] = reactivation_strength(patterns.all.reward , cond.dHPC.reward , [bins' , Spikes] , is.sws , th , 'R' , config , normalization , []); clear templates
                            reactivation.reward.dHPC = [reactivation.reward.dHPC ; R];
                            RDR = R; clear R
                        end
                    elseif strcmp(W,'V')
                        if length(ripples.vHPC.uncoordinated.all) > length(ripples.vHPC.coordinated.all)
                            tmp = [];
                            for i = 1:iterations
                                r = randperm(length(ripples.vHPC.uncoordinated.all));
                                r = ripples.vHPC.uncoordinated.all(r(1:length(ripples.vHPC.coordinated.all)),:);
                                r = sort([r(:,1) r(:,3)]);
                                is.sws.baseline = InIntervals(bins,Restrict(r,NREM.baseline));
                                is.sws.reward = InIntervals(bins,Restrict(r,NREM.reward));
                                is.sws.aversive = InIntervals(bins,Restrict(r,NREM.aversive));
                                [R] = reactivation_strength(patterns.all.reward , cond.dHPC.reward , [bins' , Spikes] , is.sws , th , 'R' , config , normalization , []); clear templates
                                tmp(:,:,i) = R;
                                clear R r
                            end
                            
                            R = [];
                            for i = 1 : size(tmp,1)
                                r = [];
                                for ii = 1 : size(tmp,2)
                                    r = [r , nanmean(squeeze(tmp(i,ii,:)))];
                                end
                                R = [R ; r]; clear r
                            end
                            reactivation.reward.dHPC = [reactivation.reward.dHPC ; R];
                            RDR = R; clear R
                        else
                            [R] = reactivation_strength(patterns.all.reward , cond.dHPC.reward , [bins' , Spikes] , is.sws , th , 'R' , config , normalization , []); clear templates
                            reactivation.reward.dHPC = [reactivation.reward.dHPC ; R];
                            RDR = R; clear R
                        end
                    end
                end
                
                % vHPC reward Assemblies
                if sum(cond.vHPC.reward)>=1
                    if and(not(strcmp(W,'V')) , not(strcmp(W,'D')))
                        [R] = reactivation_strength(patterns.all.reward , cond.vHPC.reward , [bins' , Spikes] , is.sws , th , 'R' , config , normalization , []); clear templates
                        reactivation.reward.vHPC = [reactivation.reward.vHPC ; R];
                        RVR = R; clear R
                    elseif strcmp(W,'D')
                        if length(ripples.dHPC.uncoordinated.all) > length(ripples.dHPC.coordinated.all)
                            tmp = [];
                            for i = 1:iterations
                                r = randperm(length(ripples.dHPC.uncoordinated.all));
                                r = ripples.dHPC.uncoordinated.all(r(1:length(ripples.dHPC.coordinated.all)),:);
                                r = sort([r(:,1) r(:,3)]);
                                is.sws.baseline = InIntervals(bins,Restrict(r,NREM.baseline));
                                is.sws.reward = InIntervals(bins,Restrict(r,NREM.reward));
                                is.sws.aversive = InIntervals(bins,Restrict(r,NREM.aversive));
                                [R] = reactivation_strength(patterns.all.reward , cond.vHPC.reward , [bins' , Spikes] , is.sws , th , 'R' , config , normalization , []); clear templates
                                tmp(:,:,i) = R;
                                clear R r
                            end
                            
                            R = [];
                            for i = 1 : size(tmp,1)
                                r = [];
                                for ii = 1 : size(tmp,2)
                                    r = [r , nanmean(squeeze(tmp(i,ii,:)))];
                                end
                                R = [R ; r]; clear r
                            end
                            reactivation.reward.vHPC = [reactivation.reward.vHPC ; R];
                            RVR = R; clear R
                        else
                            [R] = reactivation_strength(patterns.all.reward , cond.vHPC.reward , [bins' , Spikes] , is.sws , th , 'R' , config , normalization , []); clear templates
                            reactivation.reward.vHPC = [reactivation.reward.vHPC ; R];
                            RVR = R; clear R
                        end
                    elseif strcmp(W,'V')
                        if length(ripples.vHPC.uncoordinated.all) > length(ripples.vHPC.coordinated.all)
                            tmp = [];
                            for i = 1:iterations
                                r = randperm(length(ripples.vHPC.uncoordinated.all));
                                r = ripples.vHPC.uncoordinated.all(r(1:length(ripples.vHPC.coordinated.all)),:);
                                r = sort([r(:,1) r(:,3)]);
                                is.sws.baseline = InIntervals(bins,Restrict(r,NREM.baseline));
                                is.sws.reward = InIntervals(bins,Restrict(r,NREM.reward));
                                is.sws.aversive = InIntervals(bins,Restrict(r,NREM.aversive));
                                [R] = reactivation_strength(patterns.all.reward , cond.vHPC.reward , [bins' , Spikes] , is.sws , th , 'R' , config , normalization , []); clear templates
                                tmp(:,:,i) = R;
                                clear R r
                            end
                            
                            R = [];
                            for i = 1 : size(tmp,1)
                                r = [];
                                for ii = 1 : size(tmp,2)
                                    r = [r , nanmean(squeeze(tmp(i,ii,:)))];
                                end
                                R = [R ; r]; clear r
                            end
                            reactivation.reward.vHPC = [reactivation.reward.vHPC ; R];
                            RVR = R; clear R
                        else
                            [R] = reactivation_strength(patterns.all.reward , cond.vHPC.reward , [bins' , Spikes] , is.sws , th , 'R' , config , normalization , []); clear templates
                            reactivation.reward.vHPC = [reactivation.reward.vHPC ; R];
                            RVR = R; clear R
                        end
                    end
                end
                
                if exist('RBA')
                    save([cd,'\react_coord_ripples_aversive.mat'],'RBA')
                end
                
                if exist('RBR')
                    save([cd,'\react_coord_ripples_reward.mat'],'RBR')
                end
                clear RBA RBR
            end
            
            %% Triggered assemblies activity
            if TA
                disp('Triggered assemblies activity sourrounding cooridnated events')
                % Aversive
                %                 if RBA(:,end)>=1
                %                     condd = logical(RBA(:,end));
                if sum(cond.both.aversive) >= 1
                    if aversiveTS_run(1) < rewardTS_run(1)
                        tmp = [ripple_event.baseline(:,1)-0.05 ripple_event.baseline(:,3)+0.05];
                        [R.baseline] = triggered_activity(patterns.all.aversive , cond.both.aversive , [bins' , Spikes], 2 , ripple_event.filtered.baseline , 0 , 0 , [] , [] , [size(clusters.dHPC,1) size(clusters.vHPC,1)]); clear tmp
                        tmp = [ripple_event.aversive(:,1)-0.05 ripple_event.aversive(:,3)+0.05];
                        [R.aversive] = triggered_activity(patterns.all.aversive , cond.both.aversive , [bins' , Spikes], 2 , ripple_event.filtered.aversive , 0 , 0 , [] , [] , [size(clusters.dHPC,1) size(clusters.vHPC,1)]); clear tmp
                        
                        if isempty(gain.both.aversive.pre)
                            gain.both.aversive.pre = [gain.both.aversive.pre ; R.baseline];
                            gain.both.aversive.post = [gain.both.aversive.post ; R.aversive];
                        else
                            gain.both.aversive.pre = [gain.both.aversive.pre ; R.baseline(:,1:size(gain.both.aversive.pre,2))];
                            gain.both.aversive.post = [gain.both.aversive.post ; R.aversive(:,1:size(gain.both.aversive.pre,2))];
                        end
                        
                        clear R
                    else
                        tmp = [ripple_event.reward(:,1)-0.05 ripple_event.reward(:,3)+0.05];
                        [R.reward] = triggered_activity(patterns.all.aversive , cond.both.aversive , [bins' , Spikes], 2 , ripple_event.filtered.reward , 0 , 0 , [] , [] , [size(clusters.dHPC,1) size(clusters.vHPC,1)]); clear tmp
                        tmp = [ripple_event.aversive(:,1)-0.05 ripple_event.aversive(:,3)+0.05];
                        [R.aversive] = triggered_activity(patterns.all.aversive , cond.both.aversive , [bins' , Spikes], 2 , ripple_event.filtered.aversive , 0 , 0 , [] , [] , [size(clusters.dHPC,1) size(clusters.vHPC,1)]); clear tmp
                        
                        if isempty(gain.both.aversive.pre)
                            gain.both.aversive.pre = [gain.both.aversive.pre ; R.reward];
                            gain.both.aversive.post = [gain.both.aversive.post ; R.aversive];
                        else
                            gain.both.aversive.pre = [gain.both.aversive.pre ; R.reward(:,1:size(gain.both.aversive.pre,2))];
                            gain.both.aversive.post = [gain.both.aversive.post ; R.aversive(:,1:size(gain.both.aversive.pre,2))];
                        end
                        clear R
                    end
                    %                     end
                    %                     clear condd
                end
                
                % Reward
                %                 if RBR(:,end)>=1
                %                     condd = logical(RBR(:,end));
                if sum(cond.both.reward) >= 1
                    if aversiveTS_run(1) > rewardTS_run(1)
                        tmp = [ripple_event.baseline(:,1)-0.05 ripple_event.baseline(:,3)+0.05];
                        [R.baseline] = triggered_activity(patterns.all.reward , cond.both.reward , [bins' , Spikes], 2 , ripple_event.filtered.baseline , 0 , 0 , [] , [] , [size(clusters.dHPC,1) size(clusters.vHPC,1)]);clear tmp
                        tmp = [ripple_event.reward(:,1)-0.05 ripple_event.reward(:,3)+0.05];
                        [R.reward] = triggered_activity(patterns.all.reward , cond.both.reward , [bins' , Spikes], 2 , ripple_event.filtered.reward , 0 , 0 , [] , [] , [size(clusters.dHPC,1) size(clusters.vHPC,1)]); clear tmp
                        
                        if isempty(gain.both.reward.pre)
                            gain.both.reward.pre = [gain.both.reward.pre ; R.baseline];
                            gain.both.reward.post = [gain.both.reward.post ; R.reward];
                        else
                            gain.both.reward.pre = [gain.both.reward.pre ; R.baseline(:,1:size(gain.both.reward.pre,2))];
                            gain.both.reward.post = [gain.both.reward.post ; R.reward(:,1:size(gain.both.reward.pre,2))];
                        end
                        clear R
                    else
                        tmp = [ripple_event.aversive(:,1)-0.05 ripple_event.aversive(:,3)+0.05];
                        [R.aversive] = triggered_activity(patterns.all.reward , cond.both.reward , [bins' , Spikes], 2 , ripple_event.filtered.aversive , 0 , 0 , [] , [] , [size(clusters.dHPC,1) size(clusters.vHPC,1)]); clear tmp
                        tmp = [ripple_event.reward(:,1)-0.05 ripple_event.reward(:,3)+0.05];
                        [R.reward] = triggered_activity(patterns.all.reward , cond.both.reward , [bins' , Spikes], 2 , ripple_event.filtered.reward , 0 , 0 , [] , [] , [size(clusters.dHPC,1) size(clusters.vHPC,1)]); clear tmp
                        
                        if isempty(gain.both.reward.pre)
                            gain.both.reward.pre = [gain.both.reward.pre ; R.aversive];
                            gain.both.reward.post = [gain.both.reward.post ; R.reward];
                        else
                            gain.both.reward.pre = [gain.both.reward.pre ; R.aversive(:,1:size(gain.both.reward.pre,2))];
                            gain.both.reward.post = [gain.both.reward.post ; R.reward(:,1:size(gain.both.reward.pre,2))];
                        end
                        clear R
                    end
                end
                %                     clear cond
                %                 end
            end
            
            if Shock_Activity
                disp('Calculating mean response to the Shocks')
                baseline = SubtractIntervals(aversiveTS_run./1000,[Shocks_filt-1 Shocks_filt+2]);
                if sum(cond.both.aversive)>=1
                    [R] = triggered_activity(patterns.all.aversive , cond.both.aversive , [bins' , Spikes] , 10 , Shocks_filt, 0 , 1 , [] , baseline , []);
                    Shock_response.both.aversive = [Shock_response.both.aversive ; R]; clear R
                end
                
                if sum(cond.dHPC.aversive)>=1
                    [R] = triggered_activity(patterns.all.aversive , cond.dHPC.aversive , [bins' , Spikes] , 10 ,  Shocks_filt, 0 , 1 ,[] , baseline , []);
                    Shock_response.dHPC.aversive = [Shock_response.dHPC.aversive ; R]; clear R
                end
                
                if sum(cond.vHPC.aversive)>=1
                    [R] = triggered_activity(patterns.all.aversive , cond.vHPC.aversive , [bins' , Spikes] , 10 ,  Shocks_filt, 0 , 1 ,[] , baseline , []);
                    Shock_response.vHPC.aversive = [Shock_response.vHPC.aversive ; R]; clear R
                end
                
                if sum(cond.both.reward)>=1
                    [R] = triggered_activity(patterns.all.reward , cond.both.reward , [bins' , Spikes] , 10 , Shocks_filt, 0 , 1 , [] , baseline , []);
                    Shock_response.both.reward = [Shock_response.both.reward ; R]; clear R
                end
                
                if sum(cond.dHPC.reward)>=1
                    [R] = triggered_activity(patterns.all.reward , cond.dHPC.reward , [bins' , Spikes] , 10 ,  Shocks_filt, 0 , 1 ,[] , baseline , []);
                    Shock_response.dHPC.reward = [Shock_response.dHPC.reward ; R]; clear R
                end
                
                if sum(cond.vHPC.reward)>=1
                    [R] = triggered_activity(patterns.all.reward , cond.vHPC.reward , [bins' , Spikes] , 10 ,  Shocks_filt, 0 , 1 ,[] , baseline , []);
                    Shock_response.vHPC.reward = [Shock_response.vHPC.reward ; R]; clear R
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
        clear cooridnated_eventDV cooridnated_eventVD segments
        clear RBA RBR RDA RDR RVA RVR Shocks_filt Rewards_filt
    end
    
    Number_of_assemblies.aversive = [Number_of_assemblies.aversive ; sum(num_assembliesA)];
    Number_of_assemblies.reward = [Number_of_assemblies.reward ; sum(num_assembliesR)];
    clear num_assembliesA num_assembliesR
end

% save([cd,'\Reactivation_Strength_Data_Normalized_eliminating_withingD_5SD.mat'],'reactivation')
%% Correlations
%Reactivation and shocks
% Aversive
figure
y = reactivation.aversive.dvHPC(:,1);
x =  reactivation.aversive.dvHPC(:,end-1);
fitlm(x,y)
plot(ans),hold on
xlim([2 14])
ylim([-0.6 0.6])
% Reward
figure
y = reactivation.reward.dvHPC(:,1);
x =  reactivation.reward.dvHPC(:,end-1);
fitlm(x,y)
plot(ans)
ylim([-0.6 0.6])
xlim([0 90])

% Scatters
% Aversive
figure
y = reactivation.aversive.dvHPC(:,1);
x =  reactivation.aversive.dvHPC(:,end-1);
scatter(x,y,'filled')
ylim([-0.6 0.6])
xlim([2 14])
% Reward
figure
y = reactivation.reward.dvHPC(:,1);
x =  reactivation.reward.dvHPC(:,end-1);
scatter(x,y,'filled')
ylim([-0.6 0.6])
xlim([0 90])


%Reactivation and movment
% Aversive
figure
y = reactivation.aversive.dvHPC(:,1);
x =  reactivation.aversive.dvHPC(:,end-1);
fitlm(x,y)
plot(ans),hold on
% Reward
y = reactivation.reward.dvHPC(:,1);
x =  reactivation.reward.dvHPC(:,end-1);
fitlm(x,y)
plot(ans)
ylim([-0.6 0.6])
xlim([100 500])

% Scatters
% Aversive
figure
y = reactivation.aversive.dvHPC(:,1);
x =  reactivation.aversive.dvHPC(:,end-1);
scatter(x,y,'filled')
ylim([-0.6 0.6])
xlim([100 500])
% Reward
figure
y = reactivation.reward.dvHPC(:,1);
x =  reactivation.reward.dvHPC(:,end-1);
scatter(x,y,'filled')
ylim([-0.6 0.6])
xlim([100 500])

%% Peaks mean for joint assemblies
figure
x = reactivation.reward.dvHPC(:,1);
y = reactivation.aversive.dvHPC(:,1);

kstest(x)
kstest(y)
[h, p] = ranksum(x,y,'tail','left')  
[h, p] = signrank(y,0,'tail','right')
[h, p] = signrank(x,0,'tail','right')

subplot(131),
grps = [ones(size(x,1),1) ; ones(size(y,1),1)*2];
Y = [x;y];
scatter(grps,Y,"filled",'jitter','on', 'jitterAmount',0.1),hold on
scatter([1 2] , [nanmean(x) nanmean(y)],'filled'),xlim([0 3]),ylim([-1 1])

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
