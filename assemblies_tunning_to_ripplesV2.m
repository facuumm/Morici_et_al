clear
clc
close all

%% Parameters
path = {'E:\Rat126\Ephys\in_Pyr';'E:\Rat103\usable';'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat132\recordings\in_pyr';'E:\Rat165\in_pyr\'};%List of folders from the path

%Sleep
time_criteria = 1; % minimal time to include a NREM epoch (in min)

% What par of the code I want to run
S = logical(0);   % Reactivation Strength Calculation
MUAselection = logical(0); % to select ripples by their MUA
W = 'N'; % to select what kind of ripples I want to check
% E= all coordinated ripples, DV dRipple-vRipple, VD vRipple-dRipple
% D= uncoordinated dorsal, V= uncoordinated ventral
% CB = cooridnated bursts
% N= NREM, R= REM
TA =  logical(0); % Trigger Reactivation Strength
TPks = logical(0); %trigger CCG assemblies peaks to events
REC = logical(0); % Assemblie Recruitment during cooridnated events
SRC = logical(0); % If I want to calculate the tuning curve for shocks
C = logical(0); % if I want to calculate CumSum of peaks
SU = logical(1); %if I want to lock the activity of vHPC SUs to dorsal ripples

% for SU
criteria_fr = 0; %criteria to include or not a SU into the analysis
criteria_n = [3 3]; % minimal number of neurons from each structure [vHPC dHPC]
criteria_type = 0; %criteria for celltype (0:pyr, 1:int, 2:all)
binSize = 0.025; %for qssemblie detection qnd qxctivity strength
n_SU_V = 0;
n_SU_D = 0;

win = 60; % time window for bin construction

% Behavior
minimal_speed = 7; % minimal speed to detect quite periods
minimal_speed_time = 2; % minimal time to detect quite periods

reactivation.aversive.dvHPC = [];
reactivation.reward.dvHPC = [];
reactivation.aversive.dHPC = [];
reactivation.reward.dHPC = [];
reactivation.aversive.vHPC = [];
reactivation.reward.vHPC = [];

normalization = false; % to define if normalization over Reactivation Strength is applied or not
th = 5; % threshold for peak detection

gain.both.reward.pre = [];     gain.both.reward.post = [];
gain.both.aversive.pre = [];   gain.both.aversive.post = [];

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
        disp('Uploading session time stamps')
        x = dir([cd,'\*.cat.evt']);
        segments = readtable([cd,'\',x.name],'FileType','text');
        clear x
        % TimeStamps of begening and end of the sleep and awake trials
        % Reward and Aversive trials
        aversiveTS = [];
        aversiveTS_run = [];
        rewardTS = [];
        rewardTS_run = [];
        for y = 1 : height(segments)
            % Baseline sleep session TS detection
            if y == 1
                baselineTS(1,1) = segments.Var1(y);
            elseif y ==2
                baselineTS(1,2) = segments.Var1(y);
            end
            % Aversive sleep session TS detection
            if strcmp(segments.Var2{y},'aversive')
                if strcmp(segments.Var3{y},'End')
                    aversiveTS(1,1) = segments.Var1(y+1);
                    aversiveTS(1,2) = segments.Var1(y+2);
                    aversiveTS_run(1,1) = segments.Var1(y-1);
                    aversiveTS_run(1,2) = segments.Var1(y);
                    A = y;
                end
                % Rewarded sleep session TS detection
            elseif strcmp(segments.Var2{y},'reward')
                if strcmp(segments.Var3{y},'End')
                    rewardTS(1,1) = segments.Var1(y+1);
                    rewardTS(1,2) = segments.Var1(y+2);
                    rewardTS_run(1,1) = segments.Var1(y-1);
                    rewardTS_run(1,2) = segments.Var1(y);
                    R = y;
                end
            end
        end
        clear y A R
        
        % Defining what condition was first
        if aversiveTS_run(1) < rewardTS_run(1)
            config = 1;
        else
            config = 2;
        end
        
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
        ripplesD = table2array(readtable('ripplesD_customized2.csv'));
        ripplesV = table2array(readtable('ripplesV_customized2.csv'));
        % coordination
        coordinated = [];
        coordinatedV = [];
        coordinatedV_refined = [];
        cooridnated_event = [];
        cooridnated_eventDV = [];
        cooridnated_eventVD = [];
        coordinatedD1 = [];
        coordinatedV1 = [];
        for i = 1:length(ripplesD)
            r = ripplesD(i,:);
            tmp = sum(and(ripplesV(:,2)>= r(1,2)-0.025, ripplesV(:,2)<= r(1,2)+0.025));
            if tmp>0
                z = ripplesV(and(ripplesV(:,2)>= r(1,2)-0.025, ripplesV(:,2)<= r(1,2)+0.025),:);
                coordinatedV = [coordinatedV ; z];
                [p,indice] = min(abs(r(2)-z(:,2)));
                coordinatedV_refined = [coordinatedV_refined ; z(indice,:)];
                coordinated = [coordinated ; r];
                
                cooridnated_event = [cooridnated_event ; min([r(1) , z(indice,1)]) , min(min(z(indice,2),r(2)))+abs(z(indice,2)-r(2))/2 , max([r(3) , z(indice,3)])];
                
                if r(2)<z(indice,2)
                    cooridnated_eventDV = [cooridnated_eventDV ; min([r(1) , z(indice,1)]) , min(min(z(indice,2),r(2)))+abs(z(indice,2)-r(2))/2 , max([r(3) , z(indice,3)])];
                    coordinatedD1 = [coordinatedD1 ; r];
                    
                else
                    cooridnated_eventVD = [cooridnated_eventVD ; min([r(1) , z(indice,1)]) , min(min(z(indice,2),r(2)))+abs(z(indice,2)-r(2))/2 , max([r(3) , z(indice,3)])];
                    coordinatedV1 = [coordinatedV1 ; z(indice,:)];
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
        ripples.dHPC.coordinated.baseline = Restrict(coordinated , NREM.baseline);
        ripples.dHPC.coordinated.reward = Restrict(coordinated , NREM.reward);
        ripples.dHPC.coordinated.aversive = Restrict(coordinated , NREM.aversive);
        clear coordinated
        
        % coordinated vRipples
        ripples.vHPC.coordinated.baseline = Restrict(coordinatedV_refined , NREM.baseline);
        ripples.vHPC.coordinated.reward = Restrict(coordinatedV_refined , NREM.reward);
        ripples.vHPC.coordinated.aversive = Restrict(coordinatedV_refined , NREM.aversive);
        clear coordinatedV_refined
        
        %coordinated event
        cooridnated_event((cooridnated_event(:,3)-cooridnated_event(:,1)<0.04),:) = [];
        ripple_event.baseline = Restrict(cooridnated_event,baselineTS./1000);
        ripple_event.reward = Restrict(cooridnated_event,rewardTS./1000);
        ripple_event.aversive = Restrict(cooridnated_event,aversiveTS./1000);
        ripple_event.all = cooridnated_event; clear cooridnated_event
        
        % coordinated event when dRipple was first
        ripple_event.DV.baseline = Restrict(cooridnated_eventDV,baselineTS./1000);
        ripple_event.DV.reward = Restrict(cooridnated_eventDV,rewardTS./1000);
        ripple_event.DV.aversive = Restrict(cooridnated_eventDV,aversiveTS./1000);
        ripple_event.DV.all = cooridnated_eventDV; clear cooridnated_eventDV
        % same but keeping the timestamps from the dorsal ripple
        ripple_event.DV.unique.baseline = Restrict(coordinatedD1,baselineTS./1000);
        ripple_event.DV.unique.reward = Restrict(coordinatedD1,rewardTS./1000);
        ripple_event.DV.unique.aversive = Restrict(coordinatedD1,aversiveTS./1000);
        ripple_event.DV.unique.all = coordinatedD1; clear coordinatedD1
        
        % coordinated event when vRipple was first
        ripple_event.VD.baseline = Restrict(cooridnated_eventVD,baselineTS./1000);
        ripple_event.VD.reward = Restrict(cooridnated_eventVD,rewardTS./1000);
        ripple_event.VD.aversive = Restrict(cooridnated_eventVD,aversiveTS./1000);
        ripple_event.VD.all = cooridnated_eventVD; clear cooridnated_eventVD
        % same but keeping the timestamps from the dorsal ripple 
        ripple_event.VD.unique.baseline = Restrict(coordinatedV1,baselineTS./1000);
        ripple_event.VD.unique.reward = Restrict(coordinatedV1,rewardTS./1000);
        ripple_event.VD.unique.aversive = Restrict(coordinatedV1,aversiveTS./1000);
        ripple_event.VD.unique.all = coordinatedV1; clear coordinatedV1
        
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
        K = [K , Cell_type_classification(:,6:7)];
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
        

        if or(numberD > 2 , numberV > 2) % Assemblies detection
            %% --- Aversive ---
            disp('Lets go for the assemblies')
            if isfile('dorsalventral_assemblies_aversive.mat')
                disp('Loading Aversive template')
                load('dorsalventral_assemblies_aversive.mat')
            else
                disp('Detection of assemblies using Aversive template')
                % --- Options for assemblies detection ---
                opts.Patterns.method = 'ICA';
                opts.threshold.method= 'MarcenkoPastur';
                opts.Patterns.number_of_iterations= 500;
                opts.threshold.permutations_percentile = 0.9;
                opts.threshold.number_of_permutations = 500;
                opts.Patterns.number_of_iterations = 500;
                opts.Members.method = 'Sqrt';
                
                limits = aversiveTS_run./1000;
                events = [];
                events = movement.aversive;
                [SpksTrains.all.aversive , Bins.aversive , Cluster.all.aversive] = spike_train_construction([spks_dHPC;spks_vHPC], clusters.all, cellulartype, binSize, limits, events, false,true);
                [Th , pat] = assembly_patternsJFM([SpksTrains.all.aversive'],opts);
                save([cd,'\dorsalventral_assemblies_aversive.mat'],'Th' , 'pat' , 'criteria_fr' , 'criteria_n')
            end
            
            Thresholded.aversive.all = Th;
            patterns.all.aversive = pat;
            clear cond Th pat
            
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
                cond1 =  logical(0); %checking of dHPC SU
                cond2 =  logical(0); %checking of vHPC SU
                cond.dHPC.aversive = and(cond1 , not(cond2));
                cond.vHPC.aversive = and(cond2 , not(cond1));
                cond.both.aversive = and(cond1 , cond2); clear cond1 cond2
            end
            num_assembliesA = [num_assembliesA ; sum(cond.both.aversive) sum(cond.dHPC.aversive) sum(cond.vHPC.aversive)];
            
            %% --- Reward ---
            disp('Loading Reward template')
            if isfile('dorsalventral_assemblies_reward.mat')
                load('dorsalventral_assemblies_reward.mat')
            else
                disp('Detection of assemblies using Rewarded template')
                % --- Options for assemblies detection ---
                opts.Patterns.method = 'ICA';
                opts.threshold.method= 'MarcenkoPastur';
                opts.Patterns.number_of_iterations= 500;
                opts.threshold.permutations_percentile = 0.9;
                opts.threshold.number_of_permutations = 500;
                opts.Patterns.number_of_iterations = 500;
                opts.Members.method = 'Sqrt';
                
                limits = rewardTS_run./1000;
                events = [];
                events = movement.reward;
                [SpksTrains.all.reward , Bins.reward , Cluster.all.reward] = spike_train_construction([spks_dHPC;spks_vHPC], clusters.all, cellulartype, binSize, limits, events, false,true);
                [Th , pat] = assembly_patternsJFM([SpksTrains.all.reward'],opts);
                save([cd,'\dorsalventral_assemblies_reward.mat'],'Th' , 'pat' , 'criteria_fr' , 'criteria_n')
            end
            
            Thresholded.reward.all = Th;
            patterns.all.reward = pat;
            clear Th pat
            
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

            %% Assemblies activation in the entier recording            
            if SU
                disp('Triggered assemblies activity sourrounding cooridnated events')
                baseline = SubtractIntervals(NREM.all,[ripplesD(:,1)-0.05 , ripplesD(:,3)+0.05 ; ripplesV(:,1)-0.05 , ripplesV(:,3)+0.05]);
%                 baseline = SubtractIntervals(NREM.all,[ripplesD(:,1)-0.05 , ripplesD(:,3)]);
                
                % Upload the tag for ripple modulated SU
                load('RippleModulatedSU.mat')
                
                % Aversive
                if sum(cond.both.aversive) >= 1
                    assemblies = Thresholded.aversive.all(:,cond.both.aversive);
                    p = patterns.all.aversive(:,cond.both.aversive);
                    for i = 1:size(assemblies,2)
                        SUs = clusters.all(assemblies(:,i));
                        SUs = SUs(ismember(SUs,clusters.dHPC));
                        weights = p(assemblies(:,i));
                        q = quantile(p(:,i),0.75);
                        qq = quantile(p(:,i),0.25);
%                         q = 0;
                        qq = q;
                        criteria = pInc.dvHPC.dHPC(ismember(pInc.dvHPC.dHPC(:,1),SUs),:);
                        for ii = 1 : size(SUs,1)
                            if and(or(weights(ii)>=q ,weights(ii)<=qq)  , criteria(ii,3))
%                             if weights(ii)>q
                                tmp = spks(spks(:,1)==SUs(ii),2);
                                if aversiveTS_run(1) < rewardTS_run(1)
                                    b = Restrict(baseline,baselineTS./1000);
                                    [R.baseline , ttttt] = PHIST_Ripple_SU(ripple_event.VD.unique.baseline,tmp,b,2,0.01,1,'Gain'); clear b
                                    b = Restrict(baseline,aversiveTS./1000);
                                    [R.aversive , ttttt] = PHIST_Ripple_SU(ripple_event.VD.unique.aversive,tmp,b,2,0.01,1,'Gain'); clear b
                                    
                                    gain.both.aversive.pre = [gain.both.aversive.pre , R.baseline];
                                    gain.both.aversive.post = [gain.both.aversive.post , R.aversive];
                                    clear R
                                else
                                    b = Restrict(baseline,rewardTS./1000);
                                    [R.reward , ttttt] = PHIST_Ripple_SU(ripple_event.VD.unique.reward,tmp,b,2,0.01,1,'Gain'); clear b
                                    b = Restrict(baseline,aversiveTS./1000);
                                    [R.aversive , ttttt] = PHIST_Ripple_SU(ripple_event.VD.unique.aversive,tmp,b,2,0.01,1,'Gain'); clear b
                                    
                                    gain.both.aversive.pre = [gain.both.aversive.pre , R.reward];
                                    gain.both.aversive.post = [gain.both.aversive.post , R.aversive];
                                    clear R
                                end
                                clear tmp
                            end
                        end
                        clear SUs q weights
                    end
                    clear assemblies p
                end
                
                
                % Reward
                if sum(cond.both.reward) >= 1
                    assemblies = Thresholded.reward.all(:,cond.both.reward);
                    p = patterns.all.reward(:,cond.both.reward);
                    for i = 1:size(assemblies,2)
                        SUs = clusters.all(assemblies(:,i));
                        SUs = SUs(ismember(SUs,clusters.dHPC));
                        weights = p(assemblies(:,i));
                        q = quantile(p(:,i),0.75);
                        qq = quantile(p(:,i),0.25);
%                         q = 0;
                        qq = q;
                        criteria = pInc.dvHPC.dHPC(ismember(pInc.dvHPC.dHPC(:,1),SUs),:);
                        for ii = 1 : size(SUs,1)
                            if and(or(weights(ii)>=q , weights(ii)<=qq) , criteria(ii,3))
%                             if weights(ii)>q
                                tmp = spks(spks(:,1)==SUs(ii),2);
                                if rewardTS_run(1) < aversiveTS_run(1)
                                    b = Restrict(baseline,baselineTS./1000);
                                    [R.baseline , ttttt] = PHIST_Ripple_SU(ripple_event.VD.unique.baseline,tmp,b,2,0.01,1,'Gain'); clear b
                                    b = Restrict(baseline,rewardTS./1000);
                                    [R.reward , ttttt] = PHIST_Ripple_SU(ripple_event.VD.unique.reward,tmp,b,2,0.01,1,'Gain'); clear b
                                    
                                    gain.both.reward.pre = [gain.both.reward.pre , R.baseline];
                                    gain.both.reward.post = [gain.both.reward.post , R.reward];
                                    clear R
                                else
                                    b = Restrict(baseline,aversiveTS./1000);
                                    [R.aversive , ttttt] = PHIST_Ripple_SU(ripple_event.VD.unique.aversive,tmp,b,2,0.01,1,'Gain'); clear b
                                    b = Restrict(baseline,rewardTS./1000);
                                    [R.reward , ttttt] = PHIST_Ripple_SU(ripple_event.VD.unique.reward,tmp,b,2,0.01,1,'Gain'); clear b
                                    
                                    gain.both.reward.pre = [gain.both.reward.pre , R.aversive];
                                    gain.both.reward.post = [gain.both.reward.post , R.reward];
                                    clear R
                                end
                                clear tmp
                            end
                        end
                        clear SUs SUs q weights
                    end
                    clear assemblies
                end
                clear baseline clear pInc pDec surp
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
    end
    
    Number_of_assemblies.aversive = [Number_of_assemblies.aversive ; sum(num_assembliesA)];
    Number_of_assemblies.reward = [Number_of_assemblies.reward ; sum(num_assembliesR)];
    clear num_assembliesA num_assembliesR
    
end

save([cd,'\Reactivation_Strength_Data.mat'] , reactivation)

%% Plot Strenght Reactivation
%  for joint assemblies
figure
x = logical(reactivation.reward.dvHPC(:,6));
y = logical(reactivation.aversive.dvHPC(:,6));

% reactivation.reward.dvHPC(isnan(reactivation.reward.dvHPC(:,1)),:) = [];
% reactivation.aversive.dvHPC(isnan(reactivation.aversive.dvHPC(:,1)),:) = [];
x = reactivation.reward.dvHPC(:,1);
y = reactivation.aversive.dvHPC(:,1);

kstest(x)
kstest(y)
[h, p] = ttest2(x,y,'tail','left')  
[h, p] = ttest(y)
[h, p] = ttest(x)


xx = [1 2];
yy = [nanmean(x) nanmean(y)];
% err = [(std(x)/sqrt(length(x))) (std(y)/sqrt(length(y)))];
err = [nansem(x) nansem(y)];

subplot(131),
bar(xx,yy),hold on
er = errorbar(xx,yy,err)%;ylim([-0.07 0.07])
er.Color = [0 0 0];                            
er.LineStyle = 'none';
hold off

%  for dHPC assemblies
x = logical(reactivation.reward.dHPC(:,6));
y = logical(reactivation.aversive.dHPC(:,6));
% reactivation.reward.dHPC(isnan(reactivation.reward.dHPC(:,1)),:) = [];
% reactivation.aversive.dHPC(isnan(reactivation.aversive.dHPC(:,1)),:) = [];
x = reactivation.reward.dHPC(:,1);
y = reactivation.aversive.dHPC(:,1);
kstest(x)
kstest(y)
[h, p] = kstest2(x,y,'Tail','larger')  
[h, p] = ttest(y)
[h, p] = ttest(x)

xx = [1 2];
yy = [nanmean(x) nanmean(y)];
% err = [(std(x)/sqrt(length(x))) (std(y)/sqrt(length(y)))];
err = [nansem(x) nansem(y)];

subplot(132),
bar(xx,yy),hold on
er = errorbar(xx,yy,err)%;ylim([-0.07 0.07])
er.Color = [0 0 0];                            
er.LineStyle = 'none';
hold off

%  for vHPC assemblies
x = logical(reactivation.reward.vHPC(:,6));
y = logical(reactivation.aversive.vHPC(:,6));
% reactivation.reward.vHPC(isnan(reactivation.reward.vHPC(:,1)),:) = [];
% reactivation.aversive.vHPC(isnan(reactivation.aversive.vHPC(:,1)),:) = [];
x = reactivation.reward.vHPC(:,1);
y = reactivation.aversive.vHPC(:,1);
kstest(x)
kstest(y)
[h, p] = kstest2(x,y,'Tail','larger')  
[h, p] = ttest(y)
[h, p] = ttest(x)

xx = [1 2];
yy = [nanmean(x) nanmean(y)];
% err = [(std(x)/sqrt(length(x))) (std(y)/sqrt(length(y)))];
err = [nansem(x) nansem(y)];

subplot(133),
bar(xx,yy),hold on
er = errorbar(xx,yy,err)%;ylim([-0.07 0.07])
er.Color = [0 0 0];                            
er.LineStyle = 'none';
hold off

%% Plot cumulative distribution
%  for joint assemblies
figure
x = reactivation.reward.dvHPC(:,1);
y = reactivation.aversive.dvHPC(:,1);
[h p] = kstest2(x,y,'Tail','larger')

subplot(131),
cdfplot(x)
hold on
cdfplot(y),xlim([-0.6 0.6])

%  for dHPC assemblies
x = reactivation.reward.dHPC(:,1);
y = reactivation.aversive.dHPC(:,1);
[h p] = kstest2(x,y,'Tail','larger')

subplot(132),
cdfplot(x)
hold on
cdfplot(y),xlim([-0.6 0.6])

%  for vHPC assemblies
x = reactivation.reward.vHPC(:,1);
y = reactivation.aversive.vHPC(:,1);
[h p] = kstest2(x,y,'Tail','larger')

subplot(133),
cdfplot(x)
hold on
cdfplot(y),xlim([-0.6 0.6])


%% Plot Correlations
figure
subplot(131)
% scatter(reactivation.reward.dvHPC(:,4) , reactivation.reward.dvHPC(:,1),'filled','b'),hold on,xlim([0 60]),ylim([-1 1])
scatter(reactivation.aversive.dvHPC(:,4) , reactivation.aversive.dvHPC(:,1),'filled','r'),hold on,xlim([0 60]),ylim([-1 1])

subplot(132)
% scatter(reactivation.reward.dHPC(:,4) , reactivation.reward.dHPC(:,1),'filled','b'),hold on,xlim([0 60]),ylim([-1 1])
scatter(reactivation.aversive.dHPC(:,4) , reactivation.aversive.dHPC(:,1),'filled','r'),hold on,xlim([0 60]),ylim([-1 1])

subplot(133)
% scatter(reactivation.reward.vHPC(:,4) , reactivation.reward.vHPC(:,1),'filled','b'),hold on,xlim([0 60]),ylim([-1 1])
scatter(reactivation.aversive.vHPC(:,4) , reactivation.aversive.vHPC(:,1),'filled','r'),hold on,xlim([0 60]),ylim([-1 1])

figure
subplot(131)
fitlm(reactivation.aversive.dvHPC(:,4) , reactivation.aversive.dvHPC(:,1))
plot(ans),xlim([0 60]),ylim([-1 1]),hold on
fitlm(reactivation.reward.dvHPC(:,4) , reactivation.reward.dvHPC(:,1))
plot(ans),xlim([0 60]),ylim([-1 1]),hold on

subplot(132)
fitlm(reactivation.aversive.dHPC(:,4) , reactivation.aversive.dHPC(:,1))
plot(ans),xlim([0 60]),ylim([-1 1]),hold on
fitlm(reactivation.reward.dHPC(:,4) , reactivation.reward.dHPC(:,1))
plot(ans),xlim([0 60]),ylim([-1 1]),hold on

subplot(133)
fitlm(reactivation.aversive.vHPC(:,4) , reactivation.aversive.vHPC(:,1))
plot(ans),xlim([0 60]),ylim([-1 1]),hold on
fitlm(reactivation.reward.vHPC(:,4) , reactivation.reward.vHPC(:,1))
plot(ans),xlim([0 60]),ylim([-1 1]),hold on

%% Plot cumulative
figure,
BothA = []; BothR = []; dHPCA = []; dHPCR = []; vHPCA = []; vHPCR = [];
dur = 100;
for ttt = 1:3
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
%             disp(['-- Initiating analysis of folder #' , num2str(t) , ' from rat #',num2str(tt) , ' --'])
            session = [subFolders(t+2).folder,'\',subFolders(t+2).name];
            cd([session,'\Spikesorting'])
            if isfile('CumSumPeaks.mat')
                load('CumSumPeaks.mat')
                
                % for Joint Assemblies
                if ttt==1
                    if not(isempty(CumBR))
                        for i = 1 : size(CumBR,1)
                            tmp = nan(1,dur);
                            tmp(1:size(CumBR(i,:),2)) = CumBR(i,:);
                            BothR = [BothR ; tmp];
                            clear tmp
                        end
                    end
                    
                    if not(isempty(CumBA))
                        for i = 1 : size(CumBA,1)
                            tmp = nan(1,dur);
                            tmp(1:size(CumBA(i,:),2)) = CumBA(i,:);
                            BothA = [BothA ; tmp];
                            clear tmp
                        end
                    end
                end
                
                % for Dorsal Assemblies
                if ttt==2
                    if not(isempty(CumDR))
                        for i = 1 : size(CumDR,1)
                            tmp = nan(1,dur);
                            tmp(1:size(CumDR(i,:),2)) = CumDR(i,:);
                            dHPCR = [dHPCR ; tmp];
                            clear tmp
                        end
                    end
                    
                    if not(isempty(CumDA))
                        for i = 1 : size(CumDA,1)
                            tmp = nan(1,dur);
                            tmp(1:size(CumDA(i,:),2)) = CumDA(i,:);
                            dHPCA = [dHPCA ; tmp];
                            clear tmp
                        end
                    end
                end
                
                % for Ventral Assemblies
                if ttt==3
                   if not(isempty(CumVR))
                        for i = 1 : size(CumVR,1)
                            tmp = nan(1,dur);
                            tmp(1:size(CumVR(i,:),2)) = CumVR(i,:);
                            vHPCR = [vHPCR ; tmp];
                            clear tmp
                        end
                    end
                    
                    if not(isempty(CumVA))
                        for i = 1 : size(CumVA,1)
                            tmp = nan(1,dur);
                            tmp(1:size(CumVA(i,:),2)) = CumVA(i,:);
                            vHPCA = [vHPCA ; tmp];
                            clear tmp
                        end
                    end
                end
                clear CumBA CumDA CumVA CumBR CumDR CumVR
                
            end
        end
    end
end

subplot(311)
id = logical(reactivation.aversive.dvHPC(:,end));
plot(nanmean(BothA(id,:)),'r'),hold on

id = logical(reactivation.reward.dvHPC(:,end));
plot(nanmean(BothR(id,:)),'b'),hold on

ciplot(nanmean(BothA)-nansem(BothA) , nanmean(BothA)+nansem(BothA) , [1:1:100],'r'), alpha 0.5

ciplot(nanmean(BothR)-nansem(BothR) , nanmean(BothR)+nansem(BothR) , [1:1:100],'b'), alpha 0.5
xlim([0 60])

subplot(312)
plot(nanmean(dHPCA),'r'),hold on
plot(nanmean(dHPCR),'b'),hold on
ciplot(nanmean(dHPCA)-nansem(dHPCA) , nanmean(dHPCA)+nansem(dHPCA) , [1:1:100],'r'), alpha 0.5
ciplot(nanmean(dHPCR)-nansem(dHPCR) , nanmean(dHPCR)+nansem(dHPCR) , [1:1:100],'b'), alpha 0.5
xlim([0 60])

subplot(313)
plot(nanmean(vHPCA),'r'),hold on
plot(nanmean(vHPCR),'b'),hold on
ciplot(nanmean(vHPCA)-nansem(vHPCA) , nanmean(vHPCA)+nansem(vHPCA) , [1:1:100],'r'), alpha 0.5
ciplot(nanmean(vHPCR)-nansem(vHPCR) , nanmean(vHPCR)+nansem(vHPCR) , [1:1:100],'b'), alpha 0.5
xlim([0 60])