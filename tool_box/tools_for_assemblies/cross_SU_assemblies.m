function [cross_members cross_nonmembers time lags] = cross_SU_assemblies(path,events,type_assemblies)
% This function a CCG between SU members and non-membersof assemblies.
% It iterates within the subfolders within the path. It will do it for the
% Pre and Post Sleep.
%
% INPUTS
% path: cell, inside each cell you should put the path of each session you
%       want to analyze
%
% events: str, you define is is 'NREM', 'REM' or 'Coordinated'
%
% type_assemblies: str, it defines which type of assemblies you want to
%                  correlate. So far, I just implemanted for 'Both'. That
%                  are dHPC-vHPC joint assemblies.
%
% OUTPUT
% cross_members: structure, it contains the CCG for aversive and reward
%                joint assemblies, during the Pre and Post sleep.
%                cross_members.aversive.pre = [];    cross_members.aversive.post = [];
%                cross_members.reward.pre = [];    cross_members.reward.post = [];
%
% cross_nonmembers: structure, same as the output described above.
%
% time: column vector, it contains the time axis for the CCG plot.
%
% lags: structure, it contains the max value of each CCG for each
%       condition. SImilar structure as thoe outputs described above.
%
%
% other functions: CCG from FMA toolbox
% Morci Juan Facundo 01/2024


% Initialization of structures that wil contain the outputs
cross_members.aversive.pre = [];    cross_members.aversive.post = [];
cross_nonmembers.aversive.pre = [];    cross_nonmembers.aversive.post = [];
cross_members.reward.pre = [];    cross_members.reward.post = [];
cross_nonmembers.reward.pre = [];    cross_nonmembers.reward.post = [];

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
        
        NREM.all(NREM.all(:,2)-NREM.all(:,1)<60*1,:)=[];
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
        ripplesD = [ripplesD(:,1) ripplesD(:,3)];
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
%         cooridnated_event((cooridnated_event(:,3)-cooridnated_event(:,1)<0.04),:) = [];
        cooridnated_event = [cooridnated_event(:,1)-1 cooridnated_event(:,3)+1];
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
        cellulartype = [K(:,1) , K(:,4)];

        
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
                if or(a > 0 ,  r > 0)
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
                if or(a > 0 ,  r > 0)
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
        
        
        %% --- Aversive ---
        disp('Lets go for the assemblies')
        if isfile('dorsalventral_assemblies_aversive.mat')
            disp('Loading Aversive template')
            load('dorsalventral_assemblies_aversive.mat')
        else
            Th = [];
            pat = [];
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
        
        %% --- Reward ---
        disp('Loading Reward template')
        if isfile('dorsalventral_assemblies_reward.mat')
            load('dorsalventral_assemblies_reward.mat')
        else
            Th = [];
            pat = [];
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
        
        if strcmp(type_assemblies,'Both')
            %             if strcmp(events,'NREM')
            %                 Events = NREM;
            %             elseif strcmp(events,'REM')
            %                 Events = REM;
            %             elseif strcmp(events,'Coordinated')
            %                 Events = ripple_event;
            %             end
            
%             %% SpikeTrain
%             limits = [0 segments.Var1(end)/1000];
%             [Spikes , bins , Clusters] = spike_train_construction([spks_dHPC;spks_vHPC], clusters.all, cellulartype, 0.01, limits, [] , true, true);
%             clear limits
            
            if strcmp(events,'NREM')
                Events = NREM;
            elseif strcmp(events,'REM')
                Events = REM;
            elseif strcmp(events,'Coordinated')
                Events = ripple_event;
            end
            
            % Load tags
            load('RippleModulatedSU.mat')
            
            if sum(cond.both.aversive)>0
                %Defining the assemblies to work with
                p = Thresholded.aversive.all(:,cond.both.aversive);
                w = patterns.all.aversive(:,cond.both.aversive);
                %Definint Events to study for this assembly
                if config==1
                    pre = Events.baseline;
                    post = Events.aversive;
                else
                    pre = Events.reward;
                    post = Events.aversive;
                end
                
                % Here I will iterate across dHPC members and perform a CCG
                % between members from the vHPC.
                for i = 1 : sum(cond.both.aversive)
                    parameter = p(:,i);
                    weigth = w(:,i);
                    q = quantile(weigth,0.8);
                    for ii = 1 :size(clusters.dHPC,1)
                        criteriaD = pInc.dvHPC.dHPC(ismember(pInc.dvHPC.dHPC(:,1),clusters.dHPC(ii)),3);
                        if and(parameter(ii),criteriaD)
                            if weigth(ii)>q
                                %                             clusterD = clusters.dHPC(ii);
                                SpikesD = spks_dHPC(spks_dHPC(:,1)==clusters.dHPC(ii),2);
                                for iii = 1 :size(clusters.vHPC,1)
                                    criteriaV = pInc.dvHPC.vHPC(ismember(pInc.dvHPC.vHPC(:,1),clusters.vHPC(iii)),3);
                                    if and(parameter(size(clusters.dHPC,1)+iii) , criteriaV)
                                        if weigth(size(clusters.dHPC,1)+iii)>q
                                            SpikesV = spks_vHPC(spks_vHPC(:,1)==clusters.vHPC(iii),2);
                                            
                                            % For Pre
                                            x = Restrict(SpikesD,pre);
                                            y = Restrict(SpikesV,pre);
                                            if and(length(x)>5,length(y)>5)
                                                [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                                                [ccg,tttt] = CCG(s,ids,'binSize',0.005,'duration',2,'smooth',1,'mode','ccg'); %ccg calculation
                                                time = tttt;
                                                cross_members.aversive.pre = [cross_members.aversive.pre , ccg(:,1,2)];
                                                clear x y s ids groups ccg tttt
                                            end
                                            
                                            % for post
                                            x = Restrict(SpikesD,post);
                                            y = Restrict(SpikesV,post);
                                            if and(length(x)>5,length(y)>5)
                                                [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                                                [ccg,tttt] = CCG(s,ids,'binSize',0.005,'duration',2,'smooth',1,'mode','ccg'); %ccg calculation
                                                time = tttt;
                                                cross_members.aversive.post = [cross_members.aversive.post , ccg(:,1,2)];
                                                clear x y s ids groups ccg tttt
                                            end
                                        else
                                            SpikesV = spks_vHPC(spks_vHPC(:,1)==clusters.vHPC(iii),2);
                                            
                                            % For Pre
                                            x = Restrict(SpikesD,pre);
                                            y = Restrict(SpikesV,pre);
                                            if and(length(x)>5,length(y)>5)
                                                [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                                                [ccg,tttt] = CCG(s,ids,'binSize',0.005,'duration',2,'smooth',1,'mode','ccg'); %ccg calculation
                                                time = tttt;
                                                cross_nonmembers.aversive.pre = [cross_nonmembers.aversive.pre , ccg(:,1,2)];
                                                clear x y s ids groups ccg tttt
                                            end
                                            
                                            % for post
                                            x = Restrict(SpikesD,post);
                                            y = Restrict(SpikesV,post);
                                            if and(length(x)>10,length(y)>10)
                                                [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                                                [ccg,tttt] = CCG(s,ids,'binSize',0.005,'duration',2,'smooth',1,'mode','ccg'); %ccg calculation
                                                time = tttt;
                                                cross_nonmembers.aversive.post = [cross_nonmembers.aversive.post , ccg(:,1,2)];
                                                clear x y s ids groups ccg tttt
                                            end
                                        end
                                    end
                                    clear criteriaV
                                end
                            end
                        end
                        clear criteriaD
                    end
                    clear parameter weigth
                end
                clear p w q
            end
            
            if sum(cond.both.reward)>0
                %Defining the assemblies to work with
                p = Thresholded.reward.all(:,cond.both.reward);
                w = patterns.all.reward(:,cond.both.reward);
                %Definint Events to study for this assembly
                if config==1
                    pre = Events.aversive;
                    post = Events.reward;
                else
                    pre = Events.baseline;
                    post = Events.reward;
                end
                
                % Here I will iterate across dHPC members and perform a CCG
                % between members from the vHPC.
                for i = 1 : sum(cond.both.reward)
                    parameter = p(:,i);
                    weigth = w(:,i);
                    q = quantile(weigth,0.8);
                    for ii = 1 :size(clusters.dHPC,1)
                        criteriaD = pInc.dvHPC.dHPC(ismember(pInc.dvHPC.dHPC(:,1),clusters.dHPC(ii)),3);
                        if and(parameter(ii) , criteriaD)
                            if weigth(ii)>q
                                SpikesD = spks_dHPC(spks_dHPC(:,1)==clusters.dHPC(ii),2);
                                for iii = 1 :size(clusters.vHPC,1)
                                    criteriaV = pInc.dvHPC.vHPC(ismember(pInc.dvHPC.vHPC(:,1),clusters.vHPC(iii)),3);
                                    if and(parameter(size(clusters.dHPC,1)+iii) , criteriaV)
                                        if weigth(size(clusters.dHPC,1)+iii)>q
                                            SpikesV = spks_vHPC(spks_vHPC(:,1)==clusters.vHPC(iii),2);
                                            
                                            % For Pre
                                            x = Restrict(SpikesD,pre);
                                            y = Restrict(SpikesV,pre);
                                            if and(length(x)>5,length(y)>5)
                                                [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                                                [ccg,tttt] = CCG(s,ids,'binSize',0.005,'duration',0.4,'smooth',1,'mode','ccg'); %ccg calculation
                                                time = tttt;
                                                cross_members.reward.pre = [cross_members.reward.pre , ccg(:,1,2)];
                                                clear x y s ids groups ccg tttt
                                            end
                                            
                                            % for post
                                            x = Restrict(SpikesD,post);
                                            y = Restrict(SpikesV,post);
                                            if and(length(x)>10,length(y)>10)
                                                [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                                                [ccg,tttt] = CCG(s,ids,'binSize',0.005,'duration',2,'smooth',1,'mode','ccg'); %ccg calculation
                                                time = tttt;
                                                cross_members.reward.post = [cross_members.reward.post , ccg(:,1,2)];
                                                clear x y s ids groups ccg tttt SpikesV
                                            end
                                        else
                                            SpikesV = spks_vHPC(spks_vHPC(:,1)==clusters.vHPC(iii),2);
                                            
                                            % For Pre
                                            x = Restrict(SpikesD,pre);
                                            y = Restrict(SpikesV,pre);
                                            if and(length(x)>5,length(y)>5)
                                                [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                                                [ccg,tttt] = CCG(s,ids,'binSize',0.005,'duration',2,'smooth',1,'mode','ccg'); %ccg calculation
                                                time = tttt;
                                                cross_nonmembers.reward.pre = [cross_nonmembers.reward.pre , ccg(:,1,2)];
                                                clear x y s ids groups ccg tttt
                                            end
                                            
                                            % for post
                                            x = Restrict(SpikesD,post);
                                            y = Restrict(SpikesV,post);
                                            if and(length(x)>10,length(y)>10)
                                                [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                                                [ccg,tttt] = CCG(s,ids,'binSize',0.005,'duration',2,'smooth',1,'mode','ccg'); %ccg calculation
                                                time = tttt;
                                                cross_nonmembers.reward.post = [cross_nonmembers.reward.post , ccg(:,1,2)];
                                                clear x y s ids groups ccg tttt SpikesV
                                            end
                                        end
                                    end
                                    clear criteriaV
                                end
                            end
                        end
                        clear criteriaD
                    end
                    clear parameter weigth
                end
                clear p w q
            end
                
            clear aversiveTS aversiveTS_run baselineTS bins Cell_type_classification
            clear cellulartype clusters Events group_dHPC group_vHPC i ii iii K
            clear Kinfo NREM REM numberD numberV p parameter patterns post pre
            clear SpikesD SpikesV spks spks_dHPC spks_vHPC Thresholded WAKE
            clear Kinfo K ii iii ripples ripplesD ripplesV ripple_event
        end
    end

%Detection of lag where maximal value was detected
%Aversive
[h p] = max(cross_members.aversive.pre);
lags.aversive.pre = [];
for i = 1 : size(p,2)
    lags.aversive.pre = [lags.aversive.pre ; time(p(i))];
end

[h p] = max(cross_members.aversive.post);
lags.aversive.post = [];
for i = 1 : size(p,2)
    lags.aversive.post = [lags.aversive.post ; time(p(i))];
end

% Reward
[h p] = max(cross_members.reward.pre);
lags.reward.pre = [];
for i = 1 : size(p,2)
    lags.reward.pre = [lags.reward.pre ; time(p(i))];
end

[h p] = max(cross_members.reward.post);
lags.reward.post = [];
for i = 1 : size(p)
    lags.reward.post = [lags.reward.post ; time(p(i))];
end



end