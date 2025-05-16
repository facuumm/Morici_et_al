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

Probability.aversive.Space = [];
Probability.aversive.Shock = [];

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
        % Awake
        disp('Uploading digital imputs')
%         load('behavioral_data.mat')
        
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
        
        
        %% Decoding
        if and(and(RD,RV),and(numberD>3 , numberV>3)) %check if I have a minimun of cells and also ripples
            % SpikeTrain Construction
%             limits = [0 segments.Var1(end)/1000];
            limits = aversiveTS./1000;
            events = NREM.aversive;
            [Spikes , bins , Clusters] = spike_train_construction([spks_dHPC;spks_vHPC], clusters.all, cellulartype, 0.025, limits, events, false, false);
            clear limits events
            
            % load FiringCurves
            % SpaceMaps
            load('dHPC_all.mat')
            load('vHPC_all.mat')
            % Shock-responsivness curves
            load('dHPC_shock.mat')
            load('vHPC_shock.mat')
            
%             % PC ids
%             id.d = [];
%             for i = 1:size(dHPC,2)
%                 tmp = dHPC{i};
%                 id.d = [id.d ; tmp.id];
%                 clear tmp
%             end
%             id.v = [];
%             for i = 1:size(vHPC,2)
%                 tmp = vHPC{i};
%                 id.v = [id.v ; tmp.id];
%                 clear tmp
%             end
%             
            % RateMaps
            pc.aversive = zeros(60,numberD+numberV);            pc.reward = zeros(60,numberD+numberV);
            for i = 1:size(clusters.dHPC,1)
%                 if ismember(clusters.dHPC(i),id.d)
%                     tmp = (id.d == clusters.dHPC(i));
%                     tmp = dHPC{tmp};
                    tmp = dHPC{i};
                    pc.aversive(:,i) = tmp.frMap_ave';
                    pc.reward(:,i) = tmp.frMap_rew';
                    clear tmp
%                 end
            end
            for i = 1:size(clusters.vHPC,1)
%                 if ismember(clusters.vHPC(i),id.v)
%                     tmp = (id.v == clusters.vHPC(i));
%                     tmp = vHPC{tmp};
                    tmp = vHPC{i};
                    pc.aversive(:,i+numberD) = tmp.frMap_ave';
                    pc.reward(:,i+numberD) = tmp.frMap_rew';
                    clear tmp
%                 end
            end
            
%             clear aversiveTS aversiveTS_run baselineTS behavior Cell_type_classification
%             clear cellulartype clusters Clusters coordinated coordinatedV coordinatedV_refined
%             clear cooridnated_event cooridnated_eventDV cooridnated_eventVD dHPC group_dHPC group_vHPC
%             clear movement NREM REM Rewards_filt rewardTS rewardTS_run
%             clear ripples ripplesD ripplesV segments Shocks_filt spks spks_dHPC spks_vHPC vHPC WAKE
%             
            % Bayesian decoding
            % Space
            [Pr prMax] = placeBayes(Spikes./binSize , pc.aversive' , binSize);
%             
%             r = ripple_event.aversive(:,2);
%             tmp = [];
%             for i = 1 : size(r,1)
%                 IN = InIntervals(bins,[r(i)-2 r(i)+2]);
%                 if not(sum(IN)>=160)
%                     IN = InIntervals(bins,[r(i)-2 r(i)+2.2]);
%                 end
%                 T = max(Pr(IN,:)');
%                 tmp = [tmp , T(1:160)'];
%             end
            
            
            curvesShocks = [dHPC_shock.curve , vHPC_shock.curve];
            % Shock
            [Pr prMax] = placeBayes(Spikes./binSize , curvesShocks' , binSize);
            
            r = ripple_event.aversive(:,2);
            tmp1 = [];
            for i = 1 : size(r,1)
                IN = InIntervals(bins,[r(i)-2 r(i)+2]);
                if not(sum(IN)>=160)
                    IN = InIntervals(bins,[r(i)-2 r(i)+3]);
                    if sum(IN)<160
                        continue
                    end
                end
%                 temporal = nanmean([Pr(IN,1:40) , Pr(IN,121:160)]');
%                 T = max(Pr(IN,:)');
                T =  max(Pr(IN,80:120)'); clear temporal
%                 T(isnan(T)) = 0;
%                 T = Smooth(T,1)
                tmp1 = [tmp1 , T(1:160)'];
            end

%             Probability.aversive.Space = [ Probability.aversive.Space , nanmean(tmp')'];
            Probability.aversive.Shock = [Probability.aversive.Shock , nanmean(tmp1')'];
            
%             P = nanmean(tmp')-nanmean(nanmean(tmp'));
%             P1 = nanmean(tmp1')-nanmean(nanmean(tmp1'));
            
            
        end
        clear aversiveTS aversiveTS_run baselineTS behavior bins Cell_type_classification curvesShocks
        clear dHPC vHPC dHPC_shock vHPC_shock numberD numberV Pr prMax RD RV REM Rewards_filt rewardTS rewardTS_run
        clear ripplesD ripplesV segments spks spks_dHPC spks_vHPC Spikes Shocks_filt
        clear tmp tmp1 ripples ripple_event pc NREM K IN Kinfo group_dHPC group_vHPC
        clear config cellulartype clusters Clusters coordinated coordinatedV coordinatedV_refined
        clear cooridnated_event cooridnated_eventDV cooridnated_eventVD movement WAKE
    end
end


figure
tmp = nanmean([Probability.aversive.Space(1:40,:) ; Probability.aversive.Space(120:160,:)]);
plot(nanmean((Probability.aversive.Space./tmp)')),hold on
tmp = nanmean([Probability.aversive.Shock(1:40,:) ; Probability.aversive.Shock(120:160,:)]);
plot(nanmean((Probability.aversive.Shock)')),hold on
