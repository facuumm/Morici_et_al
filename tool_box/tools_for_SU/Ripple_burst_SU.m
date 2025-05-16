clear
clc
close all
%% PARAMETERS - run this section before run any other
path = {'\\Maryjackson\e\Rat127\Ephys\pyr';'\\Maryjackson\e\Rat128\Ephys\in_pyr\ready';'\\Maryjackson\e\Rat103\usable';'\\Maryjackson\e\Rat132\recordings\in_pyr'; '\\Maryjackson\e\Rat165\in_pyr'; ...
    '\\Maryjackson\e\Rat126\Ephys\in_Pyr'};%List of folders from the path
% path = {'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat103\usable';'E:\Rat132\recordings\in_pyr'};%List of folders from the path
% for SU
criteria_fr = 0.01; %criteria to include or not a SU into the analysis
criteria_n = 6; % minimal number of neurons from each structure
criteria_type = 0; %criteria for celltype (0:pyr, 1:int, 2:all)
n_SU_V = [];
n_SU_D = [];

binSize = 0.001; % bin size for replay events detection
bin_size = 1; % to bin pos ans spks in between/within  
% Behavior
minimal_speed = 2.5;% minimal speed to detect quite periods
minimal_speed_time = 2; % minimal time to detect quite periods
%Place cells map
%Parameters for FiringCurve.m
min_time = 0.15; %min time in a given bin to be consider
min_peak = 0.2; %
min_size = 4;% min number of bins for the place fields
Xedges = 60; %number of bins for RateMap construction - 2.5cm bin (because of the extremes cutting)
sigma = 2;%round(15/(180/Xedges)); %defined for gauss kernel of 15cm
%% MAIN LOOP ORIGINAL, to iterate across sessions  
ccg_dhpc = [];ccg_dhpc_out = [];
ccg_vhpc = [];ccg_vhpc_out = [];

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
        clear y
        %% Load REM segments 
        disp('Uploading sleep scoring')
        x = dir([cd,'\*-states.mat']);    states = load([cd,'\',x.name]);    states = states.states;
        REM.all = ToIntervals(states==5);    NREM.all = ToIntervals(states==3);
        clear x states
        baselineTS = baselineTS./1000;                  aversiveTS = aversiveTS./1000;                    rewardTS = rewardTS./1000;
        aversiveTS_run = aversiveTS_run./1000;          rewardTS_run = rewardTS_run./1000;
        NREM.baseline = Restrict(NREM.all,baselineTS);   NREM.aversive = Restrict(NREM.all,aversiveTS);   NREM.reward = Restrict(NREM.all,rewardTS);
        REM.baseline = Restrict(REM.all,baselineTS);     REM.aversive = Restrict(REM.all,aversiveTS);     REM.reward = Restrict(REM.all,rewardTS);
        
        %% Load ALL Ripples 
        if exist('ripplesD_customized2.csv','file')
            ripplesD = table2array(readtable('ripplesD_customized2.csv'));
            RD = true;
        else
            RD = false;
        end
        
        if exist('ripplesV_customized2.csv','file')
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
                    %                     if r(2)<z(2) % keep only when dorsal happen first
                    % %                         cooridnated_event = [cooridnated_event ; r];
                    %                         coordinated = [coordinated ; r];
                    %                         coordinatedV = [coordinatedV ; z];
                    %                     end
                    peak = min(min(z(indice,2),r(2)))+abs(z(indice,2)-r(2))/2;
                    low = min([r(1) , z(indice,1)]);
                    up = max([r(3) , z(indice,3)]);
                    cooridnated_event = [cooridnated_event ; low , peak , up];
                    clear tmp2 tmp1 p indice z peak low up
                 end
                clear r
             end
        end
        
        % Eliminate nan
        ripplesD(any(isnan(ripplesD), 2), :) = [];
        ripplesV(any(isnan(ripplesV), 2), :) = [];
        %% Load burst ripples info 
         try
            load('coordinated_ripple_bursts.mat');
         catch 
             disp('No burst file')
             %Jump to next session 
         end 
         
         % Select ripples outside bursts
          [status,~,~] = InIntervals(ripplesD(:,2),bursts.dHPC.events(:,[1,3]));
          ripplesD_out = ripplesD(~status,:); 
        
          
          [status,~,~] = InIntervals(ripplesV(:,2),bursts.vHPC.events(:,[1,3]));
          ripplesV_out = ripplesV(~status,:); 
       
        %% Spikes
        %Load Units
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
                n_SU_D = [n_SU_D ; length(group_dHPC)];
                spks_dHPC = spks(ismember(spks(:,1),group_dHPC(:,1)),:); %keep spks from good clusters
            else
                n_SU_V = [n_SU_V ; length(group_vHPC)];
                spks_vHPC = spks(ismember(spks(:,1),group_vHPC(:,1)),:); %keep spks from good clusters
            end
        end
        clear z
        spks_vHPC(:,2) = double(spks_vHPC(:,2))./20000;
        spks_dHPC(:,2) = double(spks_dHPC(:,2))./20000;
        % Selection of celltype to analyze
        if criteria_type == 0 %pyr
            cellulartype = [K(:,1) , K(:,4)];
        elseif criteria_type == 1 % int
            cellulartype = [K(:,1) , not(K(:,4))];
        elseif criteria_type == 2 % all
            cellulartype = [K(:,1) , ones(length(K),1)];
        end
        % Definition of the emotional transition
        if aversiveTS_run(1)<rewardTS_run(1)
            cond = 1; % 1 if the order was Aversive -> Reward
        else
            cond = 2;% 1 if the order was Reward -> Aversive
        end
        
        %% Firing maps calculation
        disp('dHPC CCG calculation')
        %Load pc ID
        try
            load('dHPC_pc_skaggs_circular');
            dhpc_pc = [];
            for d=1:size(dHPC,2)
                pc = dHPC{d}; 
                dhpc_pc = [dhpc_pc;pc.id]; 
            end
        catch 
            dhpc_pc = [];
            disp('No pc in this session'); 
        end    
            
        for ii=1:size(group_dHPC,1)
            %Load spk  
            cluster = group_dHPC(ii,1);
            spks = spks_dHPC(spks_dHPC(:,1)==cluster,2); % select tspk from cluster
            
            if any(dhpc_pc==cluster) % check if it is place cell
                id=1; 
            else 
                id=0;
            end
            
            %Calculate hist aligned to burst members ripple peak 
            times1=bursts.dHPC.members(:,1:3);  times2=spks;
            baseline=NREM.all; %sec
            d=2; b=0.01;sm=2;
            [p , t_window] = PHIST_Ripple_SU(times1,times2,baseline,d,b,sm,[]); 
            [m] = meanFR_outside_ripples(ripplesD , baseline , times2 , b);% mean outside ripple 
            p = p./m;
            ccg_dhpc =[ccg_dhpc; [p',id]];%save output
            
%             Calculate hist aligned to NO burst members ripple peak 
            times1=ripplesD_out(:,1:3); 
            [p ,  t_window] = PHIST_Ripple_SU(times1,times2,baseline,d,b,sm,[]);
            p = p./m;
            ccg_dhpc_out =[ccg_dhpc_out; [p',id]];%save output
           
            
        end
        
        
        disp('vHPC CCG calculation')
        
        %Load pc ID
        try
            load('vHPC_pc_skaggs_circular.mat');
            vhpc_pc = [];
            for d=1:size(vHPC,2)
                pc = vHPC{d}; 
                vhpc_pc = [vhpc_pc;pc.id]; 
            end
        catch
            vhpc_pc = [];
            disp('No pc in this session'); 
        end    
            
        
        for ii=1:size(group_vHPC,1)
            %Load spk  
            cluster = group_vHPC(ii,1);
            spks = spks_vHPC(spks_vHPC(:,1)==cluster,2); % select tspk from cluster
            
            if any(vhpc_pc==cluster) % check if it is place cell
                id=1; 
            else 
                id=0;
            end
            %Calculate hist aligned to burst members ripple peak 
            times1=bursts.vHPC.members(:,1:3);  times2=spks;
            baseline=NREM.all; %sec
            d=2; b=0.01;sm=2;
            [p ,  t_window] = PHIST_Ripple_SU(times1,times2,baseline,d,b,sm,[]); 
            [m] = meanFR_outside_ripples(ripplesD , baseline , times2 , b);
            p = p./m;
            ccg_vhpc =[ccg_vhpc; [p',id]];% Save output 
            
%             Calculate hist aligned to NO burst members ripple peak 
            times1=ripplesV_out(:,1:3); 
            [p ,  t_window] = PHIST_Ripple_SU(times1,times2,baseline,d,b,sm,[]);
            p = p./m;
            ccg_vhpc_out =[ccg_vhpc_out; [p',id]];%save output
           
        end
        
     
    disp(['-- Finished folder #' , num2str(t) , ' from rat #' , num2str(tt) , ' ---'])
    disp('  ')
    end
    disp(['-------- Finished rat#' , num2str(tt) , ' --------'])
    disp(' ')
end


%% PLOTS

figure(2);clf;hold on; 
title('dhpc Norm gain pc');
h=plot([1:1:size(ccg_dhpc,2)-1],nanmean(ccg_dhpc(ccg_dhpc(:,202)==1,[1:end-1]))); axis tight; set(h, 'Color', 'r');
ciplot(nanmean(ccg_dhpc(ccg_dhpc(:,202)==1,[1:end-1]))-nansem(ccg_dhpc(ccg_dhpc(:,202)==1,[1:end-1])),...
nanmean(ccg_dhpc(ccg_dhpc(:,202)==1,[1:end-1]))+nansem(ccg_dhpc(ccg_dhpc(:,202)==1,[1:end-1])),[1:1:size(ccg_dhpc,2)-1],'r'); alpha 0.1; hold on;

h=plot([1:1:size(ccg_dhpc_out,2)-1],nanmean(ccg_dhpc_out(ccg_dhpc_out(:,202)==1,[1:end-1]))); axis tight;set(h, 'Color', 'k');
ciplot(nanmean(ccg_dhpc_out(ccg_dhpc_out(:,202)==1,[1:end-1]))-nansem(ccg_dhpc_out(ccg_dhpc_out(:,202)==1,[1:end-1])),...
nanmean(ccg_dhpc_out(ccg_dhpc_out(:,202)==1,[1:end-1]))+nansem(ccg_dhpc_out(ccg_dhpc_out(:,202)==1,[1:end-1])),[1:1:size(ccg_dhpc_out,2)-1],'k'); alpha 0.1; hold on;


figure(3);clf;hold on; 
title('vhpc gain pc')
h=plot([1:1:size(ccg_vhpc,2)-1],nanmean(ccg_vhpc(ccg_vhpc(:,202)==1,[1:end-1]))); axis tight;  set(h, 'Color', 'r');
ciplot(nanmean(ccg_vhpc(ccg_vhpc(:,202)==1,[1:end-1]))-nansem(ccg_vhpc(ccg_vhpc(:,202)==1,[1:end-1])),nanmean(ccg_vhpc(ccg_vhpc(:,202)==1,[1:end-1]))+nansem(ccg_vhpc(ccg_vhpc(:,202)==1,[1:end-1])),[1:1:size(ccg_vhpc,2)-1],'r'); alpha 0.1; hold on;
 
h=plot([1:1:size(ccg_vhpc_out,2)-1],nanmean(ccg_vhpc_out(ccg_vhpc_out(:,202)==1,[1:end-1]))); axis tight;set(h, 'Color', 'k');
ciplot(nanmean(ccg_vhpc_out(ccg_vhpc_out(:,202)==1,[1:end-1]))-nansem(ccg_vhpc_out(ccg_vhpc_out(:,202)==1,[1:end-1])),nanmean(ccg_vhpc_out(ccg_vhpc_out(:,202)==1,[1:end-1]))+nansem(ccg_vhpc_out(ccg_vhpc_out(:,202)==1,[1:end-1])),[1:1:size(ccg_vhpc_out,2)-1],'k'); alpha 0.1; hold on;


figure(4);clf;hold on; 
title('dhpc gain ALL')
plot([1:1:size(ccg_dhpc,2)-1],nanmean(ccg_dhpc(:,[1:end-1]))); axis tight
ciplot(nanmean(ccg_dhpc(:,[1:end-1]))-nansem(ccg_dhpc(:,[1:end-1])),nanmean(ccg_dhpc(:,[1:end-1]))+nansem(ccg_dhpc(:,[1:end-1])),[1:1:size(ccg_dhpc,2)-1]); alpha 0.1; hold on;
 
plot([1:1:size(ccg_dhpc_out,2)-1],nanmean(ccg_dhpc_out(:,[1:end-1]))); axis tight
ciplot(nanmean(ccg_dhpc_out(:,[1:end-1]))-nansem(ccg_dhpc_out(:,[1:end-1])),nanmean(ccg_dhpc_out(:,[1:end-1]))+nansem(ccg_dhpc_out(:,[1:end-1])),[1:1:size(ccg_dhpc_out,2)-1]); alpha 0.1; hold on;


figure(5);clf;hold on; 
title('vhpc gain ALL')
plot([1:1:size(ccg_vhpc,2)-1],nanmean(ccg_vhpc(:,[1:end-1]))); axis tight
ciplot(nanmean(ccg_vhpc(:,[1:end-1]))-nansem(ccg_vhpc(:,[1:end-1])),nanmean(ccg_vhpc(:,[1:end-1]))+nansem(ccg_vhpc(:,[1:end-1])),[1:1:size(ccg_vhpc,2)-1]); alpha 0.1; hold on;
 
plot([1:1:size(ccg_vhpc_out,2)-1],nanmean(ccg_vhpc_out(:,[1:end-1]))); axis tight
ciplot(nanmean(ccg_vhpc_out(:,[1:end-1]))-nansem(ccg_vhpc_out(:,[1:end-1])),nanmean(ccg_vhpc_out(:,[1:end-1]))+nansem(ccg_vhpc_out(:,[1:end-1])),[1:1:size(ccg_vhpc_out,2)-1]); alpha 0.1; hold on;


%% MAIN LOOP: amplitud and duration of burst ripples  
r_dhpc = [];r_dhpc_out = [];
r_vhpc = [];r_vhpc_out = [];

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
        clear y
        %% Load REM segments 
%         disp('Uploading sleep scoring')
%         x = dir([cd,'\*-states.mat']);    states = load([cd,'\',x.name]);    states = states.states;
%         REM.all = ToIntervals(states==5);    NREM.all = ToIntervals(states==3);
%         clear x states
%         baselineTS = baselineTS./1000;                  aversiveTS = aversiveTS./1000;                    rewardTS = rewardTS./1000;
%         aversiveTS_run = aversiveTS_run./1000;          rewardTS_run = rewardTS_run./1000;
%         NREM.baseline = Restrict(NREM.all,baselineTS);   NREM.aversive = Restrict(NREM.all,aversiveTS);   NREM.reward = Restrict(NREM.all,rewardTS);
%         REM.baseline = Restrict(REM.all,baselineTS);     REM.aversive = Restrict(REM.all,aversiveTS);     REM.reward = Restrict(REM.all,rewardTS);
%         
        %% Load ALL Ripples 
        if exist('ripplesD_customized2.csv','file')
            ripplesD = table2array(readtable('ripplesD_customized2.csv'));
            RD = true;
        else
            RD = false;
        end
        
        if exist('ripplesV_customized2.csv','file')
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
                    %                     if r(2)<z(2) % keep only when dorsal happen first
                    % %                         cooridnated_event = [cooridnated_event ; r];
                    %                         coordinated = [coordinated ; r];
                    %                         coordinatedV = [coordinatedV ; z];
                    %                     end
                    peak = min(min(z(indice,2),r(2)))+abs(z(indice,2)-r(2))/2;
                    low = min([r(1) , z(indice,1)]);
                    up = max([r(3) , z(indice,3)]);
                    cooridnated_event = [cooridnated_event ; low , peak , up];
                    clear tmp2 tmp1 p indice z peak low up
                 end
                clear r
             end
        end
        
        % Eliminate nan
        ripplesD(any(isnan(ripplesD), 2), :) = [];
        ripplesV(any(isnan(ripplesV), 2), :) = [];
        %% Load burst ripples info 
         try
            load('coordinated_ripple_bursts.mat');
         catch 
             disp('No burst file')
             %Jump to next session 
         end 
         
         % Select ripples outside bursts
          [status,~,~] = InIntervals(ripplesD(:,2),bursts.dHPC.events(:,[1,3]));
          ripplesD_out = ripplesD(~status,:); 
        
          
          [status,~,~] = InIntervals(ripplesV(:,2),bursts.vHPC.events(:,[1,3]));
          ripplesV_out = ripplesV(~status,:); 
       
        
        %% Ripples properties calculation
        
         r_dhpc = [r_dhpc; [bursts.dHPC.members(:,4),(bursts.dHPC.members(:,3)-bursts.dHPC.members(:,1))*1000]];
         r_dhpc_out = [r_dhpc_out;[ ripplesD_out(:,4), (ripplesD_out(:,3)- ripplesD_out(:,1))*1000]];
         
         r_vhpc = [r_vhpc; [bursts.vHPC.members(:,4),(bursts.vHPC.members(:,3)-bursts.vHPC.members(:,1))*1000]];
         r_vhpc_out = [r_dhpc_out;[ ripplesV_out(:,4), (ripplesV_out(:,3)- ripplesV_out(:,1))*1000]];
     
    disp(['-- Finished folder #' , num2str(t) , ' from rat #' , num2str(tt) , ' ---'])
    disp('  ')
    end
    disp(['-------- Finished rat#' , num2str(tt) , ' --------'])
    disp(' ')
end

%%%%%%%%% Plots and stats%%%%%%%%%%%%%

%dhpc
outliers = isoutlier(r_dhpc_out(:,1));
out= r_dhpc_out(~outliers,1);
outliers = isoutlier(r_dhpc(:,1));
burst= r_dhpc(~outliers,1);
%%% HISTOGRAM
figure(1);clf;hold on;
sgtitle('dHPC')
ylabel('Probability'); xlabel('Amplitud')
histogram(out,'Normalization', 'probability','FaceColor','k','EdgeColor','none','BinWidth',2);
histogram(burst, 'Normalization', 'probability','FaceColor','g','EdgeColor','none','BinWidth',2);
%%% box plot
figure(2);clf;hold on; 
x = [ones(size(burst,1),1);ones(size(out,1),1)*2];
y= [burst;out]; c = [.3,.3,.3];
scatter(x,y,[],c,"filled",'jitter','on', 'jitterAmount',0.1); xlim([0.5 2.5]);
ylabel('Amplitud');xticks([1 2]);xticklabels({'Burst', 'No burst'});
hold on;x = [1 2];
y = [nanmedian(burst),nanmedian(out)]; 
scatter(x,y, "filled");

%Stats

[p,h,stats] = ranksum(out,burst);
 
[P,ANOVATAB,STATS] = kruskalwallis(y,x);
c = multcompare(STATS)

tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])

%%%%%%%%%%%%%%%%%%%%%
%vhpc
outliers = isoutlier(r_vhpc_out(:,1));
out= r_vhpc_out(~outliers,1);
outliers = isoutlier(r_vhpc(:,1));
burst= r_vhpc(~outliers,1);

figure(6);clf;hold on;
sgtitle('vHPC')
ylabel('Probability'); xlabel('Amplitud')
histogram(out,'Normalization', 'probability','FaceColor','k','EdgeColor','none','BinWidth',5);
histogram(burst,'Normalization', 'probability','FaceColor','b','EdgeColor','none','BinWidth',5);

x = [ones(size(burst,1),1);ones(size(out,1),1)*2];
y= [burst;out];
[P,ANOVATAB,STATS] = kruskalwallis(y,x);
c = multcompare(STATS)

tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])

%%%%%%%%%%%%%%%%
outliers = isoutlier(r_dhpc_out(:,2));
out= r_dhpc_out(~outliers,2);
outliers = isoutlier(r_dhpc(:,2));
burst= r_dhpc(~outliers,2);

figure(4);clf;hold on;
sgtitle('dHPC')
ylabel('Probability'); xlabel('Duration(ms)')
histogram(out, 'Normalization', 'probability','FaceColor','k','EdgeColor','none','BinWidth',3);
histogram(burst, 'Normalization', 'probability','FaceColor','g','EdgeColor','none','BinWidth',3);

[p,h,stats] = ranksum(out,burst);
%%%%%%%%%%%%%%%%
outliers = isoutlier(r_vhpc_out(:,2));
out= r_vhpc_out(~outliers,2);
outliers = isoutlier(r_vhpc(:,2));
burst= r_vhpc(~outliers,2);
figure(5);clf;hold on;
sgtitle('vHPC')
ylabel('Probability'); xlabel('Duration(ms)')
histogram(out,'Normalization', 'probability','FaceColor','k','EdgeColor','none','BinWidth',3);
histogram(burst, 'Normalization', 'probability','FaceColor','b','EdgeColor','none','BinWidth',3);

[p,h,stats] = ranksum(out,burst)

