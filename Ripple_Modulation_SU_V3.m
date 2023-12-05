clear
clc
close all

%% Parameters
path = {'E:\Rat103\usable';'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat132\recordings\in_pyr';'E:\Rat165\in_pyr\'};%List of folders from the path
% for SU
criteria_fr = 0; %criteria to include or not a SU into the analysis
pval = 0.001;

temp.dHPC.all = []; temp.dHPC.aversive = []; temp.dHPC.reward = []; temp.dHPC.baseline = [];
temp.vHPC.all = []; temp.vHPC.aversive = []; temp.vHPC.reward = []; temp.vHPC.baseline = [];

% Sacar el filtro que puse del FR en el counts de neuronas
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
        
        %% load sleep states
        disp('Uploading sleep scoring')
        x = dir([cd,'\*-states.mat']);    states = load([cd,'\',x.name]);    states = states.states;
        REM.all = ToIntervals(states==5);    NREM.all = ToIntervals(states==3);    WAKE.all = ToIntervals(states==1);
        REM.all(REM.all(:,2)-REM.all(:,1)<60,:) = [];
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
        
        %% load coordinated ripple bursts
        load('coordinated_ripple_bursts.mat')
        coordinated_ripple_bursts = [coordinated_ripple_bursts(:,1)  coordinated_ripple_bursts(:,3)];
        %         [coordinated_ripple_bursts] = merge_events(coordinated_ripple_bursts, 0.05);
        
        ripple_bursts.baseline = Restrict(coordinated_ripple_bursts,baselineTS./1000);
        ripple_bursts.reward = Restrict(coordinated_ripple_bursts,rewardTS./1000);
        ripple_bursts.aversive = Restrict(coordinated_ripple_bursts,aversiveTS./1000);
        ripple_bursts.all = coordinated_ripple_bursts;
        clear coordinated_ripple_bursts
        
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
%         % Selection of celltype to analyze
%         if criteria_type == 0 %pyr
%             cellulartype = [K(:,1) , K(:,3)];
%         elseif criteria_type == 1 % int
%             cellulartype = [K(:,1) , K(:,4)];
%         elseif criteria_type == 2 % all
%             cellulartype = [K(:,1) , ones(length(K),1)];
%         end
        
        %% Counting the Number f SU
        numberD = 0;
        clusters.all = [];
        clusters.dHPC = [];
        for ii=1:size(group_dHPC,1)
            cluster = group_dHPC(ii,1);
%             Cellulartype = logical(cellulartype(cellulartype(:,1) == cluster,2));
%             if Cellulartype
                a = length(Restrict(spks(spks(:,1)==cluster,2),aversiveTS_run./1000)) / ((aversiveTS_run(2)-aversiveTS_run(1))/1000);
                r = length(Restrict(spks(spks(:,1)==cluster,2),rewardTS_run./1000)) / ((rewardTS_run(2)-rewardTS_run(1))/1000);
                if or(a > criteria_fr ,  r > criteria_fr)
                    numberD = numberD+1;
                    clusters.all = [clusters.all ; cluster];
                    clusters.dHPC = [clusters.dHPC ; cluster];
                end
%             end
            clear tmp cluster Cellulartype fr1 fr2 fr3 fr4 fr5 r a
        end
        
        numberV = 0;
        clusters.vHPC = [];
        for ii=1:size(group_vHPC,1)
            cluster = group_vHPC(ii,1);
%             Cellulartype = logical(cellulartype(cellulartype(:,1) == cluster,2));
%             if Cellulartype
                a = length(Restrict(spks(spks(:,1)==cluster,2),aversiveTS_run./1000)) / ((aversiveTS_run(2)-aversiveTS_run(1))/1000);
                r = length(Restrict(spks(spks(:,1)==cluster,2),rewardTS_run./1000)) / ((rewardTS_run(2)-rewardTS_run(1))/1000);
                if or(a > criteria_fr ,  r > criteria_fr)
                    numberV = numberV+1;
                    clusters.all = [clusters.all ; cluster];
                    clusters.vHPC = [clusters.vHPC ; cluster];
                end
%             end
            clear tmp cluster Cellulartype fr1 fr2 fr3 fr4 fr5 r a
        end
        clear freq limits
        clear ejeX ejeY dX dY dX_int dY_int
        
        %%
        
        tmp.vHPC.all = []; tmp.vHPC.aversive = []; tmp.vHPC.reward = []; tmp.vHPC.baseline = [];
        if not(isempty(spks_dHPC))
            disp('PoissonTest for dHPC SUs')
            iterator = ismember(K(:,1),clusters.dHPC); %keep logical to differentiate pyr from rest
            iterator  = K(iterator,3);
            [pI pD sur] = RippleModulation(ripplesD,spks_dHPC,clusters.dHPC,NREM.all);
            ripple_modulated.dHPC.all = [clusters.dHPC iterator pI'<=pval pD'<=pval];
            clear pI pD sur
            
            [pI pD sur] = RippleModulation(ripplesD,spks_dHPC,clusters.dHPC,NREM.baseline);
            ripple_modulated.dHPC.baseline = [clusters.dHPC iterator pI'<=pval pD'<=pval];
            clear pI pD sur
            
            [pI pD sur] = RippleModulation(ripplesD,spks_dHPC,clusters.dHPC,NREM.aversive);
            ripple_modulated.dHPC.aversive = [clusters.dHPC iterator pI'<=pval pD'<=pval];
            clear pI pD sur
            
            
            [pI pD sur] = RippleModulation(ripplesD,spks_dHPC,clusters.dHPC,NREM.reward);
            ripple_modulated.dHPC.reward = [clusters.dHPC iterator pI'<=pval pD'<=pval];
            clear pI pD sur iterator   
            
            temp.dHPC.all = [temp.dHPC.all ; ripple_modulated.dHPC.all];
            temp.dHPC.aversive = [temp.dHPC.aversive ; ripple_modulated.dHPC.aversive];
            temp.dHPC.reward = [temp.dHPC.reward ; ripple_modulated.dHPC.reward];
            temp.dHPC.baseline = [temp.dHPC.baseline ; ripple_modulated.dHPC.baseline];
        end
        
        if not(isempty(spks_vHPC))
            disp('PoissonTest for vHPC SUs')
            iterator = ismember(K(:,1),clusters.vHPC);%keep logical to differentiate pyr from rest
            iterator  = K(iterator,3);
            [pI pD sur] = RippleModulation(ripplesV,spks_vHPC,clusters.vHPC,NREM.all);
            ripple_modulated.vHPC.all = [clusters.vHPC iterator pI'<=pval pD'<=pval];
            clear pI pD sur
            
            [pI pD sur] = RippleModulation(ripplesV,spks_vHPC,clusters.vHPC,NREM.aversive);
            ripple_modulated.vHPC.aversive = [clusters.vHPC iterator pI'<=pval pD'<=pval];
            clear pI pD sur            
            
            [pI pD sur] = RippleModulation(ripplesV,spks_vHPC,clusters.vHPC,NREM.reward);
            ripple_modulated.vHPC.reward = [clusters.vHPC iterator pI'<=pval pD'<=pval];
            clear pI pD sur
            
            [pI pD sur] = RippleModulation(ripplesV,spks_vHPC,clusters.vHPC,NREM.baseline);
            ripple_modulated.vHPC.baseline = [clusters.vHPC iterator pI'<=pval pD'<=pval];
            clear pI pD sur iterator
            
            temp.vHPC.all = [temp.vHPC.all ; ripple_modulated.vHPC.all];
            temp.vHPC.aversive = [temp.vHPC.aversive ; ripple_modulated.vHPC.aversive];
            temp.vHPC.reward = [temp.vHPC.reward ; ripple_modulated.vHPC.reward];
            temp.vHPC.baseline = [temp.vHPC.baseline ; ripple_modulated.vHPC.baseline];
        end
        
        save([cd,'\ripple_modulated_SU.mat'],'ripple_modulated')
        clear ripple_modulated
        
    end
    disp(['-- Finishing analysis from rat #',num2str(tt) , ' --'])
    disp('  ')
end


dHPC.pyr = temp.dHPC.all(:,2)==1;
dHPC.int = temp.dHPC.all(:,2)==0;
vHPC.pyr = temp.vHPC.all(:,2)==1;
vHPC.int = temp.vHPC.all(:,2)==0;

%% Percentage calculation dHPC pyr
% all
n = length(temp.dHPC.all(dHPC.pyr));
i = sum(temp.dHPC.all(dHPC.pyr,3)==1);
percentage.dHPC.pyr.all.up = (i/n)*100;

i = sum(temp.dHPC.all(dHPC.pyr,4)==1);
percentage.dHPC.pyr.all.down = (i/n)*100;

i = sum(and(temp.dHPC.all(dHPC.pyr,4)==0 , temp.dHPC.all(dHPC.pyr,3)==0));
percentage.dHPC.pyr.all.non = (i/n)*100;

% baseline
n = length(temp.dHPC.baseline(dHPC.pyr));
i = sum(temp.dHPC.baseline(dHPC.pyr,3)==1);
percentage.dHPC.pyr.baseline.up = (i/n)*100;

i = sum(temp.dHPC.baseline(dHPC.pyr,4)==1);
percentage.dHPC.pyr.baseline.down = (i/n)*100;

i = sum(and(temp.dHPC.baseline(dHPC.pyr,4)==0 , temp.dHPC.baseline(dHPC.pyr,3)==0));
percentage.dHPC.pyr.baseline.non = (i/n)*100;

% reward
n = length(temp.dHPC.reward(dHPC.pyr));
i = sum(temp.dHPC.reward(dHPC.pyr,3)==1);
percentage.dHPC.pyr.reward.up = (i/n)*100;

i = sum(temp.dHPC.reward(dHPC.pyr,4)==1);
percentage.dHPC.pyr.reward.down = (i/n)*100;

i = sum(and(temp.dHPC.reward(dHPC.pyr,4)==0 , temp.dHPC.reward(dHPC.pyr,3)==0));
percentage.dHPC.pyr.reward.non = (i/n)*100;

% aversive
n = length(temp.dHPC.aversive(dHPC.pyr));
i = sum(temp.dHPC.aversive(dHPC.pyr,3)==1);
percentage.dHPC.pyr.aversive.up = (i/n)*100;

i = sum(temp.dHPC.aversive(dHPC.pyr,4)==1);
percentage.dHPC.pyr.aversive.down = (i/n)*100;

i = sum(and(temp.dHPC.aversive(dHPC.pyr,4)==0 , temp.dHPC.aversive(dHPC.pyr,3)==0));
percentage.dHPC.pyr.aversive.non = (i/n)*100;

%% Percentage calculation dHPC int
% all
n = length(temp.dHPC.all(dHPC.int));
i = sum(temp.dHPC.all(dHPC.int,3)==1);
percentage.dHPC.int.all.up = (i/n)*100;

i = sum(temp.dHPC.all(dHPC.int,4)==1);
percentage.dHPC.int.all.down = (i/n)*100;

i = sum(and(temp.dHPC.all(dHPC.int,4)==0 , temp.dHPC.all(dHPC.int,3)==0));
percentage.dHPC.int.all.non = (i/n)*100;

% baseline
n = length(temp.dHPC.baseline(dHPC.int));
i = sum(temp.dHPC.baseline(dHPC.int,3)==1);
percentage.dHPC.int.baseline.up = (i/n)*100;

i = sum(temp.dHPC.baseline(dHPC.int,4)==1);
percentage.dHPC.int.baseline.down = (i/n)*100;

i = sum(and(temp.dHPC.baseline(dHPC.int,4)==0 , temp.dHPC.baseline(dHPC.int,3)==0));
percentage.dHPC.int.baseline.non = (i/n)*100;

% reward
n = length(temp.dHPC.reward(dHPC.int));
i = sum(temp.dHPC.reward(dHPC.int,3)==1);
percentage.dHPC.int.reward.up = (i/n)*100;

i = sum(temp.dHPC.reward(dHPC.int,4)==1);
percentage.dHPC.int.reward.down = (i/n)*100;

i = sum(and(temp.dHPC.reward(dHPC.int,4)==0 , temp.dHPC.reward(dHPC.int,3)==0));
percentage.dHPC.int.reward.non = (i/n)*100;

% aversive
n = length(temp.dHPC.aversive(dHPC.int));
i = sum(temp.dHPC.aversive(dHPC.int,3)==1);
percentage.dHPC.int.aversive.up = (i/n)*100;

i = sum(temp.dHPC.aversive(dHPC.int,4)==1);
percentage.dHPC.int.aversive.down = (i/n)*100;

i = sum(and(temp.dHPC.aversive(dHPC.int,4)==0 , temp.dHPC.aversive(dHPC.int,3)==0));
percentage.dHPC.int.aversive.non = (i/n)*100;

%% Percentage calculation vHPC pyr
% all
n = length(temp.vHPC.all(vHPC.pyr));
i = sum(temp.dHPC.all(vHPC.pyr,3)==1);
percentage.vHPC.pyr.all.up = (i/n)*100;

i = sum(temp.vHPC.all(vHPC.pyr,4)==1);
percentage.vHPC.pyr.all.down = (i/n)*100;

i = sum(and(temp.vHPC.all(vHPC.pyr,4)==0 , temp.vHPC.all(vHPC.pyr,3)==0));
percentage.vHPC.pyr.all.non = (i/n)*100;

% baseline
n = length(temp.vHPC.baseline(vHPC.pyr));
i = sum(temp.vHPC.baseline(vHPC.pyr,3)==1);
percentage.vHPC.pyr.baseline.up = (i/n)*100;

i = sum(temp.vHPC.baseline(vHPC.pyr,4)==1);
percentage.vHPC.pyr.baseline.down = (i/n)*100;

i = sum(and(temp.vHPC.baseline(vHPC.pyr,4)==0 , temp.vHPC.baseline(vHPC.pyr,3)==0));
percentage.vHPC.pyr.baseline.non = (i/n)*100;

% reward
n = length(temp.vHPC.reward(vHPC.pyr));
i = sum(temp.vHPC.reward(vHPC.pyr,3)==1);
percentage.vHPC.pyr.reward.up = (i/n)*100;

i = sum(temp.vHPC.reward(vHPC.pyr,4)==1);
percentage.vHPC.pyr.reward.down = (i/n)*100;

i = sum(and(temp.vHPC.reward(vHPC.pyr,4)==0 , temp.vHPC.reward(vHPC.pyr,3)==0));
percentage.vHPC.pyr.reward.non = (i/n)*100;

% aversive
n = length(temp.vHPC.aversive(vHPC.pyr));
i = sum(temp.vHPC.aversive(vHPC.pyr,3)==1);
percentage.vHPC.pyr.aversive.up = (i/n)*100;

i = sum(temp.vHPC.aversive(vHPC.pyr,4)==1);
percentage.vHPC.pyr.aversive.down = (i/n)*100;

i = sum(and(temp.vHPC.aversive(vHPC.pyr,4)==0 , temp.vHPC.aversive(vHPC.pyr,3)==0));
percentage.vHPC.pyr.aversive.non = (i/n)*100;

%% Percentage calculation vHPC int
% all
n = length(temp.vHPC.all(vHPC.int));
i = sum(temp.vHPC.all(vHPC.int,3)==1);
percentage.vHPC.int.all.up = (i/n)*100;

i = sum(temp.vHPC.all(vHPC.int,4)==1);
percentage.vHPC.int.all.down = (i/n)*100;

i = sum(and(temp.vHPC.all(vHPC.int,4)==0 , temp.vHPC.all(vHPC.int,3)==0));
percentage.vHPC.int.all.non = (i/n)*100;

% baseline
n = length(temp.vHPC.baseline(vHPC.int));
i = sum(temp.vHPC.baseline(vHPC.int,3)==1);
percentage.vHPC.int.baseline.up = (i/n)*100;

i = sum(temp.vHPC.baseline(vHPC.int,4)==1);
percentage.vHPC.int.baseline.down = (i/n)*100;

i = sum(and(temp.vHPC.baseline(vHPC.int,4)==0 , temp.vHPC.baseline(vHPC.int,3)==0));
percentage.vHPC.int.baseline.non = (i/n)*100;

% reward
n = length(temp.vHPC.reward(vHPC.int));
i = sum(temp.vHPC.reward(vHPC.int,3)==1);
percentage.vHPC.int.reward.up = (i/n)*100;

i = sum(temp.vHPC.reward(vHPC.int,4)==1);
percentage.vHPC.int.reward.down = (i/n)*100;

i = sum(and(temp.vHPC.reward(vHPC.int,4)==0 , temp.vHPC.reward(vHPC.int,3)==0));
percentage.vHPC.int.reward.non = (i/n)*100;

% aversive
n = length(temp.vHPC.aversive(vHPC.int));
i = sum(temp.vHPC.aversive(vHPC.int,3)==1);
percentage.vHPC.int.aversive.up = (i/n)*100;

i = sum(temp.vHPC.aversive(vHPC.int,4)==1);
percentage.vHPC.int.aversive.down = (i/n)*100;

i = sum(and(temp.vHPC.aversive(vHPC.int,4)==0 , temp.vHPC.aversive(vHPC.int,3)==0));
percentage.vHPC.int.aversive.non = (i/n)*100;



%% Plots

%All
figure,
x = [percentage.dHPC.pyr.all.up percentage.dHPC.pyr.all.down percentage.dHPC.pyr.all.non];
subplot(2,2,1),pie(x,{'up' , 'down' , 'non'})

x = [percentage.vHPC.pyr.all.up percentage.vHPC.pyr.all.down percentage.vHPC.pyr.all.non];
subplot(2,2,2),pie(x,{'up' , 'down' , 'non'})

x = [percentage.dHPC.int.all.up percentage.dHPC.int.all.down percentage.dHPC.int.all.non];
subplot(2,2,3),pie(x,{'up' , 'down' , 'non'})

x = [percentage.vHPC.int.all.up percentage.vHPC.int.all.down percentage.vHPC.int.all.non];
subplot(2,2,4),pie(x,{'up' , 'down' , 'non'})

%Baseline
figure,
x = [percentage.dHPC.pyr.baseline.up percentage.dHPC.pyr.baseline.down percentage.dHPC.pyr.baseline.non];
subplot(2,2,1),pie(x,{'up' , 'down' , 'non'})

x = [percentage.vHPC.pyr.baseline.up percentage.vHPC.pyr.baseline.down percentage.vHPC.pyr.baseline.non];
subplot(2,2,2),pie(x,{'up' , 'down' , 'non'})

x = [percentage.dHPC.int.baseline.up percentage.dHPC.int.baseline.down percentage.dHPC.int.baseline.non];
subplot(2,2,3),pie(x,{'up' , 'down' , 'non'})

x = [percentage.vHPC.int.baseline.up percentage.vHPC.int.baseline.down percentage.vHPC.int.baseline.non];
subplot(2,2,4),pie(x,{'up' , 'down' , 'non'})

%Reward
figure,
x = [percentage.dHPC.pyr.reward.up percentage.dHPC.pyr.reward.down percentage.dHPC.pyr.reward.non];
subplot(2,2,1),pie(x,{'up' , 'down' , 'non'})

x = [percentage.vHPC.pyr.reward.up percentage.vHPC.pyr.reward.down percentage.vHPC.pyr.reward.non];
subplot(2,2,2),pie(x,{'up' , 'down' , 'non'})

x = [percentage.dHPC.int.reward.up percentage.dHPC.int.reward.down percentage.dHPC.int.reward.non];
subplot(2,2,3),pie(x,{'up' , 'down' , 'non'})

x = [percentage.vHPC.int.reward.up percentage.vHPC.int.reward.down percentage.vHPC.int.reward.non];
subplot(2,2,4),pie(x,{'up' , 'down' , 'non'})

% Aversive
figure,
x = [percentage.dHPC.pyr.aversive.up percentage.dHPC.pyr.aversive.down percentage.dHPC.pyr.aversive.non];
subplot(2,2,1),pie(x,{'up' , 'down' , 'non'})

x = [percentage.vHPC.pyr.aversive.up percentage.vHPC.pyr.aversive.down percentage.vHPC.pyr.aversive.non];
subplot(2,2,2),pie(x,{'up' , 'down' , 'non'})

x = [percentage.dHPC.int.aversive.up percentage.dHPC.int.aversive.down percentage.dHPC.int.aversive.non];
subplot(2,2,3),pie(x,{'up' , 'down' , 'non'})

x = [percentage.vHPC.int.aversive.up percentage.vHPC.int.aversive.down percentage.vHPC.int.aversive.non];
subplot(2,2,4),pie(x,{'up' , 'down' , 'non'})