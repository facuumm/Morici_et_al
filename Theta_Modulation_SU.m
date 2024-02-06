
clear
clc
close all

%% Parameters
path = {'E:\Rat126\Ephys\in_Pyr';'E:\Rat103\usable';'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat132\recordings\in_pyr';'E:\Rat165\in_pyr\'};%List of folders from the path
% for SU
criteria_fr = 0; %criteria to include or not a SU into the analysis
pval = 0.001;

% for speed selection
minimal_speed = 7; % minimal speed to detect quite periods
minimal_speed_time = 2; % minimal time to detect quite periods

theta_modulated.dHPC.aversive = [];
theta_modulated.dHPC.reward = [];
theta_modulated.vHPC.aversive = [];
theta_modulated.vHPC.reward = [];

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
        
        %load LFP
        if isfile('lfp.mat')
            load('lfp.mat')
        elseif isfile('lfp1.mat')
            load('lfp1.mat')
        end
        
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
        
        %% Awake
        disp('Uploading digital imputs')
        % Load digitalin.mat
        load('digitalin.mat')
                
        % periods of movment during eacj condition
        camara = ((camara(:,2)-camara(:,1))/2)+camara(:,1);
        if rewardTS_run(1) < aversiveTS_run(1)
            load('laps1.mat','posx','posy');
            [camaraR,~] = find((camara(:,1)-rewardTS_run(1)/1000)>0,1,'first'); %TimeStamp of the begining of aversive
            pos = [camara(camaraR : camaraR+length(posx)-1),posx,posy];
            %interpolation of dropped frames
            ejeX = pos(~isnan(pos(:,2)),1); dX = pos(~isnan(pos(:,2)),2); dX_int = interp1(ejeX , dX , pos(:,1));
            ejeY = pos(~isnan(pos(:,3)),1); dY = pos(~isnan(pos(:,3)),3); dY_int = interp1(ejeY , dY , pos(:,1));
            pos(:,2) =dX_int; pos(:,3) =dY_int;%saving corrected pos
            behavior.pos.reward = [pos];
            behavior.speed.reward = LinearVelocity(pos,0);
            behavior.speed.reward(:,2) = smoothdata(behavior.speed.reward(:,2),'gaussian',1,'SamplePoints',behavior.speed.reward(:,1));
            behavior.quiet.reward = QuietPeriods( behavior.speed.reward , minimal_speed , minimal_speed_time);
            clear pos camaraR posx posy
            load('laps2.mat','posx','posy');
            [camaraA,~] = find((camara(:,1)-aversiveTS_run(1)/1000)>0,1,'first'); %TimeStamp of the begining of aversive
            pos = [camara(camaraA : camaraA+length(posx)-1),posx,posy];
            %interpolation of dropped frames
            ejeX = pos(~isnan(pos(:,2)),1); dX = pos(~isnan(pos(:,2)),2); dX_int = interp1(ejeX , dX , pos(:,1));
            ejeY = pos(~isnan(pos(:,3)),1); dY = pos(~isnan(pos(:,3)),3); dY_int = interp1(ejeY , dY , pos(:,1));
            pos(:,2) =dX_int; pos(:,3) =dY_int;%saving corrected pos
            behavior.pos.aversive = [pos];
            behavior.speed.aversive = LinearVelocity(pos,0);
            behavior.speed.aversive(:,2) = smoothdata(behavior.speed.aversive(:,2),'gaussian',1,'SamplePoints',behavior.speed.aversive(:,1));
            behavior.quiet.aversive = QuietPeriods(behavior.speed.aversive , minimal_speed , minimal_speed_time);
            clear pos camaraR
        else
            load('laps2.mat','posx','posy');
            [camaraR,~] = find((camara(:,1)-rewardTS_run(1)/1000)>0,1,'first'); %TimeStamp of the begining of aversive
            %         [camaraR2,~] = find((camara(:,1)-rewardTS_run(2)/1000)<0,1,'last'); %TimeStamp of the ending of aversive
            pos = [camara(camaraR : camaraR+length(posx)-1),posx,posy];
            %interpolation of dropped frames
            ejeX = pos(~isnan(pos(:,2)),1); dX = pos(~isnan(pos(:,2)),2); dX_int = interp1(ejeX , dX , pos(:,1));
            ejeY = pos(~isnan(pos(:,3)),1); dY = pos(~isnan(pos(:,3)),3); dY_int = interp1(ejeY , dY , pos(:,1));
            pos(:,2) =dX_int; pos(:,3) =dY_int;%saving corrected pos
            behavior.pos.reward = [pos];
            behavior.speed.reward = LinearVelocity(pos,0);
            behavior.speed.reward(:,2) = smoothdata(behavior.speed.reward(:,2),'gaussian',1,'SamplePoints',behavior.speed.reward(:,1));
            behavior.quiet.reward = QuietPeriods( behavior.speed.reward , minimal_speed , minimal_speed_time);
            clear pos camaraR posx posy
            load('laps1.mat','posx','posy');
            [camaraA,~] = find((camara(:,1)-aversiveTS_run(1)/1000)>0,1,'first'); %TimeStamp of the begining of aversive
            pos = [camara(camaraA : camaraA+length(posx)-1),posx,posy];
            %interpolation of dropped frames
            ejeX = pos(~isnan(pos(:,2)),1); dX = pos(~isnan(pos(:,2)),2); dX_int = interp1(ejeX , dX , pos(:,1));
            ejeY = pos(~isnan(pos(:,3)),1); dY = pos(~isnan(pos(:,3)),3); dY_int = interp1(ejeY , dY , pos(:,1));
            pos(:,2) =dX_int; pos(:,3) =dY_int; %saving corrected pos
            behavior.pos.aversive = [pos];
            behavior.speed.aversive = LinearVelocity(pos,0);
            behavior.speed.aversive(:,2) = smoothdata(behavior.speed.aversive(:,2),'gaussian',1,'SamplePoints',behavior.speed.aversive(:,1));
            behavior.quiet.aversive = QuietPeriods(behavior.speed.aversive , minimal_speed , minimal_speed_time);
            clear pos camaraA posx posy
        end
        
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
%             cellulartype = [K(:,1) , K(:,4)];
%         elseif criteria_type == 1 % int
%             cellulartype = [K(:,1) , not(K(:,4))];
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
        if not(isempty(spks_dHPC))
            lfp.raw.dHPC = dHPC;
            lfp.filtered.dHPC = FilterLFP(lfp.raw.dHPC,'passband',[6 10]);% filtering LFP
            disp('Rayleight for dHPC SUs')
            iterator = ismember(K(:,1),clusters.dHPC); %keep logical to differentiate pyr from rest
            iterator  = K(iterator,4);
            
            Periods = SubtractIntervals(aversiveTS_run./1000,behavior.quiet.aversive);
            Periods = Restrict(Periods,[behavior.pos.aversive(1,1) behavior.pos.aversive(end,1)]);
            [p,theta,rbar,delta] = ThetaModulation(spks_dHPC,clusters.dHPC,lfp.filtered.dHPC,Periods);
            theta_modulated.dHPC.aversive = [theta_modulated.dHPC.aversive ; clusters.dHPC iterator p' theta' rbar' delta'];            
            temp.dHPC.aversive = [clusters.dHPC iterator p' theta' rbar' delta'];
            clear Periods p theta rbar delta
            
            Periods = SubtractIntervals(rewardTS_run./1000,behavior.quiet.reward);
            Periods = Restrict(Periods,[behavior.pos.reward(1,1) behavior.pos.reward(end,1)]);
            [p,theta,rbar,delta] = ThetaModulation(spks_dHPC,clusters.dHPC,lfp.filtered.dHPC,Periods);
            theta_modulated.dHPC.reward = [theta_modulated.dHPC.reward ; clusters.dHPC iterator p' theta' rbar' delta'];  
            temp.dHPC.reward = [clusters.dHPC iterator p' theta' rbar' delta'];
            clear Periods p theta rbar delta iterator 
        end
        
        if not(isempty(spks_vHPC))
            PeriodsA = SubtractIntervals(aversiveTS_run./1000,behavior.quiet.aversive);
            PeriodsR = SubtractIntervals(rewardTS_run./1000,behavior.quiet.reward);
            P = sort([PeriodsA ; PeriodsR]); clear PeriodsA PeriodsR
            [spectrum1,f,s] = MTSpectrum(Restrict(vHPC1,P));
            [spectrum2,f,s] = MTSpectrum(Restrict(vHPC2,P));
            [spectrum3,f,s] = MTSpectrum(Restrict(vHPC3,P));
            [spectrum4,f,s] = MTSpectrum(Restrict(vHPC4,P));
            [spectrum5,f,s] = MTSpectrum(Restrict(vHPC5,P));
            
            [m mm] = min(abs(f-6));
            [m mmm] = min(abs(f-10));
            
            [m s] = max([mean(spectrum1(mm : mmm)) mean(spectrum2(mm : mmm)) mean(spectrum3(mm : mmm)) mean(spectrum4(mm : mmm)) mean(spectrum5(mm : mmm))]);
            
            if s == 1
                lfp.raw.vHPC = vHPC1;
            elseif s == 2
                lfp.raw.vHPC = vHPC2;
            elseif s == 3
                lfp.raw.vHPC = vHPC3;
            elseif s == 4
                lfp.raw.vHPC = vHPC4;
            elseif s == 5
                lfp.raw.vHPC = vHPC5;
            end
            clear dHPC vHPC1 vHPC2 vHPC3 vHPC4 vHPC5
            clear spectrum1 spectrum2 spectrum3 spectrum4 spectrum5 s m mm mmm
            
            lfp.filtered.vHPC = FilterLFP(lfp.raw.vHPC,'passband',[6 10]);% filtering LFP
            disp('Rayleight for vHPC SUs')
            iterator = ismember(K(:,1),clusters.vHPC); %keep logical to differentiate pyr from rest
            iterator  = K(iterator,4);
            
            Periods = SubtractIntervals(aversiveTS_run./1000,behavior.quiet.aversive);
            Periods = Restrict(Periods,[behavior.pos.aversive(1,1) behavior.pos.aversive(end,1)]);
            [p,theta,rbar,delta] = ThetaModulation(spks_vHPC,clusters.vHPC,lfp.filtered.vHPC,Periods);
            theta_modulated.vHPC.aversive = [theta_modulated.vHPC.aversive ; clusters.vHPC iterator p' theta' rbar' delta'];
            temp.vHPC.aversive = [clusters.vHPC iterator p' theta' rbar' delta'];
            clear Periods p theta rbar delta
            
            Periods = SubtractIntervals(rewardTS_run./1000,behavior.quiet.reward);
            Periods = Restrict(Periods,[behavior.pos.reward(1,1) behavior.pos.reward(end,1)]);
            [p,theta,rbar,delta] = ThetaModulation(spks_vHPC,clusters.vHPC,lfp.filtered.vHPC,Periods);
            theta_modulated.vHPC.reward = [theta_modulated.vHPC.reward ; clusters.vHPC iterator p' theta' rbar' delta'];
            temp.vHPC.reward = [clusters.vHPC iterator p' theta' rbar' delta'];
            clear Periods p theta rbar delta iterator
        end
        
        save([cd,'\theta_modulated_SU.mat'],'temp')
        clear temp
        clear aversiveTS aversiveTS_run rewardTS rewardTS_run NREM REM WAKE
        clear lfp ripplesD ripplesV spks spks_dHPC spks_vHPC leftvalve rightvalve
        clear ripple_bursts behavior baselineTS camara Cell_type_classification
        clear clusters coordinated coordinatedV cooridnated_event cooridnated_eventDV cooridnated_eventVD
        clear segments shock group_dHPC group_vHPC K Kinfo coordinatedV_refined
        
    end
    disp(['-- Finishing analysis from rat #',num2str(tt) , ' --'])
    disp('  ')
end


dHPC.pyr = theta_modulated.dHPC.aversive(:,2)==1;
dHPC.int = theta_modulated.dHPC.aversive(:,2)==0;
vHPC.pyr = theta_modulated.vHPC.aversive(:,2)==1;
vHPC.int = theta_modulated.vHPC.aversive(:,2)==0;

%% Percentage calculation dHPC pyr
% aversive
n = length(theta_modulated.dHPC.aversive(dHPC.pyr));
i = sum(theta_modulated.dHPC.aversive(dHPC.pyr,3)<0.05);
percentage.dHPC.pyr.aversive.modulated = (i/n)*100;

i = sum(isnan(theta_modulated.dHPC.aversive(dHPC.pyr,3)));
percentage.dHPC.pyr.aversive.nan = (i/n)*100;

%Reward
n = length(theta_modulated.dHPC.reward(dHPC.pyr));
i = sum(theta_modulated.dHPC.reward(dHPC.pyr,3)<0.05);
percentage.dHPC.pyr.reward.modulated = (i/n)*100;

i = sum(isnan(theta_modulated.dHPC.reward(dHPC.pyr,3)));
percentage.dHPC.pyr.reward.nan = (i/n)*100;

%% Percentage calculation dHPC int
n = length(theta_modulated.dHPC.aversive(dHPC.int));
i = sum(theta_modulated.dHPC.aversive(dHPC.int,3)<0.05);
percentage.dHPC.int.aversive.modulated = (i/n)*100;

i = sum(isnan(theta_modulated.dHPC.aversive(dHPC.int,3)));
percentage.dHPC.int.aversive.nan = (i/n)*100;

n = length(theta_modulated.dHPC.reward(dHPC.int));
i = sum(theta_modulated.dHPC.reward(dHPC.int,3)<0.05);
percentage.dHPC.int.reward.modulated = (i/n)*100;

i = sum(isnan(theta_modulated.dHPC.reward(dHPC.int,3)));
percentage.dHPC.int.reward.nan = (i/n)*100;

%% Percentage calculation vHPC pyr
n = length(theta_modulated.vHPC.aversive(vHPC.pyr));
i = sum(theta_modulated.vHPC.aversive(vHPC.pyr,3)<0.05);
percentage.vHPC.pyr.aversive.modulated = (i/n)*100;

i = sum(isnan(theta_modulated.vHPC.aversive(vHPC.pyr,3)));
percentage.vHPC.pyr.aversive.nan = (i/n)*100;


n = length(theta_modulated.vHPC.reward(vHPC.pyr));
i = sum(theta_modulated.vHPC.reward(vHPC.pyr,3)<0.05);
percentage.vHPC.pyr.reward.modulated = (i/n)*100;

i = sum(isnan(theta_modulated.vHPC.reward(vHPC.pyr,3)));
percentage.vHPC.pyr.reward.nan = (i/n)*100;

%% Percentage calculation vHPC int
n = length(theta_modulated.vHPC.aversive(vHPC.int));
i = sum(theta_modulated.vHPC.aversive(vHPC.int,3)<0.05);
percentage.vHPC.int.aversive.modulated = (i/n)*100;

i = sum(isnan(theta_modulated.vHPC.aversive(vHPC.int,3)));
percentage.vHPC.int.aversive.nan = (i/n)*100;

n = length(theta_modulated.vHPC.reward(vHPC.int));
i = sum(theta_modulated.vHPC.reward(vHPC.int,3)<0.05);
percentage.vHPC.int.reward.modulated = (i/n)*100;

i = sum(isnan(theta_modulated.vHPC.reward(vHPC.int,3)));
percentage.vHPC.int.reward.nan = (i/n)*100;

%% Plots
%All
figure,
x = [percentage.dHPC.pyr.aversive.modulated percentage.dHPC.pyr.aversive.nan 100-percentage.dHPC.pyr.aversive.modulated-percentage.dHPC.pyr.aversive.nan];
subplot(2,2,1),pie(x,{'modulated' , 'non-specified' , 'non-modulated'})

x = [percentage.dHPC.int.aversive.modulated percentage.dHPC.int.aversive.nan 100-percentage.dHPC.int.aversive.modulated-percentage.dHPC.int.aversive.nan];
subplot(2,2,2),pie(x,{'modulated' , 'non-specified' , 'non-modulated'})

x = [percentage.dHPC.pyr.reward.modulated percentage.dHPC.pyr.reward.nan 100-percentage.dHPC.pyr.reward.modulated-percentage.dHPC.pyr.reward.nan];
subplot(2,2,3),pie(x,{'modulated' , 'non-specified' , 'non-modulated'})

x = [percentage.dHPC.int.reward.modulated percentage.dHPC.int.reward.nan 100-percentage.dHPC.int.reward.modulated-percentage.dHPC.int.reward.nan];
subplot(2,2,4),pie(x,{'modulated' , 'non-specified' , 'non-modulated'})


figure,
x = [percentage.vHPC.pyr.aversive.modulated percentage.vHPC.pyr.aversive.nan 100-percentage.vHPC.pyr.aversive.modulated-percentage.vHPC.pyr.aversive.nan];
subplot(2,2,1),pie(x,{'modulated' , 'non-specified' , 'non-modulated'})

x = [percentage.vHPC.int.aversive.modulated percentage.vHPC.int.aversive.nan 100-percentage.vHPC.int.aversive.modulated-percentage.vHPC.int.aversive.nan];
subplot(2,2,2),pie(x,{'modulated' , 'non-specified' , 'non-modulated'})

x = [percentage.vHPC.pyr.reward.modulated percentage.vHPC.pyr.reward.nan 100-percentage.vHPC.pyr.reward.modulated-percentage.vHPC.pyr.reward.nan];
subplot(2,2,3),pie(x,{'modulated' , 'non-specified' , 'non-modulated'})

x = [percentage.vHPC.int.reward.modulated percentage.vHPC.int.reward.nan 100-percentage.vHPC.int.reward.modulated-percentage.vHPC.int.reward.nan];
subplot(2,2,4),pie(x,{'modulated' , 'non-specified' , 'non-modulated'})

