clear
clc
close all

%% Parameters
path = {'E:\Rat126\Ephys\in_Pyr';'E:\Rat103\usable';'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat132\recordings\in_pyr';'E:\Rat165\in_pyr\'};%List of folders from the path
% for SU
criteria_fr = 0; %criteria to include or not a SU into the analysis
pval = 0.05;
binSize = 1; % time in sec to perform correlations

% for speed selection
minimal_speed = 7; % minimal speed to detect quite periods
minimal_speed_time = 2; % minimal time to detect quite periods

Correlations.dHPC.pyr.Aversive = []; Correlations.dHPC.pyr.Reward = [];
Correlations.vHPC.pyr.Aversive = []; Correlations.vHPC.pyr.Reward = [];
Correlations.dHPC.int.Aversive = []; Correlations.dHPC.int.Reward = [];
Correlations.vHPC.int.Aversive = []; Correlations.vHPC.int.Reward = []; 

R.dHPC.pyr.Aversive = []; R.dHPC.pyr.Reward = [];
R.vHPC.pyr.Aversive = []; R.vHPC.pyr.Reward = [];
R.dHPC.int.Aversive = []; R.dHPC.int.Reward = [];
R.vHPC.int.Aversive = []; R.vHPC.int.Reward = []; 

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
        
        %         %load LFP
        %         if isfile('lfp.mat')
        %             load('lfp.mat')
        %         elseif isfile('lfp1.mat')
        %             load('lfp1.mat')
        %         end
        %
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
%                     A = y;
                end
                % Rewarded sleep session TS detection
            elseif strcmp(segments.Var2{y},'reward')
                if strcmp(segments.Var3{y},'End')
                    rewardTS(1,1) = segments.Var1(y+1);
                    rewardTS(1,2) = segments.Var1(y+2);
                    rewardTS_run(1,1) = segments.Var1(y-1);
                    rewardTS_run(1,2) = segments.Var1(y);
%                     R = y;
                end
            end
        end
%         clear y A R
        
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
        cellulartype = [K(:,1) , K(:,4)];
        %         elseif criteria_type == 1 % int
        %             cellulartype = [K(:,1) , not(K(:,4))];
        %         elseif criteria_type == 2 % all
        %             cellulartype = [K(:,1) , ones(length(K),1)];
        %         end
        
        %% Counting the Number f SU
        numberD = 0;
        clusters.all = [];
        clusters.dHPC.pyr = [];        clusters.dHPC.int = [];
        
        for ii=1:size(group_dHPC,1)
            cluster = group_dHPC(ii,1);
            Cellulartype = logical(cellulartype(cellulartype(:,1) == cluster,2));
            if Cellulartype
                a = length(Restrict(spks(spks(:,1)==cluster,2),aversiveTS_run./1000)) / ((aversiveTS_run(2)-aversiveTS_run(1))/1000);
                r = length(Restrict(spks(spks(:,1)==cluster,2),rewardTS_run./1000)) / ((rewardTS_run(2)-rewardTS_run(1))/1000);
                if or(a > criteria_fr ,  r > criteria_fr)
                    numberD = numberD+1;
                    clusters.all = [clusters.all ; cluster];
                    clusters.dHPC.pyr = [clusters.dHPC.pyr ; cluster];
                end
            else
                clusters.all = [clusters.all ; cluster];
                clusters.dHPC.int = [clusters.dHPC.int ; cluster];
            end
            clear tmp cluster Cellulartype fr1 fr2 fr3 fr4 fr5 r a
        end
        
        numberV = 0;
        clusters.vHPC.pyr = [];        clusters.vHPC.int = [];
        for ii=1:size(group_vHPC,1)
            cluster = group_vHPC(ii,1);
            Cellulartype = logical(cellulartype(cellulartype(:,1) == cluster,2));
            if Cellulartype
                a = length(Restrict(spks(spks(:,1)==cluster,2),aversiveTS_run./1000)) / ((aversiveTS_run(2)-aversiveTS_run(1))/1000);
                r = length(Restrict(spks(spks(:,1)==cluster,2),rewardTS_run./1000)) / ((rewardTS_run(2)-rewardTS_run(1))/1000);
                if or(a > criteria_fr ,  r > criteria_fr)
                    numberV = numberV+1;
                    clusters.all = [clusters.all ; cluster];
                    clusters.vHPC.pyr = [clusters.vHPC.pyr ; cluster];
                end
            else
                clusters.all = [clusters.all ; cluster];
                clusters.vHPC.int = [clusters.vHPC.int ; cluster];
            end
            clear tmp cluster Cellulartype fr1 fr2 fr3 fr4 fr5 r a
        end
        clear freq limits
        clear ejeX ejeY dX dY dX_int dY_int
        
        
        %% Metrics calculation
        % Speed Aversive
        Periods =[behavior.speed.aversive(1,1) behavior.speed.aversive(end,1)];
        TimeScale = [Periods(1) : binSize : Periods(2)];
        [N_bin,edges, bins] = histcounts(behavior.speed.aversive(:,1),TimeScale);
        MeanSpeed.Aversive = [];
        for i = 1 : size(TimeScale,2)
            iterator = bins ==i;
            MeanSpeed.Aversive = [MeanSpeed.Aversive ; TimeScale(i) nanmean(behavior.speed.aversive(iterator,2))];
            clear iterator
        end
        MeanSpeed.Aversive(:,2) = (MeanSpeed.Aversive(:,2)-nanmean(MeanSpeed.Aversive(:,2)))/nanstd(MeanSpeed.Aversive(:,2));
        clear i clear N_bin edges bins
        
        % Speed Reward
        Periods =[behavior.speed.reward(1,1) behavior.speed.reward(end,1)];
        TimeScale = [Periods(1) : binSize : Periods(2)];
        [N_bin,edges, bins] = histcounts(behavior.speed.reward(:,1),TimeScale);
        MeanSpeed.Reward = [];
        for i = 1 : size(TimeScale,2)
            iterator = bins ==i;
            MeanSpeed.Reward = [MeanSpeed.Reward ; TimeScale(i) nanmean(behavior.speed.reward(iterator,2))];
            clear iterator
        end
        MeanSpeed.Reward(:,2) = (MeanSpeed.Reward(:,2)-nanmean(MeanSpeed.Reward(:,2)))/nanstd(MeanSpeed.Reward(:,2));
        clear i clear N_bin edges bins
        
        %% SU
        % dHPC Pyr
        if not(isempty(clusters.dHPC.pyr))
            % Aversive
            Periods =[behavior.speed.aversive(1,1) behavior.speed.aversive(end,1)];
            limits = [Periods(1) Periods(2)];
            [Spikes.dHPC.Pyr.aversive , bins , Clusters] = spike_train_construction(spks, clusters.dHPC.pyr, cellulartype, binSize, limits, [] , true, false);
            clear limits bins Clusters Periods
            % Concatenation
            tmp = [];
            for i = 1 : size(Spikes.dHPC.Pyr.aversive,2)
                tmp = [tmp ; MeanSpeed.Aversive(:,2) , Spikes.dHPC.Pyr.aversive(:,i)];
                r = fitlm(MeanSpeed.Aversive(:,2) , Spikes.dHPC.Pyr.aversive(:,i));
                R.dHPC.pyr.Aversive = [R.dHPC.pyr.Aversive ; r.Coefficients.Estimate(2) r.Rsquared.Ordinary r.Coefficients.pValue(2)<pval]; clear r                
            end
            Spikes.dHPC.Pyr.aversive = tmp; clear tmp i
            
            % Reward
            Periods =[behavior.speed.reward(1,1) behavior.speed.reward(end,1)];
            limits = [Periods(1) Periods(2)];
            [Spikes.dHPC.Pyr.reward , bins , Clusters] = spike_train_construction(spks, clusters.dHPC.pyr, cellulartype, binSize, limits, [] , true, false);
            clear limits bins Clusters Periods
            % Concatenation
            tmp = [];
            for i = 1 : size(Spikes.dHPC.Pyr.reward,2)
                tmp = [tmp ; MeanSpeed.Reward(:,2) , Spikes.dHPC.Pyr.reward(:,i)];
                r = fitlm(MeanSpeed.Reward(:,2) , Spikes.dHPC.Pyr.reward(:,i));
                R.dHPC.pyr.Reward = [R.dHPC.pyr.Reward ; r.Coefficients.Estimate(2) r.Rsquared.Ordinary r.Coefficients.pValue(2)<pval]; clear r               
            end
            Spikes.dHPC.Pyr.reward = tmp; clear tmp i
        else
            Spikes.dHPC.Pyr.aversive = [];
            Spikes.dHPC.Pyr.reward = [];
        end
        
        % dHPC Int
        if not(isempty(clusters.dHPC.int))
            % Aversive
            Periods =[behavior.speed.aversive(1,1) behavior.speed.aversive(end,1)];
            limits = [Periods(1) Periods(2)];
            [Spikes.dHPC.int.aversive , bins , Clusters] = spike_train_construction(spks, clusters.dHPC.int, [cellulartype(:,1) not(cellulartype(:,2))], binSize, limits, [] , true, false);
            clear limits bins Clusters
            % Concatenation
            tmp = [];
            for i = 1 : size(Spikes.dHPC.int.aversive,2)
                tmp = [tmp ; MeanSpeed.Aversive(:,2) , Spikes.dHPC.int.aversive(:,i)];
                r = fitlm(MeanSpeed.Aversive(:,2) , Spikes.dHPC.int.aversive(:,i));
                R.dHPC.int.Aversive = [R.dHPC.int.Aversive ; r.Coefficients.Estimate(2) r.Rsquared.Ordinary r.Coefficients.pValue(2)<pval]; clear r  
            end
            Spikes.dHPC.int.aversive = tmp; clear tmp i
            
            % Reward
            Periods =[behavior.speed.reward(1,1) behavior.speed.reward(end,1)];
            limits = [Periods(1) Periods(2)];
            [Spikes.dHPC.int.reward , bins , Clusters] = spike_train_construction(spks, clusters.dHPC.int, [cellulartype(:,1) not(cellulartype(:,2))], binSize, limits, [] , true, false);
            clear limits bins Clusters    Periods
            % Concatenation
            tmp = [];
            for i = 1 : size(Spikes.dHPC.int.reward,2)
                tmp = [tmp ; MeanSpeed.Reward(:,2) , Spikes.dHPC.int.reward(:,i)];
                r = fitlm(MeanSpeed.Reward(:,2) , Spikes.dHPC.int.reward(:,i));
                R.dHPC.int.Reward = [R.dHPC.int.Reward ; r.Coefficients.Estimate(2) r.Rsquared.Ordinary r.Coefficients.pValue(2)<pval]; clear r                    
            end
            Spikes.dHPC.int.reward = tmp; clear tmp i
        else
            Spikes.dHPC.int.aversive = [];
            Spikes.dHPC.int.reward = [];
        end
        
        % vHPC Pyr
        if not(isempty(clusters.vHPC.pyr))
            % Aversive
            Periods =[behavior.speed.aversive(1,1) behavior.speed.aversive(end,1)];
            limits = [Periods(1) Periods(2)];
            [Spikes.vHPC.Pyr.aversive , bins , Clusters] = spike_train_construction(spks, clusters.vHPC.pyr, cellulartype, binSize, limits, [] , true, false);
            clear limits bins Clusters Periods
            % Concatenation
            tmp = [];
            for i = 1 : size(Spikes.vHPC.Pyr.aversive,2)
                tmp = [tmp ; MeanSpeed.Aversive(:,2) , Spikes.vHPC.Pyr.aversive(:,i)];
                r = fitlm(MeanSpeed.Aversive(:,2) , Spikes.vHPC.Pyr.aversive(:,i));
                R.vHPC.pyr.Aversive = [R.vHPC.pyr.Aversive ; r.Coefficients.Estimate(2) r.Rsquared.Ordinary r.Coefficients.pValue(2)<pval]; clear r  
            end
            Spikes.vHPC.Pyr.aversive = tmp; clear tmp i
            
            % Reward
            Periods =[behavior.speed.reward(1,1) behavior.speed.reward(end,1)];
            limits = [Periods(1) Periods(2)];
            [Spikes.vHPC.Pyr.reward , bins , Clusters] = spike_train_construction(spks, clusters.vHPC.pyr, cellulartype, binSize, limits, [] , true, false);
            clear limits bins Clusters Periods
            % Concatenation
            tmp = [];
            for i = 1 : size(Spikes.vHPC.Pyr.reward,2)
                tmp = [tmp ; MeanSpeed.Reward(:,2) , Spikes.vHPC.Pyr.reward(:,i)];
                r = fitlm(MeanSpeed.Reward(:,2) , Spikes.vHPC.Pyr.reward(:,i));
                R.vHPC.pyr.Reward = [R.vHPC.pyr.Reward ; r.Coefficients.Estimate(2) r.Rsquared.Ordinary r.Coefficients.pValue(2)<pval]; clear r                
            end
            Spikes.vHPC.Pyr.reward = tmp; clear tmp i
        else
            Spikes.vHPC.Pyr.aversive = [];
            Spikes.vHPC.Pyr.reward = [];
        end
        
        
        if not(isempty(clusters.vHPC.int))
            % Aversive
            Periods =[behavior.speed.aversive(1,1) behavior.speed.aversive(end,1)];
            limits = [Periods(1) Periods(2)];
            [Spikes.vHPC.int.aversive , bins , Clusters] = spike_train_construction(spks, clusters.vHPC.int, [cellulartype(:,1) not(cellulartype(:,2))], binSize, limits, [] , true, false);
            clear limits bins Clusters Periods
            % Concatenation
            tmp = [];
            for i = 1 : size(Spikes.vHPC.int.aversive,2)
                tmp = [tmp ; MeanSpeed.Aversive(:,2) , Spikes.vHPC.int.aversive(:,i)];
                r = fitlm(MeanSpeed.Aversive(:,2) , Spikes.vHPC.int.aversive(:,i));
                R.vHPC.int.Aversive = [R.vHPC.int.Aversive ; r.Coefficients.Estimate(2) r.Rsquared.Ordinary r.Coefficients.pValue(2)<pval]; clear r                  
            end
            Spikes.vHPC.int.aversive = tmp; clear tmp i
            
            % Reward
            Periods =[behavior.speed.reward(1,1) behavior.speed.reward(end,1)];
            limits = [Periods(1) Periods(2)];
            [Spikes.vHPC.int.reward , bins , Clusters] = spike_train_construction(spks, clusters.vHPC.int, [cellulartype(:,1) not(cellulartype(:,2))], binSize, limits, [] , true, false);
            clear limits bins Clusters Periods
            % Concatenation
            tmp = [];
            for i = 1 : size(Spikes.vHPC.int.reward,2)
                tmp = [tmp ; MeanSpeed.Reward(:,2) , Spikes.vHPC.int.reward(:,i)];
                r = fitlm(MeanSpeed.Reward(:,2) , Spikes.vHPC.int.reward(:,i));
                R.vHPC.int.Reward = [R.vHPC.int.Reward ; r.Coefficients.Estimate(2) r.Rsquared.Ordinary r.Coefficients.pValue(2)<pval]; clear r                   
            end
            Spikes.vHPC.int.reward = tmp; clear tmp i
        else
            Spikes.vHPC.int.aversive = [];
            Spikes.vHPC.int.reward = [];
        end
            
        
        %% Save data
        Correlations.dHPC.pyr.Aversive = [Correlations.dHPC.pyr.Aversive ; Spikes.dHPC.Pyr.aversive];
        Correlations.dHPC.pyr.Reward = [Correlations.dHPC.pyr.Reward ; Spikes.dHPC.Pyr.reward];
        Correlations.vHPC.pyr.Aversive = [Correlations.vHPC.pyr.Aversive ; Spikes.vHPC.Pyr.aversive];
        Correlations.vHPC.pyr.Reward = [Correlations.vHPC.pyr.Reward ; Spikes.vHPC.Pyr.reward]; 
        
        Correlations.dHPC.int.Aversive = [Correlations.dHPC.int.Aversive ; Spikes.dHPC.int.aversive];
        Correlations.dHPC.int.Reward = [Correlations.dHPC.int.Reward ; Spikes.dHPC.int.reward];
        Correlations.vHPC.int.Aversive = [Correlations.vHPC.int.Aversive ; Spikes.vHPC.int.aversive];
        Correlations.vHPC.int.Reward = [Correlations.vHPC.int.Reward ; Spikes.vHPC.int.reward];           
        
%         save([cd,'\theta_modulated_SU.mat'],'temp')
        clear temp
        clear aversiveTS aversiveTS_run rewardTS rewardTS_run NREM REM WAKE
        clear lfp ripplesD ripplesV spks spks_dHPC spks_vHPC leftvalve rightvalve
        clear ripple_bursts behavior baselineTS camara Cell_type_classification
        clear clusters coordinated coordinatedV cooridnated_event cooridnated_eventDV cooridnated_eventVD
        clear segments shock group_dHPC group_vHPC K Kinfo coordinatedV_refined Spikes numberD numberV TimeScale cellulartype
        
    end
    disp(['-- Finishing analysis from rat #',num2str(tt) , ' --'])
    disp('  ')
end


scatter(Correlations.dHPC.pyr.Aversive(:,1) , Correlations.dHPC.pyr.Aversive(:,2))

%% Plot scatter for correlations R2
figure
subplot(221),
histogram(R.dHPC.pyr.Reward(:,1),[-1:0.05:1],'Normalization','probability'),hold on
histogram(R.dHPC.pyr.Aversive(:,1),[-1:0.05:1],'Normalization','probability'),xlim([-1 1]) , ylim([0 0.35])

subplot(222),
histogram(R.vHPC.pyr.Reward(:,1),[-1:0.05:1],'Normalization','probability'),hold on
histogram(R.vHPC.pyr.Aversive(:,1),[-1:0.05:1],'Normalization','probability'),xlim([-1 1]), ylim([0 0.35])

subplot(223),
histogram(R.dHPC.int.Reward(:,1),[-1:0.05:1],'Normalization','probability'),hold on
histogram(R.dHPC.int.Aversive(:,1),[-1:0.05:1],'Normalization','probability'),xlim([-1 1]), ylim([0 0.35])

subplot(224),
histogram(R.vHPC.int.Reward(:,1),[-1:0.05:1],'Normalization','probability'),hold on
histogram(R.vHPC.int.Aversive(:,1),[-1:0.05:1],'Normalization','probability'),xlim([-1 1]), ylim([0 0.35])


%%

subplot(221)
grps = [ones(length(R.dHPC.pyr.Aversive),1) ; ones(length(R.dHPC.pyr.Reward),1)*2];
x = [R.dHPC.pyr.Aversive(:,2) ; R.dHPC.pyr.Reward(:,2)];
scatter(grps,x,"filled",'jitter','on', 'jitterAmount',0.1), xlim([0 3]),hold on
scatter([1 2],[nanmean(R.dHPC.pyr.Aversive(:,2)) , nanmean(R.dHPC.pyr.Reward(:,2))],'filled'), xlim([0 3]),ylim([0 0.7])
ttest2(R.dHPC.pyr.Aversive(:,2) , R.dHPC.pyr.Reward(:,2))

subplot(222)
grps = [ones(length(R.vHPC.pyr.Aversive),1) ; ones(length(R.vHPC.pyr.Reward),1)*2];
x = [R.vHPC.pyr.Aversive(:,2) ; R.vHPC.pyr.Reward(:,2)];
scatter(grps,x,"filled",'jitter','on', 'jitterAmount',0.1), xlim([0 3]),hold on
scatter([1 2],[nanmean(R.vHPC.pyr.Aversive(:,2)) , nanmean(R.vHPC.pyr.Reward(:,2))],'filled'), xlim([0 3]),ylim([0 0.7])
ttest2(R.vHPC.pyr.Aversive(:,2) , R.vHPC.pyr.Reward(:,2))

subplot(223)
grps = [ones(length(R.dHPC.int.Aversive(:,2)),1) ; ones(length(R.dHPC.int.Reward(:,2)),1)*2];
x = [R.dHPC.int.Aversive(:,2) ; R.dHPC.int.Reward(:,2)];
scatter(grps,x,"filled",'jitter','on', 'jitterAmount',0.1), xlim([0 3]),hold on
scatter([1 2],[nanmean(R.dHPC.int.Aversive(:,2)) , nanmean(R.dHPC.int.Reward(:,2))],'filled'), xlim([0 3]),ylim([0 0.7])

subplot(224)
grps = [ones(length(R.vHPC.int.Aversive),1) ; ones(length(R.vHPC.int.Reward),1)*2];
x = [R.vHPC.int.Aversive(:,2) ; R.vHPC.int.Reward(:,2)];
scatter(grps,x,"filled",'jitter','on', 'jitterAmount',0.1), xlim([0 3]),hold on
scatter([1 2],[nanmean(R.vHPC.int.Aversive(:,2)) , nanmean(R.vHPC.int.Reward(:,2))],'filled'), xlim([0 3]),ylim([0 0.7])



figure

subplot(221)
iterator1 = logical(R.dHPC.pyr.Aversive(:,2));
iterator2 = logical(R.dHPC.pyr.Reward(:,2));

grps = [ones(length(R.dHPC.pyr.Aversive(iterator1)),1) ; ones(length(R.dHPC.pyr.Reward(iterator2)),1)*2];
x = [R.dHPC.pyr.Aversive(iterator1,1) ; R.dHPC.pyr.Reward(iterator2,1)];
scatter(grps,x,"filled",'jitter','on', 'jitterAmount',0.1), xlim([0 3]),hold on
scatter([1 2],[nanmean(R.dHPC.pyr.Aversive(iterator1,1)) , nanmean(R.dHPC.pyr.Reward(iterator2,1))],'filled'), xlim([0 3]),ylim([0 0.7])
ttest2(R.dHPC.pyr.Aversive(iterator1,1) , R.dHPC.pyr.Reward(iterator2,1))

subplot(222)
iterator1 = logical(R.vHPC.pyr.Aversive(:,2));
iterator2 = logical(R.vHPC.pyr.Reward(:,2));

grps = [ones(length(R.vHPC.pyr.Aversive(iterator1)),1) ; ones(length(R.vHPC.pyr.Reward(iterator2)),1)*2];
x = [R.vHPC.pyr.Aversive(iterator1,1) ; R.vHPC.pyr.Reward(iterator2,1)];
scatter(grps,x,"filled",'jitter','on', 'jitterAmount',0.1), xlim([0 3]),hold on
scatter([1 2],[nanmean(R.vHPC.pyr.Aversive(iterator1,1)) , nanmean(R.vHPC.pyr.Reward(iterator2,1))],'filled'), xlim([0 3]),ylim([0 0.7])
ttest2(R.vHPC.pyr.Aversive(iterator1,1) , R.vHPC.pyr.Reward(iterator2,1))

subplot(223)
iterator1 = logical(R.dHPC.int.Aversive(:,2));
iterator2 = logical(R.dHPC.int.Reward(:,2));

grps = [ones(length(R.dHPC.int.Aversive(iterator1,1)),1) ; ones(length(R.dHPC.int.Reward(iterator2,1)),1)*2];
x = [R.dHPC.int.Aversive(iterator1,1) ; R.dHPC.int.Reward(iterator2,1)];
scatter(grps,x,"filled",'jitter','on', 'jitterAmount',0.1), xlim([0 3]),hold on
scatter([1 2],[nanmean(R.dHPC.int.Aversive(iterator1,1)) , nanmean(R.dHPC.int.Reward(iterator2,1))],'filled'), xlim([0 3]),ylim([0 0.7])

subplot(224)
iterator1 = logical(R.vHPC.int.Aversive(:,2));
iterator2 = logical(R.vHPC.int.Reward(:,2));

grps = [ones(length(R.vHPC.int.Aversive(iterator1,1)),1) ; ones(length(R.vHPC.int.Reward(iterator2,1)),1)*2];
x = [R.vHPC.int.Aversive(iterator1,1) ; R.vHPC.int.Reward(iterator2,1)];
scatter(grps,x,"filled",'jitter','on', 'jitterAmount',0.1), xlim([0 3]),hold on
scatter([1 2],[nanmean(R.vHPC.int.Aversive(iterator1,1)) , nanmean(R.vHPC.int.Reward(iterator2,1))],'filled'), xlim([0 3]),ylim([0 0.7])
