function [cross time] = cross_between_members_assemblies(path)
% This function a CCG between member-member and member-non_member.
% It iterates acrsoss the path and perform this operation in each session.
% The data is restricted to the active periods (>7cm/sec, at least 2 sec)
% that we used to detect the assemblies.
%
% INPUTS
% path: cell, inside each cell you should put the path of each session you
%       want to analyze
%
% OUTPUT
% cross_members: structure, it contains the CCG for aversive and reward
%                joint assemblies.
%                cross_members.aversive = [];  cross_members.reward = [];
%
% cross_nonmembers: structure, same as the output described above.
%
% time: column vector, it contains the time axis for the CCG plot.
%
% External functions: FMA toolbox
% Morci Juan Facundo 03/2024
%
% Reference: Oberto et al, 2022: 'Distributed cell assemblies spanning 
%            prefrontal cortex and striatum'. Current Biology.
%            doi: https://doi.org/10.1016/j.cub.2021.10.007

% Parameters to define periods of activity
minimal_speed = 7;
minimal_speed_time = 2;

% for CCG
binSize = 0.003; % bisize in seconds
smooth = 2; % smooth factor

% Initialization of structures that wil contain the outputs
cross.members.aversive = [];  cross.members.reward = [];
cross.nonmembers.aversive = [];  cross.nonmembers.reward = [];
time = [];

for tt = 1:length(path)
    %List of folders from the path
    files = dir(path{tt});
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
    clear files dirFlags
    
    for t = 1 : length(subFolders)-2
        disp(' ')
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
        
        %% Awake
        disp('Uploading digital imputs')
        % Load digitalin.mat
        load('digitalin.mat')
        
        % Behavioral calculations
        disp('Uploading DLC outputs')
        camara = ((camara(:,2)-camara(:,1))/2)+camara(:,1);
        % periods of movment during eacj condition
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
        
        % Generation of no-movements periods
        % Reward
        start = behavior.speed.reward(1,1);   stop = behavior.speed.reward(end,1);
        movement.reward = InvertIntervals(behavior.quiet.reward , start , stop); %keep only those higher than criteria
%         movement.reward(movement.reward(:,2) - movement.reward(:,1) <1,:)=[]; %eliminate 1sec segments
        clear tmp start stop
        % Aversive
        start = behavior.speed.aversive(1,1);   stop = behavior.speed.aversive(end,1);
        movement.aversive = InvertIntervals(behavior.quiet.aversive , start , stop);%keep only those higher than criteria
%         movement.aversive(movement.aversive(:,2) - movement.aversive(:,1) <1,:)=[];
        clear tmp start stop
        
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
        clusters.int.dHPC = [];
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
            else
                clusters.int.dHPC = [clusters.int.dHPC ; cluster];
            end
            clear tmp cluster Cellulartype fr1 fr2 fr3 fr4 fr5 r a
        end
        
        numberV = 0;
        clusters.vHPC = [];
        clusters.int.vHPC = [];
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
            else
                clusters.int.vHPC = [clusters.int.vHPC ; cluster];
            end
            clear tmp cluster Cellulartype fr1 fr2 fr3 fr4 fr5 r a
        end
        clear freq limits
        clear camara shock rightvalve leftvalve
        clear ejeX ejeY dX dY dX_int dY_int
        
        if or(numberD >2 , numberV > 2)
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
            
            
            %% CCG
            % Joint assemblies
            % Aversive assemblies
            if sum(cond.both.aversive)>0
                T = Thresholded.aversive.all(:,cond.both.aversive);
                C = [clusters.dHPC ; clusters.vHPC];
                
                % Between members
                for i = 1 : size(T,2)
                    for ii = 1 :size(clusters.dHPC,1)
                        clustD = clusters.dHPC(ii);
                        if T(ii,i)
                            x = spks(spks(:,1)==clustD,2);
                            x = Restrict(x,movement.aversive);
                            for iii = 1:size(clusters.vHPC,1)
                                clustV = clusters.vHPC(iii);
                                if T(size(clusters.dHPC,1)+iii,i)
                                    y = spks(spks(:,1)==clustV,2);
                                    y = Restrict(y,movement.aversive);
                                    if and(length(x) >= 5 , length(y) >= 5)
                                        [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                                        [ccg,time] = CCG(s,ids,'binSize',binSize,'duration',0.6,'smooth',smooth,'mode','ccg');
%                                         if not(sum(isnan(ccg(:,1,2)))>0)
                                            cross.members.aversive = [cross.members.aversive , zscore(ccg(:,1,2))];
%                                         end
                                        clear y s ids groups ccg
                                    end
                                end
                                clear clustV
                            end
                        end
                        clear clustD
                    end
                end
                
                T = sum(T,2)>0; % to select only SU that were not members in any assemblie
                % Between nonmember-nonmember
                for i = 1 : size(T,2)
                    for ii = 1 :size(clusters.dHPC,1)
                        clustD = clusters.dHPC(ii);
                        if not(T(ii,i))
                            x = spks(spks(:,1)==clustD,2);
                            x = Restrict(x,movement.aversive);
                            for iii = 1:size(clusters.vHPC,1)
                                clustV = clusters.vHPC(iii);
                                if not(T(size(clusters.dHPC,1)+iii,i))
                                    y = spks(spks(:,1)==clustV,2);
                                    y = Restrict(y,movement.aversive);
                                    if and(length(x) >= 5 , length(y) >= 5)
                                        [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                                        [ccg,time] = CCG(s,ids,'binSize',binSize,'duration',0.6,'smooth',smooth,'mode','ccg');
%                                         if not(sum(isnan(ccg(:,1,2)))>0)
                                            cross.nonmembers.aversive = [cross.nonmembers.aversive , zscore(ccg(:,1,2))];
%                                         end
                                        clear y s ids groups ccg
                                    end
                                end
                                clear clustV
                            end
                        end
                        clear clustD
                    end
                end
                clear T C
            end
            
            % Reward assemblies
            if sum(cond.both.reward)>0
                T = Thresholded.reward.all(:,cond.both.reward);
                C = [clusters.dHPC ; clusters.vHPC];
                
                % Between members
                for i = 1 : size(T,2)
                    for ii = 1 :size(clusters.dHPC,1)
                        clustD = clusters.dHPC(ii);
                        if T(ii,i)
                            x = spks(spks(:,1)==clustD,2);
                            x = Restrict(x,movement.reward);
                            for iii = 1:size(clusters.vHPC,1)
                                clustV = clusters.vHPC(iii);
                                if T(size(clusters.dHPC,1)+iii,i)
                                    y = spks(spks(:,1)==clustV,2);
                                    y = Restrict(y,movement.reward);
                                    if and(length(x) >= 5 , length(y) >= 5)
                                        [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                                        [ccg,time] = CCG(s,ids,'binSize',binSize,'duration',0.6,'smooth',smooth,'mode','ccg');
%                                         if not(sum(isnan(ccg(:,1,2)))>0)
                                            cross.members.reward = [cross.members.reward , zscore(ccg(:,1,2))];
%                                         end
                                        clear y s ids groups ccg
                                    end
                                end
                                clear clustV
                            end
                        end
                        clear clustD
                    end
                end
                
                T = sum(T,2)>0;% to select only SU that were not members in any assemblie
                % Between nonmembers
                for i = 1 : size(T,2)
                    for ii = 1 :size(clusters.dHPC,1)
                        clustD = clusters.dHPC(ii);
                        if not(T(ii,i))
                            x = spks(spks(:,1)==clustD,2);
                            x = Restrict(x,movement.reward);
                            for iii = 1:size(clusters.vHPC,1)
                                clustV = clusters.vHPC(iii);
                                if not(T(size(clusters.dHPC,1)+iii,i))
                                    y = spks(spks(:,1)==clustV,2);
                                    y = Restrict(y,movement.reward);
                                    if and(length(x) >= 5 , length(y) >= 5)
                                        [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                                        [ccg,time] = CCG(s,ids,'binSize',binSize,'duration',0.6,'smooth',smooth,'mode','ccg');
%                                         if not(sum(isnan(ccg(:,1,2)))>0)
                                            cross.nonmembers.reward = [cross.nonmembers.reward , zscore(ccg(:,1,2))];
%                                         end
                                        clear y s ids groups ccg
                                    end
                                end
                                clear clustV
                            end
                        end
                        clear clustD
                    end
                end    
                clear T C
            end            
            
            clear aversiveTS aversiveTS_run baselineTS bins Cell_type_classification
            clear cellulartype clusters Events group_dHPC group_vHPC i ii iii K x
            clear Kinfo NREM REM numberD numberV p parameter patterns post pre
            clear SpikesD SpikesV spks spks_dHPC spks_vHPC Thresholded WAKE
            clear Kinfo K ii iii ripples ripplesD ripplesV ripple_event Spikes
            clear behaviior rewardTS rewardTS_run aversiveTS aversiveTS_run
        end
    end
    disp(' ')
    disp(' ')
    disp(' ')
end

%% ColorMaps
figure,
subplot(221)
cross.members.aversive(:,sum(isnan(cross.members.aversive))>0) = [];
[i ii] = max(cross.members.aversive);
[i ii] = sort(ii);
imagesc(time,[1:1:size(cross.members.aversive,2)],cross.members.aversive(:,ii)'),caxis([-3 3]), colormap 'jet'

cross.nonmembers.aversive(:,sum(isnan(cross.nonmembers.aversive))>0) = [];
x = randperm(size(cross.nonmembers.aversive,2));
tmp = cross.nonmembers.aversive(:,x);
tmp = tmp(:,1:size(cross.members.aversive,2));
subplot(222)
[i ii] = max(tmp);
[i ii] = sort(ii);
imagesc(time,[1:1:size(tmp,2)],tmp(:,ii)'),caxis([-3 3]), colormap 'jet'

subplot(223)
cross.members.reward(:,sum(isnan(cross.members.reward))>0) = [];
[i ii] = max(cross.members.reward);
[i ii] = sort(ii);
imagesc(time,[1:1:size(cross.members.reward,2)],cross.members.reward(:,ii)'),caxis([-3 3]), colormap 'jet'

cross.nonmembers.reward(:,sum(isnan(cross.nonmembers.reward))>0) = [];
x = randperm(size(cross.nonmembers.reward,2));
tmp1 = cross.nonmembers.reward(:,x);
tmp1 = tmp1(:,1:size(cross.members.reward,2));
subplot(224)
[i ii] = max(tmp1);
[i ii] = sort(ii);
imagesc(time,[1:1:size(tmp1,2)],tmp1(:,ii)'),caxis([-3 3]), colormap 'jet'

%% Average
figure
subplot(121)
plot(time,nanmean(cross.members.aversive,2),'r'), hold on
ciplot(nanmean(cross.members.aversive,2) - nansem(cross.members.aversive')' , nanmean(cross.members.aversive,2) + nansem(cross.members.aversive')' , time , 'r'), alpha 0.5
plot(time,nanmean(tmp,2)','k'),ylim([-0.3 0.5]) , xlim([-0.3 0.3])
ciplot(nanmean(tmp,2) - nansem(tmp')' , nanmean(tmp,2) + nansem(tmp')' , time , 'k'), alpha 0.5
xline(-0.025,'--')
xline(0.025,'--')


subplot(122)
plot(time,nanmean(cross.members.reward,2),'b'), hold on
ciplot(nanmean(cross.members.reward,2) - nansem(cross.members.reward')' , nanmean(cross.members.reward,2) + nansem(cross.members.reward')' , time , 'b'), alpha 0.5
plot(time,nanmean(tmp1,2)','k'),ylim([-0.3 0.5]) , xlim([-0.3 0.3])
ciplot(nanmean(tmp1,2) - nansem(tmp1')' , nanmean(tmp1,2) + nansem(tmp1')' , time , 'k'), alpha 0.5
xline(-0.025,'--')
xline(0.025,'--')
end