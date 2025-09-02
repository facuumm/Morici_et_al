function [Response , time ,  modulated] = AssemblyResponse_to_valve(path)
% This function calculate the average assemblies activity sourrounding the 
% shock delivery.
%
% INPUTS
% path: cell, inside each cell you should put the path of each session you
%       want to analyze
%
% OUTPUT
% Response: structure, it contains the triggered average response to the shocks.
%           Response.dHPC.aversive      Response.dHPC.reward
%           Response.vHPC.aversive      Response.vHPC.reward
%           Response.joint.aversive     Response.joint.reward
%
% time: column vector, it contains the time axis for the plot.
%
% modulated: structure, it contains 1 if the average activity is above the
%            90th percentile of a surrogated distribution constructed by
%            shuffling the shocks timestamps 100 times.
%            Same organization than Response
%
% Response: structure, contains the speed sourrounding the reward delivery
%
% other functions: CCG from FMA toolbox
% Morci Juan Facundo 07/2025


% Initialization of structures that wil contain the outputs
Response.dHPC.aversive = [];   Response.vHPC.aversive = [];  Response.joint.aversive = [];
Response.dHPC.reward = [];   Response.vHPC.reward = [];  Response.joint.reward = [];

modulated.dHPC.aversive = [];   modulated.vHPC.aversive = [];  modulated.joint.aversive = [];
modulated.dHPC.reward = [];   modulated.vHPC.reward = [];  modulated.joint.reward = [];

Response.speed.data = [];          Response.speed.id = [];

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
        
        %% Awake
        disp('Uploading digital imputs')
        % Load digitalin.mat
        load('digitalin.mat')
        load('behavioral_dataVF.mat')
        %Shocks selection
        Shocks_filt = Restrict(shock,aversiveTS_run ./1000);
        % Keep only the first shock of each TTL (first from 20)
        count = 1;
        deff = [];
        for i = 1:length(Shocks_filt)
            if count == 1
                deff = [deff ; Shocks_filt(i,1)];
            end
            if count ==20
                count = 0;
            end
            count = count + 1;
        end
        Shocks_filt = deff;
        clear count deff shock i
        
        %Rewards selection
        Rewards_filt = sort([leftvalve ; rightvalve]);
        Rewards_filt = Restrict(Rewards_filt, [behavior.speed.reward(1,1) behavior.speed.reward(end,1)]);
        
        %% speed sourrounding the Valves opening
        means = meanInGroups(behavior.speed.reward(:,2), 3);
        downsampled_t = downsampleTimeVector(behavior.speed.reward(:,1), 1/30, 0.1);
        
        if size(means,2) < size(downsampled_t,2)
            means = [downsampled_t(1:size(means,2))' , means']; clear downsampled_t
        elseif size(means,2) > size(downsampled_t,2)
            means = [downsampled_t' , means(1:size(downsampled_t,2))']; clear downsampled_t
        else
            means = [downsampled_t' , means']; clear downsampled_t
        end
        
        tmp = [];
        for i = 1 : length(Rewards_filt(:,1))
            [~ , ii] = min(abs(means(:,1)-Rewards_filt(i,1)));
            if ii+50 < length(means)
                if ii-50>0
                    tmp = [tmp , means(ii-50 : ii+50 , 2)];
                end
            end
        end
        Mean = nanmean(tmp'); clear tmp means
        Response.speed.data = [Response.speed.data , Mean']; clear Mean
        Response.speed.id = [Response.speed.id ; tt t length(Rewards_filt(:,1))];        
        
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
        
        if or(numberD >3 , numberV > 3)
            %% --- Aversive ---
            disp('Lets go for the assemblies')
            if isfile('dorsalventral_assemblies_aversiveVF.mat')
                disp('Loading Aversive template')
                load('dorsalventral_assemblies_aversiveVF.mat')
            else
                Th = [];
                pat = [];
            end
            
            Thresholded.aversive.all = Th;
            patterns.all.aversive = pat.*Th;
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
            if isfile('dorsalventral_assemblies_rewardVF.mat')
                load('dorsalventral_assemblies_rewardVF.mat')
            else
                Th = [];
                pat = [];
            end
            Thresholded.reward.all = Th;
            patterns.all.reward = pat.*Th;
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
            
            
            %% SpikeTrain
            limits = [0 segments.Var1(end)/1000];
            cellulartype = [K(:,1) , K(:,4)];
            [Spikes , bins , Clusters] = spike_train_construction([spks_dHPC;spks_vHPC], clusters.all, cellulartype, 0.025, limits, [] , false, true);
            clear limits
            
            % Joint assemblies
            if sum(cond.both.aversive)>0
                [curve , time , responsive] = Assemblies_responsivness_fixed(patterns.all.aversive , cond.both.aversive , [bins' , Spikes] , Rewards_filt(:,1) , [0 1] ,10 , 0.1 , 0 ,'zscore' , 2);
                Response.joint.aversive = [Response.joint.aversive , curve]; clear curve
                
                tmp = ones(length(responsive'),2);
                tmp = [tmp(:,1)*tt , tmp(:,2)*t];
                
                modulated.joint.aversive = [modulated.joint.aversive ; responsive' , tmp]; clear surrogate tmp
            end
            
            if sum(cond.both.reward)>0
                [curve , time , responsive] = Assemblies_responsivness_fixed(patterns.all.reward , cond.both.reward , [bins' , Spikes] , Rewards_filt(:,1) , [0 1] ,10 , 0.1 , 0 ,'zscore' , 2);
                Response.joint.reward = [Response.joint.reward , curve]; clear curve
                
                tmp = ones(length(responsive'),2);
                tmp = [tmp(:,1)*tt , tmp(:,2)*t];
                
                modulated.joint.reward = [modulated.joint.reward ; responsive' tmp]; clear surrogate tmp                
            end
            
            % dHPC assemblies
            if sum(cond.dHPC.aversive)>0
                [curve , time , responsive] = Assemblies_responsivness_fixed(patterns.all.aversive , cond.dHPC.aversive , [bins' , Spikes] , Rewards_filt(:,1) , [0 1] ,10 , 0.1 , 0 ,'zscore' , 2);
                Response.dHPC.aversive = [Response.dHPC.aversive , curve]; clear curve
                
                tmp = ones(length(responsive'),2);
                tmp = [tmp(:,1)*tt , tmp(:,2)*t];
                
                modulated.dHPC.aversive = [modulated.dHPC.aversive ; responsive' tmp]; clear surrogate tmp
            end
            
            if sum(cond.dHPC.reward)>0
                [curve , time , responsive] = Assemblies_responsivness_fixed(patterns.all.reward , cond.dHPC.reward , [bins' , Spikes] , Rewards_filt(:,1) , [0 1] ,10 , 0.1 , 0 ,'zscore' , 2);
                Response.dHPC.reward = [Response.dHPC.reward , curve]; clear curve
                
                tmp = ones(length(responsive'),2);
                tmp = [tmp(:,1)*tt , tmp(:,2)*t];                
                
                modulated.dHPC.reward = [modulated.dHPC.reward ; responsive' tmp]; clear surrogate tmp                
            end
            
            
            % vHPC assemblies
            if sum(cond.vHPC.aversive)>0
                [curve , time , responsive] = Assemblies_responsivness_fixed(patterns.all.aversive , cond.vHPC.aversive , [bins' , Spikes] , Rewards_filt(:,1) , [0 1] ,10 , 0.1 , 0 ,'zscore' , 2);
                Response.vHPC.aversive = [Response.vHPC.aversive , curve]; clear curve
                
                tmp = ones(length(responsive'),2);
                tmp = [tmp(:,1)*tt , tmp(:,2)*t];
                
                modulated.vHPC.aversive = [modulated.vHPC.aversive ; responsive' tmp]; clear surrogate tmp
            end
            
            if sum(cond.vHPC.reward)>0
                [curve , time , responsive] = Assemblies_responsivness_fixed(patterns.all.reward , cond.vHPC.reward , [bins' , Spikes] , Rewards_filt(:,1) , [0 1] ,10 , 0.1 , 0 ,'zscore' , 2);
                Response.vHPC.reward = [Response.vHPC.reward , curve]; clear curve
                
                tmp = ones(length(responsive'),2);
                tmp = [tmp(:,1)*tt , tmp(:,2)*t];                
                
                modulated.vHPC.reward = [modulated.vHPC.reward ; responsive' tmp]; clear surrogate tmp                
            end
            
            clear aversiveTS aversiveTS_run baselineTS bins Cell_type_classification
            clear cellulartype clusters Events group_dHPC group_vHPC i ii iii K
            clear Kinfo NREM REM numberD numberV p parameter patterns post pre
            clear SpikesD SpikesV spks spks_dHPC spks_vHPC Thresholded WAKE
            clear Kinfo K ii iii ripples ripplesD ripplesV ripple_event Spikes
            clear Shocks_filt pos posx posy rewardTS rewardTS_run 
            clear start stop Rewards_filt segments Clusters behavior camaraR cond 
        end
    end
end

s = 0;
%% Reward
% Joint
% Mean Response Reward Increased response
% joint
figure
x = Response.joint.reward(:,modulated.joint.reward(:,1)==1);
plot(time,nanmean(x,2),'r'),hold on
ciplot(nanmean(x,2)-nansem(x')' , nanmean(x,2)+nansem(x')' , time,'r'),alpha 0.5
hold on

x = Response.joint.reward(:,not(modulated.joint.reward(:,1)==1));
plot(time,nanmean(x,2),'k'),hold on
ciplot(nanmean(x,2)-nansem(x')' , nanmean(x,2)+nansem(x')' , time,'k'),alpha 0.5

xline(0)
xline(1)

ylim([-1 3])

% dHPC
figure
x = Response.dHPC.reward(:,modulated.dHPC.reward(:,1)==1);
plot(time,nanmean(x,2),'r'),hold on
ciplot(nanmean(x,2)-nansem(x')' , nanmean(x,2)+nansem(x')' , time,'r'),alpha 0.5
hold on

x = Response.dHPC.reward(:,not(modulated.dHPC.reward(:,1)==1));
plot(time,nanmean(x,2),'k'),hold on
ciplot(nanmean(x,2)-nansem(x')' , nanmean(x,2)+nansem(x')' , time,'k'),alpha 0.5

xline(0)
xline(1)

ylim([-1 3])

% vHPC
figure
x = Response.vHPC.reward(:,modulated.vHPC.reward(:,1)==1);
plot(time,nanmean(x,2),'r'),hold on
ciplot(nanmean(x,2)-nansem(x')' , nanmean(x,2)+nansem(x')' , time,'r'),alpha 0.5
hold on

x = Response.vHPC.reward(:,not(modulated.vHPC.reward(:,1)==1));
plot(time,nanmean(x,2),'k'),hold on
ciplot(nanmean(x,2)-nansem(x')' , nanmean(x,2)+nansem(x')' , time,'k'),alpha 0.5

xline(0)
xline(1)

ylim([-1 3])

% Percentages
figure,
p1 = (sum(modulated.joint.reward(:,1)==1)/length(modulated.joint.reward(:,1)))*100;
p2 = (sum(not(modulated.joint.reward(:,1)==1))/length(modulated.joint.reward(:,1)))*100;

subplot(131),pie([p1 , p2] , {'Responsive' , 'none'})

p1 = (sum(modulated.dHPC.reward(:,1)==1)/length(modulated.dHPC.reward(:,1)))*100;
p2 = (sum(not(modulated.dHPC.reward(:,1)==1))/length(modulated.dHPC.reward(:,1)))*100;

subplot(132),pie([p1 , p2] , {'increased' ,'none'})

p1 = (sum(modulated.vHPC.reward(:,1)==1)/length(modulated.vHPC.reward(:,1)))*100;
p2 = (sum(not(modulated.vHPC.reward(:,1)==1))/length(modulated.vHPC.reward(:,1)))*100;

subplot(133),pie([p1 , p2] , {'increased' ,'none'})


%% Aversive
% Joint
% Mean Response Shock Increased response
% joint
figure
x = Response.joint.aversive(:,modulated.joint.aversive(:,1)==1);
plot(time,nanmean(x,2),'r'),hold on
ciplot(nanmean(x,2)-nansem(x')' , nanmean(x,2)+nansem(x')' , time,'r'),alpha 0.5
hold on

x = Response.joint.aversive(:,not(modulated.joint.aversive(:,1)==1));
plot(time,nanmean(x,2),'k'),hold on
ciplot(nanmean(x,2)-nansem(x')' , nanmean(x,2)+nansem(x')' , time,'k'),alpha 0.5

xline(0)
xline(1)

ylim([-1 3])

% dHPC
figure
x = Response.dHPC.aversive(:,modulated.dHPC.aversive(:,1)==1);
plot(time,nanmean(x,2),'r'),hold on
ciplot(nanmean(x,2)-nansem(x')' , nanmean(x,2)+nansem(x')' , time,'r'),alpha 0.5
hold on

x = Response.dHPC.aversive(:,not(modulated.dHPC.aversive(:,1)==1));
plot(time,nanmean(x,2),'k'),hold on
ciplot(nanmean(x,2)-nansem(x')' , nanmean(x,2)+nansem(x')' , time,'k'),alpha 0.5

xline(0)
xline(1)

ylim([-1 3])

% vHPC
figure
x = Response.vHPC.aversive(:,modulated.vHPC.aversive(:,1)==1);
plot(time,nanmean(x,2),'r'),hold on
ciplot(nanmean(x,2)-nansem(x')' , nanmean(x,2)+nansem(x')' , time,'r'),alpha 0.5
hold on

x = Response.vHPC.aversive(:,not(modulated.vHPC.aversive(:,1)==1));
plot(time,nanmean(x,2),'k'),hold on
ciplot(nanmean(x,2)-nansem(x')' , nanmean(x,2)+nansem(x')' , time,'k'),alpha 0.5

xline(0)
xline(1)

ylim([-1 3])

% Percentages
figure,
p1 = (sum(modulated.joint.aversive(:,1)==1)/length(modulated.joint.aversive(:,1)))*100;
p2 = (sum(not(modulated.joint.aversive(:,1)==1))/length(modulated.joint.aversive(:,1)))*100;

subplot(131),pie([p1 , p2] , {'Responsive' , 'none'})

p1 = (sum(modulated.dHPC.aversive(:,1)==1)/length(modulated.dHPC.aversive(:,1)))*100;
p2 = (sum(not(modulated.dHPC.aversive(:,1)==1))/length(modulated.dHPC.aversive(:,1)))*100;

subplot(132),pie([p1 , p2] , {'increased' ,'none'})

p1 = (sum(modulated.vHPC.aversive(:,1)==1)/length(modulated.vHPC.aversive(:,1)))*100;
p2 = (sum(not(modulated.vHPC.aversive(:,1)==1))/length(modulated.vHPC.aversive(:,1)))*100;

subplot(133),pie([p1 , p2] , {'increased' ,'none'})



end