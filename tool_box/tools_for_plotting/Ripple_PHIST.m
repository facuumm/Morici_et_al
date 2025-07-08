function output = Ripple_PHIST(path,Type,events)
% This function construct a CCG using SU activity lock the occurrence
% of dorsal or ventral ripples. It will iterate in the subfolder from the
% paths you introduce.
%
% Syntax: [output] = Ripple_PHIST(path,type)
%
% INPUTS
% path: cell, contains the paths of the sessions you want to analyze.
%
% type: int, put 0 if you want putative Pyr cells put 1 if you want putative
%       itnernerons
%
% events: strg, if 'Coordinated' will compute compute Gain in coordinated
%         and not coordinated events within structure.
%
%
% OUTPUT
% output.dHPC: PHIST of dorsal neurons
% output.vHPC: PHIST of ventral neurons
% output.up.dHPC/vHPC: PHIST of dorsal/ventral ripple up-modulated neurons
% output.down.dHPC/vHPC: PHIST of dorsal/ventral ripple down-modulated neurons
% output.non.dHPC/vHPC: PHIST of dorsal/ventral ripple non-modulated neurons
% output.time: time vector for plotting
%
% other functions: PHIST
%                  CCG (FMA toolbox)
% Morci Juan Facundo 12/2023

% Storage creation depending on the ripple type
if strcmp(events,'Coordinated')
    output.dHPC.coordinated = [];
    output.dHPC.uncoordinated = [];
    output.vHPC.coordinated = [];
    output.vHPC.uncoordinated = [];
    output.time = [];
else
    output.dHPC = [];
    output.vHPC = [];
    output.time = [];
    output.up.dHPC = [];
    output.up.vHPC = [];
    output.down.dHPC = [];
    output.down.vHPC = [];
    output.non.dHPC = [];
    output.non.vHPC = [];
end


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
        
        %% Segments of the session
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
        
        %% Sleep states
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
        
        %% Load ripples
        if isfile('ripplesD_customized2.csv')
            ripplesD = table2array(readtable('ripplesD_customized2.csv'));
            RD = true;
        else
            RD = false;
        end
        
        if isfile('ripplesV_customized2.csv')
            ripplesV = table2array(readtable('ripplesV_customized2.csv'));
            RV = true;
        else
            RV = False;
        end
        
        %% Spikes
        % Load Units
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
        if Type == 0 %pyr
            cellulartype = [K(:,1) , K(:,4)];
        elseif Type == 1 % int
            cellulartype = [K(:,1) , not(K(:,4))];
        end
        
        
        %% load file to separate data deppending on the type of response to the ripple
        load('ripple_modulated_SU.mat')
        
        %% Counting the Number f SU
        clusters.dHPC = [];
        for ii=1:size(group_dHPC,1)
            cluster = group_dHPC(ii,1);
            Cellulartype = logical(cellulartype(cellulartype(:,1) == cluster,2));
            if Cellulartype
                clusters.dHPC = [clusters.dHPC ; cluster];
            end
            clear tmp cluster Cellulartype fr1 fr2 fr3 fr4 fr5 r a
        end
        
        clusters.vHPC = [];
        for ii=1:size(group_vHPC,1)
            cluster = group_vHPC(ii,1);
            Cellulartype = logical(cellulartype(cellulartype(:,1) == cluster,2));
            if Cellulartype
                clusters.vHPC = [clusters.vHPC ; cluster];
            end
            clear tmp cluster Cellulartype fr1 fr2 fr3 fr4 fr5 r a
        end
        
        
        if strcmp(events,'Coordinated')
            if and(RV , RD)
                % coordinated events
                coordinatedD = [];
                coordinatedV = [];
                for i = 1:length(ripplesD)
                    r = ripplesD(i,:);
                    tmp = sum(and(ripplesV(:,2)>= r(1,2)-0.1, ripplesV(:,2)<= r(1,2)+0.1));
                    if tmp>0
                        z = ripplesV(and(ripplesV(:,2)>= r(1,2)-0.1, ripplesV(:,2)<= r(1,2)+0.1),:);
                        [p,indice] = min(abs(r(2)-z(:,2)));
                        coordinatedV = [coordinatedV ; z(indice,:)];
                        coordinatedD = [coordinatedD ; r];
                        clear tmp2 tmp1 p indice z
                    end
                    clear r
                end
                clear x tmp i
                [C,IA,IC] = unique(coordinatedV(:,1));
                coordinatedV  = coordinatedV(IA,:); clear C IA IC
                
                uncoordinatedD = ripplesD(~ismember(ripplesD(:,1),coordinatedD(:,1)),:);
                uncoordinatedV = ripplesV(~ismember(ripplesV(:,1),coordinatedV(:,1)),:);
                
                % CCG dHPC
                if not(isempty(spks_dHPC))
                    baseline = SubtractIntervals(NREM.all,[ripplesD(:,1)-0.05 ripplesD(:,3)+0.05]);
                    for i = 1:size(clusters.dHPC,1)
                        [p , time] = PHIST_Ripple_SU(coordinatedD,spks_dHPC(spks_dHPC(:,1)==clusters.dHPC(i),2),baseline,2,0.01,1,'Gain');
                        output.dHPC.coordinated = [output.dHPC.coordinated , p];
                        clear p
                        
                        [p , time] = PHIST_Ripple_SU(uncoordinatedD,spks_dHPC(spks_dHPC(:,1)==clusters.dHPC(i),2),baseline,2,0.01,1,'Gain');
                        output.dHPC.uncoordinated = [output.dHPC.uncoordinated , p];
                        clear p
                    end
                    clear i
                end
                
                if not(isempty(spks_vHPC))
                    baseline = SubtractIntervals(NREM.all,[ripplesV(:,1)-0.05 ripplesV(:,3)+0.05]);
                    for i = 1:size(clusters.vHPC,1)
                        [p , time] = PHIST_Ripple_SU(coordinatedV,spks_vHPC(spks_vHPC(:,1)==clusters.vHPC(i),2),baseline,2,0.01,1,'Gain');
                        output.vHPC.coordinated = [output.vHPC.coordinated , p];
                        clear p
                        
                        [p , time] = PHIST_Ripple_SU(uncoordinatedV,spks_vHPC(spks_vHPC(:,1)==clusters.vHPC(i),2),baseline,2,0.01,1,'Gain');
                        output.vHPC.uncoordinated = [output.vHPC.uncoordinated , p];
                        clear p
                    end
                    clear i
                end
                %         else
                %             %% CCG
                %             if not(isempty(spks_dHPC))
                %                 baseline = SubtractIntervals(NREM.all,[ripplesD(:,1)-0.05 ripplesD(:,3)+0.05]);
                %                 for i = 1:size(clusters.dHPC,1)
                %                     [p , time] = PHIST(ripplesD(:,2),spks_dHPC(spks_dHPC(:,1)==clusters.dHPC(i),2),baseline,1,0.01,0,'NormGain');
                %                     output.dHPC = [output.dHPC , p];
                %
                %                     c = ripple_modulated.dHPC.all(clusters.dHPC(i) == ripple_modulated.dHPC.all(:,1),:);
                %                     [p , time] = PHIST(ripplesD(:,2),spks_dHPC(spks_dHPC(:,1)==clusters.dHPC(i),2),baseline,1,0.01,0,'Gain');
                %                     if logical(c(1,3))
                %                         output.up.dHPC = [output.up.dHPC , p];
                %                     elseif logical(c(1,4))
                %                         output.down.dHPC = [output.down.dHPC , p];
                %                     elseif and(not(logical(c(1,3))) , not(logical(c(1,4))))
                %                         output.non.dHPC = [output.non.dHPC , p];
                %                     end
                %                     clear p c
                %                 end
                %             end
                %
                %             if not(isempty(spks_vHPC))
                %                 baseline = SubtractIntervals(NREM.all,[ripplesV(:,1)-0.05 ripplesV(:,3)+0.05]);
                %                 for i = 1:size(clusters.vHPC,1)
                %                     [p , time] = PHIST(ripplesV(:,2),spks_vHPC(spks_vHPC(:,1)==clusters.vHPC(i),2),baseline,1,0.01,0,'NormGain');
                %                     output.vHPC = [output.vHPC , p];
                %
                %                     c = ripple_modulated.vHPC.all(clusters.vHPC(i) == ripple_modulated.vHPC.all(:,1),:);
                %                     [p , time] = PHIST(ripplesV(:,2),spks_vHPC(spks_vHPC(:,1)==clusters.vHPC(i),2),baseline,1,0.01,0,'Gain');
                %                     if logical(c(1,3))
                %                         output.up.vHPC = [output.up.vHPC , p];
                %                     elseif logical(c(1,4))
                %                         output.down.vHPC = [output.down.vHPC , p];
                %                     elseif and(not(logical(c(1,3))) , not(logical(c(1,4))))
                %                         output.non.vHPC = [output.non.vHPC , p];
                %                     end
                %                     clear p c
                %                 end
                %             end
            end
        end
        disp('-- Finishing --')
    end
    disp(['-- Finishing analysis from rat #',num2str(tt) , ' --'])
    disp('  ')
    
end


% if strcmp(events,'Coordinated')
%     output.time = time;
%     [n i] = min(abs(output.time - (-0.05)));
%     [n ii] = min(abs(output.time - 0.05)); clear n
%     
%     % sort dHPC
%     % [i ii] = max(output.dHPC);
%     n = nanmean(output.dHPC.coordinated(i:ii,:));
%     m = nanmean(output.dHPC.coordinated(1:i,:));
%     m = mean([m;nanmean(output.dHPC.coordinated(ii:end,:))]);
%     n = m-n;
%     [n m] = sort(n,'ascend'); clear n
%     output.dHPC.coordinated = output.dHPC.coordinated(:,m);
%     output.dHPC.uncoordinated = output.dHPC.uncoordinated(:,m);
%     clear m
%     
%     % sort vHPC
%     % [i ii] = max(output.vHPC);
%     % [i ii] = max(output.dHPC);
%     n = nanmean(output.vHPC.coordinated(i:ii,:));
%     m = nanmean(output.vHPC.coordinated(1:i,:));
%     m = mean([m;nanmean(output.vHPC.coordinated(ii:end,:))]);
%     n = m-n;
%     [n m] = sort(n,'ascend'); clear n
%     output.vHPC.coordinated = output.vHPC.coordinated(:,m);
%     output.vHPC.uncoordinated = output.vHPC.uncoordinated(:,m);
%     clear m
% else
%     output.time = time;
%     [n i] = min(abs(output.time - (-0.05)));
%     [n ii] = min(abs(output.time - 0.05)); clear n
%     
%     % sort dHPC
%     % [i ii] = max(output.dHPC);
%     n = nanmean(output.dHPC(i:ii,:));
%     m = nanmean(output.dHPC(1:i,:));
%     m = mean([m;nanmean(output.dHPC(ii:end,:))]);
%     n = m-n;
%     [n m] = sort(n,'ascend'); clear n
%     output.dHPC = output.dHPC(:,m); clear m
%     
%     % sort vHPC
%     % [i ii] = max(output.vHPC);
%     n = nanmean(output.vHPC(i:ii,:));
%     m = nanmean(output.vHPC(1:i,:));
%     m = mean([m;nanmean(output.vHPC(ii:end,:))]);
%     n = m-n;
%     [n m] = sort(n,'ascend'); clear n
%     output.vHPC = output.vHPC(:,m); clear m i ii
% end

%% ---> For plotting data for Extended data requested by R2 and R3
figure
T = [-1 : 0.01 : 1];

subplot(121)
s = nansem((output.dHPC.uncoordinated)');
m = nanmean((output.dHPC.uncoordinated)');
plot(T,m,'k'),hold on
ciplot(m-s,m+s,T,'k'), alpha 0.5
s = nansem((output.dHPC.coordinated)');
m = nanmean((output.dHPC.coordinated)');
plot(T,m,'r'),hold on
ciplot(m-s,m+s,T,'r'), alpha 0.5
ylim([1 12])
xlim([-0.5 0.5])

subplot(122)
s = nansem((output.vHPC.uncoordinated)');
m = nanmean((output.vHPC.uncoordinated)');
plot(T,m,'k'),hold on
ciplot(m-s,m+s,T,'k'), alpha 0.5
s = nansem((output.vHPC.coordinated)');
m = nanmean((output.vHPC.coordinated)');
plot(T,m,'b'),hold on
ciplot(m-s,m+s,T,'b'), alpha 0.5
ylim([1 6])
xlim([-0.5 0.5])

% To compare differences
[~ , ff] = min(abs(T-0.05));
[~ , f] = min(abs(T-(-0.05)));
m1 = nanmean(output.dHPC.uncoordinated(f:ff,:));
m2 = nanmean(output.dHPC.coordinated(f:ff,:));
[~ , p] = lillietest(m1)
[~ , p] = lillietest(m2)
[~ , p] = ttest2(m1,m2)

m3 = nanmean(output.vHPC.uncoordinated(f:ff,:));
m4 = nanmean(output.vHPC.coordinated(f:ff,:));
[~ , p] = lillietest(m3)
[~ , p] = lillietest(m4)
[~ , p] = ttest2(m3,m4)

end
