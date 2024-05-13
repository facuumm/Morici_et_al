function [Pre , Post , Run] = correlation_between_assemblies(path)
% This function calculates a CCG between dorsal and ventral assemblies
% in Pre, Post sleep, and Run. Also, it defines which assemblies are
% correlated (mean activation > 99percentile from a surrogate)
%
% INPUTS
% path: cell, inside each cell you should put the path of each session you
%       want to analyze
%
% OUTPUT
% Pre and Post: Structure, it store the CCG Assemblies between assemblies.
%
% Run: Structure, it store the CCG between assemblies during run.
%
% Iterator: Structure, it contains a logical array (1: correlated / 0: not)
%
% T: time vector.
%
%               Architecture of each output:
%                   Pre.aversive / reward
%                   Post.aversive / reward
%                   Run.aversive / reward
%                   Iterator.aversive / reward
%
% Morci Juan Facundo 04/2024

%variables to use in the script
criteria_fr = 0;

%variables for CCG construction
th = 5;           % threshold for peak detection
dur = 1;       % window for train of assemblies
b = 0.005;        % binsize

% storage variables
Pre.aversive = [];       Post.aversive = [];       Run. aversive = [];
Pre.reward = [];         Post.reward = [];         Run.reward = [];

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
        load('behavioral_data.mat')
        
        %% load sleep states
        disp('Uploading sleep scoring')
        x = dir([cd,'\*-states.mat']);    states = load([cd,'\',x.name]);    states = states.states;
        REM.all = ToIntervals(states==5);    NREM.all = ToIntervals(states==3);    WAKE.all = ToIntervals(states==1);
        clear x states
        
        NREM.baseline = Restrict(NREM.all,baselineTS./1000);
        NREM.aversive = Restrict(NREM.all,aversiveTS./1000);
        NREM.reward = Restrict(NREM.all,rewardTS./1000);
        
        REM.baseline = Restrict(REM.all,baselineTS./1000);
        REM.aversive = Restrict(REM.all,aversiveTS./1000);
        REM.reward = Restrict(REM.all,rewardTS./1000);
        
        
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
        % pyr
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
        
        %% Assemblies detection
        if and(numberD > 3 , numberV > 3)
            % --- Aversive ---
            disp('Lets go for the assemblies')
            if isfile('separated_assemblies.mat')
                load('separated_assemblies.mat')
            end
            
            % Definition of limits for assemblie activity strength
            limits = [0 segments.Var1(end)/1000];
            events = [];
            
            % Spiketrains construction
            [SpksTrains.dHPC , Binss , Cluster] = spike_train_construction(spks_dHPC, clusters.dHPC, cellulartype, b, limits, events, false,true);
            [SpksTrains.vHPC , Binss , Cluster] = spike_train_construction(spks_vHPC, clusters.vHPC, cellulartype, b, limits, events, false,true);
            
            if isfield(patterns,'aversive')
                disp('Aversive assemblies')
                if aversiveTS_run(1)<rewardTS_run(1)
                    TS.pre = Restrict(NREM.baseline,[NREM.baseline(end,2)-1800 NREM.baseline(end,2)]);
                    TS.post = Restrict(NREM.aversive,[NREM.aversive(1,1) NREM.aversive(1,1)+1800]);
                    TS.run = movement.aversive;
                else
                    TS.pre = Restrict(NREM.reward,[NREM.reward(end,2)-1800 NREM.reward(end,2)]);
                    TS.post = Restrict(NREM.aversive,[NREM.aversive(1,1) NREM.aversive(1,1)+1800]);
                    TS.run = movement.aversive;
                end
                
                % Peaks detection
                pks.vHPC.aversive = assemblies_peaks([Binss' SpksTrains.vHPC] , patterns.aversive.vHPC , th);
                pks.dHPC.aversive = assemblies_peaks([Binss' SpksTrains.dHPC] , patterns.aversive.dHPC , th);
                
                Trains.vHPC = [];
                for i = 1:size(pks.vHPC.aversive,1)
                    [tmp,bins]=binspikes(pks.vHPC.aversive{i}(:,1),1/dur,limits);
                    Trains.vHPC = [Trains.vHPC , zscore(tmp)]; clear tmp bins
                end
                
                Trains.dHPC = [];
                for i = 1:size(pks.dHPC.aversive,1)
                    [tmp,bins]=binspikes(pks.dHPC.aversive{i}(:,1),1/dur,limits);
                    Trains.dHPC = [Trains.dHPC , zscore(tmp)]; clear tmp
                end
                
                
                % Correlation RUN
                for i = 1 : size(Trains.dHPC,2)
                    tmp = Restrict([bins',Trains.dHPC(:,i)] , TS.run);
                    for ii = 1 : size(Trains.vHPC,2)
                        tmp1 = Restrict([bins',Trains.vHPC(:,ii)] , TS.run);
                        c = fitlm(tmp(:,2),tmp1(:,2));
%                         [c p] = corr(tmp(:,2),tmp1(:,2),'rows','pairwise');
%                         Run.aversive = [Run.aversive ; c p]; clear c tmp1
                        Run.aversive = [Run.aversive ; c.Rsquared.Ordinary c.Coefficients.pValue(end)]; clear c tmp1
                    end
                    clear tmp
                end
                
                % Correlation PRE
                for i = 1 : size(Trains.dHPC,2)
                    tmp = Restrict([bins',Trains.dHPC(:,i)] , TS.pre);
                    for ii = 1 : size(Trains.vHPC,2)
                        tmp1 = Restrict([bins',Trains.vHPC(:,ii)] , TS.pre);
                        c = fitlm(tmp(:,2),tmp1(:,2));
%                         [c p] = corr(tmp(:,2),tmp1(:,2),'rows','pairwise');
%                         Pre.aversive = [Pre.aversive ; c p]; clear c tmp1
                        Pre.aversive = [Pre.aversive ; c.Rsquared.Ordinary c.Coefficients.pValue(end)]; clear c tmp1
                    end
                    clear tmp
                end
                
                % Correlation POST
                for i = 1 : size(Trains.dHPC,2)
                    tmp = Restrict([bins',Trains.dHPC(:,i)] , TS.post);
                    for ii = 1 : size(Trains.vHPC,2)
                        tmp1 = Restrict([bins',Trains.vHPC(:,ii)] , TS.post);
                        c = fitlm(tmp(:,2),tmp1(:,2));
%                         [c p] = corr(tmp(:,2),tmp1(:,2),'rows','pairwise');
%                         Post.aversive = [Post.aversive ; c p]; clear c tmp1
                        Post.aversive = [Post.aversive ; c.Rsquared.Ordinary c.Coefficients.pValue(end)]; clear c tmp1
                    end
                    clear tmp
                end
                
                %                 Bins = [0:win:segments.Var1(end)/1000];
                %
                %                 if and(size(Trains.dHPC,2)>1 , size(Trains.vHPC,2)>1)
                %                     Corr  = [];
                %                     for i = 1 : size(Bins,2)-1
                %                         In = InIntervals(bins,[Bins(i) Bins(i+1)]);
                %                         % Correlation Matrix Calculation
                %                         x = Trains.dHPC(In,:);
                %                         y = Trains.vHPC(In,:);
                %                         [S1 , p] = corr(x,y);
                %                         tmp = [];
                %                         for ii = 1 : size(Bins,2)-1
                %                             In = InIntervals(bins,[Bins(ii) Bins(ii+1)]);
                %                             % Correlation Matrix Calculation
                %                             x = Trains.dHPC(In,:);
                %                             y = Trains.vHPC(In,:);
                %                             [S2 , p] = corr(x,y);
                %                             S3 = corrcoef(S1,S2,'rows','complete');
                %                             if size(S3,2)>1
                %                                 tmp = [tmp , S3(1,2)];
                %                             else
                %                                 tmp = [tmp , nan];
                %                             end
                %                         end
                %                         Corr  = [Corr ; tmp]; clear tmp
                %                     end
                %                     Corr = Corr - diag(diag(Corr));
                %                     %                 imagesc(Bins,Bins,Corr'),hold on
                %                     %                 yline(aversiveTS_run(1)/1000,'LineWidth',3)
                %                     %                 yline(aversiveTS_run(2)/1000,'LineWidth',3)
                %                     %                 xline(aversiveTS_run(1)/1000,'LineWidth',3)
                %                     %                 xline(aversiveTS_run(2)/1000,'LineWidth',3)
                %                     %                 yline(rewardTS_run(1)/1000,'LineWidth',1)
                %                     %                 yline(rewardTS_run(2)/1000,'LineWidth',1)
                %                     %                 xline(rewardTS_run(1)/1000,'LineWidth',1)
                %                     %                 xline(rewardTS_run(2)/1000,'LineWidth',1)
                %                     clear Trains
                %
                %                     x = InIntervals(Bins,baselineTS./1000); x = x*x';
                %                     y = InIntervals(Bins,aversiveTS./1000); y = y*y';
                %                     z = InIntervals(Bins,rewardTS./1000); z = z*z';
                %                     a = InIntervals(Bins,aversiveTS_run./1000); a = a*a';
                %                     r = InIntervals(Bins,rewardTS_run./1000); r = r*r';
                %                     template = logical(x+y+z+a+r); clear x y z a r
                %
                %                     x = Corr.*template(2:end,2:end);
                %                     x = x - triu(x);
                %                     x(x==0) = nan;
                %                     x = nanmean(x);
                %                     %                figure,plot(Bins(2:end),nanmean(x))
                %                     %                xline(aversiveTS_run(1)/1000,'LineWidth',3)
                %                     %                xline(aversiveTS_run(2)/1000,'LineWidth',3)
                %                     y =  nanmean(x(InIntervals(Bins(2:end),TS.pre)));
                %                     z =  nanmean(x(InIntervals(Bins(2:end),TS.post)));
                %
                %                     Pre.aversive = [Pre.aversive ; y];
                %                     Post.aversive = [Post.aversive ; z]; clear y z x Corr TS pks
                %                 end
            end
            
            if isfield(patterns,'reward')
                disp('Reward assemblies')
                if aversiveTS_run(1)<rewardTS_run(1)
                    TS.pre = Restrict(NREM.aversive,[NREM.aversive(end,2)-1800 NREM.aversive(end,2)]);
                    TS.post = Restrict(NREM.reward,[NREM.reward(1,1) NREM.reward(1,1)+1800]);
                    TS.run = movement.reward;
                else
                    TS.pre = Restrict(NREM.baseline,[NREM.baseline(end,2)-1800 NREM.baseline(end,2)]);
                    TS.post = Restrict(NREM.reward,[NREM.reward(1,1) NREM.reward(1,1)+1800]);
                    TS.run = movement.reward;
                end
                
                % Peaks detection
                pks.vHPC.reward = assemblies_peaks([Binss' SpksTrains.vHPC] , patterns.reward.vHPC , th);
                pks.dHPC.reward = assemblies_peaks([Binss' SpksTrains.dHPC] , patterns.reward.dHPC , th);
                clear SpksTrains Bins Cluster
                
                
                Trains.vHPC = [];
                for i = 1:size(pks.vHPC.reward,1)
                    [tmp,bins]=binspikes(pks.vHPC.reward{i}(:,1),1/dur,limits);
                    Trains.vHPC = [Trains.vHPC , zscore(tmp)]; clear tmp bins
                end
                
                Trains.dHPC = [];
                for i = 1:size(pks.dHPC.reward,1)
                    [tmp,bins]=binspikes(pks.dHPC.reward{i}(:,1),1/dur,limits);
                    Trains.dHPC = [Trains.dHPC , zscore(tmp)]; clear tmp
                end
                
                % Correlation RUN
                for i = 1 : size(Trains.dHPC,2)
                    tmp = Restrict([bins',Trains.dHPC(:,i)] , TS.run);
                    for ii = 1 : size(Trains.vHPC,2)
                        tmp1 = Restrict([bins',Trains.vHPC(:,ii)] , TS.run);
                        c = fitlm(tmp(:,2),tmp1(:,2));
%                         [c p] = corr(tmp(:,2),tmp1(:,2),'rows','pairwise');
%                         Run.reward = [Run.reward ; c p]; clear c tmp1
                        Run.reward = [Run.reward ; c.Rsquared.Ordinary c.Coefficients.pValue(end)]; clear c tmp1
                    end
                    clear tmp
                end
                
                % Correlation PRE
                for i = 1 : size(Trains.dHPC,2)
                    tmp = Restrict([bins',Trains.dHPC(:,i)] , TS.pre);
                    for ii = 1 : size(Trains.vHPC,2)
                        tmp1 = Restrict([bins',Trains.vHPC(:,ii)] , TS.pre);
                        c = fitlm(tmp(:,2),tmp1(:,2));
%                         [c p] = corr(tmp(:,2),tmp1(:,2),'rows','pairwise');
%                         Pre.reward = [Pre.reward ; c p]; clear c tmp1
                        Pre.reward = [Pre.reward ; c.Rsquared.Ordinary c.Coefficients.pValue(end)]; clear c tmp1
                    end
                    clear tmp
                end
                
                % Correlation POST
                for i = 1 : size(Trains.dHPC,2)
                    tmp = Restrict([bins',Trains.dHPC(:,i)] , TS.post);
                    for ii = 1 : size(Trains.vHPC,2)
                        tmp1 = Restrict([bins',Trains.vHPC(:,ii)] , TS.post);
                        c = fitlm(tmp(:,2),tmp1(:,2));
%                         [c p] = corr(tmp(:,2),tmp1(:,2),'rows','pairwise');
%                         Post.reward = [Post.reward ; c p]; clear c tmp1
                        Post.reward = [Post.reward ; c.Rsquared.Ordinary c.Coefficients.pValue(end)]; clear c tmp1
                    end
                    clear tmp
                end
                
                %                 Bins = [0:win:segments.Var1(end)/1000];
                %                 if and(size(Trains.dHPC,2)>1 , size(Trains.vHPC,2)>1)
                %                     Corr  = [];
                %                     for i = 1 : size(Bins,2)-1
                %                         In = InIntervals(bins,[Bins(i) Bins(i+1)]);
                %                         % Correlation Matrix Calculation
                %                         x = Trains.dHPC(In,:);
                %                         y = Trains.vHPC(In,:);
                %                         [S1 , p] = corr(x,y);
                %                         tmp = [];
                %                         for ii = 1 : size(Bins,2)-1
                %                             In = InIntervals(bins,[Bins(ii) Bins(ii+1)]);
                %                             % Correlation Matrix Calculation
                %                             x = Trains.dHPC(In,:);
                %                             y = Trains.vHPC(In,:);
                %                             [S2 , p] = corr(x,y);
                %                             S3 = corrcoef(S1,S2,'rows','complete');
                %                             if size(S3,2)>1
                %                                 tmp = [tmp , S3(1,2)];
                %                             else
                %                                 tmp = [tmp , nan];
                %                             end
                %                         end
                %                         Corr  = [Corr ; tmp]; clear tmp
                %                     end
                %                     Corr = Corr - diag(diag(Corr));
                %                     %                 imagesc(Bins,Bins,Corr'),hold on
                %                     %                 yline(aversiveTS_run(1)/1000,'LineWidth',3)
                %                     %                 yline(aversiveTS_run(2)/1000,'LineWidth',3)
                %                     %                 xline(aversiveTS_run(1)/1000,'LineWidth',3)
                %                     %                 xline(aversiveTS_run(2)/1000,'LineWidth',3)
                %                     %                 yline(rewardTS_run(1)/1000,'LineWidth',1)
                %                     %                 yline(rewardTS_run(2)/1000,'LineWidth',1)
                %                     %                 xline(rewardTS_run(1)/1000,'LineWidth',1)
                %                     %                 xline(rewardTS_run(2)/1000,'LineWidth',1)
                %                     clear Trains
                %
                %                     x = InIntervals(Bins,baselineTS./1000); x = x*x';
                %                     y = InIntervals(Bins,aversiveTS./1000); y = y*y';
                %                     z = InIntervals(Bins,rewardTS./1000); z = z*z';
                %                     a = InIntervals(Bins,aversiveTS_run./1000); a = a*a';
                %                     r = InIntervals(Bins,rewardTS_run./1000); r = r*r';
                %                     template = logical(x+y+z+a+r); clear x y z a r
                %
                %                     x = Corr.*template(2:end,2:end);
                %                     x = x - triu(x);
                %                     x(x==0) = nan;
                %                     x = nanmean(x);
                %                     %                figure,plot(Bins(2:end),nanmean(x))
                %                     %                xline(aversiveTS_run(1)/1000,'LineWidth',3)
                %                     %                xline(aversiveTS_run(2)/1000,'LineWidth',3)
                %                     y =  nanmean(x(InIntervals(Bins(2:end),TS.pre)));
                %                     z =  nanmean(x(InIntervals(Bins(2:end),TS.post)));
                %
                %                     Pre.reward = [Pre.reward ; y];
                %                     Post.reward = [Post.reward ; z]; clear y z x Corr TS pks
                %                 end
            end
            
            clear patterns aversiveTS aversiveTS_run rewardTS rewardTS_run baselineTS
            clear behavior movement Cell_type_classification cellulartype Cellulartype
            clear Bins Cluster clusters config group_dHPC group_vHPC K Kinfo limits
            clear numberD numberV NREM REM WAKE Rewards_filt Shocks_filt segments
            clear spks spks_dHPC spks_vHPC SpksTrains Thresholded
        end
    end
end


end



%% for plotting
x = Run.aversive(:,2)<0.05;
y = Run.reward(:,2)<0.05;

figure
% subplot(121),boxplot([Pre.aversive(x,1) , Post.aversive(x,1)])
[h p] = signrank(Pre.aversive(:,1) , Post.aversive(:,1))
% subplot(122),boxplot([Pre.reward(y,1) , Post.reward(y,1)])
[h p] = signrank(Pre.reward(:,1) , Post.reward(:,1))


figure
subplot(121)
for i = 1:length(Pre.aversive(:,1))
    x = 1+rand*0.1;
    scatter(x,Pre.aversive(i,1),'filled','k'),hold on
    y = 2+rand*0.1;
    scatter(y,Post.aversive(i,1),'filled','r'),hold on
    plot([x,y],[Pre.aversive(i,1) , Post.aversive(i,1)],'k-')
end
xlim([0 3])
ylim([-0.02 0.07])
scatter([1 2],[nanmedian(Pre.aversive(:,1)) nanmedian(Post.aversive(:,1))],'filled')

subplot(122)
for i = 1:length(Pre.reward(:,1))
    x = 1+rand*0.1;
    scatter(x,Pre.reward(i,1),'filled','k'),hold on
    y = 2+rand*0.1;
    scatter(y,Post.reward(i,1),'filled','b'),hold on
    plot([x,y],[Pre.reward(i,1) , Post.reward(i,1)],'k-')
end
xlim([0 3])
ylim([-0.02 0.07])
scatter([1 2],[nanmedian(Pre.reward(:,1)) nanmedian(Post.reward(:,1))],'filled')
