clear
clc
% close all

%% Parameters
path = {'E:\Rat126\Ephys\in_Pyr';'E:\Rat103\usable';'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat132\recordings\in_pyr';'E:\Rat165\in_pyr\'};%List of folders from the path



%% For revisions
% dHPC
Amplitude.dHPC.aversive.pre.mean = [];
Amplitude.dHPC.aversive.post.mean = [];
Amplitude.dHPC.reward.pre.mean = [];
Amplitude.dHPC.reward.post.mean = [];

Amplitude.dHPC.aversive.pre.all = [];
Amplitude.dHPC.aversive.post.all = [];
Amplitude.dHPC.reward.pre.all = [];
Amplitude.dHPC.reward.post.all = [];

Duration.dHPC.aversive.pre.mean = [];
Duration.dHPC.aversive.post.mean = [];
Duration.dHPC.reward.pre.mean = [];
Duration.dHPC.reward.post.mean = [];

Duration.dHPC.aversive.pre.all = [];
Duration.dHPC.aversive.post.all = [];
Duration.dHPC.reward.pre.all = [];
Duration.dHPC.reward.post.all = [];

% vHPC
Amplitude.vHPC.aversive.pre.mean = [];
Amplitude.vHPC.aversive.post.mean = [];
Amplitude.vHPC.reward.pre.mean = [];
Amplitude.vHPC.reward.post.mean = [];

Amplitude.vHPC.aversive.pre.all = [];
Amplitude.vHPC.aversive.post.all = [];
Amplitude.vHPC.reward.pre.all = [];
Amplitude.vHPC.reward.post.all = [];

Duration.vHPC.aversive.pre.mean = [];
Duration.vHPC.aversive.post.mean = [];
Duration.vHPC.reward.pre.mean = [];
Duration.vHPC.reward.post.mean = [];

Duration.vHPC.aversive.pre.all = [];
Duration.vHPC.aversive.post.all = [];
Duration.vHPC.reward.pre.all = [];
Duration.vHPC.reward.post.all = [];

%% Dorsal-ventrla ripple coordination
Proportions.DorsalVentral.aversive.pre = [];
Proportions.DorsalVentral.aversive.post = [];
Proportions.DorsalVentral.reward.pre = [];
Proportions.DorsalVentral.reward.post = [];

Proportions.VentralDorsal.aversive.pre = [];
Proportions.VentralDorsal.aversive.post = [];
Proportions.VentralDorsal.reward.pre = [];
Proportions.VentralDorsal.reward.post = [];


%% Main tloop to iterate across sessions
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
        
        disp('Uploading session time stamps')
        %Loading TS of the sessions
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
        
        %% Sleep
        disp('Uploading sleep scoring')
        x = dir([cd,'\*-states.mat']);    states = load([cd,'\',x.name]);    states = states.states;
        
        REM.all = ToIntervals(states==5);    NREM.all = ToIntervals(states==3);    WAKE.all = ToIntervals(states==1);
        clear x states
        
        %keep only WAKE in HomeCage
        %         WAKE = Restrict(WAKE, [aversiveTS ; rewardTS ; baselineTS] ./1000);
        
        % NREM events restriction according conditions
        NREM.B = Restrict(NREM.all,baselineTS./1000);
        NREM.R = Restrict(NREM.all,rewardTS./1000);
        NREM.A = Restrict(NREM.all,aversiveTS./1000);
        % REM events restriction according conditions
        REM.B = Restrict(REM.all,baselineTS./1000);
        REM.R = Restrict(REM.all,rewardTS./1000);
        REM.A = Restrict(REM.all,aversiveTS./1000);
        
        %% Upload ripples
        if isfile('ripplesD_customized2.csv')
            disp('Uploading dRipples')
            ripplesD = table2array(readtable('ripplesD_customized2.csv'));
            D = true;
        else
            disp('No data related to dRipples')
            D = false;
        end
        
        if isfile('ripplesV_customized2.csv')
            disp('Uploading vRipples')
            ripplesV = table2array(readtable('ripplesV_customized2.csv'));
            V = true;
        else
            disp('No data related to vRipples')
            V = false;
        end
        
        %% ---- dRipples ----
        if D 
            disp('Lets analyze dRipples')
            if aversiveTS(1) < rewardTS(1)
                x = Restrict(ripplesD,NREM.B);
                y = Restrict(ripplesD,NREM.A);
                Amplitude.dHPC.aversive.pre.mean = [Amplitude.dHPC.aversive.pre.mean ; nanmean(x(:,4))];
                Amplitude.dHPC.aversive.post.mean = [Amplitude.dHPC.aversive.post.mean ; nanmean(y(:,4))];
                Amplitude.dHPC.aversive.pre.all = [Amplitude.dHPC.aversive.pre.all ; x(:,4)];
                Amplitude.dHPC.aversive.post.all = [Amplitude.dHPC.aversive.post.all ; y(:,4)];
                Duration.dHPC.aversive.pre.mean = [Duration.dHPC.aversive.pre.mean ; nanmean([x(:,3) - x(:,1)])];
                Duration.dHPC.aversive.post.mean = [Duration.dHPC.aversive.post.mean ; nanmean([y(:,3) - y(:,1)])];
                Duration.dHPC.aversive.pre.all = [Duration.dHPC.aversive.pre.all ; [x(:,3) - x(:,1)]];
                Duration.dHPC.aversive.post.all = [Duration.dHPC.aversive.post.all ; [y(:,3) - y(:,1)]];
                clear x y
                
                x = Restrict(ripplesD,NREM.A);
                y = Restrict(ripplesD,NREM.R);
                Amplitude.dHPC.reward.pre.mean = [Amplitude.dHPC.reward.pre.mean ; nanmean(x(:,4))];
                Amplitude.dHPC.reward.post.mean = [Amplitude.dHPC.reward.post.mean ; nanmean(y(:,4))];
                Amplitude.dHPC.reward.pre.all = [Amplitude.dHPC.reward.pre.all ; x(:,4)];
                Amplitude.dHPC.reward.post.all = [Amplitude.dHPC.reward.post.all ; y(:,4)];
                Duration.dHPC.reward.pre.mean = [Duration.dHPC.reward.pre.mean ; nanmean([x(:,3) - x(:,1)])];
                Duration.dHPC.reward.post.mean = [Duration.dHPC.reward.post.mean ; nanmean([y(:,3) - y(:,1)])];
                Duration.dHPC.reward.pre.all = [Duration.dHPC.reward.pre.all ; [x(:,3) - x(:,1)]];
                Duration.dHPC.reward.post.all = [Duration.dHPC.reward.post.all ; [y(:,3) - y(:,1)]];
                clear x y
            else
                x = Restrict(ripplesD,NREM.R);
                y = Restrict(ripplesD,NREM.A);
                Amplitude.dHPC.aversive.pre.mean = [Amplitude.dHPC.aversive.pre.mean ; nanmean(x(:,4))];
                Amplitude.dHPC.aversive.post.mean = [Amplitude.dHPC.aversive.post.mean ; nanmean(y(:,4))];
                Amplitude.dHPC.aversive.pre.all = [Amplitude.dHPC.aversive.pre.all ; x(:,4)];
                Amplitude.dHPC.aversive.post.all = [Amplitude.dHPC.aversive.post.all ; y(:,4)];
                Duration.dHPC.aversive.pre.mean = [Duration.dHPC.aversive.pre.mean ; nanmean([x(:,3) - x(:,1)])];
                Duration.dHPC.aversive.post.mean = [Duration.dHPC.aversive.post.mean ; nanmean([y(:,3) - y(:,1)])];
                Duration.dHPC.aversive.pre.all = [Duration.dHPC.aversive.pre.all ; [x(:,3) - x(:,1)]];
                Duration.dHPC.aversive.post.all = [Duration.dHPC.aversive.post.all ; [y(:,3) - y(:,1)]];
                clear x y
                
                x = Restrict(ripplesD,NREM.B);
                y = Restrict(ripplesD,NREM.R);
                Amplitude.dHPC.reward.pre.mean = [Amplitude.dHPC.reward.pre.mean ; nanmean(x(:,4))];
                Amplitude.dHPC.reward.post.mean = [Amplitude.dHPC.reward.post.mean ; nanmean(y(:,4))];
                Amplitude.dHPC.reward.pre.all = [Amplitude.dHPC.reward.pre.all ; x(:,4)];
                Amplitude.dHPC.reward.post.all = [Amplitude.dHPC.reward.post.all ; y(:,4)];
                Duration.dHPC.reward.pre.mean = [Duration.dHPC.reward.pre.mean ; nanmean([x(:,3) - x(:,1)])];
                Duration.dHPC.reward.post.mean = [Duration.dHPC.reward.post.mean ; nanmean([y(:,3) - y(:,1)])];
                Duration.dHPC.reward.pre.all = [Duration.dHPC.reward.pre.all ; [x(:,3) - x(:,1)]];
                Duration.dHPC.reward.post.all = [Duration.dHPC.reward.post.all ; [y(:,3) - y(:,1)]];
                clear x y
            end
        end
        
        %% ---- vRipples ----
        if V
            disp('Lets analyze vRipples')
            if aversiveTS(1) < rewardTS(1)
                x = Restrict(ripplesV,NREM.B);
                y = Restrict(ripplesV,NREM.A);
                Amplitude.vHPC.aversive.pre.mean = [Amplitude.vHPC.aversive.pre.mean ; nanmean(x(:,4))];
                Amplitude.vHPC.aversive.post.mean = [Amplitude.vHPC.aversive.post.mean ; nanmean(y(:,4))];
                Amplitude.vHPC.aversive.pre.all = [Amplitude.vHPC.aversive.pre.all ; x(:,4)];
                Amplitude.vHPC.aversive.post.all = [Amplitude.vHPC.aversive.post.all ; y(:,4)];
                Duration.vHPC.aversive.pre.mean = [Duration.vHPC.aversive.pre.mean ; nanmean([x(:,3) - x(:,1)])];
                Duration.vHPC.aversive.post.mean = [Duration.vHPC.aversive.post.mean ; nanmean([y(:,3) - y(:,1)])];
                Duration.vHPC.aversive.pre.all = [Duration.vHPC.aversive.pre.all ; [x(:,3) - x(:,1)]];
                Duration.vHPC.aversive.post.all = [Duration.vHPC.aversive.post.all ; [y(:,3) - y(:,1)]];
                clear x y
                
                x = Restrict(ripplesV,NREM.A);
                y = Restrict(ripplesV,NREM.R);
                Amplitude.vHPC.reward.pre.mean = [Amplitude.vHPC.reward.pre.mean ; nanmean(x(:,4))];
                Amplitude.vHPC.reward.post.mean = [Amplitude.vHPC.reward.post.mean ; nanmean(y(:,4))];
                Amplitude.vHPC.reward.pre.all = [Amplitude.vHPC.reward.pre.all ; x(:,4)];
                Amplitude.vHPC.reward.post.all = [Amplitude.vHPC.reward.post.all ; y(:,4)];
                Duration.vHPC.reward.pre.mean = [Duration.vHPC.reward.pre.mean ; nanmean([x(:,3) - x(:,1)])];
                Duration.vHPC.reward.post.mean = [Duration.vHPC.reward.post.mean ; nanmean([y(:,3) - y(:,1)])];
                Duration.vHPC.reward.pre.all = [Duration.vHPC.reward.pre.all ; [x(:,3) - x(:,1)]];
                Duration.vHPC.reward.post.all = [Duration.vHPC.reward.post.all ; [y(:,3) - y(:,1)]];
                clear x y
            else
                x = Restrict(ripplesV,NREM.R);
                y = Restrict(ripplesV,NREM.A);
                Amplitude.vHPC.aversive.pre.mean = [Amplitude.vHPC.aversive.pre.mean ; nanmean(x(:,4))];
                Amplitude.vHPC.aversive.post.mean = [Amplitude.vHPC.aversive.post.mean ; nanmean(y(:,4))];
                Amplitude.vHPC.aversive.pre.all = [Amplitude.vHPC.aversive.pre.all ; x(:,4)];
                Amplitude.vHPC.aversive.post.all = [Amplitude.vHPC.aversive.post.all ; y(:,4)];
                Duration.vHPC.aversive.pre.mean = [Duration.vHPC.aversive.pre.mean ; nanmean([x(:,3) - x(:,1)])];
                Duration.vHPC.aversive.post.mean = [Duration.vHPC.aversive.post.mean ; nanmean([y(:,3) - y(:,1)])];
                Duration.vHPC.aversive.pre.all = [Duration.vHPC.aversive.pre.all ; [x(:,3) - x(:,1)]];
                Duration.vHPC.aversive.post.all = [Duration.vHPC.aversive.post.all ; [y(:,3) - y(:,1)]];
                clear x y
                
                x = Restrict(ripplesV,NREM.B);
                y = Restrict(ripplesV,NREM.R);
                Amplitude.vHPC.reward.pre.mean = [Amplitude.vHPC.reward.pre.mean ; nanmean(x(:,4))];
                Amplitude.vHPC.reward.post.mean = [Amplitude.vHPC.reward.post.mean ; nanmean(y(:,4))];
                Amplitude.vHPC.reward.pre.all = [Amplitude.vHPC.reward.pre.all ; x(:,4)];
                Amplitude.vHPC.reward.post.all = [Amplitude.vHPC.reward.post.all ; y(:,4)];
                Duration.vHPC.reward.pre.mean = [Duration.vHPC.reward.pre.mean ; nanmean([x(:,3) - x(:,1)])];
                Duration.vHPC.reward.post.mean = [Duration.vHPC.reward.post.mean ; nanmean([y(:,3) - y(:,1)])];
                Duration.vHPC.reward.pre.all = [Duration.vHPC.reward.pre.all ; [x(:,3) - x(:,1)]];
                Duration.vHPC.reward.post.all = [Duration.vHPC.reward.post.all ; [y(:,3) - y(:,1)]];
                clear x y
            end            
            
            
        end
        
        %% Coordinated dHPC ripples
        if and(D,V)
            coordinatedD = [];
            coordinatedV = [];
            DorsalVentral = [];
            VentralDorsal = [];
            Total = [];
            for i = 1:length(ripplesD)
                r = ripplesD(i,:);
                tmp = sum(and(ripplesV(:,2)>= r(1,2)-0.1, ripplesV(:,2)<= r(1,2)+0.1));
                if tmp>0
                    z = ripplesV(and(ripplesV(:,2)>= r(1,2)-0.1, ripplesV(:,2)<= r(1,2)+0.1),:);
                    coordinatedV = [coordinatedV ; z];
                    [p,indice] = min(abs(r(2)-z(:,2)));
                    coordinatedD = [coordinatedD ; r];

                    if z(indice,2) < r(2)
                        DorsalVentral = [DorsalVentral ; r(2)];
                    else
                        VentralDorsal = [VentralDorsal ; z(indice,2)];
                    end
                    Total = [Total ; z(indice,2)];
                    clear tmp2 tmp1 p indice z
                end
                clear r
            end
            
            if aversiveTS(1,1) < rewardTS(1,1)
                % --- Dorsal-Ventral ---
                % Aversive
                P1 = Restrict(DorsalVentral,NREM.B);
                P1 = size(P1,1)/size(Restrict(Total,NREM.B),1);
                
                P2 = Restrict(DorsalVentral,NREM.A);
                P2 = size(P2,1)/size(Restrict(Total,NREM.A),1);
                
                Proportions.DorsalVentral.aversive.pre = [Proportions.DorsalVentral.aversive.pre ; P1];
                Proportions.DorsalVentral.aversive.post = [Proportions.DorsalVentral.aversive.post ; P2]; clear P1 P2
                
                % Reward
                P1 = Restrict(DorsalVentral,NREM.A);
                P1 = size(P1,1)/size(Restrict(Total,NREM.A),1);
                
                P2 = Restrict(DorsalVentral,NREM.R);
                P2 = size(P2,1)/size(Restrict(Total,NREM.R),1);
                                
                Proportions.DorsalVentral.reward.pre = [Proportions.DorsalVentral.reward.pre ; P1];
                Proportions.DorsalVentral.reward.post = [Proportions.DorsalVentral.reward.post ; P2]; clear P1 P2
            
                % --- Ventral-Dorsal ----
                % Aversive
                P1 = Restrict(VentralDorsal,NREM.B);
                P1 = size(P1,1)/size(Restrict(Total,NREM.B),1);
                
                P2 = Restrict(VentralDorsal,NREM.A);
                P2 = size(P2,1)/size(Restrict(Total,NREM.A),1);
                
                Proportions.VentralDorsal.aversive.pre = [Proportions.VentralDorsal.aversive.pre ; P1];
                Proportions.VentralDorsal.aversive.post = [Proportions.VentralDorsal.aversive.post ; P2]; clear P1 P2
                
                % Reward
                P1 = Restrict(VentralDorsal,NREM.A);
                P1 = size(P1,1)/size(Restrict(Total,NREM.A),1);
                
                P2 = Restrict(VentralDorsal,NREM.R);
                P2 = size(P2,1)/size(Restrict(Total,NREM.R),1);
                                
                Proportions.VentralDorsal.reward.pre = [Proportions.VentralDorsal.reward.pre ; P1];
                Proportions.VentralDorsal.reward.post = [Proportions.VentralDorsal.reward.post ; P2]; clear P1 P2
            else
                % --- Dorsal-Ventral ---
                % Aversive
                P1 = Restrict(DorsalVentral,NREM.R);
                P1 = size(P1,1)/size(Restrict(Total,NREM.R),1);
                
                P2 = Restrict(DorsalVentral,NREM.A);
                P2 = size(P2,1)/size(Restrict(Total,NREM.A),1);
                
                Proportions.DorsalVentral.aversive.pre = [Proportions.DorsalVentral.aversive.pre ; P1];
                Proportions.DorsalVentral.aversive.post = [Proportions.DorsalVentral.aversive.post ; P2]; clear P1 P2
                
                % Reward
                P1 = Restrict(DorsalVentral,NREM.B);
                P1 = size(P1,1)/size(Restrict(Total,NREM.B),1);
                
                P2 = Restrict(DorsalVentral,NREM.R);
                P2 = size(P2,1)/size(Restrict(Total,NREM.R),1);
                                
                Proportions.DorsalVentral.reward.pre = [Proportions.DorsalVentral.reward.pre ; P1];
                Proportions.DorsalVentral.reward.post = [Proportions.DorsalVentral.reward.post ; P2]; clear P1 P2
            
                % --- Ventral-Dorsal ----
                % Aversive
                P1 = Restrict(VentralDorsal,NREM.R);
                P1 = size(P1,1)/size(Restrict(Total,NREM.R),1);
                
                P2 = Restrict(VentralDorsal,NREM.A);
                P2 = size(P2,1)/size(Restrict(Total,NREM.A),1);
                
                Proportions.VentralDorsal.aversive.pre = [Proportions.VentralDorsal.aversive.pre ; P1];
                Proportions.VentralDorsal.aversive.post = [Proportions.VentralDorsal.aversive.post ; P2]; clear P1 P2
                
                % Reward
                P1 = Restrict(VentralDorsal,NREM.B);
                P1 = size(P1,1)/size(Restrict(Total,NREM.B),1);
                
                P2 = Restrict(VentralDorsal,NREM.R);
                P2 = size(P2,1)/size(Restrict(Total,NREM.R),1);
                                
                Proportions.VentralDorsal.reward.pre = [Proportions.VentralDorsal.reward.pre ; P1];
                Proportions.VentralDorsal.reward.post = [Proportions.VentralDorsal.reward.post ; P2]; clear P1 P2                
            end
            
        end
        clear ripplesD ripplesV
        clear aversiveTS aversiveTS_run rewardTS rewardTS_run baselineTS
        clear coordinated coordinatedA coordinatedB coordinatedR
        clear coordinatedV coordinatedA_V coordinatedB_V coordinatedR_V
        clear coordinatedA_V_non_refined coordinatedB_V_non_refined coordinatedR_V_non_refined
        clear uncoordinated uncoordinatedA uncoordinatedA_V uncoordinatedB uncoordinatedB_V
        clear uncoordinatedR uncoordinatedR_V uncoordinatedV
        clear REM REM_A REM_B REM_R NREM NREM_A NREM_B NREM_R WAKE D V A
        clear burstD burstV ccg coordinated_ripple_bursts q R tmp tmpA_D tmpA_V tmpB_D tmpB_V tmpR_D tmpR_V
    end
end

%% Mean
% --- Amplitude ---
% dHPC
figure,
% Aversive
m1 = Amplitude.dHPC.aversive.pre.mean;
m2 = Amplitude.dHPC.aversive.post.mean;
[h p] = ranksum(m1,m2)
x = [ones(size(m1)) ; ones(size(m2))*2];
y = [m1 ; m2];

subplot(221),
hold on

for i = 1:length(m1)
    plot([1 2], [m1(i) m2(i)], 'k-', 'Color', [0.5 0.5 0.5 0.6])  % líneas grises semi-transparentes
end
% scatter(x, y, 50, 'b', 'filled', 'jitter','on', 'jitterAmount',0.1)
% scatter([1 2] , [nanmean(m1) nanmean(m2)],'filled')
boxplot(y,x)
xlim([0 3])
ylim([0 60])

% Reward
m1 = Amplitude.dHPC.reward.pre.mean;
m2 = Amplitude.dHPC.reward.post.mean;
[h p] = ranksum(m1,m2)
x = [ones(size(m1)) ; ones(size(m2))*2];
y = [m1 ; m2];

subplot(222)
hold on
for i = 1:length(m1)
    plot([1 2], [m1(i) m2(i)], 'k-', 'Color', [0.5 0.5 0.5 0.6])  % líneas grises semi-transparentes
end
% scatter(x, y, 50, 'b', 'filled', 'jitter','on', 'jitterAmount',0.1)
% scatter([1 2] , [nanmean(m1) nanmean(m2)],'filled')
boxplot(y,x)
xlim([0 3])
ylim([0 60])

% vHPC
% Aversive
m1 = Amplitude.vHPC.aversive.pre.mean;
m2 = Amplitude.vHPC.aversive.post.mean;
[h p] = ranksum(m1,m2)
x = [ones(size(m1)) ; ones(size(m2))*2];
y = [m1 ; m2];

subplot(223)
hold on
for i = 1:length(m1)
    plot([1 2], [m1(i) m2(i)], 'k-', 'Color', [0.5 0.5 0.5 0.6])  % líneas grises semi-transparentes
end
% scatter(x, y, 50, 'b', 'filled', 'jitter','on', 'jitterAmount',0.1)
% scatter([1 2] , [nanmean(m1) nanmean(m2)],'filled')
boxplot(y,x)
xlim([0 3])
ylim([5 30])

% Reward
m1 = Amplitude.vHPC.reward.pre.mean;
m2 = Amplitude.vHPC.reward.post.mean;
[h p] = ranksum(m1,m2)
x = [ones(size(m1)) ; ones(size(m2))*2];
y = [m1 ; m2];

subplot(224)
hold on
for i = 1:length(m1)
    plot([1 2], [m1(i) m2(i)], 'k-', 'Color', [0.5 0.5 0.5 0.6])  % líneas grises semi-transparentes
end
% scatter(x, y, 50, 'b', 'filled', 'jitter','on', 'jitterAmount',0.1)
% scatter([1 2] , [nanmean(m1) nanmean(m2)],'filled')
boxplot(y,x)
xlim([0 3])
ylim([5 30])

% --- Duration ---
% dHPC
figure,
% Aversive
m1 = Duration.dHPC.aversive.pre.mean;
m2 = Duration.dHPC.aversive.post.mean;
[h p] = ranksum(m1,m2)
x = [ones(size(m1)) ; ones(size(m2))*2];
y = [m1 ; m2];

subplot(221),
hold on

for i = 1:length(m1)
    plot([1 2], [m1(i) m2(i)], 'k-', 'Color', [0.5 0.5 0.5 0.6])  % líneas grises semi-transparentes
end
% scatter(x, y, 50, 'b', 'filled', 'jitter','on', 'jitterAmount',0.1)
% scatter([1 2] , [nanmean(m1) nanmean(m2)],'filled')
boxplot(y,x)
xlim([0 3])
ylim([0.04 0.055])

% Reward
m1 = Duration.dHPC.reward.pre.mean;
m2 = Duration.dHPC.reward.post.mean;
[h p] = ranksum(m1,m2)
x = [ones(size(m1)) ; ones(size(m2))*2];
y = [m1 ; m2];

subplot(222)
hold on
for i = 1:length(m1)
    plot([1 2], [m1(i) m2(i)], 'k-', 'Color', [0.5 0.5 0.5 0.6])  % líneas grises semi-transparentes
end
% scatter(x, y, 50, 'b', 'filled', 'jitter','on', 'jitterAmount',0.1)
% scatter([1 2] , [nanmean(m1) nanmean(m2)],'filled')
boxplot(y,x)
xlim([0 3])
ylim([0.04 0.055])

% vHPC
% Aversive
m1 = Duration.vHPC.aversive.pre.mean;
m2 = Duration.vHPC.aversive.post.mean;
[h p] = ranksum(m1,m2)
x = [ones(size(m1)) ; ones(size(m2))*2];
y = [m1 ; m2];

subplot(223)
hold on
for i = 1:length(m1)
    plot([1 2], [m1(i) m2(i)], 'k-', 'Color', [0.5 0.5 0.5 0.6])  % líneas grises semi-transparentes
end
% scatter(x, y, 50, 'b', 'filled', 'jitter','on', 'jitterAmount',0.1)
% scatter([1 2] , [nanmean(m1) nanmean(m2)],'filled')
boxplot(y,x)
xlim([0 3])
ylim([0.03 0.06])

% Reward
m1 = Duration.vHPC.reward.pre.mean;
m2 = Duration.vHPC.reward.post.mean;
[h p] = ranksum(m1,m2)
x = [ones(size(m1)) ; ones(size(m2))*2];
y = [m1 ; m2];

subplot(224)
hold on
for i = 1:length(m1)
    plot([1 2], [m1(i) m2(i)], 'k-', 'Color', [0.5 0.5 0.5 0.6])  % líneas grises semi-transparentes
end
% scatter(x, y, 50, 'b', 'filled', 'jitter','on', 'jitterAmount',0.1)
% scatter([1 2] , [nanmean(m1) nanmean(m2)],'filled')
boxplot(y,x)
xlim([0 3])
ylim([0.03 0.06])

%% Dorsal-Ventral Ripple Coordination
m1 = Proportions.DorsalVentral.aversive.pre;
m2 = Proportions.DorsalVentral.aversive.post;
m3 = m2./m1;

m4 = Proportions.DorsalVentral.reward.pre;
m5 = Proportions.DorsalVentral.reward.post;
m6 = m5./m4;

m7 = Proportions.VentralDorsal.aversive.pre;
m8 = Proportions.VentralDorsal.aversive.post;
m9 = m8./m7;


m10 = Proportions.VentralDorsal.reward.pre;
m11 = Proportions.VentralDorsal.reward.post;
m12 = m11./m10;


% --- Analisis ---
% DorsalVentral
[h p] = ranksum(m3 , m6)
[p h] = signrank(m3 ,1)
[p h] = signrank(m6 ,1)
% VentralDorsal
[h p] = ranksum(m9 , m12)
[p h] = signrank(m9 ,1)
[p h] = signrank(m12 ,1)

x = [ones(size(m3)) ; ones(size(m6))*2 ; ones(size(m9))*3 ; ones(size(m12))*4];
y = [m3 ; m6 ; m9 ; m12];

figure,
scatter(x , y , 'filled' , 'jitter' , 0.1),hold on
scatter([1 2 3 4] , [nanmean(m3) nanmean(m6) nanmean(m9) nanmean(m12)] , 'filled')
xlim([0 5]),
ylim([])
yline(1,'--')



%% Plot raw data
% 1 and 3 are pre-sleep, 2 and 4 are post-sleep
% 1 and 2 are dorsal-ventral events, while 3 and 4 are ventral-dorsal
figure,
% Aversive
subplot(121)
x = [ones(size(m1)) ; ones(size(m2))*2 ; ones(size(m7))*3 ; ones(size(m8))*4];
y = [m1 ; m2 ; m7 ; m8];

scatter(x , y , 'filled' , 'jitter' , 0.1),hold on
scatter([1 2 3 4] , [nanmean(m1) nanmean(m2) nanmean(m7) nanmean(m8)] , 'filled')
xlim([0 5]),
ylim([0.30 0.70]), hold on
ylabel('Proportions of events (%)')
title('Aversive')

% Reward
subplot(122)
x = [ones(size(m4)) ; ones(size(m5))*2 ; ones(size(m10))*3 ; ones(size(m11))*4];
y = [m4 ; m5 ; m10 ; m11];

scatter(x , y , 'filled' , 'jitter' , 0.1),hold on
scatter([1 2 3 4] , [nanmean(m4) nanmean(m5) nanmean(m10) nanmean(m11)] , 'filled')
xlim([0 5]),
ylim([0.30 0.70]), hold on
ylabel('Proportions of events (%)')
title('Reward')