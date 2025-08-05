%% With RankOrder

clear
clc
close all

% Parameters
path = 'E:\output_replay';%List of folders from the path
files = dir(path);

path2 = 'Z:\Facundo\light_dataset_FM';
template = dir(path2);

folders = []; for x = 3 : length(template), folders = [folders ; template(x).name]; end

Replay.aversive.dHPC.cooridnated = [];
Replay.aversive.dHPC.uncooridnated = [];
Replay.aversive.dHPC.all = [];

Replay.reward.dHPC.cooridnated = [];
Replay.reward.dHPC.uncooridnated = [];
Replay.reward.dHPC.all = [];

for t = 3:length(files)
    f = [files(t).folder,'\',files(t).name];
    load(f)
    
    template = [path2,'\',files(t).name(1:6),'\',files(t).name(1:end-4)];
    load([template,'\Session_Structure.mat'])
    
    %% Aversive
    % --- Coordinated ---
    q = quantile(abs(replay.coordinated.dHPC.ave.rankOrder.rankOrd_shuf),0.94);
    % Pre-sleep
    if TimeStamps.Sleep.Aversive(1) < TimeStamps.Sleep.Reward(1)
        p = TimeStamps.Sleep.Baseline;
    else
        p = TimeStamps.Sleep.Reward;
    end
    
    y = replay.coordinated.dHPC.timestamps(:,2);
    x = and(y > p(1) , y < p(2));
    
    % Check Replay with more than 3 neurons active
    y = replay.coordinated.dHPC.ave.rankOrder.nCells > 2;
    total = and(x , y');
    CoordPre = sum(total);
    % Check Replay higher than shuffle
    y = abs(replay.coordinated.dHPC.ave.rankOrder.rankOrd)>=q;
    x = and(total , y'); 
    Suma = sum(x);
    
    %store
    Pre = sum(x);
    Pre = Pre/sum(total); clear total

    
    % Post-sleep
    y = replay.coordinated.dHPC.timestamps(:,2);
    x = and(y > TimeStamps.Sleep.Aversive(1) , y < TimeStamps.Sleep.Aversive(2));
    % Check Replay with more than 3 neurons active
    y = replay.coordinated.dHPC.ave.rankOrder.nCells > 2;
    total = and(x , y');
    CoordPost = sum(total);
    % Check Replay higher than shuffle
    y = abs(replay.coordinated.dHPC.ave.rankOrder.rankOrd)>=q;
    x = and(total , y');    
    Suma1 = sum(x);
    %store
    Post = sum(x);
    Post = Post /sum(total);
    
    Replay.aversive.dHPC.cooridnated = [Replay.aversive.dHPC.cooridnated ; Pre Post Suma Suma1]; clear Pre Post total q
    
    % --- Uncooridnated ---
    q = quantile(abs(replay.uncoordinated.dHPC.ave.rankOrder.rankOrd_shuf),0.94);
    % Pre-sleep
    if TimeStamps.Sleep.Aversive(1) < TimeStamps.Sleep.Reward(1)
        p = TimeStamps.Sleep.Baseline;
    else
        p = TimeStamps.Sleep.Reward;
    end
    
    y = replay.uncoordinated.dHPC.timestamps(:,2);
    x = and(y > p(1) , y < p(2));
    
    % Check Replay with more than 3 neurons active
    y = replay.uncoordinated.dHPC.ave.rankOrder.nCells > 2;
    total = and(x , y');
    
       % Subsampling
    if sum(total) > CoordPost
        Pre = [];
        for i = 1 : 1000
            % Subsampling if needed
            idx_1 = find(total == 1);      % x = 1
            idx_0 = find(total == 0);      % x = 0
            keep = randsample(idx_1, CoordPre);
            Total = zeros(size(total));
            Total(keep) = 1;
            
            % Check Replay higher than shuffle
            y = abs(replay.uncoordinated.dHPC.ave.rankOrder.rankOrd)>=q;
            x = and(Total , y');
            %store
            x = sum(x);
            Pre = [Pre ; x /sum(Total)]; clear Total
        end
        Pre = nanmean(Pre);
    else
        % Check Replay higher than shuffle
        y = abs(replay.uncoordinated.dHPC.ave.rankOrder.rankOrd)>=q;
        x = and(total , y');
        %store
        x = sum(x);
        Pre =  x /sum(total);
    end
    
    % Post-sleep
    y = replay.uncoordinated.dHPC.timestamps(:,2);
    x = and(y > TimeStamps.Sleep.Aversive(1) , y < TimeStamps.Sleep.Aversive(2));
    % Check Replay with more than 3 neurons active
    y = replay.uncoordinated.dHPC.ave.rankOrder.nCells > 2;
    total = and(x , y');
    
    % Subsampling
    if sum(total) > CoordPost
        Post = [];
        for i = 1 : 1000
            % Subsampling if needed
            idx_1 = find(total == 1);      % x = 1
            idx_0 = find(total == 0);      % x = 0
            keep = randsample(idx_1, CoordPre);
            Total = zeros(size(total));
            Total(keep) = 1;
            
            % Check Replay higher than shuffle
            y = abs(replay.uncoordinated.dHPC.ave.rankOrder.rankOrd)>=q;
            x = and(Total , y');
            %store
            x = sum(x);
            Post = [Post ; x /sum(Total)]; clear Total
        end
        Post = nanmean(Post);
    else
        % Check Replay higher than shuffle
        y = abs(replay.uncoordinated.dHPC.ave.rankOrder.rankOrd)>=q;
        x = and(total , y');
        %store
        x = sum(x);
        Post =  x /sum(total);
    end   
    
    Replay.aversive.dHPC.uncooridnated = [Replay.aversive.dHPC.uncooridnated ; Pre Post]; clear Pre Post total
    
    %% Reward
    % --- Coordinated ---
    q = quantile(abs(replay.coordinated.dHPC.rew.rankOrder.rankOrd_shuf),0.94);
    % Pre-sleep
    if TimeStamps.Sleep.Reward(1) < TimeStamps.Sleep.Aversive(1)
        p = TimeStamps.Sleep.Baseline;
    else
        p = TimeStamps.Sleep.Aversive;
    end
    
    y = replay.coordinated.dHPC.timestamps(:,2);
    x = and(y > p(1) , y < p(2));
    
    % Check Replay with more than 3 neurons active
    y = replay.coordinated.dHPC.rew.rankOrder.nCells > 2;
    total = and(x , y');
    CoordPre = sum(total);
    % Check Replay higher than shuffle
    y = abs(replay.coordinated.dHPC.rew.rankOrder.rankOrd)>=q;
    x = and(total , y'); 
    Suma = sum(x);

    %store
    Pre = sum(x);
    Pre = Pre/sum(total); clear total

    
    % Post-sleep
    y = replay.coordinated.dHPC.timestamps(:,2);
    x = and(y > TimeStamps.Sleep.Reward(1) , y < TimeStamps.Sleep.Reward(2));
    % Check Replay with more than 3 neurons active
    y = replay.coordinated.dHPC.rew.rankOrder.nCells > 2;
    total = and(x , y');
    CoordPost = sum(total);
    % Check Replay higher than shuffle
    y = abs(replay.coordinated.dHPC.rew.rankOrder.rankOrd)>=q;
    x = and(total , y');  
    Suma1 = sum(x);
    
    %store
    Post = sum(x);
    Post = Post /sum(total);
    
    Replay.reward.dHPC.cooridnated = [Replay.reward.dHPC.cooridnated ; Pre Post Suma Suma1]; clear Pre Post total q Suma Suma1
    
    % --- Uncooridnated ---
    q = quantile(abs(replay.uncoordinated.dHPC.rew.rankOrder.rankOrd_shuf),0.94);
    % Pre-sleep
    if TimeStamps.Sleep.Reward(1) < TimeStamps.Sleep.Aversive(1)
        p = TimeStamps.Sleep.Baseline;
    else
        p = TimeStamps.Sleep.Aversive;
    end
    
    y = replay.uncoordinated.dHPC.timestamps(:,2);
    x = and(y > p(1) , y < p(2));
    
    % Check Replay with more than 3 neurons active
    y = replay.uncoordinated.dHPC.rew.rankOrder.nCells > 2;
    total = and(x , y');
    
       % Subsampling
    if sum(total) > CoordPost
        Pre = [];
        for i = 1 : 1000
            % Subsampling if needed
            idx_1 = find(total == 1);      % x = 1
            idx_0 = find(total == 0);      % x = 0
            keep = randsample(idx_1, CoordPre);
            Total = zeros(size(total));
            Total(keep) = 1;
            
            % Check Replay higher than shuffle
            y = abs(replay.uncoordinated.dHPC.rew.rankOrder.rankOrd)>=q;
            x = and(Total , y');
            %store
            x = sum(x);
            Pre = [Pre ; x /sum(Total)]; clear Total
        end
        Pre = nanmean(Pre);
    else
        % Check Replay higher than shuffle
        y = abs(replay.uncoordinated.dHPC.rew.rankOrder.rankOrd)>=q;
        x = and(total , y');
        %store
        x = sum(x);
        Pre =  x /sum(total);
    end
    
    % Post-sleep
    y = replay.uncoordinated.dHPC.timestamps(:,2);
    x = and(y > TimeStamps.Sleep.Reward(1) , y < TimeStamps.Sleep.Reward(2));
    % Check Replay with more than 3 neurons active
    y = replay.uncoordinated.dHPC.rew.rankOrder.nCells > 2;
    total = and(x , y');
    
    % Subsampling
    if sum(total) > CoordPost
        Post = [];
        for i = 1 : 1000
            % Subsampling if needed
            idx_1 = find(total == 1);      % x = 1
            idx_0 = find(total == 0);      % x = 0
            keep = randsample(idx_1, CoordPre);
            Total = zeros(size(total));
            Total(keep) = 1;
            
            % Check Replay higher than shuffle
            y = abs(replay.uncoordinated.dHPC.rew.rankOrder.rankOrd)>=q;
            x = and(Total , y');
            %store
            x = sum(x);
            Post = [Post ; x /sum(Total)]; clear Total
        end
        Post = nanmean(Post);
    else
        % Check Replay higher than shuffle
        y = abs(replay.uncoordinated.dHPC.rew.rankOrder.rankOrd)>=q;
        x = and(total , y');
        %store
        x = sum(x);
        Post =  x /sum(total);
    end   
    
    Replay.reward.dHPC.uncooridnated = [Replay.reward.dHPC.uncooridnated ; Pre Post]; clear Pre Post total
    
end

%% Calcular relaciones Aversive
R1 = Replay.aversive.dHPC.uncooridnated(:,2) ./ Replay.aversive.dHPC.uncooridnated(:,1);
R2 = Replay.aversive.dHPC.cooridnated(:,2) ./ Replay.aversive.dHPC.cooridnated(:,1);
% 
% Cleaning Data with more than 10 events
r = not(Replay.aversive.dHPC.cooridnated(:,4)>1);
R1(r) = []; R2(r) = []; 

% Cleaning Nans
r = or(isnan(R1) , isnan(R2));
R1(r) = []; R2(r) = []; 
% Cleaning Inf
r = or(R1 == inf , R2 == inf);
R1(r) = []; R2(r) = []; 


% Generar valores con jitter manual
n = numel(R1);
jitterAmount = 0.1;
x1 = 1 + (rand(n,1)-0.5)*jitterAmount;
x2 = 2 + (rand(n,1)-0.5)*jitterAmount;

% Subplot
subplot(1,2,1)
hold on

% Dibujar líneas entre puntos pareados
for i = 1:n
    plot([x1(i), x2(i)], [R1(i), R2(i)], 'k-', 'LineWidth', 0.5, 'Color', [0.6 0.6 0.6])
end

% Dibujar los puntos individuales con jitter
scatter(x1, R1, 30, 'b', 'filled')
scatter(x2, R2, 30, 'r', 'filled')

% Dibujar los promedios
scatter([1 2], [nanmean(R1) nanmean(R2)], 100, 'k', 'filled', 'MarkerEdgeColor','k')

% Ajustes de ejes
xlim([0.5 2.5])
ylim([0 4])
xticks([1 2])
xticklabels({'Uncoordinated', 'Coordinated'})
ylabel('Post / Pre (%Replay Events)')
title('Replay Comparison')

% Estadística
[h, p] = ttest(R1 , R2);
text(1.5, 3.3, sprintf('p = %.3f', p), 'HorizontalAlignment', 'center')

[h p] = ttest(R1,1)
[h p] = ttest(R2,1)


%% Calcular relaciones Reward
R1 = Replay.reward.dHPC.uncooridnated(:,2) ./ Replay.reward.dHPC.uncooridnated(:,1);
R2 = Replay.reward.dHPC.cooridnated(:,2) ./ Replay.reward.dHPC.cooridnated(:,1);

% % Cleaning Data with more than 10 events
% r = not(Replay.reward.dHPC.cooridnated(:,4)>0);
% R1(r) = []; R2(r) = []; 

% Cleaning Nans
r = or(isnan(R1) , isnan(R2));
R1(r) = []; R2(r) = []; 
% Cleaning Inf
r = or(R1 == inf , R2 == inf);
R1(r) = []; R2(r) = []; 

% Generar valores con jitter manual
n = numel(R1);
jitterAmount = 0.1;
x1 = 1 + (rand(n,1)-0.5)*jitterAmount;
x2 = 2 + (rand(n,1)-0.5)*jitterAmount;

% Subplot
subplot(1,2,2)
hold on

% Dibujar líneas entre puntos pareados
for i = 1:n
    plot([x1(i), x2(i)], [R1(i), R2(i)], 'k-', 'LineWidth', 0.5, 'Color', [0.6 0.6 0.6])
end

% Dibujar los puntos individuales con jitter
scatter(x1, R1, 30, 'b', 'filled')
scatter(x2, R2, 30, 'r', 'filled')

% Dibujar los promedios
scatter([1 2], [nanmean(R1) nanmean(R2)], 100, 'k', 'filled', 'MarkerEdgeColor','k')

% Ajustes de ejes
xlim([0.5 2.5])
ylim([0 4])
xticks([1 2])
xticklabels({'Uncoordinated', 'Coordinated'})
ylabel('Post / Pre (%Replay Events)')
title('Replay Comparison')

% Estadística
[h, p] = ttest(R1 , R2);
text(1.5, 3.3, sprintf('p = %.3f', p), 'HorizontalAlignment', 'center')

[h p] = ttest(R1,1)
[h p] = ttest(R2,1)

%% All Riples Stats
figure
R1 = nanmean([Replay.aversive.dHPC.cooridnated(:,1:2) , Replay.aversive.dHPC.uncooridnated]').*100;
R2 = nanmean([Replay.reward.dHPC.cooridnated(:,1:2) , Replay.reward.dHPC.uncooridnated]').*100;

% Cleaning Nans
r = or(isnan(R1) , isnan(R2));
R1(r) = []; R2(r) = []; 
% Cleaning Inf
r = or(R1 == inf , R2 == inf);
R1(r) = []; R2(r) = []; 

% Generar valores con jitter manual
n = numel(R1);
jitterAmount = 0.1;
x1 = 1 + (rand(n,1)-0.5)*jitterAmount;
x2 = 2 + (rand(n,1)-0.5)*jitterAmount;

% Subplot
hold on

% Dibujar líneas entre puntos pareados
for i = 1:n
    plot([x1(i), x2(i)], [R1(i), R2(i)], 'k-', 'LineWidth', 0.5, 'Color', [0.6 0.6 0.6])
end

% Dibujar los puntos individuales con jitter
scatter(x1, R1, 30, 'b', 'filled')
scatter(x2, R2, 30, 'r', 'filled')

% Dibujar los promedios
scatter([1 2], [nanmean(R1) nanmean(R2)], 100, 'k', 'filled', 'MarkerEdgeColor','k')

% Ajustes de ejes
xlim([0.5 2.5])
% ylim([0 7])
xticks([1 2])
xticklabels({'Aversive', 'Reward'})
ylabel('Post / Pre (%Replay Events)')
title('Replay Comparison')

% Estadística
[h, p] = ttest(R1 , R2);
text(1.5, 3.3, sprintf('p = %.3f', p), 'HorizontalAlignment', 'center')

