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
    % Pre-sleep
    if TimeStamps.Sleep.Aversive(1) < TimeStamps.Sleep.Reward(1)
        p = TimeStamps.Sleep.Baseline;
    else
        p = TimeStamps.Sleep.Reward;
    end
    
    y = replay.coordinated.dHPC.timestamps(:,2);
    x = and(y > p(1) , y < p(2));
    
    % Check Replay with more than 3 neurons active
    y = replay.coordinated.dHPC.ave.rankOrder.nCells >= 3;
    total = and(x , y');
    % Check Replay higher than shuffle
    y = replay.coordinated.dHPC.ave.rankOrder.rankOrd_shuf < replay.coordinated.dHPC.ave.rankOrder.rankOrd;
    x = and(total , y'); 
    
    %store
    Pre = sum(x);
    Pre = Pre/sum(total); clear total

    
    % Post-sleep
    y = replay.coordinated.dHPC.timestamps(:,2);
    x = and(y > TimeStamps.Sleep.Aversive(1) , y < TimeStamps.Sleep.Aversive(2));
    % Check Replay with more than 3 neurons active
    y = replay.coordinated.dHPC.ave.rankOrder.nCells > 3;
    total = and(x , y');
    % Check Replay higher than shuffle
    y = replay.coordinated.dHPC.ave.rankOrder.rankOrd_shuf < replay.coordinated.dHPC.ave.rankOrder.rankOrd;
    x = and(total , y');    
    %store
    Post = sum(x);
    Post = Post /sum(total);
    
    Replay.aversive.dHPC.cooridnated = [Replay.aversive.dHPC.cooridnated ; Pre Post]; clear Pre Post total
    
    % --- Uncooridnated ---
    % Pre-sleep
    if TimeStamps.Sleep.Aversive(1) < TimeStamps.Sleep.Reward(1)
        p = TimeStamps.Sleep.Baseline;
    else
        p = TimeStamps.Sleep.Reward;
    end
    
    y = replay.uncoordinated.dHPC.timestamps(:,2);
    x = and(y > p(1) , y < p(2));
    
    % Check Replay with more than 3 neurons active
    y = replay.uncoordinated.dHPC.ave.rankOrder.nCells > 3;
    total = and(x , y');
    % Check Replay higher than shuffle
    y = replay.uncoordinated.dHPC.ave.rankOrder.rankOrd_shuf < replay.uncoordinated.dHPC.ave.rankOrder.rankOrd;
    x = and(total , y');    
    %store
    Pre = sum(x);
    Pre = Pre /sum(total); clear total
    
    % Post-sleep
    y = replay.uncoordinated.dHPC.timestamps(:,2);
    x = and(y > TimeStamps.Sleep.Aversive(1) , y < TimeStamps.Sleep.Aversive(2));
    % Check Replay with more than 3 neurons active
    y = replay.uncoordinated.dHPC.ave.rankOrder.nCells > 3;
    total = and(x , y');
    % Check Replay higher than shuffle
    y = replay.uncoordinated.dHPC.ave.rankOrder.rankOrd_shuf < replay.uncoordinated.dHPC.ave.rankOrder.rankOrd;
    x = and(total , y');    
    %store
    Post = sum(x);
    Post = Post /sum(total);    
    
    Replay.aversive.dHPC.uncooridnated = [Replay.aversive.dHPC.uncooridnated ; Pre Post]; clear Pre Post total
    
    %% Reward
    % --- Coordinated ---
    % Check replay within Pre-sleep
    if TimeStamps.Sleep.Aversive(1) < TimeStamps.Sleep.Reward(1)
        p = TimeStamps.Sleep.Aversive;
    else
        p = TimeStamps.Sleep.Baseline;
    end
    
    y = replay.coordinated.dHPC.timestamps(:,2);
    x = and(y > p(1) , y < p(2));
    
    % Check Replay with more than 3 neurons active
    y = replay.coordinated.dHPC.rew.rankOrder.nCells > 3;
    total = and(x , y');
    % Check Replay higher than shuffle
    y = replay.coordinated.dHPC.rew.rankOrder.rankOrd_shuf < replay.coordinated.dHPC.rew.rankOrder.rankOrd;
    x = and(total , y');    
    %store
%     Pre = nanmean(replay.coordinated.dHPC.ave.rankOrder.rankOrd(x));
    Pre = sum(x);
    Pre = Pre /sum(total); clear total
    
    % Check replay within post-sleep
    y = replay.coordinated.dHPC.timestamps(:,2);
    x = and(y > TimeStamps.Sleep.Reward(1) , y < TimeStamps.Sleep.Reward(2));
    % Check Replay with more than 3 neurons active
    y = replay.coordinated.dHPC.rew.rankOrder.nCells > 3;
    total = and(x , y');
    % Check Replay higher than shuffle
    y = replay.coordinated.dHPC.rew.rankOrder.rankOrd_shuf < replay.coordinated.dHPC.rew.rankOrder.rankOrd;
    x = and(total , y');    
    %store
    Post = sum(x);
    Post = Post/sum(total);    
    
    Replay.reward.dHPC.cooridnated = [Replay.reward.dHPC.cooridnated ; Pre Post]; clear Pre Post total
    
    % --- Uncooridnated ---
    % Check replay within Pre-sleep
    if TimeStamps.Sleep.Aversive(1) < TimeStamps.Sleep.Reward(1)
        p = TimeStamps.Sleep.Aversive;
    else
        p = TimeStamps.Sleep.Baseline;
    end
    
    y = replay.uncoordinated.dHPC.timestamps(:,2);
    x = and(y > p(1) , y < p(2));
    
    % Check Replay with more than 3 neurons active
    y = replay.uncoordinated.dHPC.rew.rankOrder.nCells > 3;
    total = and(x , y');
    % Check Replay higher than shuffle
    y = replay.uncoordinated.dHPC.rew.rankOrder.rankOrd_shuf < replay.uncoordinated.dHPC.rew.rankOrder.rankOrd;
    x = and(total , y');    
    %store
    Pre = sum(x);
    Pre = Pre /sum(total); clear total    

    % Check replay within post-sleep
    y = replay.uncoordinated.dHPC.timestamps(:,2);
    x = and(y > TimeStamps.Sleep.Reward(1) , y < TimeStamps.Sleep.Reward(2));
    % Check Replay with more than 3 neurons active
    y = replay.uncoordinated.dHPC.rew.rankOrder.nCells > 3;
    total = and(x , y');
    % Check Replay higher than shuffle
    y = replay.uncoordinated.dHPC.rew.rankOrder.rankOrd_shuf < replay.uncoordinated.dHPC.rew.rankOrder.rankOrd;
    x = and(total , y');    
    %store
    Post = sum(x);
    Post = Post /sum(total);        
    
    Replay.reward.dHPC.uncooridnated = [Replay.reward.dHPC.uncooridnated ; Pre Post]; clear Pre Post total
    
end

%% Calcular relaciones Aversive
R1 = Replay.aversive.dHPC.uncooridnated(:,1) ./ Replay.aversive.dHPC.uncooridnated(:,1);
R2 = Replay.aversive.dHPC.cooridnated(:,1) ./ Replay.aversive.dHPC.cooridnated(:,1);

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
ylim([0.5 3.5])
xticks([1 2])
xticklabels({'Uncoordinated', 'Coordinated'})
ylabel('Post / Pre (%Replay Events)')
title('Replay Comparison')

% Estadística
[h, p] = ttest(R1 , R2);
text(1.5, 3.3, sprintf('p = %.3f', p), 'HorizontalAlignment', 'center')



%% Calcular relaciones Reward
R1 = Replay.reward.dHPC.uncooridnated(:,2) ./ Replay.reward.dHPC.uncooridnated(:,1);
R2 = Replay.reward.dHPC.cooridnated(:,2) ./ Replay.reward.dHPC.cooridnated(:,1);

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
ylim([0.5 3.5])
xticks([1 2])
xticklabels({'Uncoordinated', 'Coordinated'})
ylabel('Post / Pre (%Replay Events)')
title('Replay Comparison')

% Estadística
[h, p] = ttest(R1 , R2);
text(1.5, 3.3, sprintf('p = %.3f', p), 'HorizontalAlignment', 'center')

%% All Riples Stats
A_Pre = Replay.aversive.dHPC.all(:,1)*100;
A_Post = Replay.aversive.dHPC.all(:,1)*100;

R_Pre = Replay.reward.dHPC.all(:,1)*100;
R_Post = Replay.reward.dHPC.all(:,1)*100;
% Agrupar todos los datos
allData = [A_Pre, A_Post, R_Pre, R_Post];

% Número de sujetos
n = size(allData,1);

% Jitter
jitterAmount = 0.1;
x = repmat(1:4, n, 1) + (rand(n,4)-0.5)*jitterAmount;

% Subplot
figure,hold on

% Dibujar líneas entre pares Pre/Post para Aversive y Reward
for i = 1:n
    plot(x(i,[1 2]), allData(i,[1 2]), '-', 'Color', [0.5 0.2 0.2 0.3])  % Aversive
    plot(x(i,[3 4]), allData(i,[3 4]), '-', 'Color', [0.2 0.2 0.5 0.3])  % Reward
end

% Dibujar los puntos individuales
scatter(x(:,1), allData(:,1), 30, [0.8 0.2 0.2], 'filled') % Aversive Pre
scatter(x(:,2), allData(:,2), 30, [1.0 0.4 0.4], 'filled') % Aversive Post
scatter(x(:,3), allData(:,3), 30, [0.2 0.2 0.8], 'filled') % Reward Pre
scatter(x(:,4), allData(:,4), 30, [0.4 0.4 1.0], 'filled') % Reward Post

% Dibujar medias
means = nanmean(allData);
scatter(1:4, means, 100, 'k', 'filled', 'MarkerEdgeColor','k')

% Estética
xlim([0.5 4.5])
ylim([min(allData(:))*0.9, max(allData(:))*1.1])
xticks(1:4)
xticklabels({'Aversive Pre', 'Aversive Post', 'Reward Pre', 'Reward Post'})
ylabel('Replay Count')
title('Replay Events Across Conditions')
set(gca,'FontSize',10)

% Estadísticas (opcional)
[hA,pA] = ttest(A_Pre, A_Post);
[hR,pR] = ttest(R_Pre, R_Post);
text(1.5, max(means)*1.05, sprintf('Aversive p = %.3f', pA), 'HorizontalAlignment','center')
text(3.5, max(means)*1.05, sprintf('Reward p = %.3f', pR), 'HorizontalAlignment','center')