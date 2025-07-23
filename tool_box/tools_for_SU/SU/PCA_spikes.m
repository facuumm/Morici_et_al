function [dHPC vHPC] = PCA_spikes(path)
% This function calculates the Pre and Post sleep assemblies rate.
%
% Morci Juan Facundo 08/2024

% variables to use in the script
criteria_fr = 0;
criteria_type = 0; % criteria for celltype (0:pyr, 1:int, 2:all)
sm = 0;
bin = 1;
count=0;

% storage variables
dHPC = [];      vHPC = [];      dvHPC = [];


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
        count = count+1;
        disp(['-- Initiating analysis of folder #' , num2str(t) , ' from rat #',num2str(tt) , ' --'])
        session = [subFolders(t+2).folder,'\',subFolders(t+2).name];
        cd(session)
        
        %Loading TS of the sessions
        disp('Uploading session time stamps')
        load('session_organization.mat')
        load('behavioral_data.mat','movement')
        
        %% load sleep states
        disp('Uploading sleep scoring')
        x = dir([cd,'\*-states.mat']);    states = load([cd,'\',x.name]);    states = states.states;
        REM.all = ToIntervals(states==5);    NREM.all = ToIntervals(states==3);
        clear x states
        
        baselineTS = baselineTS./1000;                  aversiveTS = aversiveTS./1000;                    rewardTS = rewardTS./1000;
        aversiveTS_run = aversiveTS_run./1000;          rewardTS_run = rewardTS_run./1000;
        
        %% Spikes
        % Load Units
        disp('Uploading Spiking activity')
        cd 'Spikesorting'
        [clusters , numberD , numberV , spks , spks_dHPC , spks_vHPC , cellulartype] = load_SU_FM(cd,criteria_type,criteria_fr,aversiveTS_run,rewardTS_run);
        if numberD > 2
            limits = [0 segments.Var1(end)/1000];
            [SpksTrains , Bins , C] = spike_train_construction(spks, clusters.dHPC, cellulartype , bin , limits , [] , true , false);
            
            i = InIntervals(Bins,movement.aversive);
            ii = InIntervals(Bins,movement.reward);
            [x] = pca([SpksTrains(i,:) ; SpksTrains(ii,:)]');
            
            centroide1 = mean(x(1:sum(i),1:2));
            centroide2 = mean(x(sum(i)+1:end,1:2));
            
            if and(norm(centroide1 - centroide2)>0.035 , norm(centroide1 - centroide2)<0.04)
%                 y = [zeros(sum(i),1) ; ones(sum(ii),1)];
                figure,
                subplot(121),scatter(x(1:sum(i),1),x(1:sum(i),2),20,'filled'),ylim([-0.1 0.2]),xlim([-0.05 0.2])
                subplot(122),scatter(x(sum(i)+1:end,1),x(sum(i)+1:end,2),20,'filled'),ylim([-0.1 0.2]),xlim([-0.05 0.2])
                title('dHPC')
            end            
            
            dHPC = [dHPC ; norm(centroide1 - centroide2)];
            clear centroide1 centroide2 x i ii SpksTrains Bins C limits y
        end
        
        if numberV > 2
            limits = [0 segments.Var1(end)/1000];
            [SpksTrains , Bins , C] = spike_train_construction(spks, clusters.vHPC, cellulartype , bin , limits , [] , true , false);
            
            i = InIntervals(Bins,movement.aversive);
            ii = InIntervals(Bins,movement.reward);
            [x] = pca([SpksTrains(i,:) ; SpksTrains(ii,:)]');

            
            centroide1 = mean(x(1:sum(i),1:2));
            centroide2 = mean(x(sum(i)+1:end,1:2));
            
            if norm(centroide1 - centroide2)>0.06
%                 y = [zeros(sum(i),1) ; ones(sum(ii),1)];
                figure,
                subplot(121),scatter(x(1:sum(i),1),x(1:sum(i),2),20,'filled'),ylim([-0.1 0.2]),xlim([-0.05 0.2])
                subplot(122),scatter(x(sum(i)+1:end,1),x(sum(i)+1:end,2),20,'filled'),ylim([-0.1 0.2]),xlim([-0.05 0.2])
                title('vHPC')
            end
            
             vHPC = [vHPC ; norm(centroide1 - centroide2)];
             clear centroide1 centroide2 x i ii SpksTrains Bins C limits y
        end
        
        clear aversiveTS aversiveTS_run baselineTS behaviorclusters config 
        clear minimal_speed minimal_speed_time movement REM NREM Shocks_filt
        clear spks spks_dHPC spks_vHPC SpksTrains numberD numberV Rewards_filt
        clear rewardTS rewardTS_run segments behavior ans
    end
end
x = [ones(length(dHPC),1) ; ones(length(vHPC),1)*2];
figure,boxplot([dHPC ; vHPC],x)
x1 = dHPC; x1(dHPC>0.045) = [];
x2 = vHPC;
[h p] = ranksum(x1,x2)
end