function [dHPC vHPC] = FiringRate_overlap(path)
% This function calculates the Firing Rate overlap between reward and
% aversive. 

% variables to use in the script
criteria_fr = 0;
criteria_type = 0; % criteria for celltype (0:pyr, 1:int, 2:all)
count=0;

% storage variables
vHPC.baseline = [];     vHPC.aversive = [];     vHPC.reward = [];       vHPC.overlap = [];
dHPC.baseline = [];     dHPC.aversive = [];     dHPC.reward = [];       dHPC.overlap = [];


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
        load('behavioral_data.mat')
        
        %% load sleep states
        disp('Uploading sleep scoring')
        x = dir([cd,'\*-states.mat']);    states = load([cd,'\',x.name]);    states = states.states;
        REM.all = ToIntervals(states==5);    NREM.all = ToIntervals(states==3);     WAKE.all = ToIntervals(states==1);
        
        clear x states
        
        baselineTS = baselineTS./1000;                  aversiveTS = aversiveTS./1000;                    rewardTS = rewardTS./1000;
        aversiveTS_run = aversiveTS_run./1000;          rewardTS_run = rewardTS_run./1000;
        
        WAKE.control = WAKE.all(InIntervals(WAKE.all(:,2),sort([aversiveTS_run(1)-3600 aversiveTS_run(1) ; rewardTS_run(1)-3600 rewardTS_run(1)])),:);
        
        %% Spikes
        % Load Units
        disp('Uploading Spiking activity')
        cd 'Spikesorting'
        [clusters , numberD , numberV , spks , spks_dHPC , spks_vHPC , cellulartype] = load_SU_FM(cd,criteria_type,criteria_fr,aversiveTS_run,rewardTS_run);
        if numberD > 2
            for i = 1 : size(clusters.dHPC,1)
                
                FRB = Restrict(spks_dHPC(spks_dHPC(:,1)==clusters.dHPC(i),2) , WAKE.control);
                FRB = size(FRB,1)/sum([WAKE.control(:,2) - WAKE.control(:,1)]);
                
                FRA = Restrict(spks_dHPC(spks_dHPC(:,1)==clusters.dHPC(i),2) , movement.aversive);
                FRA = size(FRA,1)/sum([movement.aversive(:,2) - movement.aversive(:,1)]);
                
                FRR = Restrict(spks_dHPC(spks_dHPC(:,1)==clusters.dHPC(i),2) , movement.reward);
                FRR = size(FRR,1)/sum([movement.reward(:,2) - movement.reward(:,1)]); 
                
%                 o = (max([FRA ; FRR]) - min([FRA ; FRR]))/nanmean([FRA ; FRR]);
                o = (FRA - FRR)/((FRA + FRR)/2);

                if and(FRA>0 , FRR>0)
                    dHPC.baseline = [dHPC.baseline , FRB];
                    dHPC.aversive = [dHPC.aversive ; FRA];
                    dHPC.reward   = [dHPC.reward ; FRR];
                    dHPC.overlap  = [dHPC.overlap ; o];
                end
                clear FRB FRA FRR o
            end
        end
        
        if numberV > 2
            for i = 1 : size(clusters.vHPC,1)
                
                FRB = Restrict(spks_vHPC(spks_vHPC(:,1)==clusters.vHPC(i),2) , WAKE.control);
                FRB = size(FRB,1)/sum([WAKE.control(:,2) - WAKE.control(:,1)]);                
                
                FRA = Restrict(spks_vHPC(spks_vHPC(:,1)==clusters.vHPC(i),2) , movement.aversive);
                FRA = size(FRA,1)/sum([movement.aversive(:,2) - movement.aversive(:,1)]);
                
                FRR = Restrict(spks_vHPC(spks_vHPC(:,1)==clusters.vHPC(i),2) , movement.reward);
                FRR = size(FRR,1)/sum([movement.reward(:,2) - movement.reward(:,1)]); 
                
%                 o = (max([FRA ; FRR]) - min([FRA ; FRR]))/nanmean([FRA ; FRR]);
                o = (FRA - FRR)/((FRA + FRR)/2);

                if and(FRA>0 , FRR>0)
                    vHPC.baseline = [vHPC.baseline , FRB];
                    vHPC.aversive = [vHPC.aversive ; FRA];
                    vHPC.reward   = [vHPC.reward ; FRR];
                    vHPC.overlap  = [vHPC.overlap ; o];
                end
                clear FRB FRA FRR o
            end
        end
        
        clear aversiveTS aversiveTS_run baselineTS behaviorclusters config 
        clear minimal_speed minimal_speed_time movement REM NREM Shocks_filt
        clear spks spks_dHPC spks_vHPC SpksTrains numberD numberV Rewards_filt
        clear rewardTS rewardTS_run segments behavior ans
    end
end

% dHPC
grps = [ones(size(dHPC.reward,1),1) ; ones(size(dHPC.aversive,1),1)*2];
y = [dHPC.reward./dHPC.baseline' ; dHPC.aversive./dHPC.baseline'];

subplot(131),
scatter(grps,y,[],"filled",'jitter','on', 'jitterAmount',0.1),hold on
scatter([1 2] , [nanmean(dHPC.reward./dHPC.baseline') nanmean(dHPC.aversive./dHPC.baseline')],'filled'),xlim([0 3]),set(gca, 'YScale', 'log')
[p h] = ttest2(dHPC.reward./dHPC.baseline' , dHPC.aversive./dHPC.baseline')

% vHPC
grps = [ones(size(vHPC.reward,1),1) ; ones(size(vHPC.aversive,1),1)*2];
y = [vHPC.reward./vHPC.baseline' ; vHPC.aversive./vHPC.baseline'];

subplot(132),
scatter(grps,y,[],"filled",'jitter','on', 'jitterAmount',0.1),hold on
scatter([1 2] , [nanmean(vHPC.reward./vHPC.baseline') nanmean(vHPC.aversive./vHPC.baseline')],'filled'),xlim([0 3]),set(gca, 'YScale', 'log')
[p h] = ttest2(vHPC.reward./vHPC.baseline' , vHPC.aversive./vHPC.baseline')


% overlap
grps = [ones(size(dHPC.overlap,1),1) ; ones(size(vHPC.overlap,1),1)*2];
y = [dHPC.overlap ; vHPC.overlap];

subplot(133),
vHPC.overlap(vHPC.overlap == Inf) = [];
dHPC.overlap(dHPC.overlap == Inf) = [];

scatter(grps,y,[],"filled",'jitter','on', 'jitterAmount',0.1),hold on
scatter([1 2] , [nanmean(dHPC.overlap) nanmean(vHPC.overlap)],'filled'),xlim([0 3])
[p h] = ttest2(dHPC.overlap , vHPC.overlap)

[h p] = ttest(dHPC.overlap,0)
[h p] = ttest(vHPC.overlap,0)

figure
subplot(121)
scatter(dHPC.aversive , dHPC.reward,'filled'),hold on
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
plot([10^-3 10^2] , [10^-3 10^2])
[p h] = signrank(dHPC.reward , dHPC.aversive)

subplot(122)
scatter(vHPC.aversive , vHPC.reward,'filled'),hold on
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
plot([10^-3 10^2] , [10^-3 10^2])
[p h] = signrank(vHPC.reward , vHPC.aversive)
end