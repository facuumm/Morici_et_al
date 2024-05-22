function [dHPC vHPC Joint] = members_count(path)
% This function counts the number of members detected in each type of
% assemblies.
%
% Syntax: [dHPC vHPC Joint] = members_count(path)
%
% INPUTS
% path: cell, inside each cell you should put the path of each session you
%       want to analyze
%
% OUTPUT
% dHPC, vHPC, Joint: Structure, it store the number of detected members.
%
%               Architecture of each output:
%                   dHPC.aversive
%                       .reward
%
%                   vHPC.aversive
%                       .reward
%
%                   Joint.aversive.dHPC/vHPC
%                        .reward.dHPC/vHPC
%
%
% Morci Juan Facundo 05/2024

% variables to use in the script
criteria_fr = 0;
criteria_type = 0; % criteria for celltype (0:pyr, 1:int, 2:all)

% Storage
dHPC.aversive = [];         dHPC.reward = [];   
vHPC.aversive = [];         vHPC.reward = [];   
Joint.aversive.dHPC = [];   Joint.reward.dHPC = [];   
Joint.aversive.vHPC = [];   Joint.reward.vHPC = [];   


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
        REM.all = ToIntervals(states==5);    NREM.all = ToIntervals(states==3);
        clear x states
        
        baselineTS = baselineTS./1000;
        aversiveTS = aversiveTS./1000;
        rewardTS = rewardTS./1000;
        aversiveTS_run = aversiveTS_run./1000;
        rewardTS_run = rewardTS_run./1000;
        
        NREM.baseline = Restrict(NREM.all,baselineTS);
        NREM.aversive = Restrict(NREM.all,aversiveTS);
        NREM.reward = Restrict(NREM.all,rewardTS);
        
        REM.baseline = Restrict(REM.all,baselineTS);
        REM.aversive = Restrict(REM.all,aversiveTS);
        REM.reward = Restrict(REM.all,rewardTS);

        
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
        if criteria_type == 0 %pyr
            cellulartype = [K(:,1) , K(:,4)];
        elseif criteria_type == 1 % int
            cellulartype = [K(:,1) , not(K(:,4))];
        elseif criteria_type == 2 % all
            cellulartype = [K(:,1) , ones(length(K),1)];
        end
        
        %% Counting the Number f SU
        numberD = 0;
        clusters.all = [];
        clusters.dHPC = [];
        for ii=1:size(group_dHPC,1)
            cluster = group_dHPC(ii,1);
            Cellulartype = logical(cellulartype(cellulartype(:,1) == cluster,2));
            if Cellulartype
                a = length(Restrict(spks(spks(:,1)==cluster,2),aversiveTS_run)) / ((aversiveTS_run(2)-aversiveTS_run(1)));
                r = length(Restrict(spks(spks(:,1)==cluster,2),rewardTS_run)) / ((rewardTS_run(2)-rewardTS_run(1)));
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
                a = length(Restrict(spks(spks(:,1)==cluster,2),aversiveTS_run)) / ((aversiveTS_run(2)-aversiveTS_run(1)));
                r = length(Restrict(spks(spks(:,1)==cluster,2),rewardTS_run)) / ((rewardTS_run(2)-rewardTS_run(1)));
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
        if or(numberD > 3 , numberV > 3)
            
            % --- Aversive ---
            disp('Lets go for the assemblies')
            if isfile('dorsalventral_assemblies_aversiveVF.mat')
                disp('Loading Aversive template')
                load('dorsalventral_assemblies_aversiveVF.mat')
                
                Thresholded.aversive = Th;
                patterns.aversive = pat;
                clear cond Th pat
                
                % Detection of members
                if not(isempty(patterns.aversive))
                    if numberD>0
                        cond1 =  sum(Thresholded.aversive(1:size(clusters.dHPC,1),:),1)>0; %checking of dHPC SU
                        cond2 =  sum(Thresholded.aversive(size(clusters.dHPC,1)+1:end,:),1)>0; %checking of vHPC SU
                        cond.dHPC = and(cond1 , not(cond2));
                        cond.vHPC = and(cond2 , not(cond1));
                        cond.both = and(cond1 , cond2); clear cond1 cond2
                    else
                        cond1 =  logical(zeros(1,size(Thresholded.aversive,2))); %checking of dHPC SU
                        cond2 =  logical(ones(1,size(Thresholded.aversive,2))); %checking of vHPC SU
                        cond.dHPC = and(cond1 , not(cond2));
                        cond.vHPC = and(cond2 , not(cond1));
                        cond.both = and(cond1 , cond2); clear cond1 cond2
                    end
                else
                    cond1 =  false; %checking of dHPC SU
                    cond2 =  logical(0); %checking of vHPC SU
                    cond.dHPC = and(cond1 , not(cond2));
                    cond.vHPC = and(cond2 , not(cond1));
                    cond.both = and(cond1 , cond2); clear cond1 cond2
                end
                
                if sum(cond.both)>0
                    tmp = Thresholded.aversive(:,cond.both);
                    tmpD = [sum(tmp(1:numberD,:))./numberD];
                    tmpV = [sum(tmp(numberD+1:end,:))./numberV];
                    Joint.aversive.dHPC = [Joint.aversive.dHPC ; tmpD'];
                    Joint.aversive.vHPC = [Joint.aversive.vHPC ; tmpV'];
                    clear tmp tmpD tmpV 
                end
                
                if sum(cond.dHPC)>0
                    tmp = Thresholded.aversive(:,cond.dHPC);
                    tmpD = [sum(tmp(1:numberD,:))./numberD];
                    dHPC.aversive = [dHPC.aversive ; tmpD'];
                    clear tmp tmpD tmpV 
                end
                
                if sum(cond.vHPC)>0
                    tmp = Thresholded.aversive(:,cond.vHPC);
                    tmpV = [sum(tmp(numberD+1:end,:))./numberV];
                    vHPC.aversive = [vHPC.aversive ; tmpV'];
                    clear tmp tmpD tmpV 
                end    
                
                clear cond Thresholded patterns
                
             end
            
            
            if isfile('dorsalventral_assemblies_rewardVF.mat')
                disp('Loading Reward template')
                load('dorsalventral_assemblies_rewardVF.mat')
                
                Thresholded.reward = Th;
                patterns.reward = pat;
                clear cond Th pat
                
                % Detection of members
                if not(isempty(patterns.reward))
                    if numberD>0
                        cond1 =  sum(Thresholded.reward(1:size(clusters.dHPC,1),:),1)>0; %checking of dHPC SU
                        cond2 =  sum(Thresholded.reward(size(clusters.dHPC,1)+1:end,:),1)>0; %checking of vHPC SU
                        cond.dHPC = and(cond1 , not(cond2));
                        cond.vHPC = and(cond2 , not(cond1));
                        cond.both = and(cond1 , cond2); clear cond1 cond2
                    else
                        cond1 =  logical(zeros(1,size(Thresholded.reward,2))); %checking of dHPC SU
                        cond2 =  logical(ones(1,size(Thresholded.reward,2))); %checking of vHPC SU
                        cond.dHPC = and(cond1 , not(cond2));
                        cond.vHPC = and(cond2 , not(cond1));
                        cond.both = and(cond1 , cond2); clear cond1 cond2
                    end
                else
                    cond1 =  false; %checking of dHPC SU
                    cond2 =  logical(0); %checking of vHPC SU
                    cond.dHPC = and(cond1 , not(cond2));
                    cond.vHPC = and(cond2 , not(cond1));
                    cond.both = and(cond1 , cond2); clear cond1 cond2
                end
                
                if sum(cond.both)>0
                    tmp = Thresholded.reward(:,cond.both);
                    tmpD = [sum(tmp(1:numberD,:))./numberD];
                    tmpV = [sum(tmp(numberD+1:end,:))./numberV];
                    Joint.reward.dHPC = [Joint.reward.dHPC ; tmpD'];
                    Joint.reward.vHPC = [Joint.reward.vHPC ; tmpV'];
                    clear tmp tmpD tmpV 
                end
                
                if sum(cond.dHPC)>0
                    tmp = Thresholded.reward(:,cond.dHPC);
                    tmpD = [sum(tmp(1:numberD,:))./numberD];
                    dHPC.reward = [dHPC.reward ; tmpD'];
                    clear tmp tmpD tmpV 
                end
                
                if sum(cond.vHPC)>0
                    tmp = Thresholded.reward(:,cond.vHPC);
                    tmpV = [sum(tmp(numberD+1:end,:))./numberV];
                    vHPC.reward = [vHPC.reward ; tmpV'];
                    clear tmp tmpD tmpV 
                end    
                
                clear cond Thresholded patterns                
 
            end
            
            
        end
        disp(' ')
        clear A aversiveTS aversiveTS_run baselineTS rewardTS rewardTS_run
        clear behavior bins Cell_type_classification cellulartype cond
        clear is K Kinfo group_dHPC group_vHPC matrixC matrixD matrixV
        clear NREM REM WAKE segmentation tmp cond config
        clear spiketrains_dHPC spiketrains_vHPC opts MUA
        clear patterns Thresholded i  ii numberD numberV movement cross crossN
        clear Spikes bins Clusters Shocks_filt Rewards_filt config n_SU_D n_SU_V
        clear clusters coordinated coordinated_ripple_bursts coordinatedV
        clear cooridnated_event coordinatedV_refined coordinatedV_refined
        clear ripple_bursts ripple_event ripplesD ripplesV
        clear spks spks_dHPC spks_vHPC ripples cooridnated_event
        clear cooridnated_eventDV cooridnated_eventVD segments movement
    end
    clear num_assembliesA num_assembliesR
    
end

figure
subplot(221)
grps = [ones(length(dHPC.reward),1) ; ones(length(dHPC.aversive),1)*2];
x = [dHPC.reward ; dHPC.aversive];
scatter(grps,x,'filled','Jitter',0.1), hold on
scatter([1 2] , [nanmedian(dHPC.reward) , nanmedian(dHPC.aversive)],'filled'),xlim([0 3]),ylim([0 1])
[h p] = ranksum(dHPC.reward , dHPC.aversive)

subplot(222)
grps = [ones(length(vHPC.reward),1) ; ones(length(vHPC.aversive),1)*2];
x = [vHPC.reward ; vHPC.aversive];
scatter(grps,x,'filled','Jitter',0.1), hold on
scatter([1 2] , [nanmedian(vHPC.reward) , nanmedian(vHPC.aversive)],'filled'),xlim([0 3]),ylim([0 1])
[h p] = ranksum(vHPC.reward , vHPC.aversive)


subplot(223)
grps = [ones(length(Joint.reward.dHPC),1) ; ones(length(Joint.aversive.dHPC),1)*2];
x = [Joint.reward.dHPC ; Joint.aversive.dHPC];
scatter(grps,x,'filled','Jitter',0.1), hold on
scatter([1 2] , [nanmedian(Joint.reward.dHPC) , nanmedian(Joint.aversive.dHPC)],'filled'),xlim([0 3]),ylim([0 1])
[h p] = ranksum(Joint.reward.dHPC , Joint.aversive.dHPC)


subplot(224)
grps = [ones(length(Joint.reward.vHPC),1) ; ones(length(Joint.aversive.vHPC),1)*2];
x = [Joint.reward.vHPC ; Joint.aversive.vHPC];
scatter(grps,x,'filled','Jitter',0.1), hold on
scatter([1 2] , [nanmedian(Joint.reward.vHPC) , nanmedian(Joint.aversive.vHPC)],'filled'),xlim([0 3]),ylim([0 1])
[h p] = ranksum(Joint.reward.vHPC , Joint.aversive.vHPC)

end