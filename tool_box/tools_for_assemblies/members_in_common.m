function [Aversive Reward] = members_in_common(path)
% This function count the nomber of Joint members in both regions and check
% howmany are also members in other atype of assemblies.
%
% Syntax: [Aversive Reward] = members_in_common(path)
%
% INPUTS
% path: cell, inside each cell you should put the path of each session you
%       want to analyze
%
% OUTPUT
% Aversive, Reward: struct,contain the total number of members and the
%                   number of shared with other assemblies.
%
% Morci Juan Facundo 07/2025

% variables to use in the script
criteria_fr = 0;
criteria_type = 0; % criteria for celltype (0:pyr, 1:int, 2:all)

% Storage
Aversive.dHPC.shared = 0;
Aversive.dHPC.all = 0;
Aversive.vHPC.shared = 0;
Aversive.vHPC.all = 0;

Reward.dHPC.shared = 0;
Reward.dHPC.all = 0;
Reward.vHPC.shared = 0;
Reward.vHPC.all = 0;


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
        load('behavioral_data.mat', 'movement')

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
        aversiveTS_run = aversiveTS_run./1000;
        rewardTS_run = rewardTS_run./1000;
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
                    tmp = sum(tmp,2)>=1;
                    
                    Aversive.dHPC.all = Aversive.dHPC.all + sum(sum(tmp(1:numberD,:)));
                    Aversive.vHPC.all = Aversive.vHPC.all + sum(sum(tmp(numberD+1:end,:)));
                    
                    tmp = Thresholded.aversive(:,not(cond.both));
                    tmp = sum(tmp,2);
                    
                    tmp1 = sum(Thresholded.aversive(:,cond.both),2);
                    
                    tmp2 = sum([tmp>=1, tmp1>=1],2)>1;
                    
                    Aversive.dHPC.shared = Aversive.dHPC.shared + sum(sum(tmp2(1:numberD,:)));
                    Aversive.vHPC.shared = Aversive.vHPC.shared + sum(sum(tmp2(numberD+1:end,:)));
                    
                    clear tmp   
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
                    tmp = sum(tmp,2)>=1;
                    
                    Reward.dHPC.all = Reward.dHPC.all + sum(sum(tmp(1:numberD,:)));
                    Reward.vHPC.all = Reward.vHPC.all + sum(sum(tmp(numberD+1:end,:)));
                    
                    tmp = Thresholded.reward(:,not(cond.both));
                    tmp = sum(tmp,2);
                    
                    tmp1 = sum(Thresholded.reward(:,cond.both),2);
                    
                    tmp2 = sum([tmp>=1, tmp1>=1],2)>1;
                    
                    Reward.dHPC.shared = Reward.dHPC.shared + sum(sum(tmp2(1:numberD,:)));
                    Reward.vHPC.shared = Reward.vHPC.shared + sum(sum(tmp2(numberD+1:end,:)));
                    
                    clear tmp
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
subplot(221),
p1 = (Aversive.dHPC.shared / Aversive.dHPC.all)*100;
p2 = 100 - p1;
pie([p1 , p2] , {'Shared' , 'Rest'})

subplot(222),
p3 = (Aversive.vHPC.shared / Aversive.vHPC.all)*100;
p4 = 100 - p3;
pie([p3 , p4] , {'Shared' , 'Rest'})

subplot(223),
p5 = (Reward.dHPC.shared / Reward.dHPC.all)*100;
p6 = 100 - p5;
pie([p5 , p6] , {'Shared' , 'Rest'})

subplot(224),
p7 = (Reward.vHPC.shared / Reward.vHPC.all)*100;
p8 = 100 - p7;
pie([p7 , p8] , {'Shared' , 'Rest'})


end