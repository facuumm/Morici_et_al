function [Aversive Reward] = countAssemblies_shock_responsive_member(path) 
% This function calculates the Pre and Post sleep assemblies rate.
%
% INPUTS
% path: cell, inside each cell you should put the path of each session you
%       want to analyze
%
% OUTPUT
% Pre, Post: Structure, it store the SU response sourrounding the ripple.
%
% T: vector, it contains time vector for visualization
%
% I: Structure, it contains the tag to determine if the SU is or not
%    shock responssive
%
% curve: Structure, it contains the response of the neuron sourrounding the
%        shock. It is zscored.
%
% Morci Juan Facundo 08/2024

% variables to use in the script
criteria_fr = 0;
criteria_type = 0; % criteria for celltype (0:pyr, 1:int, 2:all)
sm = 0;
dur = 1;
bin = 0.01;

% storage variables
Aversive.dHPC = [];           Reward.dHPC = [];
Aversive.vHPC = [];           Reward.vHPC = [];
Aversive.joint = [];          Reward.joint = [];

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
        disp('Uploading session data')
        load('session_organization.mat')
        baselineTS = baselineTS./1000;
        aversiveTS = aversiveTS./1000;      aversiveTS_run = aversiveTS_run./1000;
        rewardTS = rewardTS./1000;          rewardTS_run = rewardTS_run./1000;
        
        %% Spikes
        % Load Units
        disp('Uploading Spiking activity')
        cd 'Spikesorting'
        [clusters , numberD , numberV , spks , spks_dHPC , spks_vHPC , cellulartype] = load_SU_FM(cd,criteria_type,criteria_fr,aversiveTS_run,rewardTS_run);
        
        %% Load Shock info
        if exist('dHPC_responsivness_all.mat')>=1
            load('dHPC_responsivness_all.mat')
            dHPC = [dHPC_resp.id , dHPC_resp.resp_ave];
            dHPC = dHPC(dHPC(:,2)==1);
        end
        
        if exist('vHPC_responsivness_all.mat')>=1
            load('vHPC_responsivness_all.mat')
            vHPC = [vHPC_resp.id , vHPC_resp.resp_ave];
            vHPC = vHPC(vHPC(:,2)==1);
        end        
        
        
        %% Load assemblies info
        if isfile('dorsalventral_assemblies_aversiveVF.mat')
            disp('Loading Aversive template')
            [patterns , cond , Thresholded] = load_assemblies(cd , 'dorsalventral_assemblies_aversiveVF.mat', clusters, numberD , 'aversive');
            
            cond = sum([cond.dHPC*1 ; cond.vHPC*2 ; cond.both*3]); %1: dHPC, 2:vHPC, 3:both
            
            for i = 1:size(Thresholded.aversive,2)
                ii = clusters.all(Thresholded.aversive(:,i));
                if cond(i)==1
                    c = sum(ismember(ii,dHPC)>=1);
                    Aversive.dHPC = [Aversive.dHPC ; c numberD sum(Thresholded.aversive(:,i))]; clear c
                elseif cond(i)==2
                    c = sum(ismember(ii,vHPC)>=1);
                    Aversive.vHPC = [Aversive.vHPC ; c numberV sum(Thresholded.aversive(:,i))]; clear c
                elseif cond(i)==3
                    c1 = sum(ismember(ii,dHPC)>=1);
                    c2 = sum(ismember(ii,vHPC)>=1);
                    Aversive.joint = [Aversive.joint ; c1 c2 numberD numberV sum(Thresholded.aversive(1:numberD,i)) sum(Thresholded.aversive(numberD+1:end,i))]; clear c1 c2
                end
                
            end
        end
        

        if isfile('dorsalventral_assemblies_rewardVF.mat')
            disp('Loading Reward template')
            [patterns , cond , Thresholded] = load_assemblies(cd , 'dorsalventral_assemblies_rewardVF.mat', clusters, numberD , 'reward');
            
            cond = sum([cond.dHPC*1 ; cond.vHPC*2 ; cond.both*3]); %1: dHPC, 2:vHPC, 3:both
            
            for i = 1:size(Thresholded.reward,2)
                ii = clusters.all(Thresholded.reward(:,i));
                if cond(i)==1
                    c = sum(ismember(ii,dHPC)>=1);
                    Reward.dHPC = [Reward.dHPC ; c numberD sum(Thresholded.reward(:,i))]; clear c
                elseif cond(i)==2
                    c = sum(ismember(ii,vHPC)>=1);
                    Reward.vHPC = [Reward.vHPC ; c numberV sum(Thresholded.reward(:,i))]; clear c
                elseif cond(i)==3
                    c1 = sum(ismember(ii,dHPC)>=1);
                    c2 = sum(ismember(ii,vHPC)>=1);
                    Reward.joint = [Reward.joint ; c1 c2 numberD numberV sum(Thresholded.reward(1:numberD,i)) sum(Thresholded.reward(numberD+1:end,i))]; clear c1 c2
                end
                
            end
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
    clear cooridnated_eventDV cooridnated_eventVD segments movement RD RV RB
end
end
% % To plot data as it is in the paper
% p1 = (sum(Aversive.joint(:,1)>=1)/length(Aversive.joint(:,1)))*100;
% p2 = (sum(Aversive.joint(:,2)>=1)/length(Aversive.joint(:,1)))*100;
% 
% p3 = (sum(Reward.joint(:,1)>=1)/length(Reward.joint(:,1)))*100;
% p4 = (sum(Reward.joint(:,2)>=1)/length(Reward.joint(:,1)))*100;
% 
% figure,bar([1 2],[p1 100-p1]),title('Joint Aversive - dHPC'),ylim([0 100])
% figure,bar([1 2],[p2 100-p2]),title('Joint Aversive - vHPC'),ylim([0 100])
% 
% figure,bar([1 2],[p3 100-p1]),title('Joint Reward - dHPC'),ylim([0 100])
% figure,bar([1 2],[p4 100-p2]),title('Joint Reward - vHPC'),ylim([0 100])
% 
% p1 = (sum(Aversive.dHPC(:,1)>=1)/length(Aversive.dHPC(:,1)))*100;
% figure,bar([1 2],[p1 100-p1]),title('dHPC Aversive'),ylim([0 100])
% p2 = (sum(Reward.dHPC(:,1)>=1)/length(Reward.dHPC(:,1)))*100;
% figure,bar([1 2],[p2 100-p2]),title('dHPC Reward'),ylim([0 100])
% 
% 
% p1 = (sum(Aversive.vHPC(:,1)>=1)/length(Aversive.vHPC(:,1)))*100;
% figure,bar([1 2],[p1 100-p1]),title('vHPC Aversive'),ylim([0 100])
% p2 = (sum(Reward.vHPC(:,1)>=1)/length(Reward.vHPC(:,1)))*100;
% figure,bar([1 2],[p2 100-p2]),title('vHPC Reward'),ylim([0 100])
% 
% 
% percentage = sum(Aversive.joint);